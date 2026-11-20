//! Low level writer: sections of a KFF file
//!
//! [Writer] is the equivalent of `Kff_file` opened in "w" mode in the C++ API, and the section
//! writers are the equivalent of `Section_GV`, `Section_Index`, `Section_Raw` and
//! `Section_Minimizer`. The files written are byte for byte the files written by the C++ API with
//! the same calls.

/* std use */
use std::collections::BTreeMap;
use std::io::{Seek, SeekFrom, Write};

/* project use */
use crate::block::{check_buffer, BlockFormat};
use crate::error;
use crate::nucleotides::{bytes_for_nucleotides, remove_nucleotides};
use crate::section::{index_section_size, store_big_endian, SectionType, Values};
use crate::Encoding;

/// Size of the write buffer: the file is written by chunks of this size
const WRITE_BUFFER_CAPACITY: usize = 1 << 20;

/// Size of the footer written by the writer: values section with "first_index" and "footer_size"
const FOOTER_SIZE: u64 = 1 + 8 + (12 + 8) + (12 + 8);

type Reopen<W> = Box<dyn FnMut() -> std::io::Result<W> + Send>;

/// Low level writer of a KFF file.
///
/// The header is written at the creation (signature, version, default encoding A=0 C=1 G=3 T=2,
/// flags set to false). The encoding and the flags can be changed at any time, the metadata must be
/// written before the first section. The sections are then written in the order of the calls.
///
/// By default, the file is indexed as with the C++ API: on [Writer::close], an index section that
/// lists all the sections and a footer (values section with "first_index" and "footer_size") are
/// written before the final signature. Use [Writer::set_indexation] to disable it.
///
/// ```no_run
/// # fn main() -> kff::error::Result<()> {
/// use kff::Writer;
///
/// let mut writer = Writer::create("file.kff")?;
/// writer.set_uniqueness(true)?;
/// writer.write_metadata(b"my comment")?;
///
/// let mut values = writer.values_section()?;
/// values.write_var("k", 5)?;
/// values.write_var("max", 10)?;
/// values.write_var("data_size", 1)?;
/// values.close()?;
///
/// let encoding = writer.encoding();
/// let mut raw = writer.raw_section()?;
/// raw.write_compacted_sequence(&encoding.encode(b"ACGTACG")?, 7, &[1, 2, 3])?;
/// raw.close()?;
///
/// writer.close()?;
/// # Ok(())
/// # }
/// ```
pub struct Writer<W>
where
    W: Write + Seek,
{
    /// None when the file is released
    inner: Option<W>,
    /// Reopen the file after a release
    reopen: Option<Reopen<W>>,
    /// Path of the file when created from a path
    path: Option<std::path::PathBuf>,
    /// Position of the beginning of the file in `inner`
    origin: u64,
    /// Bytes not written in `inner` yet
    buffer: Vec<u8>,
    /// Number of bytes written in `inner`
    flushed: u64,

    header_over: bool,
    encoding: Encoding,
    uniqueness: bool,
    canonicity: bool,

    indexed: bool,
    section_positions: BTreeMap<u64, SectionType>,
    global_vars: Values,
    /// Working space to remove the minimizers
    scratch: Vec<u8>,
    closed: bool,
}

impl Writer<std::fs::File> {
    /// Create a KFF file at `path` (an existing file is overwritten)
    pub fn create<P>(path: P) -> error::Result<Self>
    where
        P: AsRef<std::path::Path>,
    {
        let path = path.as_ref().to_path_buf();
        let mut writer = Self::new(std::fs::File::create(&path)?)?;

        let reopen_path = path.clone();
        writer.reopen = Some(Box::new(move || {
            std::fs::OpenOptions::new().write(true).open(&reopen_path)
        }));
        writer.path = Some(path);

        Ok(writer)
    }
}

impl<W> Writer<W>
where
    W: Write + Seek,
{
    /// Write a KFF file in a seekable stream (a file, a `std::io::Cursor`...), from its current
    /// position.
    pub fn new(mut inner: W) -> error::Result<Self> {
        let origin = inner.stream_position()?;
        let encoding = Encoding::default();

        let mut writer = Self {
            inner: Some(inner),
            reopen: None,
            path: None,
            origin,
            buffer: Vec::new(),
            flushed: 0,
            header_over: false,
            encoding,
            uniqueness: false,
            canonicity: false,
            indexed: true,
            section_positions: BTreeMap::new(),
            global_vars: Values::new(),
            scratch: Vec::new(),
            closed: false,
        };

        writer.write_bytes(&[
            b'K',
            b'F',
            b'F',
            crate::KFF_FORMAT_VERSION_MAJOR,
            crate::KFF_FORMAT_VERSION_MINOR,
            encoding.as_byte(),
            0, // uniqueness
            0, // canonicity
        ])?;

        Ok(writer)
    }

    // ----- Raw access to the file -----

    fn inner(&mut self) -> error::Result<&mut W> {
        if self.inner.is_none() {
            let reopen = self.reopen.as_mut().ok_or(error::Error::CannotRelease)?;
            let mut inner = reopen()?;
            inner.seek(SeekFrom::Start(self.origin + self.flushed))?;
            self.inner = Some(inner);
        }
        Ok(self.inner.as_mut().unwrap())
    }

    fn flush_buffer(&mut self) -> error::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let buffer = std::mem::take(&mut self.buffer);
        let result = self.inner().and_then(|inner| Ok(inner.write_all(&buffer)?));
        self.flushed += buffer.len() as u64;
        self.buffer = buffer;
        self.buffer.clear();
        result
    }

    /// Write bytes at the end of the file. They are kept in a buffer and written by chunks.
    pub fn write_bytes(&mut self, bytes: &[u8]) -> error::Result<()> {
        if bytes.len() >= WRITE_BUFFER_CAPACITY {
            self.flush_buffer()?;
            self.inner()?.write_all(bytes)?;
            self.flushed += bytes.len() as u64;
            return Ok(());
        }
        self.buffer.extend_from_slice(bytes);
        if self.buffer.len() >= WRITE_BUFFER_CAPACITY {
            self.flush_buffer()?;
        }
        Ok(())
    }

    /// Overwrite bytes at a position of the file. The bytes after the end of the file (if any)
    /// are added to the file.
    pub fn write_at(&mut self, bytes: &[u8], position: u64) -> error::Result<()> {
        if position > self.tell() {
            return Err(error::Error::OutOfFile(position));
        }
        let mut bytes = bytes;
        let mut position = position;

        // Part already written in the stream
        if position < self.flushed && !bytes.is_empty() {
            let nb = std::cmp::min(bytes.len() as u64, self.flushed - position) as usize;
            let (origin, flushed) = (self.origin, self.flushed);
            let inner = self.inner()?;
            inner.seek(SeekFrom::Start(origin + position))?;
            inner.write_all(&bytes[..nb])?;
            inner.seek(SeekFrom::Start(origin + flushed))?;
            bytes = &bytes[nb..];
            position += nb as u64;
        }

        // Part in the buffer
        if !bytes.is_empty() {
            let offset = (position - self.flushed) as usize;
            let nb = std::cmp::min(bytes.len(), self.buffer.len() - offset);
            self.buffer[offset..offset + nb].copy_from_slice(&bytes[..nb]);
            if nb < bytes.len() {
                self.write_bytes(&bytes[nb..])?;
            }
        }

        Ok(())
    }

    /// Current position (size of the file written so far)
    pub fn tell(&self) -> u64 {
        self.flushed + self.buffer.len() as u64
    }

    /// Release the file: the buffer is written and the file is closed, it is reopened at the next
    /// write. Useful when a lot of files are written simultaneously (`tmp_close` in the C++ API).
    ///
    /// Only for the writers created with [Writer::create].
    pub fn release(&mut self) -> error::Result<()> {
        if self.reopen.is_none() {
            return Err(error::Error::CannotRelease);
        }
        self.flush_buffer()?;
        if let Some(mut inner) = self.inner.take() {
            inner.flush()?;
        }
        Ok(())
    }

    /// True if the file is currently released
    pub fn is_released(&self) -> bool {
        self.inner.is_none()
    }

    // ----- Header -----

    /// Set the encoding of the nucleotides
    pub fn set_encoding(&mut self, encoding: Encoding) -> error::Result<()> {
        self.encoding = encoding;
        self.write_at(&[encoding.as_byte()], 5)
    }

    /// Encoding of the nucleotides
    pub fn encoding(&self) -> Encoding {
        self.encoding
    }

    /// Set the uniqueness flag: true means that no kmer is present twice in the file
    pub fn set_uniqueness(&mut self, uniqueness: bool) -> error::Result<()> {
        self.uniqueness = uniqueness;
        self.write_at(&[uniqueness as u8], 6)
    }

    /// Uniqueness flag
    pub fn uniqueness(&self) -> bool {
        self.uniqueness
    }

    /// Set the canonicity flag: true means that if a kmer is in the file, its reverse complement
    /// is absent
    pub fn set_canonicity(&mut self, canonicity: bool) -> error::Result<()> {
        self.canonicity = canonicity;
        self.write_at(&[canonicity as u8], 7)
    }

    /// Canonicity flag
    pub fn canonicity(&self) -> bool {
        self.canonicity
    }

    /// Write the metadata (free block) of the header. Must be called before the first section;
    /// without call, an empty metadata is written.
    pub fn write_metadata(&mut self, metadata: &[u8]) -> error::Result<()> {
        if self.header_over {
            return Err(error::Error::MetadataAfterSection);
        }
        let size = u32::try_from(metadata.len()).map_err(|_| error::Error::BufferTooSmall {
            needed: metadata.len(),
            available: u32::MAX as usize,
        })?;
        self.write_bytes(&size.to_be_bytes())?;
        self.write_bytes(metadata)?;
        self.header_over = true;
        Ok(())
    }

    /// Complete the header before the first section
    fn complete_header(&mut self) -> error::Result<()> {
        if !self.header_over {
            self.write_metadata(&[])?;
        }
        Ok(())
    }

    // ----- Index -----

    /// Enable or disable the index and footer written at the closing (enabled by default)
    pub fn set_indexation(&mut self, indexed: bool) {
        self.indexed = indexed;
    }

    /// True if an index and a footer are written at the closing
    pub fn is_indexed(&self) -> bool {
        self.indexed
    }

    /// Register a section at the current position in the index written at the closing (the
    /// section writers do it automatically)
    pub fn register_position(&mut self, section_type: SectionType) {
        if self.indexed {
            let position = self.tell();
            self.section_positions.insert(position, section_type);
        }
    }

    /// Sections registered for the index: absolute position and type
    pub fn section_positions(&self) -> &BTreeMap<u64, SectionType> {
        &self.section_positions
    }

    /// Variables of the last values section (the current scope)
    pub fn global_vars(&self) -> &Values {
        &self.global_vars
    }

    /// Value of a variable of the current scope
    pub fn var(&self, name: &str) -> Option<u64> {
        self.global_vars.get(name).copied()
    }

    // ----- Sections -----

    /// Open a values section ('v'). It starts a new scope: the variables of the previous values
    /// section are forgotten.
    pub fn values_section(&mut self) -> error::Result<ValuesSectionWriter<'_, W>> {
        self.complete_header()?;
        self.global_vars.clear();
        self.register_position(SectionType::Values);
        let beginning = self.tell();

        Ok(ValuesSectionWriter {
            writer: self,
            beginning,
            vars: Values::new(),
            closed: false,
        })
    }

    /// Write a values section containing `values`
    pub fn write_values(&mut self, values: &Values) -> error::Result<()> {
        let mut section = self.values_section()?;
        for (name, value) in values {
            section.write_var(name, *value)?;
        }
        section.close()
    }

    /// Open an index section ('i')
    pub fn index_section(&mut self) -> error::Result<IndexSectionWriter<'_, W>> {
        self.complete_header()?;
        let beginning = self.tell();

        Ok(IndexSectionWriter {
            writer: self,
            beginning,
            entries: BTreeMap::new(),
            next_index: 0,
            closed: false,
        })
    }

    /// Open a raw section ('r'). The variables k, max and data_size must be defined in the current
    /// values section.
    pub fn raw_section(&mut self) -> error::Result<RawSectionWriter<'_, W>> {
        let format = BlockFormat::raw(&self.global_vars)?;
        self.complete_header()?;
        self.register_position(SectionType::Raw);
        let beginning = self.tell();
        self.write_bytes(&[b'r', 0, 0, 0, 0, 0, 0, 0, 0])?;

        Ok(RawSectionWriter {
            writer: self,
            format,
            beginning,
            nb_blocks: 0,
            closed: false,
        })
    }

    /// Open a minimizer section ('m'). The variables k, m, max and data_size must be defined in
    /// the current values section.
    pub fn minimizer_section(&mut self) -> error::Result<MinimizerSectionWriter<'_, W>> {
        let format = BlockFormat::minimizer(&self.global_vars)?;
        self.complete_header()?;
        self.register_position(SectionType::Minimizer);
        let beginning = self.tell();
        self.write_bytes(b"m")?;
        self.write_bytes(&vec![0; format.mini_bytes])?;
        self.write_bytes(&[0; 8])?;

        Ok(MinimizerSectionWriter {
            writer: self,
            format,
            beginning,
            nb_blocks: 0,
            closed: false,
        })
    }

    // ----- Closing -----

    fn write_footer(&mut self) -> error::Result<()> {
        let positions = std::mem::take(&mut self.section_positions);

        let mut index = self.index_section()?;
        let end = index.beginning() + index_section_size(positions.len() as u64);
        for (position, section_type) in positions {
            index.register_section(section_type, position as i64 - end as i64);
        }
        let first_index = index.beginning();
        index.close()?;

        let mut footer = self.values_section()?;
        footer.write_var("first_index", first_index)?;
        footer.write_var("footer_size", FOOTER_SIZE)?;
        footer.close()
    }

    fn finish(&mut self) -> error::Result<()> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;

        self.complete_header()?;
        if self.indexed {
            self.write_footer()?;
        }
        self.write_bytes(b"KFF")?;
        self.flush_buffer()?;
        self.inner()?.flush()?;
        Ok(())
    }

    /// Write the index and the footer (if indexed) and the final signature, then return the
    /// underlying stream
    pub fn close(mut self) -> error::Result<W> {
        self.finish()?;
        Ok(self.inner.take().unwrap())
    }

    /// Stop the writing without completing the file. A file created with [Writer::create] is
    /// deleted (`close(false)` in the C++ API).
    pub fn discard(mut self) -> error::Result<()> {
        self.closed = true;
        self.inner = None;
        if let Some(path) = self.path.take() {
            std::fs::remove_file(path)?;
        }
        Ok(())
    }
}

impl<W> Drop for Writer<W>
where
    W: Write + Seek,
{
    /// A writer that is not closed is closed when dropped, errors are ignored
    fn drop(&mut self) {
        let _ = self.finish();
    }
}

impl<W> std::fmt::Debug for Writer<W>
where
    W: Write + Seek,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Writer")
            .field("path", &self.path)
            .field("position", &self.tell())
            .field("encoding", &self.encoding)
            .field("indexed", &self.indexed)
            .field("global_vars", &self.global_vars)
            .finish()
    }
}

/// Writer of a values section ('v'), equivalent of `Section_GV` in the C++ API.
///
/// The variables are written when the section is closed (or dropped), sorted by name.
pub struct ValuesSectionWriter<'a, W>
where
    W: Write + Seek,
{
    writer: &'a mut Writer<W>,
    beginning: u64,
    vars: Values,
    closed: bool,
}

impl<'a, W> ValuesSectionWriter<'a, W>
where
    W: Write + Seek,
{
    /// Position of the section in the file
    pub fn beginning(&self) -> u64 {
        self.beginning
    }

    /// Variables of the section
    pub fn vars(&self) -> &Values {
        &self.vars
    }

    /// Declare a variable. A variable declared twice is written once, with its last value.
    pub fn write_var(&mut self, name: &str, value: u64) -> error::Result<()> {
        if name.as_bytes().contains(&0) {
            return Err(error::Error::InvalidVariableName(name.to_string()));
        }
        self.vars.insert(name.to_string(), value);
        self.writer.global_vars.insert(name.to_string(), value);
        Ok(())
    }

    fn finish(&mut self) -> error::Result<()> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;

        let mut bytes = vec![b'v'];
        bytes.extend_from_slice(&(self.vars.len() as u64).to_be_bytes());
        for (name, value) in &self.vars {
            bytes.extend_from_slice(name.as_bytes());
            bytes.push(0);
            bytes.extend_from_slice(&value.to_be_bytes());
        }
        self.writer.write_bytes(&bytes)
    }

    /// Write the section
    pub fn close(mut self) -> error::Result<()> {
        self.finish()
    }
}

impl<'a, W> Drop for ValuesSectionWriter<'a, W>
where
    W: Write + Seek,
{
    fn drop(&mut self) {
        let _ = self.finish();
    }
}

/// Writer of an index section ('i'), equivalent of `Section_Index` in the C++ API.
///
/// The section is written when it is closed (or dropped).
pub struct IndexSectionWriter<'a, W>
where
    W: Write + Seek,
{
    writer: &'a mut Writer<W>,
    beginning: u64,
    entries: BTreeMap<i64, SectionType>,
    next_index: i64,
    closed: bool,
}

impl<'a, W> IndexSectionWriter<'a, W>
where
    W: Write + Seek,
{
    /// Position of the section in the file
    pub fn beginning(&self) -> u64 {
        self.beginning
    }

    /// Register a section at a position relative to the end of this index section (a section
    /// registered twice at the same position keeps its last type)
    pub fn register_section(&mut self, section_type: SectionType, relative_position: i64) {
        self.entries.insert(relative_position, section_type);
    }

    /// Set the position of the next index section, relative to the end of this one (0 for none)
    pub fn set_next_index(&mut self, relative_position: i64) {
        self.next_index = relative_position;
    }

    fn finish(&mut self) -> error::Result<()> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;

        let mut bytes = Vec::with_capacity(index_section_size(self.entries.len() as u64) as usize);
        bytes.push(b'i');
        bytes.extend_from_slice(&(self.entries.len() as u64).to_be_bytes());
        for (position, section_type) in &self.entries {
            bytes.push(section_type.as_byte());
            bytes.extend_from_slice(&position.to_be_bytes());
        }
        bytes.extend_from_slice(&self.next_index.to_be_bytes());
        self.writer.write_bytes(&bytes)
    }

    /// Write the section
    pub fn close(mut self) -> error::Result<()> {
        self.finish()
    }
}

impl<'a, W> Drop for IndexSectionWriter<'a, W>
where
    W: Write + Seek,
{
    fn drop(&mut self) {
        let _ = self.finish();
    }
}

/// Writer of a raw section ('r'), equivalent of `Section_Raw` in the C++ API.
///
/// The number of blocks is written in the section when it is closed (or dropped).
pub struct RawSectionWriter<'a, W>
where
    W: Write + Seek,
{
    writer: &'a mut Writer<W>,
    format: BlockFormat,
    beginning: u64,
    nb_blocks: u64,
    closed: bool,
}

impl<'a, W> RawSectionWriter<'a, W>
where
    W: Write + Seek,
{
    /// Position of the section in the file
    pub fn beginning(&self) -> u64 {
        self.beginning
    }

    /// Number of blocks written
    pub fn nb_blocks(&self) -> u64 {
        self.nb_blocks
    }

    /// Size of the kmers
    pub fn k(&self) -> u64 {
        self.format.k
    }

    /// Maximal number of kmers per block
    pub fn max(&self) -> u64 {
        self.format.max
    }

    /// Size of the data of a kmer in bytes
    pub fn data_size(&self) -> u64 {
        self.format.data_size
    }

    /// Write a block: a sequence of `seq_size` nucleotides (right aligned 2-bit nucleotides, see
    /// [crate::encoding]) and the data of its `seq_size - k + 1` kmers.
    pub fn write_compacted_sequence(
        &mut self,
        seq: &[u8],
        seq_size: u64,
        data: &[u8],
    ) -> error::Result<()> {
        if seq_size < self.format.k {
            return Err(error::Error::InvalidBlock(format!(
                "a sequence of {} nucleotides contains no {}-mer",
                seq_size, self.format.k
            )));
        }
        let nb_kmers = seq_size - self.format.k + 1;
        self.format.check_nb_kmers(nb_kmers)?;
        let seq_bytes = bytes_for_nucleotides(seq_size);
        let data_bytes = self.format.data_bytes(nb_kmers);
        check_buffer(seq_bytes, seq.len())?;
        check_buffer(data_bytes, data.len())?;

        let nb_kmers_bytes = self.format.nb_kmers_bytes;
        self.writer
            .write_bytes(&store_big_endian(nb_kmers, nb_kmers_bytes)[..nb_kmers_bytes])?;
        self.writer.write_bytes(&seq[..seq_bytes])?;
        self.writer.write_bytes(&data[..data_bytes])?;
        self.nb_blocks += 1;
        Ok(())
    }

    fn finish(&mut self) -> error::Result<()> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;
        self.writer
            .write_at(&self.nb_blocks.to_be_bytes(), self.beginning + 1)
    }

    /// Write the number of blocks of the section
    pub fn close(mut self) -> error::Result<()> {
        self.finish()
    }
}

impl<'a, W> Drop for RawSectionWriter<'a, W>
where
    W: Write + Seek,
{
    fn drop(&mut self) {
        let _ = self.finish();
    }
}

/// Writer of a minimizer section ('m'), equivalent of `Section_Minimizer` in the C++ API.
///
/// The number of blocks is written in the section when it is closed (or dropped).
pub struct MinimizerSectionWriter<'a, W>
where
    W: Write + Seek,
{
    writer: &'a mut Writer<W>,
    format: BlockFormat,
    beginning: u64,
    nb_blocks: u64,
    closed: bool,
}

impl<'a, W> MinimizerSectionWriter<'a, W>
where
    W: Write + Seek,
{
    /// Position of the section in the file
    pub fn beginning(&self) -> u64 {
        self.beginning
    }

    /// Number of blocks written
    pub fn nb_blocks(&self) -> u64 {
        self.nb_blocks
    }

    /// Size of the kmers
    pub fn k(&self) -> u64 {
        self.format.k
    }

    /// Size of the minimizer
    pub fn m(&self) -> u64 {
        self.format.m
    }

    /// Maximal number of kmers per block
    pub fn max(&self) -> u64 {
        self.format.max
    }

    /// Size of the data of a kmer in bytes
    pub fn data_size(&self) -> u64 {
        self.format.data_size
    }

    /// Write the minimizer of the section (`m` right aligned 2-bit nucleotides). Can be called at
    /// any time before the closing.
    pub fn write_minimizer(&mut self, minimizer: &[u8]) -> error::Result<()> {
        check_buffer(self.format.mini_bytes, minimizer.len())?;
        self.writer
            .write_at(&minimizer[..self.format.mini_bytes], self.beginning + 1)
    }

    /// Write a block whose sequence is given without the minimizer: `seq` contains the
    /// `seq_size` nucleotides that are not part of the minimizer, and the minimizer is at position
    /// `mini_pos` of the full sequence.
    pub fn write_compacted_sequence_without_mini(
        &mut self,
        seq: &[u8],
        seq_size: u64,
        mini_pos: u64,
        data: &[u8],
    ) -> error::Result<()> {
        if seq_size + self.format.m < self.format.k {
            return Err(error::Error::InvalidBlock(format!(
                "a sequence of {} nucleotides contains no {}-mer",
                seq_size + self.format.m,
                self.format.k
            )));
        }
        let nb_kmers = seq_size + self.format.m - self.format.k + 1;
        self.format.check_nb_kmers(nb_kmers)?;
        self.format.check_mini_pos(nb_kmers, mini_pos)?;
        let seq_bytes = bytes_for_nucleotides(seq_size);
        let data_bytes = self.format.data_bytes(nb_kmers);
        check_buffer(seq_bytes, seq.len())?;
        check_buffer(data_bytes, data.len())?;

        let nb_kmers_bytes = self.format.nb_kmers_bytes;
        let mini_pos_bytes = self.format.mini_pos_bytes;
        self.writer
            .write_bytes(&store_big_endian(nb_kmers, nb_kmers_bytes)[..nb_kmers_bytes])?;
        self.writer
            .write_bytes(&store_big_endian(mini_pos, mini_pos_bytes)[..mini_pos_bytes])?;
        self.writer.write_bytes(&seq[..seq_bytes])?;
        self.writer.write_bytes(&data[..data_bytes])?;
        self.nb_blocks += 1;
        Ok(())
    }

    /// Write a block from its full sequence of `seq_size` nucleotides, where the minimizer starts
    /// at position `mini_pos`. The minimizer is removed from the sequence before the writing.
    pub fn write_compacted_sequence(
        &mut self,
        seq: &[u8],
        seq_size: u64,
        mini_pos: u64,
        data: &[u8],
    ) -> error::Result<()> {
        let m = self.format.m;
        if seq_size < self.format.k || mini_pos + m > seq_size {
            return Err(error::Error::InvalidBlock(format!(
                "minimizer at position {} in a sequence of {} nucleotides (k = {}, m = {})",
                mini_pos, seq_size, self.format.k, m
            )));
        }
        check_buffer(bytes_for_nucleotides(seq_size), seq.len())?;

        let mut scratch = std::mem::take(&mut self.writer.scratch);
        scratch.clear();
        scratch.resize(bytes_for_nucleotides(seq_size - m), 0);
        remove_nucleotides(&mut scratch, seq, seq_size, mini_pos, m);
        let result =
            self.write_compacted_sequence_without_mini(&scratch, seq_size - m, mini_pos, data);
        self.writer.scratch = scratch;
        result
    }

    fn finish(&mut self) -> error::Result<()> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;
        self.writer.write_at(
            &self.nb_blocks.to_be_bytes(),
            self.beginning + 1 + self.format.mini_bytes as u64,
        )
    }

    /// Write the number of blocks of the section
    pub fn close(mut self) -> error::Result<()> {
        self.finish()
    }
}

impl<'a, W> Drop for MinimizerSectionWriter<'a, W>
where
    W: Write + Seek,
{
    fn drop(&mut self) {
        let _ = self.finish();
    }
}
