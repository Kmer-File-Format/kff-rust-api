//! Low level reader: sections of a KFF file
//!
//! [Reader] is the equivalent of `Kff_file` opened in "r" mode in the C++ API, and the section
//! readers are the equivalent of `Section_GV`, `Section_Index`, `Section_Raw` and
//! `Section_Minimizer`.

/* std use */
use std::io::{Read, Seek, SeekFrom};

/* project use */
use crate::block::{check_buffer, Block, BlockFormat};
use crate::error;
use crate::nucleotides::{bytes_for_nucleotides, insert_nucleotides};
use crate::section::{load_big_endian, IndexSection, SectionType, Values};
use crate::{Encoding, Header, Writer};

/// Size of the read cache. Moderate because a reader can be created for each input file of a merge.
const READ_CACHE_CAPACITY: usize = 1 << 16;

/// Low level reader of a KFF file.
///
/// At the opening, the header is read, the footer and the index sections (if any) are loaded and
/// the reader is placed on the first section. The sections are then read in the order of the file
/// with [Reader::read_section_type] and the functions that read or open each type of section.
///
/// ```no_run
/// # fn main() -> kff::error::Result<()> {
/// use kff::{BlockReader, Reader, SectionType};
///
/// let mut reader = Reader::open("file.kff")?;
/// while let Some(section_type) = reader.read_section_type()? {
///     match section_type {
///         SectionType::Values => {
///             let values = reader.read_values()?;
///             println!("{:?}", values);
///         }
///         SectionType::Raw | SectionType::Minimizer => {
///             let mut section = reader.open_block_section()?.unwrap();
///             while let Some(block) = section.read_block()? {
///                 println!("{} kmers", block.nb_kmers);
///             }
///         }
///         _ => {
///             reader.skip_section()?;
///         }
///     }
/// }
/// # Ok(())
/// # }
/// ```
pub struct Reader<R> {
    inner: R,
    /// Position of the cursor of `inner`
    inner_position: u64,
    /// Cache of the file
    cache: Vec<u8>,
    cache_start: u64,
    cache_len: usize,
    /// Current position in the file
    position: u64,
    file_size: u64,

    header: Header,
    first_section: u64,
    global_vars: Values,
    footer: Option<Values>,
    footer_position: Option<u64>,
    index: Vec<IndexSection>,
    /// Working space to insert the minimizers
    scratch: Vec<u8>,
}

impl Reader<std::fs::File> {
    /// Open the KFF file at `path`
    pub fn open<P>(path: P) -> error::Result<Self>
    where
        P: AsRef<std::path::Path>,
    {
        Self::new(std::fs::File::open(path)?)
    }
}

impl<R> Reader<R>
where
    R: Read + Seek,
{
    /// Read a KFF file from a seekable source (a file, a `std::io::Cursor`...).
    ///
    /// The header is read and checked (signatures, version, encoding), then the footer and the
    /// index sections are discovered.
    pub fn new(mut inner: R) -> error::Result<Self> {
        let file_size = inner.seek(SeekFrom::End(0))?;

        let mut reader = Self {
            inner,
            inner_position: file_size,
            cache: Vec::new(),
            cache_start: 0,
            cache_len: 0,
            position: 0,
            file_size,
            header: Header::default(),
            first_section: 0,
            global_vars: Values::new(),
            footer: None,
            footer_position: None,
            index: Vec::new(),
            scratch: Vec::new(),
        };

        reader.read_header()?;
        reader.discover_footer()?;
        reader.discover_index()?;
        reader.position = reader.first_section;

        Ok(reader)
    }

    fn read_header(&mut self) -> error::Result<()> {
        if self.file_size < 3 || &self.read_array::<3>()? != b"KFF" {
            return Err(error::Error::MissingSignature("beginning"));
        }

        let [major_version, minor_version, encoding, uniqueness, canonicity] =
            self.read_array::<5>()?;
        if crate::KFF_FORMAT_VERSION_MAJOR < major_version
            || (crate::KFF_FORMAT_VERSION_MAJOR == major_version
                && crate::KFF_FORMAT_VERSION_MINOR < minor_version)
        {
            return Err(error::Error::UnsupportedVersion {
                major: major_version,
                minor: minor_version,
            });
        }
        let encoding = Encoding::from_byte(encoding)?;

        let metadata_size = self.read_u64_be(4)?;
        if metadata_size > self.file_size - self.position {
            return Err(error::Error::MissingSignature("end"));
        }
        let mut metadata = vec![0; metadata_size as usize];
        self.read_bytes(&mut metadata)?;

        self.header = Header {
            major_version,
            minor_version,
            encoding,
            uniqueness: uniqueness != 0,
            canonicity: canonicity != 0,
            metadata,
        };
        self.first_section = self.position;

        // Signature at the end of the file
        if self.file_size < self.first_section + 3 {
            return Err(error::Error::MissingSignature("end"));
        }
        self.position = self.file_size - 3;
        if &self.read_array::<3>()? != b"KFF" {
            return Err(error::Error::MissingSignature("end"));
        }

        Ok(())
    }

    /// The footer is a values section that ends with the variable "footer_size"
    fn discover_footer(&mut self) -> error::Result<()> {
        if self.file_size < 23 || self.file_size - 23 < self.first_section {
            return Ok(());
        }

        self.position = self.file_size - 23;
        if &self.read_array::<11>()? != b"footer_size" {
            return Ok(());
        }
        self.position += 1; // '\0'
        let size = self.read_u64_be(8)?;

        let start = (self.file_size - 3).checked_sub(size);
        match start {
            Some(start) if start >= self.first_section && size > 0 => {
                self.position = start;
                let footer = self.read_values_section()?;
                self.footer = Some(footer);
                self.footer_position = Some(start);
                Ok(())
            }
            _ => Err(error::Error::InvalidFooter(format!(
                "footer_size {} is out of the file",
                size
            ))),
        }
    }

    /// The first index section is given by the footer ("first_index") or is the first section
    fn discover_index(&mut self) -> error::Result<()> {
        let first_index = match self.footer.as_ref().and_then(|f| f.get("first_index")) {
            Some(position) => Some(*position),
            None => {
                self.position = self.first_section;
                if self.peek_byte()? == Some(b'i') {
                    Some(self.first_section)
                } else {
                    None
                }
            }
        };

        let mut position = match first_index {
            Some(position) => position,
            None => return Ok(()),
        };
        let mut visited = std::collections::HashSet::new();
        loop {
            if position < self.first_section || position >= self.end_position() {
                return Err(error::Error::InvalidFooter(format!(
                    "index section at position {} out of the sections of the file",
                    position
                )));
            }
            if !visited.insert(position) {
                return Err(error::Error::InvalidFooter(
                    "the chain of index sections contains a loop".to_string(),
                ));
            }

            self.position = position;
            let index = self.read_index_section()?;
            let next = index.next_index;
            let end = index.end();
            self.index.push(index);

            if next == 0 {
                return Ok(());
            }
            position = (end as i64 + next) as u64;
        }
    }

    // ----- Raw access to the file -----

    /// Read `bytes.len()` bytes at the current position
    pub fn read_bytes(&mut self, bytes: &mut [u8]) -> error::Result<()> {
        let size = bytes.len();
        if self.position + size as u64 > self.file_size {
            return Err(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                format!("read out of the file, byte {}", self.file_size),
            )
            .into());
        }

        // Fast path: all the bytes are in the cache
        if self.position >= self.cache_start
            && self.position + size as u64 <= self.cache_start + self.cache_len as u64
        {
            let offset = (self.position - self.cache_start) as usize;
            bytes.copy_from_slice(&self.cache[offset..offset + size]);
            self.position += size as u64;
            return Ok(());
        }

        let mut done = 0;
        while done < size {
            let cache_end = self.cache_start + self.cache_len as u64;
            if self.position >= self.cache_start && self.position < cache_end {
                // Copy what is available in the cache
                let offset = (self.position - self.cache_start) as usize;
                let nb = std::cmp::min(size - done, self.cache_len - offset);
                bytes[done..done + nb].copy_from_slice(&self.cache[offset..offset + nb]);
                done += nb;
                self.position += nb as u64;
            } else if size - done >= READ_CACHE_CAPACITY {
                // Large reads bypass the cache
                Self::disk_read(
                    &mut self.inner,
                    &mut self.inner_position,
                    &mut bytes[done..],
                    self.position,
                )?;
                self.position += (size - done) as u64;
                done = size;
            } else {
                self.fill_cache()?;
            }
        }

        Ok(())
    }

    fn disk_read(
        inner: &mut R,
        inner_position: &mut u64,
        bytes: &mut [u8],
        position: u64,
    ) -> error::Result<()> {
        if *inner_position != position {
            inner.seek(SeekFrom::Start(position))?;
        }
        // The position of the stream is unknown if the read fails
        *inner_position = u64::MAX;
        inner.read_exact(bytes)?;
        *inner_position = position + bytes.len() as u64;
        Ok(())
    }

    /// Load the cache from the current position
    fn fill_cache(&mut self) -> error::Result<()> {
        if self.cache.len() < READ_CACHE_CAPACITY {
            self.cache.resize(READ_CACHE_CAPACITY, 0);
        }
        let nb = std::cmp::min(READ_CACHE_CAPACITY as u64, self.file_size - self.position) as usize;
        self.cache_len = 0;
        Self::disk_read(
            &mut self.inner,
            &mut self.inner_position,
            &mut self.cache[..nb],
            self.position,
        )?;
        self.cache_start = self.position;
        self.cache_len = nb;
        Ok(())
    }

    /// Byte at the current position without moving, None at the end of the file
    fn peek_byte(&mut self) -> error::Result<Option<u8>> {
        if self.position >= self.file_size {
            return Ok(None);
        }
        if self.position < self.cache_start
            || self.position >= self.cache_start + self.cache_len as u64
        {
            self.fill_cache()?;
        }
        Ok(Some(
            self.cache[(self.position - self.cache_start) as usize],
        ))
    }

    fn read_array<const N: usize>(&mut self) -> error::Result<[u8; N]> {
        let mut bytes = [0; N];
        self.read_bytes(&mut bytes)?;
        Ok(bytes)
    }

    fn read_u8(&mut self) -> error::Result<u8> {
        Ok(self.read_array::<1>()?[0])
    }

    /// Read a big endian value of `nb_bytes` bytes (at most 8)
    fn read_u64_be(&mut self, nb_bytes: usize) -> error::Result<u64> {
        let mut bytes = [0; 8];
        self.read_bytes(&mut bytes[..nb_bytes])?;
        Ok(load_big_endian(&bytes[..nb_bytes]))
    }

    /// Current position in the file
    pub fn tell(&self) -> u64 {
        self.position
    }

    /// Relative jump (can be negative)
    pub fn jump(&mut self, offset: i64) -> error::Result<()> {
        let position = self.position as i64 + offset;
        if position < 0 {
            return Err(error::Error::OutOfFile(0));
        }
        self.jump_to(position as u64)
    }

    /// Jump to a position of the file
    pub fn jump_to(&mut self, position: u64) -> error::Result<()> {
        if position > self.file_size {
            return Err(error::Error::OutOfFile(position));
        }
        self.position = position;
        Ok(())
    }

    // ----- File properties -----

    /// Header of the file
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Encoding of the nucleotides
    pub fn encoding(&self) -> Encoding {
        self.header.encoding
    }

    /// True if each kmer is present at most once in the file
    pub fn uniqueness(&self) -> bool {
        self.header.uniqueness
    }

    /// True if a kmer and its reverse complement are never both in the file
    pub fn canonicity(&self) -> bool {
        self.header.canonicity
    }

    /// Metadata (free block) of the header
    pub fn metadata(&self) -> &[u8] {
        &self.header.metadata
    }

    /// Size of the file in bytes
    pub fn file_size(&self) -> u64 {
        self.file_size
    }

    /// Position of the first section (just after the header)
    pub fn first_section_position(&self) -> u64 {
        self.first_section
    }

    /// Position of the final signature: there is no section after it
    pub fn end_position(&self) -> u64 {
        self.file_size - 3
    }

    /// Variables of the footer, if the file has one
    pub fn footer(&self) -> Option<&Values> {
        self.footer.as_ref()
    }

    /// Position of the footer, if the file has one
    pub fn footer_position(&self) -> Option<u64> {
        self.footer_position
    }

    /// Index sections of the file, in the order of their chain
    pub fn index(&self) -> &[IndexSection] {
        &self.index
    }

    /// True if the file contains index sections
    pub fn is_indexed(&self) -> bool {
        !self.index.is_empty()
    }

    /// Absolute positions and types of all the sections registered in the index sections
    pub fn indexed_sections(&self) -> std::collections::BTreeMap<u64, SectionType> {
        self.index
            .iter()
            .flat_map(|index| index.absolute_entries())
            .collect()
    }

    /// Variables of the last values section read (the current scope)
    pub fn global_vars(&self) -> &Values {
        &self.global_vars
    }

    /// Value of a variable of the current scope
    pub fn var(&self, name: &str) -> Option<u64> {
        self.global_vars.get(name).copied()
    }

    /// Release the underlying reader
    pub fn into_inner(self) -> R {
        self.inner
    }

    // ----- Sections -----

    /// Type of the next section, without moving. None at the end of the file (final signature).
    pub fn read_section_type(&mut self) -> error::Result<Option<SectionType>> {
        if self.position >= self.end_position() {
            return Ok(None);
        }
        Ok(self.peek_byte()?.map(SectionType::from_byte))
    }

    fn expect_section(&mut self, expected: SectionType) -> error::Result<()> {
        let found = self.read_u8()?;
        if found != expected.as_byte() {
            self.position -= 1;
            return Err(error::Error::UnexpectedSection {
                expected: expected.as_byte() as char,
                found,
            });
        }
        Ok(())
    }

    fn read_values_section(&mut self) -> error::Result<Values> {
        self.expect_section(SectionType::Values)?;
        let nb_vars = self.read_u64_be(8)?;

        let mut values = Values::new();
        let mut name = Vec::new();
        for _ in 0..nb_vars {
            if self.position >= self.end_position() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    "end of the file reached before the end of the values section",
                )
                .into());
            }
            name.clear();
            loop {
                let c = self.read_u8()?;
                if c == 0 {
                    break;
                }
                name.push(c);
            }
            let name = String::from_utf8(name.clone()).map_err(|_| {
                error::Error::InvalidVariableName(String::from_utf8_lossy(&name).into_owned())
            })?;
            let value = self.read_u64_be(8)?;
            values.insert(name, value);
        }

        Ok(values)
    }

    /// Read a values section ('v'). Its variables become the current scope.
    pub fn read_values(&mut self) -> error::Result<Values> {
        let values = self.read_values_section()?;
        self.global_vars = values.clone();
        Ok(values)
    }

    fn read_index_section(&mut self) -> error::Result<IndexSection> {
        let beginning = self.position;
        self.expect_section(SectionType::Index)?;
        let nb_entries = self.read_u64_be(8)?;

        let mut index = IndexSection {
            beginning,
            ..Default::default()
        };
        for _ in 0..nb_entries {
            let entry = self.read_array::<9>()?;
            let position = load_big_endian(&entry[1..]) as i64;
            if index
                .entries
                .insert(position, SectionType::from_byte(entry[0]))
                .is_some()
            {
                return Err(error::Error::IndexCollision);
            }
        }
        index.next_index = self.read_u64_be(8)? as i64;

        Ok(index)
    }

    /// Read an index section ('i')
    pub fn read_index(&mut self) -> error::Result<IndexSection> {
        self.read_index_section()
    }

    /// Read the beginning of a raw or minimizer section: returns the format of its blocks, its
    /// minimizer (empty for a raw section) and its number of blocks
    pub(crate) fn read_block_section_header(
        &mut self,
        minimizer: bool,
    ) -> error::Result<(BlockFormat, Vec<u8>, u64)> {
        let format = if minimizer {
            BlockFormat::minimizer(&self.global_vars)?
        } else {
            BlockFormat::raw(&self.global_vars)?
        };
        self.atomic(|reader| {
            reader.expect_section(if minimizer {
                SectionType::Minimizer
            } else {
                SectionType::Raw
            })?;
            let mut mini = vec![0; format.mini_bytes];
            reader.read_bytes(&mut mini)?;
            let nb_blocks = reader.read_u64_be(8)?;
            Ok((format, mini, nb_blocks))
        })
    }

    /// Open a raw section ('r') with the variables of the current scope (k, max and data_size)
    pub fn open_raw(&mut self) -> error::Result<RawSectionReader<'_, R>> {
        let beginning = self.position;
        let (format, _, nb_blocks) = self.read_block_section_header(false)?;

        Ok(RawSectionReader {
            reader: self,
            format,
            beginning,
            nb_blocks,
            remaining_blocks: nb_blocks,
        })
    }

    /// Open a minimizer section ('m') with the variables of the current scope (k, m, max and
    /// data_size)
    pub fn open_minimizer(&mut self) -> error::Result<MinimizerSectionReader<'_, R>> {
        let beginning = self.position;
        let (format, minimizer, nb_blocks) = self.read_block_section_header(true)?;

        Ok(MinimizerSectionReader {
            reader: self,
            format,
            beginning,
            nb_blocks,
            remaining_blocks: nb_blocks,
            minimizer,
        })
    }

    /// Open the next section if it contains kmers (raw or minimizer section), None otherwise.
    ///
    /// Equivalent of `Block_section_reader::construct_section` in the C++ API.
    pub fn open_block_section(&mut self) -> error::Result<Option<BlockSectionReader<'_, R>>> {
        match self.read_section_type()? {
            Some(SectionType::Raw) => Ok(Some(BlockSectionReader::Raw(self.open_raw()?))),
            Some(SectionType::Minimizer) => {
                Ok(Some(BlockSectionReader::Minimizer(self.open_minimizer()?)))
            }
            _ => Ok(None),
        }
    }

    /// Read or open the next section, whatever its type. None at the end of the file.
    ///
    /// Equivalent of `SectionBuilder::build` in the C++ API.
    pub fn read_section(&mut self) -> error::Result<Option<Section<'_, R>>> {
        match self.read_section_type()? {
            None => Ok(None),
            Some(SectionType::Values) => Ok(Some(Section::Values(self.read_values()?))),
            Some(SectionType::Index) => Ok(Some(Section::Index(self.read_index()?))),
            Some(SectionType::Raw) => Ok(Some(Section::Raw(self.open_raw()?))),
            Some(SectionType::Minimizer) => Ok(Some(Section::Minimizer(self.open_minimizer()?))),
            Some(SectionType::Other(byte)) => Err(error::Error::UnknownSection(byte)),
        }
    }

    /// Jump over the next section if it contains kmers (raw or minimizer section).
    ///
    /// Returns false if the next section is not a raw or minimizer section, or at the end of the
    /// file. As in the C++ API, values sections are not skipped (they change the scope).
    pub fn jump_next_section(&mut self) -> error::Result<bool> {
        match self.open_block_section()? {
            Some(mut section) => {
                section.jump_section()?;
                Ok(true)
            }
            None => Ok(false),
        }
    }

    /// Go over the next section, whatever its type: the values sections are read (they change the
    /// scope), the other sections are jumped. Returns the type of the section, None at the end of
    /// the file.
    pub fn skip_section(&mut self) -> error::Result<Option<SectionType>> {
        let section_type = self.read_section_type()?;
        match section_type {
            None => {}
            Some(SectionType::Values) => {
                self.read_values()?;
            }
            Some(SectionType::Index) => {
                self.position += 1;
                let nb_entries = self.read_u64_be(8)?;
                self.jump((9 * nb_entries + 8) as i64)?;
            }
            Some(SectionType::Raw) | Some(SectionType::Minimizer) => {
                self.jump_next_section()?;
            }
            Some(SectionType::Other(byte)) => return Err(error::Error::UnknownSection(byte)),
        }
        Ok(section_type)
    }

    // ----- Blocks -----

    /// Read the number of kmers and the position of the minimizer of a block
    fn read_block_head(&mut self, format: &BlockFormat) -> error::Result<(u64, u64)> {
        let nb_kmers = if format.nb_kmers_bytes == 0 {
            1
        } else {
            self.read_u64_be(format.nb_kmers_bytes)?
        };
        format.check_nb_kmers(nb_kmers)?;

        let mini_pos = if format.minimizer {
            let mini_pos = self.read_u64_be(format.mini_pos_bytes)?;
            format.check_mini_pos(nb_kmers, mini_pos)?;
            mini_pos
        } else {
            0
        };

        Ok((nb_kmers, mini_pos))
    }

    /// Run a block reading function, the position is restored if it fails
    fn atomic<T, F>(&mut self, read: F) -> error::Result<T>
    where
        F: FnOnce(&mut Self) -> error::Result<T>,
    {
        let position = self.position;
        let result = read(self);
        if result.is_err() {
            self.position = position;
        }
        result
    }

    /// Read a block as stored in the file: the sequence without the minimizer for the minimizer
    /// sections. Returns the number of kmers and the position of the minimizer.
    fn read_stored_block(
        &mut self,
        format: &BlockFormat,
        seq: &mut [u8],
        data: &mut [u8],
    ) -> error::Result<(u64, u64)> {
        self.atomic(|reader| {
            let (nb_kmers, mini_pos) = reader.read_block_head(format)?;
            let seq_bytes = bytes_for_nucleotides(format.stored_size(nb_kmers));
            let data_bytes = format.data_bytes(nb_kmers);
            check_buffer(seq_bytes, seq.len())?;
            check_buffer(data_bytes, data.len())?;
            reader.read_bytes(&mut seq[..seq_bytes])?;
            reader.read_bytes(&mut data[..data_bytes])?;
            Ok((nb_kmers, mini_pos))
        })
    }

    /// Same as `read_stored_block` in vectors resized to the block
    fn read_stored_block_vec(
        &mut self,
        format: &BlockFormat,
        seq: &mut Vec<u8>,
        data: &mut Vec<u8>,
    ) -> error::Result<(u64, u64)> {
        self.atomic(|reader| {
            let (nb_kmers, mini_pos) = reader.read_block_head(format)?;
            seq.resize(bytes_for_nucleotides(format.stored_size(nb_kmers)), 0);
            data.resize(format.data_bytes(nb_kmers), 0);
            reader.read_bytes(seq)?;
            reader.read_bytes(data)?;
            Ok((nb_kmers, mini_pos))
        })
    }

    /// Same as `read_stored_block` with the data just after the sequence
    fn read_stored_block_joined(
        &mut self,
        format: &BlockFormat,
        seq_data: &mut [u8],
    ) -> error::Result<(u64, u64)> {
        self.atomic(|reader| {
            let (nb_kmers, mini_pos) = reader.read_block_head(format)?;
            let size =
                bytes_for_nucleotides(format.stored_size(nb_kmers)) + format.data_bytes(nb_kmers);
            check_buffer(size, seq_data.len())?;
            reader.read_bytes(&mut seq_data[..size])?;
            Ok((nb_kmers, mini_pos))
        })
    }

    /// Read a block with its full sequence (the minimizer is inserted back)
    fn read_full_block(
        &mut self,
        format: &BlockFormat,
        minimizer: &[u8],
        seq: &mut [u8],
        data: &mut [u8],
    ) -> error::Result<u64> {
        if !format.minimizer {
            return Ok(self.read_stored_block(format, seq, data)?.0);
        }

        let mut scratch = std::mem::take(&mut self.scratch);
        let result = self.atomic(|reader| {
            let (nb_kmers, mini_pos) = reader.read_block_head(format)?;
            let stored = format.stored_size(nb_kmers);
            let seq_bytes = bytes_for_nucleotides(format.sequence_size(nb_kmers));
            let data_bytes = format.data_bytes(nb_kmers);
            check_buffer(seq_bytes, seq.len())?;
            check_buffer(data_bytes, data.len())?;

            scratch.clear();
            scratch.resize(bytes_for_nucleotides(stored), 0);
            reader.read_bytes(&mut scratch)?;
            reader.read_bytes(&mut data[..data_bytes])?;

            seq[..seq_bytes].fill(0);
            insert_nucleotides(seq, &scratch, stored, mini_pos, minimizer, format.m);
            Ok(nb_kmers)
        });
        self.scratch = scratch;
        result
    }

    /// Same as `read_full_block` with the data just after the sequence
    fn read_full_block_joined(
        &mut self,
        format: &BlockFormat,
        minimizer: &[u8],
        seq_data: &mut [u8],
    ) -> error::Result<u64> {
        if !format.minimizer {
            return Ok(self.read_stored_block_joined(format, seq_data)?.0);
        }

        let mut scratch = std::mem::take(&mut self.scratch);
        let result = self.atomic(|reader| {
            let (nb_kmers, mini_pos) = reader.read_block_head(format)?;
            let stored = format.stored_size(nb_kmers);
            let seq_bytes = bytes_for_nucleotides(format.sequence_size(nb_kmers));
            let data_bytes = format.data_bytes(nb_kmers);
            check_buffer(seq_bytes + data_bytes, seq_data.len())?;

            let stored_bytes = bytes_for_nucleotides(stored);
            scratch.clear();
            scratch.resize(stored_bytes, 0);
            reader.read_bytes(&mut scratch)?;
            reader.read_bytes(&mut seq_data[seq_bytes..seq_bytes + data_bytes])?;

            seq_data[..seq_bytes].fill(0);
            insert_nucleotides(seq_data, &scratch, stored, mini_pos, minimizer, format.m);
            Ok(nb_kmers)
        });
        self.scratch = scratch;
        result
    }

    /// Read a block into an owned block, its buffers are reused
    pub(crate) fn read_block_into(
        &mut self,
        format: &BlockFormat,
        minimizer: &[u8],
        block: &mut Block,
    ) -> error::Result<()> {
        let mut scratch = std::mem::take(&mut self.scratch);
        let result = self.atomic(|reader| {
            let (nb_kmers, mini_pos) = reader.read_block_head(format)?;
            let stored = format.stored_size(nb_kmers);
            let seq_bytes = bytes_for_nucleotides(format.sequence_size(nb_kmers));

            block.k = format.k;
            block.nb_kmers = nb_kmers;
            block.data_size = format.data_size;
            block.sequence.clear();
            block.sequence.resize(seq_bytes, 0);
            block.data.clear();
            block.data.resize(format.data_bytes(nb_kmers), 0);

            if format.minimizer {
                scratch.clear();
                scratch.resize(bytes_for_nucleotides(stored), 0);
                reader.read_bytes(&mut scratch)?;
                insert_nucleotides(
                    &mut block.sequence,
                    &scratch,
                    stored,
                    mini_pos,
                    minimizer,
                    format.m,
                );
            } else {
                reader.read_bytes(&mut block.sequence)?;
            }
            reader.read_bytes(&mut block.data)?;
            Ok(())
        });
        self.scratch = scratch;
        result
    }

    fn jump_block(&mut self, format: &BlockFormat) -> error::Result<()> {
        self.atomic(|reader| {
            let nb_kmers = if format.nb_kmers_bytes == 0 {
                1
            } else {
                reader.read_u64_be(format.nb_kmers_bytes)?
            };
            format.check_nb_kmers(nb_kmers)?;
            let size = format.mini_pos_bytes
                + bytes_for_nucleotides(format.stored_size(nb_kmers))
                + format.data_bytes(nb_kmers);
            reader.jump(size as i64)
        })
    }
}

impl<T> std::fmt::Debug for Reader<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Reader")
            .field("header", &self.header)
            .field("position", &self.position)
            .field("file_size", &self.file_size)
            .field("global_vars", &self.global_vars)
            .field("footer", &self.footer)
            .finish()
    }
}

/// A section returned by [Reader::read_section]
pub enum Section<'a, R>
where
    R: Read + Seek,
{
    /// Variables of a values section
    Values(Values),
    /// Content of an index section
    Index(IndexSection),
    /// An opened raw section
    Raw(RawSectionReader<'a, R>),
    /// An opened minimizer section
    Minimizer(MinimizerSectionReader<'a, R>),
}

/// Reading of the blocks of a raw or minimizer section.
///
/// Equivalent of `Block_section_reader` in the C++ API. The sequences are right aligned 2-bit
/// nucleotides (see [crate::encoding]).
pub trait BlockReader {
    /// Type of the section
    fn section_type(&self) -> SectionType;
    /// Position of the section in the file
    fn beginning(&self) -> u64;
    /// Size of the kmers
    fn k(&self) -> u64;
    /// Maximal number of kmers in a block
    fn max(&self) -> u64;
    /// Size of the data of a kmer in bytes
    fn data_size(&self) -> u64;
    /// Number of blocks of the section
    fn nb_blocks(&self) -> u64;
    /// Number of blocks not read yet
    fn remaining_blocks(&self) -> u64;

    /// Size of a buffer that can contain the sequence of any block of the section
    fn max_sequence_bytes(&self) -> usize {
        bytes_for_nucleotides(self.k() + self.max() - 1)
    }

    /// Size of a buffer that can contain the data of any block of the section
    fn max_data_bytes(&self) -> usize {
        (self.max() * self.data_size()) as usize
    }

    /// Read the next block: its sequence in `seq` and the data of its kmers in `data`. Returns
    /// the number of kmers of the block.
    ///
    /// The buffers must be large enough for the block ([BlockReader::max_sequence_bytes] and
    /// [BlockReader::max_data_bytes] are always enough), otherwise
    /// [error::Error::BufferTooSmall] is returned and the block is not consumed.
    fn read_compacted_sequence(&mut self, seq: &mut [u8], data: &mut [u8]) -> error::Result<u64>;

    /// Read the next block in one buffer: the sequence immediately followed by the data. Returns
    /// the number of kmers of the block.
    fn read_compacted_sequence_joined(&mut self, seq_data: &mut [u8]) -> error::Result<u64>;

    /// Read the next block into `block`, whose buffers are reused. Returns false if the section
    /// has no remaining block.
    fn read_block_into(&mut self, block: &mut Block) -> error::Result<bool>;

    /// Read the next block, None if the section has no remaining block
    fn read_block(&mut self) -> error::Result<Option<Block>> {
        let mut block = Block::default();
        Ok(if self.read_block_into(&mut block)? {
            Some(block)
        } else {
            None
        })
    }

    /// Jump over the next block
    fn jump_sequence(&mut self) -> error::Result<()>;

    /// Jump over all the remaining blocks of the section
    fn jump_section(&mut self) -> error::Result<()> {
        while self.remaining_blocks() > 0 {
            self.jump_sequence()?;
        }
        Ok(())
    }
}

fn no_remaining_block() -> error::Error {
    error::Error::InvalidBlock("the section has no remaining block".to_string())
}

/// Reader of a raw section ('r'), equivalent of `Section_Raw` in the C++ API.
///
/// When the reader is closed or dropped, the remaining blocks are jumped: the [Reader] is then
/// at the beginning of the next section.
pub struct RawSectionReader<'a, R>
where
    R: Read + Seek,
{
    reader: &'a mut Reader<R>,
    format: BlockFormat,
    beginning: u64,
    nb_blocks: u64,
    remaining_blocks: u64,
}

impl<'a, R> RawSectionReader<'a, R>
where
    R: Read + Seek,
{
    /// Copy the remaining blocks into a new raw section of `writer`. The variables of the
    /// current values section of the writer are used for the new section.
    pub fn copy_to<W>(&mut self, writer: &mut Writer<W>) -> error::Result<()>
    where
        W: std::io::Write + Seek,
    {
        let mut section = writer.raw_section()?;
        let mut block = Block::default();
        while self.read_block_into(&mut block)? {
            section.write_compacted_sequence(
                &block.sequence,
                block.sequence_size(),
                &block.data,
            )?;
        }
        section.close()
    }

    /// Jump over the remaining blocks
    pub fn close(mut self) -> error::Result<()> {
        self.jump_section()
    }
}

impl<'a, R> BlockReader for RawSectionReader<'a, R>
where
    R: Read + Seek,
{
    fn section_type(&self) -> SectionType {
        SectionType::Raw
    }
    fn beginning(&self) -> u64 {
        self.beginning
    }
    fn k(&self) -> u64 {
        self.format.k
    }
    fn max(&self) -> u64 {
        self.format.max
    }
    fn data_size(&self) -> u64 {
        self.format.data_size
    }
    fn nb_blocks(&self) -> u64 {
        self.nb_blocks
    }
    fn remaining_blocks(&self) -> u64 {
        self.remaining_blocks
    }

    fn read_compacted_sequence(&mut self, seq: &mut [u8], data: &mut [u8]) -> error::Result<u64> {
        if self.remaining_blocks == 0 {
            return Err(no_remaining_block());
        }
        let nb_kmers = self.reader.read_stored_block(&self.format, seq, data)?.0;
        self.remaining_blocks -= 1;
        Ok(nb_kmers)
    }

    fn read_compacted_sequence_joined(&mut self, seq_data: &mut [u8]) -> error::Result<u64> {
        if self.remaining_blocks == 0 {
            return Err(no_remaining_block());
        }
        let nb_kmers = self
            .reader
            .read_stored_block_joined(&self.format, seq_data)?
            .0;
        self.remaining_blocks -= 1;
        Ok(nb_kmers)
    }

    fn read_block_into(&mut self, block: &mut Block) -> error::Result<bool> {
        if self.remaining_blocks == 0 {
            return Ok(false);
        }
        self.reader.read_block_into(&self.format, &[], block)?;
        self.remaining_blocks -= 1;
        Ok(true)
    }

    fn jump_sequence(&mut self) -> error::Result<()> {
        if self.remaining_blocks == 0 {
            return Err(no_remaining_block());
        }
        self.reader.jump_block(&self.format)?;
        self.remaining_blocks -= 1;
        Ok(())
    }
}

impl<'a, R> Drop for RawSectionReader<'a, R>
where
    R: Read + Seek,
{
    fn drop(&mut self) {
        while self.remaining_blocks > 0 && self.jump_sequence().is_ok() {}
    }
}

/// Reader of a minimizer section ('m'), equivalent of `Section_Minimizer` in the C++ API.
///
/// The blocks can be read with their full sequence ([BlockReader] functions) or as they are
/// stored, without the minimizer
/// ([MinimizerSectionReader::read_compacted_sequence_without_mini]).
///
/// When the reader is closed or dropped, the remaining blocks are jumped: the [Reader] is then
/// at the beginning of the next section.
pub struct MinimizerSectionReader<'a, R>
where
    R: Read + Seek,
{
    reader: &'a mut Reader<R>,
    format: BlockFormat,
    beginning: u64,
    nb_blocks: u64,
    remaining_blocks: u64,
    minimizer: Vec<u8>,
}

impl<'a, R> MinimizerSectionReader<'a, R>
where
    R: Read + Seek,
{
    /// Size of the minimizer
    pub fn m(&self) -> u64 {
        self.format.m
    }

    /// Minimizer of the section (right aligned 2-bit nucleotides)
    pub fn minimizer(&self) -> &[u8] {
        &self.minimizer
    }

    /// Read the next block as stored in the file: its sequence without the minimizer in `seq` and
    /// the data in `data`. Returns the number of kmers and the position of the minimizer in the
    /// sequence.
    pub fn read_compacted_sequence_without_mini(
        &mut self,
        seq: &mut [u8],
        data: &mut [u8],
    ) -> error::Result<(u64, u64)> {
        if self.remaining_blocks == 0 {
            return Err(no_remaining_block());
        }
        let result = self.reader.read_stored_block(&self.format, seq, data)?;
        self.remaining_blocks -= 1;
        Ok(result)
    }

    /// Same as [MinimizerSectionReader::read_compacted_sequence_without_mini] with the data just
    /// after the sequence in the same buffer
    pub fn read_compacted_sequence_without_mini_joined(
        &mut self,
        seq_data: &mut [u8],
    ) -> error::Result<(u64, u64)> {
        if self.remaining_blocks == 0 {
            return Err(no_remaining_block());
        }
        let result = self
            .reader
            .read_stored_block_joined(&self.format, seq_data)?;
        self.remaining_blocks -= 1;
        Ok(result)
    }

    /// Insert the minimizer of the section in a sequence read without it.
    ///
    /// `seq` contains the `nb_kmers + k - 1 - m` nucleotides of a block of `nb_kmers` kmers and
    /// must have room for `nb_kmers + k - 1` nucleotides. The result has null padding bits.
    pub fn add_minimizer(
        &mut self,
        nb_kmers: u64,
        seq: &mut [u8],
        mini_pos: u64,
    ) -> error::Result<()> {
        self.format.check_nb_kmers(nb_kmers)?;
        self.format.check_mini_pos(nb_kmers, mini_pos)?;
        let stored = self.format.stored_size(nb_kmers);
        let stored_bytes = bytes_for_nucleotides(stored);
        let seq_bytes = bytes_for_nucleotides(self.format.sequence_size(nb_kmers));
        check_buffer(seq_bytes, seq.len())?;

        let without = seq[..stored_bytes].to_vec();
        seq[..seq_bytes].fill(0);
        insert_nucleotides(
            seq,
            &without,
            stored,
            mini_pos,
            &self.minimizer,
            self.format.m,
        );
        Ok(())
    }

    /// Copy the minimizer and the remaining blocks into a new minimizer section of `writer`. The
    /// variables of the current values section of the writer are used for the new section.
    pub fn copy_to<W>(&mut self, writer: &mut Writer<W>) -> error::Result<()>
    where
        W: std::io::Write + Seek,
    {
        let mut section = writer.minimizer_section()?;
        section.write_minimizer(&self.minimizer)?;

        let mut seq = Vec::new();
        let mut data = Vec::new();
        while self.remaining_blocks > 0 {
            let (nb_kmers, mini_pos) =
                self.reader
                    .read_stored_block_vec(&self.format, &mut seq, &mut data)?;
            self.remaining_blocks -= 1;
            section.write_compacted_sequence_without_mini(
                &seq,
                self.format.stored_size(nb_kmers),
                mini_pos,
                &data,
            )?;
        }
        section.close()
    }

    /// Jump over the remaining blocks
    pub fn close(mut self) -> error::Result<()> {
        self.jump_section()
    }
}

impl<'a, R> BlockReader for MinimizerSectionReader<'a, R>
where
    R: Read + Seek,
{
    fn section_type(&self) -> SectionType {
        SectionType::Minimizer
    }
    fn beginning(&self) -> u64 {
        self.beginning
    }
    fn k(&self) -> u64 {
        self.format.k
    }
    fn max(&self) -> u64 {
        self.format.max
    }
    fn data_size(&self) -> u64 {
        self.format.data_size
    }
    fn nb_blocks(&self) -> u64 {
        self.nb_blocks
    }
    fn remaining_blocks(&self) -> u64 {
        self.remaining_blocks
    }

    fn read_compacted_sequence(&mut self, seq: &mut [u8], data: &mut [u8]) -> error::Result<u64> {
        if self.remaining_blocks == 0 {
            return Err(no_remaining_block());
        }
        let nb_kmers = self
            .reader
            .read_full_block(&self.format, &self.minimizer, seq, data)?;
        self.remaining_blocks -= 1;
        Ok(nb_kmers)
    }

    fn read_compacted_sequence_joined(&mut self, seq_data: &mut [u8]) -> error::Result<u64> {
        if self.remaining_blocks == 0 {
            return Err(no_remaining_block());
        }
        let nb_kmers =
            self.reader
                .read_full_block_joined(&self.format, &self.minimizer, seq_data)?;
        self.remaining_blocks -= 1;
        Ok(nb_kmers)
    }

    fn read_block_into(&mut self, block: &mut Block) -> error::Result<bool> {
        if self.remaining_blocks == 0 {
            return Ok(false);
        }
        self.reader
            .read_block_into(&self.format, &self.minimizer, block)?;
        self.remaining_blocks -= 1;
        Ok(true)
    }

    fn jump_sequence(&mut self) -> error::Result<()> {
        if self.remaining_blocks == 0 {
            return Err(no_remaining_block());
        }
        self.reader.jump_block(&self.format)?;
        self.remaining_blocks -= 1;
        Ok(())
    }
}

impl<'a, R> Drop for MinimizerSectionReader<'a, R>
where
    R: Read + Seek,
{
    fn drop(&mut self) {
        while self.remaining_blocks > 0 && self.jump_sequence().is_ok() {}
    }
}

/// A raw or a minimizer section reader, returned by [Reader::open_block_section]
pub enum BlockSectionReader<'a, R>
where
    R: Read + Seek,
{
    /// Raw section
    Raw(RawSectionReader<'a, R>),
    /// Minimizer section
    Minimizer(MinimizerSectionReader<'a, R>),
}

impl<'a, R> BlockSectionReader<'a, R>
where
    R: Read + Seek,
{
    fn inner(&self) -> &dyn BlockReader {
        match self {
            Self::Raw(section) => section,
            Self::Minimizer(section) => section,
        }
    }

    fn inner_mut(&mut self) -> &mut dyn BlockReader {
        match self {
            Self::Raw(section) => section,
            Self::Minimizer(section) => section,
        }
    }

    /// Copy the remaining blocks into a new section of the same type in `writer`
    pub fn copy_to<W>(&mut self, writer: &mut Writer<W>) -> error::Result<()>
    where
        W: std::io::Write + Seek,
    {
        match self {
            Self::Raw(section) => section.copy_to(writer),
            Self::Minimizer(section) => section.copy_to(writer),
        }
    }

    /// Jump over the remaining blocks
    pub fn close(mut self) -> error::Result<()> {
        self.jump_section()
    }
}

impl<'a, R> BlockReader for BlockSectionReader<'a, R>
where
    R: Read + Seek,
{
    fn section_type(&self) -> SectionType {
        self.inner().section_type()
    }
    fn beginning(&self) -> u64 {
        self.inner().beginning()
    }
    fn k(&self) -> u64 {
        self.inner().k()
    }
    fn max(&self) -> u64 {
        self.inner().max()
    }
    fn data_size(&self) -> u64 {
        self.inner().data_size()
    }
    fn nb_blocks(&self) -> u64 {
        self.inner().nb_blocks()
    }
    fn remaining_blocks(&self) -> u64 {
        self.inner().remaining_blocks()
    }
    fn read_compacted_sequence(&mut self, seq: &mut [u8], data: &mut [u8]) -> error::Result<u64> {
        self.inner_mut().read_compacted_sequence(seq, data)
    }
    fn read_compacted_sequence_joined(&mut self, seq_data: &mut [u8]) -> error::Result<u64> {
        self.inner_mut().read_compacted_sequence_joined(seq_data)
    }
    fn read_block_into(&mut self, block: &mut Block) -> error::Result<bool> {
        self.inner_mut().read_block_into(block)
    }
    fn jump_sequence(&mut self) -> error::Result<()> {
        self.inner_mut().jump_sequence()
    }
}
