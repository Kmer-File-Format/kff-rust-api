//! High level reader: enumerate the kmers of a KFF file
//!
//! [KmerReader] is the equivalent of `Kff_reader` in the C++ API: it hides the sections and gives
//! the kmers (or the blocks) of the file with their data.

/* std use */
use std::io::{Read, Seek};

/* project use */
use crate::block::{Block, BlockFormat};
use crate::error;
use crate::nucleotides::{bytes_for_nucleotides, padding};
use crate::{Encoding, Reader, SectionType};

/// A kmer with its data, returned by [KmerReader::kmers]
#[derive(Clone, std::fmt::Debug, PartialEq, Eq, Hash)]
pub struct Kmer {
    /// Size of the kmer
    pub k: u64,
    /// Nucleotides of the kmer: right aligned 2-bit nucleotides with null padding bits
    pub sequence: Vec<u8>,
    /// Data of the kmer
    pub data: Vec<u8>,
    /// Encoding of the nucleotides
    pub encoding: Encoding,
}

impl Kmer {
    /// Build a kmer from its nucleotides (A, C, G, T)
    pub fn from_ascii(sequence: &[u8], data: Vec<u8>, encoding: Encoding) -> error::Result<Self> {
        Ok(Self {
            k: sequence.len() as u64,
            sequence: encoding.encode(sequence)?,
            data,
            encoding,
        })
    }

    /// Nucleotides of the kmer (A, C, G, T)
    pub fn to_ascii(&self) -> Vec<u8> {
        self.encoding.decode(&self.sequence, self.k as usize)
    }
}

impl std::fmt::Display for Kmer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(&self.to_ascii()))
    }
}

/// A kmer returned by [KmerReader::next_kmer], borrowed from the reader
#[derive(Clone, Copy, std::fmt::Debug, PartialEq, Eq)]
pub struct KmerView<'a> {
    /// Size of the kmer
    pub k: u64,
    /// Nucleotides of the kmer: right aligned 2-bit nucleotides with null padding bits
    pub sequence: &'a [u8],
    /// Data of the kmer
    pub data: &'a [u8],
    /// Encoding of the nucleotides
    pub encoding: Encoding,
}

impl<'a> KmerView<'a> {
    /// Nucleotides of the kmer (A, C, G, T)
    pub fn to_ascii(&self) -> Vec<u8> {
        self.encoding.decode(self.sequence, self.k as usize)
    }

    /// Owned copy of the kmer
    pub fn to_kmer(&self) -> Kmer {
        Kmer {
            k: self.k,
            sequence: self.sequence.to_vec(),
            data: self.data.to_vec(),
            encoding: self.encoding,
        }
    }
}

/// Raw or minimizer section being read
struct CurrentSection {
    format: BlockFormat,
    minimizer: Vec<u8>,
    remaining_blocks: u64,
}

/// High level reader of a KFF file: enumerates the kmers or the blocks of all the sections.
///
/// ```no_run
/// # fn main() -> kff::error::Result<()> {
/// let mut reader = kff::KmerReader::open("file.kff")?;
/// while let Some(kmer) = reader.next_kmer()? {
///     println!("{} {:?}", String::from_utf8_lossy(&kmer.to_ascii()), kmer.data);
/// }
/// # Ok(())
/// # }
/// ```
pub struct KmerReader<R>
where
    R: Read + Seek,
{
    reader: Reader<R>,
    section: Option<CurrentSection>,
    block: Block,
    /// Number of kmers of the current block already given
    kmers_read: u64,
    /// Sequence of the current block shifted by 1, 2 and 3 nucleotides (index 0 unused): each
    /// kmer is then a copy of whole bytes in one of the 4 alignments (as in the C++ API)
    shifts: [Vec<u8>; 4],
    /// Bit s set when the shift s of the current block is computed
    shifts_ready: u8,
    kmer: Vec<u8>,
    k: u64,
    max: u64,
    data_size: u64,
}

impl KmerReader<std::fs::File> {
    /// Open the KFF file at `path`
    pub fn open<P>(path: P) -> error::Result<Self>
    where
        P: AsRef<std::path::Path>,
    {
        Ok(Self::new(Reader::open(path)?))
    }
}

impl<R> KmerReader<R>
where
    R: Read + Seek,
{
    /// Enumerate the kmers of a file opened with the low level reader. The enumeration starts at
    /// the current position of the reader.
    pub fn new(reader: Reader<R>) -> Self {
        let vars = reader.global_vars();
        let k = vars.get("k").copied().unwrap_or(0);
        let max = vars.get("max").copied().unwrap_or(0);
        let data_size = vars.get("data_size").copied().unwrap_or(0);

        Self {
            reader,
            section: None,
            block: Block::default(),
            kmers_read: 0,
            shifts: Default::default(),
            shifts_ready: 0,
            kmer: Vec::new(),
            k,
            max,
            data_size,
        }
    }

    /// Encoding of the nucleotides
    pub fn encoding(&self) -> Encoding {
        self.reader.encoding()
    }

    /// Size of the kmers of the current section
    pub fn k(&self) -> u64 {
        self.k
    }

    /// Maximal number of kmers per block of the current section
    pub fn max(&self) -> u64 {
        self.max
    }

    /// Size of the data of a kmer of the current section
    pub fn data_size(&self) -> u64 {
        self.data_size
    }

    /// Value of a variable of the current values section (`get_var` in the C++ API)
    pub fn var(&self, name: &str) -> Option<u64> {
        self.reader.var(name)
    }

    /// The low level reader (header, footer, index...)
    pub fn reader(&self) -> &Reader<R> {
        &self.reader
    }

    /// Release the low level reader
    pub fn into_inner(self) -> Reader<R> {
        self.reader
    }

    fn block_has_kmers(&self) -> bool {
        self.kmers_read < self.block.nb_kmers
    }

    /// Go to the next section that contains blocks. Returns false at the end of the file.
    fn next_section(&mut self) -> error::Result<bool> {
        loop {
            match self.reader.read_section_type()? {
                None => return Ok(false),
                Some(SectionType::Values) => {
                    let values = self.reader.read_values()?;
                    // As in the C++ API, the properties are updated by the sections that define them
                    if let Some(k) = values.get("k") {
                        self.k = *k;
                    }
                    if let Some(max) = values.get("max") {
                        self.max = *max;
                    }
                    if let Some(data_size) = values.get("data_size") {
                        self.data_size = *data_size;
                    }
                }
                Some(SectionType::Index) => {
                    self.reader.skip_section()?;
                }
                Some(section_type @ (SectionType::Raw | SectionType::Minimizer)) => {
                    let (format, minimizer, nb_blocks) = self
                        .reader
                        .read_block_section_header(section_type == SectionType::Minimizer)?;
                    self.k = format.k;
                    self.max = format.max;
                    self.data_size = format.data_size;
                    // Sections without block are skipped
                    if nb_blocks > 0 {
                        self.section = Some(CurrentSection {
                            format,
                            minimizer,
                            remaining_blocks: nb_blocks,
                        });
                        return Ok(true);
                    }
                }
                Some(SectionType::Other(byte)) => return Err(error::Error::UnknownSection(byte)),
            }
        }
    }

    /// True if there is at least one more kmer to read
    pub fn has_next(&mut self) -> error::Result<bool> {
        if self.block_has_kmers() {
            return Ok(true);
        }
        if let Some(section) = &self.section {
            if section.remaining_blocks > 0 {
                return Ok(true);
            }
        }
        self.section = None;
        self.next_section()
    }

    /// Read the next block of the current section (has_next must be true)
    fn read_next_block(&mut self) -> error::Result<()> {
        let section = self
            .section
            .as_mut()
            .expect("no section with remaining blocks");
        self.reader
            .read_block_into(&section.format, &section.minimizer, &mut self.block)?;
        section.remaining_blocks -= 1;
        self.kmers_read = 0;
        self.shifts_ready = 1;
        Ok(())
    }

    /// Next block of the file, None at the end of the file. The kmers of the current block not
    /// given by [KmerReader::next_kmer] are skipped.
    pub fn next_block(&mut self) -> error::Result<Option<&Block>> {
        // The rest of a block partially read is skipped
        self.kmers_read = self.block.nb_kmers;
        if !self.has_next()? {
            return Ok(None);
        }
        self.read_next_block()?;
        self.kmers_read = self.block.nb_kmers;
        Ok(Some(&self.block))
    }

    /// Next kmer of the file with its data, None at the end of the file.
    ///
    /// The kmer is given as right aligned 2-bit nucleotides with null padding bits
    /// (`ceil(k / 4)` bytes), see [crate::encoding].
    pub fn next_kmer(&mut self) -> error::Result<Option<KmerView<'_>>> {
        if !self.block_has_kmers() {
            if !self.has_next()? {
                return Ok(None);
            }
            self.read_next_block()?;
        }

        let i = self.kmers_read as usize;
        let k = self.block.k;
        let kmer_bytes = bytes_for_nucleotides(k);
        let kmer_padding = padding(k);
        let seq_padding = padding(self.block.sequence_size());

        // Shift s puts the first nucleotide of the kmer at the position kmer_padding of a byte:
        // the kmer is then made of whole bytes and ends at the end of a byte
        let shift = (kmer_padding + 4 * 4 - seq_padding - i % 4) % 4;
        if self.shifts_ready & (1 << shift) == 0 {
            // Right shift of 2 * shift bits. The padding bits of the sequence stay before its
            // first nucleotide: they are removed by the mask of the first byte of the kmers.
            let sequence = &self.block.sequence;
            let shifted = &mut self.shifts[shift];
            shifted.clear();
            shifted.reserve(sequence.len() + 1);
            let bits = 2 * shift as u32;
            let mut previous = 0u8;
            for byte in sequence {
                shifted.push((previous << (8 - bits)) | (byte >> bits));
                previous = *byte;
            }
            shifted.push(previous << (8 - bits));
            self.shifts_ready |= 1 << shift;
        }
        let source = if shift == 0 {
            &self.block.sequence
        } else {
            &self.shifts[shift]
        };
        let start = (seq_padding + shift + i - kmer_padding) / 4;
        self.kmer.resize(kmer_bytes, 0);
        self.kmer
            .copy_from_slice(&source[start..start + kmer_bytes]);
        // Nucleotides before the kmer in its first byte
        self.kmer[0] &= 0xFF >> (2 * kmer_padding);
        let data_size = self.block.data_size as usize;
        let data = &self.block.data[i * data_size..(i + 1) * data_size];
        self.kmers_read += 1;

        Ok(Some(KmerView {
            k,
            sequence: &self.kmer,
            data,
            encoding: self.reader.encoding(),
        }))
    }

    /// Iterator over the kmers of the file (owned kmers)
    pub fn kmers(self) -> Kmers<R> {
        Kmers { reader: self }
    }
}

/// Iterator over the kmers of a file, see [KmerReader::kmers]
pub struct Kmers<R>
where
    R: Read + Seek,
{
    reader: KmerReader<R>,
}

impl<R> Iterator for Kmers<R>
where
    R: Read + Seek,
{
    type Item = error::Result<Kmer>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader
            .next_kmer()
            .transpose()
            .map(|kmer| kmer.map(|kmer| kmer.to_kmer()))
    }
}
