//! Header of a KFF file

/* project use */
use crate::Encoding;

/// Header of a KFF file
#[derive(Clone, std::fmt::Debug, PartialEq, Eq)]
pub struct Header {
    /// Major version of the format of the file
    pub major_version: u8,

    /// Minor version of the format of the file
    pub minor_version: u8,

    /// Encoding of the nucleotides
    pub encoding: Encoding,

    /// True if each kmer is present at most once in the file
    pub uniqueness: bool,

    /// True if a kmer and its reverse complement are never both in the file
    pub canonicity: bool,

    /// Content of the free block (metadata) of the header
    pub metadata: Vec<u8>,
}

impl Header {
    /// Size of the header in bytes (signature included)
    pub fn size(&self) -> u64 {
        12 + self.metadata.len() as u64
    }
}

impl std::default::Default for Header {
    /// Header of a new file: current version of the format, default encoding, no flag and no
    /// metadata
    fn default() -> Self {
        Self {
            major_version: crate::KFF_FORMAT_VERSION_MAJOR,
            minor_version: crate::KFF_FORMAT_VERSION_MINOR,
            encoding: Encoding::default(),
            uniqueness: false,
            canonicity: false,
            metadata: Vec::new(),
        }
    }
}
