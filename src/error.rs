//! Errors of the crate

/// Errors raised when reading or writing a KFF file
#[derive(std::fmt::Debug, thiserror::Error)]
pub enum Error {
    /// Error of the underlying reader or writer
    #[error(transparent)]
    Io(#[from] std::io::Error),

    /// The "KFF" signature is absent at the beginning or at the end of the file
    #[error(
        "absent KFF signature at the {0} of the file, please check that the file is not corrupted"
    )]
    MissingSignature(&'static str),

    /// The file is written in a version of the format that this library cannot read
    #[error("this library reads the KFF format up to version {}.{}: it can't read files written in version {major}.{minor}", crate::KFF_FORMAT_VERSION_MAJOR, crate::KFF_FORMAT_VERSION_MINOR)]
    UnsupportedVersion {
        /// Major version of the file
        major: u8,
        /// Minor version of the file
        minor: u8,
    },

    /// The 4 nucleotides of an encoding must have different 2-bit values
    #[error("wrong encoding {0:#010b}: the 4 2-bit values must be different")]
    InvalidEncoding(u8),

    /// A character of a sequence is not a nucleotide (A, C, G or T)
    #[error("{0:?} is not a nucleotide")]
    InvalidNucleotide(char),

    /// A variable needed by a section is not defined in the current values section
    #[error("variable '{0}' is not defined in the current values section")]
    MissingVariable(&'static str),

    /// A variable has a value that can't be used
    #[error("invalid value for variable '{name}': {reason}")]
    InvalidVariable {
        /// Name of the variable
        name: &'static str,
        /// Why the value is rejected
        reason: String,
    },

    /// The name of a variable can't be written in a values section
    #[error("invalid variable name {0:?}: names can't contain a null byte")]
    InvalidVariableName(String),

    /// A section of a given type is expected but the file contains another byte
    #[error("expected a '{expected}' section, found byte {found:#04x}")]
    UnexpectedSection {
        /// Type of the expected section
        expected: char,
        /// Byte found in the file
        found: u8,
    },

    /// The byte at the beginning of a section is not a known section type
    #[error("unknown section type {0:#04x}")]
    UnknownSection(u8),

    /// Two entries of an index section have the same position
    #[error("index collision: two entries of an index section have the same position")]
    IndexCollision,

    /// The footer or the index chain of the file is not consistent
    #[error("invalid footer or index: {0}")]
    InvalidFooter(String),

    /// A block read or written is not consistent with the variables of its section
    #[error("invalid block: {0}")]
    InvalidBlock(String),

    /// A buffer given to the API is too small for the data to read or write
    #[error("buffer too small: {needed} bytes needed, {available} available")]
    BufferTooSmall {
        /// Number of bytes needed
        needed: usize,
        /// Size of the buffer
        available: usize,
    },

    /// The metadata of the header must be written before the first section
    #[error("the metadata have to be written prior to other content")]
    MetadataAfterSection,

    /// A jump or a write goes out of the file
    #[error("position {0} is out of the file")]
    OutOfFile(u64),

    /// A writer created on a stream can't release and reopen its file
    #[error("only a writer created from a path can release its file")]
    CannotRelease,
}

/// Result type of the crate
pub type Result<T> = core::result::Result<T, Error>;
