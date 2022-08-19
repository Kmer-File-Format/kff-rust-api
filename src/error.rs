//! Error struct of project kff

/* crate use */
use thiserror;

/// Enum to manage error
#[derive(std::fmt::Debug, thiserror::Error)]
pub enum Error {
    /// Kff error
    #[error(transparent)]
    Kff(#[from] Kff),

    /// Standard io error
    #[error(transparent)]
    Io(#[from] std::io::Error),

    /// Standard from Utf8 error
    #[error(transparent)]
    FromUtf8(#[from] std::string::FromUtf8Error),
}

/// Kff specific error
#[derive(std::fmt::Debug, thiserror::Error)]
pub enum Kff {
    /// Missing magic number at begin or end
    #[error("Missing magic number {0}")]
    MissingMagic(String),

    /// Major version number is upper than support
    #[error("Major version number {0} is upper than support (1 or lower)")]
    HighMajorVersionNumber(u8),

    /// Minor version number is upper than support
    #[error("Minor version number {0} is upper than support (0 or lower)")]
    HighMinorVersionNumber(u8),

    /// Encoding isn't valid each pair of bits must be different
    #[error("Encoding {0:#b} isn't a valid encoding, each pair of bits must be different")]
    BadEncoding(u8),

    /// Value with name 'name' isn't present in Values this kff file seems not correct
    #[error("Value with name '{0}' isn't present in Values this kff file seems not correct")]
    FieldIsMissing(String),

    /// Value max seems to be too large
    #[error("Value max `{0}` seems to be too large")]
    MaxValueIsTooLarge(u64),

    /// Footer size isn't correct or file not respect footer good practices
    #[error("Footer size isn't correct of file not respect footer good practices")]
    FooterSizeNotCorrect,

    /// Not a valid Kff section prefix
    #[error("'{0}' isn't a valid Kff section prefix")]
    NotASectionPrefix(u8),

    /// Not an index
    #[error("Section at position isn't an index section")]
    NotAnIndex,

    /// No 'first_index' in footer
    #[error("Variable 'first_index' not present in footer it's seems not be an indexed Kff file")]
    NoFirstIndex,

    /// No global index
    #[error("To read a kmers of a section, KFF file should be fully indexed")]
    NoIndex,

    /// No Value section before target section
    #[error("No value section before target section")]
    NoValueSectionBeforeTarget,

    /// Not a kmer section
    #[error("Not a kmer section, you try to read a section isn't ")]
    NotAKmerSection,
}

/// Alias of result
pub type Result<T> = core::result::Result<T, Error>;
