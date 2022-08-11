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
}

/// Alias of result
pub type Result<T> = core::result::Result<T, Error>;
