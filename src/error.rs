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
}

/// Kff specific error
#[derive(std::fmt::Debug, thiserror::Error)]
pub enum Kff {
    /// Missing magic number at begin or end
    #[error("Missing magic number {0}")]
    MissingMagic(String),
}

/// Alias of result
pub type Result<T> = core::result::Result<T, Error>;
