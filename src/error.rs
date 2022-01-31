//! Define error produce by kff crates

/* std use */

/* crate use */

/* project use */

/* mod declaration */

/// Enum of all error can be produce by KFF
#[derive(thiserror::Error, std::fmt::Debug)]
pub enum Error {
    /// std::io error
    #[error(transparent)]
    IO(#[from] std::io::Error),

    /// KFF header section start by string 'KFF'
    #[error("KFF format start by string 'KFF'")]
    HeaderMissingMarker,

    /// KFF version isn't support
    #[error("KFF version support is 1.0 or lower")]
    VersionNotSupport,

    /// Invalid encoding
    #[error("Encoding isn't valid")]
    EncodingNotValid,
}

/// A type alias to simplify usage of Result
pub type Result<T> = std::result::Result<T, Error>;
