use thiserror::Error;

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(Error, Debug)]
pub enum Header {
    #[error("Encoding isn't a valid")]
    BadEncoding,
}

#[derive(Error, Debug)]
pub enum Data {
    #[error("k global variable isn't set")]
    KMissing,

    #[error("m global variable isn't set")]
    MMissing,

    #[error("max global variable isn't set")]
    MaxMissing,

    #[error("data_size global variable isn't set")]
    DataSizeMissing,

    #[error("kmer in block is upper than max kmer in block")]
    NUpperThanMax,

    #[error("number of kmer and number of data is differente")]
    NbKmerNbDataDiff,

    #[error("to many block")]
    ToManyBlock,
}

#[derive(Error, Debug)]
pub enum Minimizer {
    #[error("minimizer size and m global variable is differente")]
    MinimizerSizeMDiff,
}

#[derive(Error, Debug)]
pub enum Kff {
    #[error("section type isn't valid")]
    UnknowSectionType,

    #[error("wrong version number, support version must be equal or upper 1.0")]
    NotSupportVersionNumber,

    #[error("comment provide is to large")]
    CommentToLarge,

    #[error("a block is already open you can't open ")]
    ABlockIsAlreadyOpen,
}
