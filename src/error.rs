use thiserror::Error;

#[derive(Error, Debug)]
pub enum Header {
    #[error("Encoding isn't a valid")]
    BadEncoding,
}

#[derive(Error, Debug)]
pub enum Variables {
    #[error("k global variable isn't set")]
    KMissing,

    #[error("m global variable isn't set")]
    MMissing,

    #[error("max global variable isn't set")]
    MaxMissing,

    #[error("data_size global variable isn't set")]
    DataSizeMissing,
}

#[derive(Error, Debug)]
pub enum Data {
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

#[derive(Error, Debug)]
pub enum Error {
    #[error(transparent)]
    Header(Header),

    #[error(transparent)]
    Variables(Variables),

    #[error(transparent)]
    Data(Data),

    #[error(transparent)]
    Minimizer(Minimizer),

    #[error(transparent)]
    Kff(Kff),

    #[error(transparent)]
    OtherError(Box<dyn std::error::Error>),
}

pub trait LocalResult<T> {
    fn map_local(self) -> Result<T, Error>;
}

impl<T, E> LocalResult<T> for std::result::Result<T, E>
where
    E: std::error::Error + 'static,
{
    fn map_local(self) -> Result<T, Error> {
        match self {
            Ok(o) => Ok(o),
            Err(e) => Err(crate::error::Error::OtherError(Box::new(e))),
        }
    }
}

pub type Result<T, E = Error> = std::result::Result<T, E>;
