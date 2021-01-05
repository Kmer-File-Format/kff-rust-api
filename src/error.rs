//! Declaration of error type

use thiserror::Error;

/// Error can be generate in Header section
#[derive(Error, Debug)]
pub enum Header {
    /// Value associate to each nucleotide must be differente
    #[error("Encoding isn't a valid")]
    BadEncoding,
}

/// Error can be generate in Variable section
#[derive(Error, Debug)]
pub enum Variables {
    /// To read Raw and Minimizer section, variable k must be set
    #[error("k global variable isn't set")]
    KMissing,

    /// To read Minimizer section, variable m must be set
    #[error("m global variable isn't set")]
    MMissing,

    /// To read Raw and Minimizer section, variable max must be set
    #[error("max global variable isn't set")]
    MaxMissing,

    /// To read Raw and Minimizer section variable data_size must be set
    #[error("data_size global variable isn't set")]
    DataSizeMissing,
}

/// Error can be generate in Raw or Minimizer section 
#[derive(Error, Debug)]
pub enum Data {
    /// Number of kmer in block can't be larger than max
    #[error("kmer in block is upper than max kmer in block")]
    NUpperThanMax,

    /// Number of kmer and number of data must be equal
    #[error("number of kmer and number of data is differente")]
    NbKmerNbDataDiff,

    /// We can't have more than u32::max() block in section
    #[error("to many block")]
    ToManyBlock,
}

/// Error can be generate in Minimizer section
#[derive(Error, Debug)]
pub enum Minimizer {
    /// Size of minimizer sequence is differente than m variables
    #[error("minimizer size and m global variable is differente")]
    MinimizerSizeMDiff,
}

/// Error can be generate at any moment durring kff reading or writing
#[derive(Error, Debug)]
pub enum Kff {
    /// This api support only Variables, Raw and Minimizer section
    #[error("section type isn't valid")]
    UnknowSectionType,

    /// Not suport version of KFF file
    #[error("wrong version number, support version is major 1 and minor upper than 0")]
    NotSupportVersionNumber,
}

/// Generale error type use in this parser
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

    /// Field used to store Error design outside of local crate
    #[error(transparent)]
    OtherError(Box<dyn std::error::Error>),
}

/// Trait to convert classic Result in Kff::Result
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

/// Just some syntaxic sugar
pub type Result<T, E = Error> = std::result::Result<T, E>;
