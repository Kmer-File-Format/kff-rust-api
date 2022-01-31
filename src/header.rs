//! Module to manage kff header

/* std use */

/* crate use */
use byteorder::ReadBytesExt;

/* project use */
use crate::*;

/// Struct store header information
struct Header {
    marker: [u8; 3],
    major: u8,
    minor: u8,
    encoding: encoding::Encoding,
    uniq_kmer: bool,
    canonical_kmer: bool,
    free_block: Vec<u8>,
}

impl Header {
    /// Create a new header
    pub fn new(
        major: u8,
        minor: u8,
        encoding: encoding::Encoding,
        uniq_kmer: bool,
        canonical_kmer: bool,
        free_block: Vec<u8>,
    ) -> Header {
        Self {
            marker: [b'K', b'F', b'F'],
            major,
            minor,
            encoding,
            uniq_kmer,
            canonical_kmer,
            free_block,
        }
    }

    /// Create a new header from reading stream
    pub fn from_reader<R>(input: &mut R) -> error::Result<Header>
    where
        R: std::io::Read,
    {
        let marker = io::read_fix::<R, 3>(input)?;

        if &marker != b"KFF" {
            return Err(error::Error::HeaderMissingMarker);
        }

        let major = input.read_u8()?;
        let minor = input.read_u8()?;

        if major > 1 && minor > 0 {
            return Err(error::Error::VersionNotSupport);
        }

        let encoding = encoding::Encoding::new(input.read_u8()?)?;

        let uniq_kmer = input.read_u8()? > 0;

        let canonical_kmer = input.read_u8()? > 0;

        let free_block_size = input.read_u32::<Endianess>()?;
        let free_block = io::read_dynamic(input, free_block_size as usize)?;

        Ok(Self {
            marker,
            major,
            minor,
            encoding,
            uniq_kmer,
            canonical_kmer,
            free_block,
        })
    }
}
