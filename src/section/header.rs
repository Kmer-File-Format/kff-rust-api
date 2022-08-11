//! Parse, manage and write Header information

/* std use */

/* crate use */

/* project use */
use crate::error;

/// Struct to parse, manage and write Header information
#[derive(std::fmt::Debug, getset::Getters, getset::Setters, getset::MutGetters)]
#[getset(get = "pub")]
pub struct Header {
    /// Major version number
    major_version: u8,

    /// Minor version number
    minor_version: u8,

    /// Encoding schema
    encoding: u8,

    /// This file contains only uniq kmer
    #[getset(set = "pub", get_mut = "pub")]
    uniq_kmer: bool,

    /// This file contains only canonical kmer
    #[getset(set = "pub", get_mut = "pub")]
    canonical_kmer: bool,

    /// Comment link to this file
    #[getset(set = "pub", get_mut = "pub")]
    free_block: Vec<u8>,
}

impl Header {
    /// Constructor of header
    pub fn new(
        major_version: u8,
        minor_version: u8,
        encoding: u8,
        uniq_kmer: bool,
        canonical_kmer: bool,
        free_block: Vec<u8>,
    ) -> error::Result<Self> {
        let obj = Self {
            major_version,
            minor_version,
            encoding,
            uniq_kmer,
            canonical_kmer,
            free_block,
        };

        obj.check()?;

        Ok(obj)
    }

    /// Read a readable to create a new header
    pub fn read<R>(inner: &mut R) -> error::Result<Self>
    where
        R: std::io::Read + crate::KffRead,
    {
        let mut obj = Self {
            major_version: 0,
            minor_version: 0,
            encoding: 0b00101110,
            uniq_kmer: false,
            canonical_kmer: false,
            free_block: Vec::new(),
        };

        let magic_number = inner.read_n_bytes::<3>()?;
        if &magic_number != b"KFF" {
            return Err(error::Kff::MissingMagic("start".to_string()).into());
        }

        obj.major_version = inner.read_u8()?;
        obj.minor_version = inner.read_u8()?;
        obj.encoding = inner.read_u8()?;
        obj.uniq_kmer = inner.read_bool()?;
        obj.canonical_kmer = inner.read_bool()?;

        let free_block_size = inner.read_u32()? as usize;

        obj.free_block = inner.read_n_bytes_dyn(free_block_size)?;

        obj.check()?;

        Ok(obj)
    }

    /// Write this Header in KFF format
    pub fn write<W>(&self, outer: &mut W) -> error::Result<()>
    where
        W: std::io::Write + crate::KffWrite,
    {
        outer.write_bytes(b"KFF")?; // Write magic number
        outer.write_u8(&self.major_version)?; // Major version
        outer.write_u8(&self.minor_version)?; // Minor version
        outer.write_u8(&self.encoding)?; // Encoding
        outer.write_bool(&self.uniq_kmer)?; // Uniq kmer
        outer.write_bool(&self.canonical_kmer)?; // Canonical kmer
        outer.write_u32(&(self.free_block.len() as u32))?; // Size of free block
        outer.write_ascii(&self.free_block)?; // Free block

        Ok(())
    }

    /// Set major version
    pub fn set_major_version(&mut self, val: u8) -> error::Result<&mut Self> {
        self.major_version = val;

        self.check_version()?;

        Ok(self)
    }

    /// Set minor version
    pub fn set_minor_version(&mut self, val: u8) -> error::Result<&mut Self> {
        self.minor_version = val;

        self.check_version()?;

        Ok(self)
    }

    /// Set encoding
    pub fn set_encoding(&mut self, val: u8) -> error::Result<&mut Self> {
        self.encoding = val;

        self.check_encoding()?;

        Ok(self)
    }

    /// Function run after construction of header to check value
    fn check(&self) -> error::Result<&Self> {
        self.check_version()?.check_encoding()
    }

    /// Function check if version number is support
    fn check_version(&self) -> error::Result<&Self> {
        if self.major_version > 1 {
            return Err(error::Kff::HighMajorVersionNumber(self.major_version).into());
        }

        if self.minor_version > 0 {
            return Err(error::Kff::HighMinorVersionNumber(self.minor_version).into());
        }

        Ok(self)
    }

    /// Function check encoding is a valid one
    fn check_encoding(&self) -> error::Result<&Self> {
        let a = self.encoding >> 6;
        let c = (self.encoding >> 4) & 0b11;
        let t = (self.encoding >> 2) & 0b11;
        let g = self.encoding & 0b11;

        if a != c && a != t && a != g && c != t && t != g {
            Ok(self)
        } else {
            Err(error::Kff::BadEncoding(self.encoding).into())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const VALID: &[u8] = &[
        b'K', b'F', b'F', 1, 0, 0b00101110, 1, 0, 0, 0, 0, 4, b't', b'e', b's', b't', 0,
    ];

    const BAD_MAGIC_NUMBER: &[u8] = b"KKF";

    #[test]
    fn new() -> error::Result<()> {
        assert!(Header::new(1, 0, 0b00101110, true, false, b"test".to_vec()).is_ok());

        assert!(Header::new(2, 0, 0b00101110, true, false, b"test".to_vec()).is_err());
        assert!(Header::new(1, 1, 0b00101110, true, false, b"test".to_vec()).is_err());
        assert!(Header::new(1, 0, 0b11111111, true, false, b"test".to_vec()).is_err());

        Ok(())
    }

    #[test]
    fn read() -> error::Result<()> {
        let mut reader = std::io::Cursor::new(VALID);

        assert!(Header::read(&mut reader).is_ok());

        let mut reader = std::io::Cursor::new(BAD_MAGIC_NUMBER);

        assert!(Header::read(&mut reader).is_err());

        Ok(())
    }

    #[test]
    fn write() -> error::Result<()> {
        let header = Header::new(1, 0, 0b00101110, true, false, b"test".to_vec())?;

        let mut writer = std::io::Cursor::new(Vec::new());

        assert!(header.write(&mut writer).is_ok());

        assert_eq!(VALID, writer.into_inner());

        Ok(())
    }

    #[test]
    fn setter() -> error::Result<()> {
        let mut reader = std::io::Cursor::new(VALID);

        let mut header = Header::read(&mut reader)?;

        assert!(header.set_major_version(0).is_ok());
        assert!(header.set_major_version(1).is_ok());
        assert!(header.set_major_version(2).is_err());
        assert!(header.set_major_version(1).is_ok());

        assert!(header.set_minor_version(0).is_ok());
        assert!(header.set_minor_version(1).is_err());
        assert!(header.set_minor_version(2).is_err());
        assert!(header.set_minor_version(0).is_ok());

        assert!(header.set_encoding(0b00101101).is_ok());
        assert!(header.set_encoding(1).is_err());
        assert!(header.set_encoding(2).is_err());

        Ok(())
    }
}
