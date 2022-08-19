//! Kmer File Format Rust parser

/* std use */

/* crate use */
use std::io::Seek;

/* project use */
use crate::error;
use crate::section;
use crate::utils;
use crate::GlobalIndex;
use crate::KffRead;
use crate::Kmer;
use crate::KmerIterator;

use crate::section::values::AbcValues as _;

/// Struct to read a kff file
#[derive(getset::Getters, getset::Setters, getset::MutGetters)]
#[getset(get = "pub")]
pub struct Kff<R>
where
    R: std::io::Read,
{
    /// Inner read source
    inner: R,

    /// Header extract from `inner`
    #[getset(set = "pub", get_mut = "pub")]
    header: section::Header,

    /// Current Values extract from `inner`
    #[getset(set = "pub", get_mut = "pub")]
    values: section::Values,

    /// GlobalIndex present only if inner is seekable and first section is index or footer contains first_index
    index: Option<utils::GlobalIndex>,
}

impl<R> Kff<R>
where
    R: std::io::Read + std::io::BufRead + crate::KffRead,
{
    /// Create a new Kff reader by accept mutable reference on [std::io::Read]
    pub fn new(mut inner: R) -> error::Result<Self> {
        let header = section::Header::read(&mut inner)?;
        let values = section::Values::default();

        Ok(Self {
            inner,
            header,
            values,
            index: None,
        })
    }

    /// Consume Kff object to create a KmerIterator
    pub fn kmers(self) -> KmerIterator<R> {
        KmerIterator::new(self)
    }

    /// Read Kff until last kmer section
    pub fn next_kmer_section(&mut self) -> std::option::Option<error::Result<Vec<Kmer>>> {
        loop {
            match self.inner.read_u8() {
                Ok(b'v') => {
                    self.values = {
                        match section::Values::read(&mut self.inner) {
                            Err(e) => return Some(Err(e)),
                            Ok(v) => v,
                        }
                    }
                }
                Ok(b'r') => match section::Raw::new(&self.values) {
                    Ok(section) => return Some(section.read(&mut self.inner)),
                    Err(e) => return Some(Err(e)),
                },
                Ok(b'm') => match section::Minimizer::new(&self.values) {
                    Ok(section) => return Some(section.read(&mut self.inner)),
                    Err(e) => return Some(Err(e)),
                },
                Ok(b'K') => return None, // It's the begin of last signature stop reading
                Ok(b'i') => match section::Index::skip(&mut self.inner) {
                    Err(e) => return Some(Err(e)),
                    Ok(_) => continue,
                },
                Ok(e) => return Some(Err(error::Kff::NotASectionPrefix(e).into())), // Any other value is an error
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

impl Kff<std::io::BufReader<std::fs::File>> {
    /// Create a new Kff by read file match with path
    pub fn open<P>(path: P) -> error::Result<Self>
    where
        P: std::convert::AsRef<std::path::Path>,
    {
        std::fs::File::open(&path)
            .map(std::io::BufReader::new)
            .map(Kff::new)?
    }

    /// Create a Kff and generate a global index
    pub fn with_index<P>(path: P) -> error::Result<Self>
    where
        P: std::convert::AsRef<std::path::Path>,
    {
        let mut inner = std::fs::File::open(&path).map(std::io::BufReader::new)?;

        let header = section::Header::read(&mut inner)?;
        let values = section::Values::default();

        let pos_first_section = inner.seek(std::io::SeekFrom::Current(0))?;
        let index = match utils::GlobalIndex::new(&mut inner, pos_first_section) {
            Ok(index) => Some(index),
            Err(error::Error::Kff(error::Kff::NotAnIndex)) => {
                let value = Kff::load_footer(&mut inner)?;

                Some(GlobalIndex::new(
                    &mut inner,
                    *value.get("first_index").ok_or(error::Kff::NoFirstIndex)?,
                )?)
            }
            Err(e) => return Err(e),
        };

        Ok(Self {
            inner,
            header,
            values,
            index,
        })
    }
}

impl<R> Kff<R>
where
    R: std::io::Read + std::io::Seek + KffRead,
{
    /// Check readable match with a KFF file
    pub fn check(&mut self) -> error::Result<bool> {
        self.inner.seek(std::io::SeekFrom::Start(0))?;
        let magic_number = self.inner.read_n_bytes::<3>()?;
        if &magic_number != b"KFF" {
            return Err(error::Kff::MissingMagic("start".to_string()).into());
        }

        self.inner.seek(std::io::SeekFrom::End(-3))?;
        let magic_number = self.inner.read_n_bytes::<3>()?;
        if &magic_number != b"KFF" {
            return Err(error::Kff::MissingMagic("end".to_string()).into());
        }

        self.inner.seek(std::io::SeekFrom::Start(0))?;

        Ok(true)
    }

    /// Load footer, assume last section is a value and last value of this section is footer_size
    fn load_footer(inner: &mut R) -> error::Result<section::Values> {
        inner.seek(std::io::SeekFrom::End(-11))?;
        let footer_size = inner.read_u64()?;

        inner.seek(std::io::SeekFrom::End(-(footer_size as i64 + 3)))?;

        let v = inner.read_u8()?;
        if v != b'v' {
            Err(error::Kff::FooterSizeNotCorrect.into())
        } else {
            section::Values::read(inner)
        }
    }

    /// Get kmer of nth section in index.
    ///
    /// If index isn't set return an Error
    /// If we didn't found section value before target section return an Error
    /// If section isn't a kmer section return an Error
    pub fn kmer_of_section(&mut self, n: usize) -> error::Result<Vec<Kmer>> {
        let index = self
            .index
            .as_ref()
            .ok_or(error::Error::Kff(error::Kff::NoIndex))?;

        self.values = match index.pair()[..n].iter().rev().find(|x| x.0 == b'r') {
            Some((_t, p)) => {
                self.inner.seek(std::io::SeekFrom::Start(p + 1))?;
                section::Values::read(&mut self.inner)?
            }
            None => return Err(error::Kff::NoValueSectionBeforeTarget.into()),
        };

        self.inner
            .seek(std::io::SeekFrom::Start(index.pair()[n].1))?;
        match self.inner.read_u8()? {
            b'r' => match section::Raw::new(&self.values) {
                Ok(section) => section.read(&mut self.inner),
                Err(e) => Err(e),
            },
            b'm' => match section::Minimizer::new(&self.values) {
                Ok(section) => section.read(&mut self.inner),
                Err(e) => Err(e),
            },
            _ => Err(error::Kff::NotAKmerSection.into()),
        }
    }
}

impl<T> std::io::Seek for Kff<T>
where
    T: std::io::Read + std::io::Seek,
{
    fn seek(&mut self, pos: std::io::SeekFrom) -> std::io::Result<u64> {
        self.inner.seek(pos)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Seek;

    const KFF_FILE: &[u8] = &[
        b'K', b'F', b'F',       // Magic number
        1,          // Major version number
        0,          // Minor version number
        0b00101110, // Encoding
        1,          // Uniq kmer
        0,          // Canonical kmer
        0, 0, 0, 4, // Free space size
        b't', b'e', b's', b't', // Free space
        b'v', // Footer
        0, 0, 0, 0, 0, 0, 0, 1, // Footer nb variables
        b'f', b'o', b'o', b't', b'e', b'r', b'_', b's', b'i', b'z', b'e',
        0, // name of variable footer_size
        0, 0, 0, 0, 0, 0, 0, 29, // value of variable footer_size
        b'K', b'F', b'F', // Magic number
    ];

    use std::io::Write;

    #[test]
    fn create_kff_reader() -> error::Result<()> {
        let inner = std::io::Cursor::new(KFF_FILE.to_vec());

        let mut reader = Kff::new(inner)?;
        assert!(reader.check()?);

        let mut tmpfile = tempfile::NamedTempFile::new()?;
        tmpfile.write_all(KFF_FILE)?;
        let mut reader: Kff<std::io::BufReader<std::fs::File>> =
            Kff::<std::io::BufReader<std::fs::File>>::open(tmpfile.path())?;

        assert!(reader.check()?);

        Ok(())
    }

    #[test]
    fn check_header() -> error::Result<()> {
        let inner = std::io::Cursor::new(KFF_FILE.to_vec());

        let reader = Kff::new(inner)?;

        assert_eq!(reader.header().major_version(), &1);
        assert_eq!(reader.header().minor_version(), &0);
        assert_eq!(reader.header().encoding(), &0b00101110);
        assert_eq!(reader.header().uniq_kmer(), &true);
        assert_eq!(reader.header().canonical_kmer(), &false);
        assert_eq!(reader.header().free_block(), b"test");

        Ok(())
    }

    #[test]
    fn check() -> error::Result<()> {
        let inner_len = KFF_FILE.len();
        let mut readable = std::io::Cursor::new(KFF_FILE.to_vec());

        let mut file = Kff::new(readable.clone())?;

        assert!(file.check()?); // Header init and check work

        readable.get_mut()[1] = b'K';
        let file = Kff::new(readable.clone());
        assert!(file.is_err()); // Header init failled

        readable.get_mut()[1] = b'F';
        readable.get_mut()[inner_len - 1] = b'K';
        let mut file = Kff::new(readable.clone())?;

        assert!(file.check().is_err()); // Header init work but check failled

        Ok(())
    }

    #[test]
    fn load_footer() -> error::Result<()> {
        let inner_len = KFF_FILE.len();
        let mut inner = std::io::Cursor::new(KFF_FILE.to_vec());

        let mut truth = section::Values::new();
        truth.insert("footer_size".to_string(), 29);

        assert_eq!(Kff::load_footer(&mut inner)?, truth);

        inner.get_mut()[inner_len - 32] = b'f';

        assert!(Kff::load_footer(&mut inner).is_err());

        Ok(())
    }

    #[test]
    fn seek() -> error::Result<()> {
        let inner = std::io::Cursor::new(KFF_FILE.to_vec());
        let mut reader = Kff::new(inner)?;

        assert_eq!(reader.seek(std::io::SeekFrom::Start(2))?, 2);

        Ok(())
    }
}
