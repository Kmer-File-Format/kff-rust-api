//! Kmer File Format Rust parser

/* std use */

/* crate use */

/* project use */
use crate::error;
use crate::section;

/// Struct to read a kff file
pub struct Kff<R>
where
    R: std::io::Read,
{
    inner: R,
    header: section::Header,
    // values: section::Values,
    // index: section::Index,
}

impl<R> Kff<R>
where
    R: std::io::Read + std::io::BufRead + crate::KffRead,
{
    /// Create a new Kff reader by accept mutable reference on [std::io::Read]
    pub fn new(mut inner: R) -> error::Result<Self> {
        let header = section::Header::read(&mut inner)?;

        Ok(Self { inner, header })
    }

    /// Create a new Kff by read file match with path in parameter
    pub fn open<P>(path: P) -> error::Result<Kff<std::io::BufReader<std::fs::File>>>
    where
        P: std::convert::AsRef<std::path::Path>,
    {
        std::fs::File::open(&path)
            .map(std::io::BufReader::new)
            .map(Kff::new)?
    }
}

impl<R> Kff<R>
where
    R: std::io::Read + std::io::Seek + crate::KffRead,
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
}

#[cfg(test)]
mod tests {
    use super::*;

    const KFF_FILE: &[u8] = b"KFF test KFF";

    use std::io::Write;

    #[test]
    #[ignore]
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
    #[ignore]
    fn check() -> error::Result<()> {
        let inner_len = KFF_FILE.len();
        let mut readable = std::io::Cursor::new(KFF_FILE.to_vec());

        let mut file = Kff::new(readable.clone())?;

        assert!(file.check()?);

        readable.get_mut()[1] = b'K';
        let mut file = Kff::new(readable.clone())?;

        assert!(file.check().is_err());

        readable.get_mut()[1] = b'F';
        readable.get_mut()[inner_len - 1] = b'K';
        let mut file = Kff::new(readable.clone())?;

        assert!(file.check().is_err());

        Ok(())
    }
}
