//! Kmer File Format Rust parser

/* std use */

/* crate use */

/* project use */
use crate::error;
use crate::section;

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

    const KFF_FILE: &[u8] = &[
        b'K', b'F', b'F',       // Magic number
        1,          // Major version number
        0,          // Minor version number
        0b00101110, // Encoding
        1,          // Uniq kmer
        0,          // Canonical kmer
        0, 0, 0, 4, // Free space size
        b't', b'e', b's', b't', // Free space
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
}
