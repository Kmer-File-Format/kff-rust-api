//! Kmer File Format Rust parser

/* std use */

/* crate use */

/* project use */
use crate::error;

/// Struct to read a kff file
pub struct Kff<'a, R>
where
    R: std::io::Read,
{
    inner: &'a mut R,
}

impl<R> Kff<'_, R>
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

        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const KFF_FILE: &[u8] = b"KFF test KFF";

    #[test]
    fn check() -> error::Result<()> {
        let mut inner = KFF_FILE.to_vec();
        let inner_len = inner.len();
        let mut readable = std::io::Cursor::new(&mut inner);

        let mut file = Kff {
            inner: &mut readable,
        };

        assert!(file.check()?);

        readable.get_mut()[1] = b'K';
        let mut file = Kff {
            inner: &mut readable,
        };

        assert!(file.check().is_err());

        readable.get_mut()[1] = b'F';
        readable.get_mut()[inner_len - 1] = b'K';
        let mut file = Kff {
            inner: &mut readable,
        };

        assert!(file.check().is_err());

        Ok(())
    }
}
