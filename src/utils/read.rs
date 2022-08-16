//! Utils function to read KFF

/* std use */

/* crate use */

/* project use */
use crate::error;

/// Define trait containts utils function to parsing kff
pub trait KffRead {
    /// Function read N bytes (N define at compile time) in a readable
    fn read_n_bytes<const N: usize>(&mut self) -> error::Result<[u8; N]>;

    /// Function read N bytes (N define at run time) in a readable
    fn read_n_bytes_dyn(&mut self, n: usize) -> error::Result<Vec<u8>>;

    /// Function read a Kff 'ascii'
    fn read_ascii(&mut self) -> error::Result<Vec<u8>>;

    /// Function some base in 2bits representation
    fn read_2bits(
        &mut self,
        k: usize,
    ) -> error::Result<bitvec::vec::BitVec<u8, bitvec::order::Msb0>>;

    /// Function that read one bit and convert it as bool
    fn read_bool(&mut self) -> error::Result<bool> {
        self.read_u8().map(|x| x != 0)
    }

    /// Function that read u8
    fn read_u8(&mut self) -> error::Result<u8> {
        self.read_n_bytes::<1>()
            .map(|x| unsafe { *x.get_unchecked(0) })
    }

    /// Function that read u16
    fn read_u16(&mut self) -> error::Result<u16> {
        self.read_n_bytes::<2>().map(u16::from_be_bytes)
    }

    /// Function that read u32
    fn read_u32(&mut self) -> error::Result<u32> {
        self.read_n_bytes::<4>().map(u32::from_be_bytes)
    }

    /// Function that read u64
    fn read_u64(&mut self) -> error::Result<u64> {
        self.read_n_bytes::<8>().map(u64::from_be_bytes)
    }

    /// Function that read i64
    fn read_i64(&mut self) -> error::Result<i64> {
        self.read_n_bytes::<8>().map(i64::from_be_bytes)
    }
}

impl<T> KffRead for T
where
    T: std::io::BufRead,
{
    fn read_n_bytes<const N: usize>(&mut self) -> error::Result<[u8; N]> {
        let mut values = [0; N];

        self.read_exact(&mut values)?;

        Ok(values)
    }

    fn read_n_bytes_dyn(&mut self, n: usize) -> error::Result<Vec<u8>> {
        let mut values = vec![0; n];

        self.read_exact(&mut values)?;

        Ok(values)
    }

    fn read_ascii(&mut self) -> error::Result<Vec<u8>> {
        let mut values = Vec::with_capacity(50);

        self.read_until(0, &mut values)?;

        if let Some(0) = values.last() {
            values.pop();
        }

        Ok(values)
    }

    fn read_2bits(
        &mut self,
        k: usize,
    ) -> error::Result<bitvec::vec::BitVec<u8, bitvec::order::Msb0>> {
        let mut values = bitvec::vec::BitVec::from_slice(
            &self.read_n_bytes_dyn(crate::bytes2store_k(k as u64) as usize)?,
        );

        values.resize(k * 2, false);

        Ok(values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Seek;

    const LOREM: &[u8] = b"Lorem ipsum dolor\0sit amet, consectetur adipiscing elit.";

    #[test]
    fn read_n_bytes() -> error::Result<()> {
        let mut reader = std::io::Cursor::new(LOREM);

        let values = reader.read_n_bytes::<11>()?;

        assert_eq!(&values, b"Lorem ipsum");

        let values = reader.read_n_bytes::<11>()?;

        assert_eq!(&values, b" dolor\0sit ");

        let values = reader.read_n_bytes::<400>();

        assert!(values.is_err());

        Ok(())
    }

    #[test]
    fn read_n_bytes_dyn() -> error::Result<()> {
        let mut reader = std::io::Cursor::new(LOREM);

        let values = reader.read_n_bytes_dyn(11)?;

        assert_eq!(&values, b"Lorem ipsum");

        let values = reader.read_n_bytes_dyn(11)?;

        assert_eq!(&values, b" dolor\0sit ");

        let values = reader.read_n_bytes_dyn(400);

        assert!(values.is_err());

        Ok(())
    }

    #[test]
    fn read_ascii() -> error::Result<()> {
        let mut reader = std::io::Cursor::new(LOREM);

        let values = reader.read_ascii()?;

        assert_eq!(&values, b"Lorem ipsum dolor");

        reader.seek(std::io::SeekFrom::Start(values.len() as u64 + 1))?; // Move after first \0
        let values = reader.read_ascii()?;

        assert_eq!(&values, b"sit amet, consectetur adipiscing elit.");

        reader.seek(std::io::SeekFrom::End(0))?;
        let values = reader.read_ascii()?;

        assert_eq!(&values, b"");

        Ok(())
    }

    #[test]
    fn read_2bits() -> error::Result<()> {
        let mut reader = std::io::Cursor::new([0b11101110, 0b00010001]);

        let kmer = reader.read_2bits(5)?;

        assert_eq!(
            kmer,
            bitvec::bitvec![u8, bitvec::order::Msb0; 1, 1, 1, 0, 1, 1, 1, 0, 0, 0]
        );

        Ok(())
    }

    #[test]
    fn read_bool() -> error::Result<()> {
        let mut reader = std::io::Cursor::new(LOREM);

        assert!(reader.read_bool()?);

        let _ = reader.read_n_bytes::<16>()?;

        assert!(!reader.read_bool()?);

        Ok(())
    }

    #[test]
    fn read_u8() -> error::Result<()> {
        let mut reader = std::io::Cursor::new(LOREM);

        assert_eq!(reader.read_u8()?, b'L');
        assert_eq!(reader.read_u8()?, b'o');

        Ok(())
    }

    #[test]
    fn read_u16() -> error::Result<()> {
        let mut reader = std::io::Cursor::new(LOREM);

        assert_eq!(reader.read_u16()?, 19567);
        assert_eq!(reader.read_u16()?, 29285);

        Ok(())
    }

    #[test]
    fn read_u32() -> error::Result<()> {
        let mut reader = std::io::Cursor::new(LOREM);

        assert_eq!(reader.read_u32()?, 1282372197);
        assert_eq!(reader.read_u32()?, 1830840688);

        Ok(())
    }

    #[test]
    fn read_u64() -> error::Result<()> {
        let mut reader = std::io::Cursor::new(LOREM);

        assert_eq!(reader.read_u64()?, 5507746649245510000);
        assert_eq!(reader.read_u64()?, 8319675872528264303);

        Ok(())
    }
}
