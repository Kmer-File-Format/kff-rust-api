//! Representation of a KFF kmer block

/* std use */

/* crate use */

/* project use */
use crate::error;

/// Represent a kmer
pub type Kmer = bitvec::boxed::BitBox<u8, bitvec::order::Msb0>;

/// Represent data associate to a kmer
pub type Data = Vec<u8>;

/// Struct to represente a KFF block
#[derive(getset::Getters, std::fmt::Debug)]
#[getset(get = "pub")]
pub struct Block {
    /// size of kmer
    k: usize,

    /// size of data associate (in bytes) to each kmer
    data_size: u64,

    /// number of kmer in block
    nb_kmer: u64,

    /// Bit field store all kmer of this block
    kmer: Kmer,

    /// Array store data associate to kmer of this block
    data: Data,

    /// Actual position of next kmer
    #[getset(skip)]
    offset: usize,
}

impl Block {
    /// Create a new block
    pub fn new(k: usize, data_size: u64, nb_kmer: u64, kmer: Kmer, data: Data) -> Self {
        Self {
            k,
            data_size,
            nb_kmer,
            kmer,
            data,
            offset: 0,
        }
    }

    /// Create a new block by read a raw block
    pub fn from_raw<R>(inner: &mut R, k: usize, max: u64, data_size: u64) -> error::Result<Self>
    where
        R: std::io::Read + crate::KffRead,
    {
        let nb_kmer = if max <= 1 {
            1
        } else {
            read_nb_kmer(inner, max)?
        };

        let kmer = Kmer::from_boxed_slice(
            inner
                .read_2bits((nb_kmer as usize + k - 1) as usize)?
                .into_boxed_slice(),
        );

        let data = inner.read_n_bytes_dyn((nb_kmer * data_size) as usize)?;

        Ok(Self {
            k,
            data_size,
            nb_kmer,
            kmer,
            data,
            offset: 0,
        })
    }

    /// Write content of the block
    pub fn to_raw<W>(&self, outer: &mut W, max: u64) -> error::Result<()>
    where
        W: std::io::Write + crate::KffWrite,
    {
        if max > 1 {
            write_nb_kmer(outer, &self.nb_kmer)?;
        }
        outer.write_bytes(self.kmer.as_raw_slice())?;
        outer.write_bytes(self.data.as_slice())?;

        Ok(())
    }
}

impl std::iter::Iterator for Block {
    type Item = (Kmer, Data);

    fn next(&mut self) -> Option<Self::Item> {
        if self.offset == self.nb_kmer as usize {
            None
        } else {
            let k_range = self.offset * 2..(self.offset + self.k) * 2;
            let d_range =
                self.offset * self.data_size as usize..(self.offset + 1) * self.data_size as usize;

            self.offset += 1;
            Some((
                bitvec::boxed::BitBox::from_bitslice(&self.kmer[k_range]),
                self.data[d_range].to_vec(),
            ))
        }
    }
}

pub(crate) fn read_nb_kmer<R>(inner: &mut R, max: u64) -> error::Result<u64>
where
    R: std::io::Read + crate::KffRead,
{
    match max.leading_zeros() {
        0..=31 => Ok(inner.read_u64()? as u64),
        32..=47 => Ok(inner.read_u32()? as u64),
        48..=55 => Ok(inner.read_u16()? as u64),
        56..=64 => Ok(inner.read_u8()? as u64),
        _ => unreachable!("You can't have more than 64 leading_zeros() with an u64"),
    }
}

pub(crate) fn write_nb_kmer<W>(outer: &mut W, max: &u64) -> error::Result<()>
where
    W: std::io::Write + crate::KffWrite,
{
    println!("{} {:064b} {}", max, max, max.leading_zeros());
    match max.leading_zeros() {
        0..=31 => Ok(outer.write_u64(max)?),
        32..=47 => Ok(outer.write_u32(&(*max as u32))?),
        48..=55 => Ok(outer.write_u16(&(*max as u16))?),
        56..=64 => Ok(outer.write_u8(&(*max as u8))?),
        _ => unreachable!("You can't have more than 64 leading_zeros() with an u64"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod raw {
        use super::*;

        #[test]
        fn full() -> error::Result<()> {
            let mut readable: &[u8] = &[3, 0b00011011, 0b11110100, 1, 2, 3];

            let block = Block::from_raw(&mut readable, 5, 255, 1)?;

            let mut kmers: Vec<Kmer> = Vec::new();
            let mut datas: Vec<Data> = Vec::new();
            for (kmer, data) in block {
                kmers.push(kmer);
                datas.push(data);
            }

            assert_eq!(
                &kmers[..],
                &[
                    bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1],
                    bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1, 1, 1, 1, 1],
                    bitvec::bitbox![u8, bitvec::order::Msb0; 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                ]
            );
            assert_eq!(&datas[..], &[vec![1], vec![2], vec![3],]);

            Ok(())
        }

        #[test]
        fn no_data() -> error::Result<()> {
            let mut readable: &[u8] = &[3, 0b00011011, 0b11110100];

            let block = Block::from_raw(&mut readable, 5, 255, 0)?;

            let mut kmers: Vec<Kmer> = Vec::new();
            let mut datas: Vec<Data> = Vec::new();
            for (kmer, data) in block {
                kmers.push(kmer);
                datas.push(data);
            }

            assert_eq!(
                &kmers[..],
                &[
                    bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1],
                    bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1, 1, 1, 1, 1],
                    bitvec::bitbox![u8, bitvec::order::Msb0; 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                ]
            );
            assert_eq!(&datas[..], &[vec![], vec![], vec![],]);

            Ok(())
        }

        #[test]
        fn max_one_kmer() -> error::Result<()> {
            let mut readable: &[u8] = &[0b00011011, 0b11000000, 1];

            let block = Block::from_raw(&mut readable, 5, 1, 1)?;

            let mut kmers: Vec<Kmer> = Vec::new();
            let mut datas: Vec<Data> = Vec::new();
            for (kmer, data) in block {
                kmers.push(kmer);
                datas.push(data);
            }

            assert_eq!(
                &kmers[..],
                &[bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1],]
            );
            assert_eq!(&datas[..], &[vec![1]]);

            Ok(())
        }

        #[test]
        fn write() -> error::Result<()> {
            let block = Block::new(
                5,
                1,
                3,
                bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0],
                vec![1, 2, 3],
            );

            let mut writable = Vec::new();

            block.to_raw(&mut writable, 255)?;

            assert_eq!(writable, vec![3, 0b00011011, 0b11110100, 1, 2, 3]);

            Ok(())
        }
    }

    #[test]
    fn max_value_read() -> error::Result<()> {
        let readable: &[u8] = &[1, 2, 3, 4, 5, 6, 7, 8];
        let value = read_nb_kmer(&mut readable.clone(), u8::MAX.into())?;
        assert_eq!(value, 1);

        let value = read_nb_kmer(&mut readable.clone(), u8::MAX as u64 + 1)?;
        assert_eq!(value, 258);

        let value = read_nb_kmer(&mut readable.clone(), u16::MAX.into())?;
        assert_eq!(value, 258);

        let value = read_nb_kmer(&mut readable.clone(), u16::MAX as u64 + 1)?;
        assert_eq!(value, 16909060);

        let value = read_nb_kmer(&mut readable.clone(), u32::MAX.into())?;
        assert_eq!(value, 16909060);

        let value = read_nb_kmer(&mut readable.clone(), u32::MAX as u64 + 1)?;
        assert_eq!(value, 72623859790382856);

        let value = read_nb_kmer(&mut readable.clone(), u64::MAX)?;
        assert_eq!(value, 72623859790382856);

        Ok(())
    }

    #[test]
    fn max_value_write() -> error::Result<()> {
        let mut writable = Vec::new();

        write_nb_kmer(&mut writable, &(u8::MAX as u64))?;
        assert_eq!(writable, vec![255]);
        writable.clear();

        write_nb_kmer(&mut writable, &((u8::MAX as u64) + 1))?;
        assert_eq!(writable, vec![1, 0]);
        writable.clear();

        write_nb_kmer(&mut writable, &(u16::MAX as u64))?;
        assert_eq!(writable, vec![255, 255]);
        writable.clear();

        write_nb_kmer(&mut writable, &((u16::MAX as u64) + 1))?;
        assert_eq!(writable, vec![0, 1, 0, 0]);
        writable.clear();

        write_nb_kmer(&mut writable, &(u32::MAX as u64))?;
        assert_eq!(writable, vec![255, 255, 255, 255]);
        writable.clear();

        write_nb_kmer(&mut writable, &((u32::MAX as u64) + 1))?;
        assert_eq!(writable, vec![0, 0, 0, 1, 0, 0, 0, 0]);
        writable.clear();

        write_nb_kmer(&mut writable, &(u64::MAX as u64))?;
        assert_eq!(writable, vec![255, 255, 255, 255, 255, 255, 255, 255]);
        writable.clear();

        Ok(())
    }
}
