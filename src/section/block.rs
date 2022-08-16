//! Representation of a KFF kmer block

/* std use */

/* crate use */

/* project use */
use crate::error;

/// Represent a kmer
pub type Kmer = bitvec::boxed::BitBox<u8, bitvec::order::Msb0>;

/// Represent data associate to a kmer
pub type Data = Vec<u8>;

/// Struct to data present in KFF Raw or Minimizer block
#[derive(getset::Getters, std::fmt::Debug)]
#[getset(get = "pub")]
pub struct Block {
    /// Size of kmer
    pub(crate) k: u64,

    /// Size of data associate (in bytes) to each kmer
    pub(crate) data_size: usize,

    /// Bit field store all kmer of this block
    pub(crate) kmer: Kmer,

    /// Array store data associate to kmer of this block
    pub(crate) data: Data,

    /// Minimizer offset
    pub(crate) minimizer_offset: usize,

    /// Actual position of next kmer
    #[getset(skip)]
    pub(crate) offset: usize,
}

impl Block {
    /// Read raw block
    pub fn read_raw<R>(inner: &mut R, k: u64, data_size: usize, max: u64) -> error::Result<Self>
    where
        R: std::io::Read + crate::KffRead,
    {
        let nb_kmer = if max <= 1 {
            1
        } else {
            read_nb_kmer(inner, max)? as usize
        };

        let kmer = inner
            .read_2bits(nb_kmer + k as usize - 1)?
            .into_boxed_bitslice();

        let data = inner.read_n_bytes_dyn(nb_kmer * data_size)?;

        Ok(Self {
            k,
            data_size,
            kmer,
            data,
            minimizer_offset: 0,
            offset: 0,
        })
    }

    /// Write raw block
    pub fn write_raw<W>(&self, outer: &mut W, max: u64) -> error::Result<()>
    where
        W: std::io::Write + crate::KffWrite,
    {
        if max > 1 {
            write_nb_kmer(
                outer,
                max,
                (self.kmer.len() / 2 - self.k as usize + 1) as u64,
            )?;
        }
        outer.write_bytes(self.kmer.as_raw_slice())?;
        outer.write_bytes(self.data.as_slice())?;

        Ok(())
    }

    /// Read minimizer block
    pub fn read_minimizer<R>(
        inner: &mut R,
        k: u64,
        m: u64,
        data_size: usize,
        max: u64,
        minimizer: &Kmer,
    ) -> error::Result<Self>
    where
        R: std::io::Read + crate::KffRead,
    {
        let nb_kmer = if max <= 1 {
            1
        } else {
            read_nb_kmer(inner, max)? as usize
        };
        println!("nb_kmer {}", nb_kmer);

        let minimizer_offset = read_nb_kmer(inner, std::cmp::min(k + max - 1, u64::MAX))? as usize;
        println!("{}, {}", k + max - 1, std::cmp::min(k + max - 1, u64::MAX));
        println!("{}", minimizer_offset);

        let kmer_without_minimizer = inner.read_2bits(nb_kmer + k as usize - 1 - m as usize)?;
        println!("{}", kmer_without_minimizer);

        let mut kmer = bitvec::vec::BitVec::from_bitslice(
            &kmer_without_minimizer[..(minimizer_offset as usize * 2)],
        );
        println!("{:?}", kmer);
        kmer.extend_from_bitslice(&minimizer);
        println!("{:?}", kmer);
        kmer.extend_from_bitslice(&kmer_without_minimizer[(minimizer_offset as usize * 2)..]);
        println!("{:?}", kmer);

        let data = inner.read_n_bytes_dyn((nb_kmer * data_size) as usize)?;

        Ok(Self {
            k,
            data_size,
            kmer: kmer.into_boxed_bitslice(),
            data,
            minimizer_offset,
            offset: 0,
        })
    }

    /// Write minimizer block
    pub fn write_minimizer<W>(&self, outer: &mut W, m: usize, max: u64) -> error::Result<()>
    where
        W: std::io::Write + crate::KffWrite,
    {
        if max > 1 {
            write_nb_kmer(
                outer,
                max,
                (self.kmer.len() / 2 - self.k as usize + 1) as u64,
            )?;
        }
        write_nb_kmer(
            outer,
            std::cmp::min(self.k + max - 1, u64::MAX),
            self.minimizer_offset as u64,
        )?;

        let mut kmer =
            bitvec::vec::BitVec::from_bitslice(&self.kmer[..(self.minimizer_offset as usize * 2)]);

        println!("{:?}", kmer);

        kmer.extend_from_bitslice(&self.kmer[((self.minimizer_offset + m * 2) as usize + 1)..]);

        println!("{:?}", kmer);

        kmer.resize((self.minimizer_offset + m * 2) as usize, false);

        println!("{:?}", kmer);

        outer.write_bytes(kmer.as_raw_slice())?;
        outer.write_bytes(self.data.as_slice())?;

        Ok(())
    }

    /// Get the next kmer of the block
    pub fn next_kmer(&mut self) -> std::option::Option<(Kmer, Data)> {
        println!("{} {} {}", self.offset, self.offset * 2, self.kmer.len());
        if (self.offset + self.k as usize) * 2 > self.kmer.len() {
            None
        } else {
            let k_range = self.offset * 2..(self.offset + self.k as usize) * 2;
            let d_range =
                self.offset * self.data_size as usize..(self.offset + 1) * self.data_size as usize;

            self.offset += 1;
            Some((
                bitvec::boxed::BitBox::from_bitslice(&self.kmer()[k_range]),
                self.data()[d_range].to_vec(),
            ))
        }
    }
}

impl std::iter::Iterator for Block {
    type Item = (Kmer, Data);

    fn next(&mut self) -> Option<Self::Item> {
        self.next_kmer()
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

pub(crate) fn write_nb_kmer<W>(outer: &mut W, max: u64, value: u64) -> error::Result<()>
where
    W: std::io::Write + crate::KffWrite,
{
    println!("{} {:064b} {}", max, max, max.leading_zeros());
    match max.leading_zeros() {
        0..=31 => Ok(outer.write_u64(&value)?),
        32..=47 => Ok(outer.write_u32(&(value as u32))?),
        48..=55 => Ok(outer.write_u16(&(value as u16))?),
        56..=64 => Ok(outer.write_u8(&(value as u8))?),
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

            let block = Block::read_raw(&mut readable, 5, 1, 255)?;

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

            let block = Block::read_raw(&mut readable, 5, 0, 255)?;

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

            let block = Block::read_raw(&mut readable, 5, 1, 1)?;

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
            let block = Block {
                k: 5,
                data_size: 1,
                kmer: bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1,],
                data: vec![1, 2, 3],
                minimizer_offset: 0,
                offset: 0,
            };

            let mut writable = Vec::new();

            block.write_raw(&mut writable, 255)?;

            assert_eq!(writable, vec![3, 0b00011011, 0b11110100, 1, 2, 3]);

            Ok(())
        }
    }

    mod minimizer {
        use super::*;

        #[test]
        fn full() -> error::Result<()> {
            let mut readable: &[u8] = &[3, 1, 0b00111101, 1, 2, 3];

            let minimizer = bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1];
            let block = Block::read_minimizer(&mut readable, 5, 3, 1, 200, &minimizer)?;

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
            let mut readable: &[u8] = &[3, 1, 0b00111101];

            let minimizer = bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1];
            let block = Block::read_minimizer(&mut readable, 5, 3, 0, 100, &minimizer)?;

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
            let mut readable: &[u8] = &[1, 0b00111101, 1];

            let minimizer = bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1];
            let block = Block::read_minimizer(&mut readable, 5, 3, 1, 1, &minimizer)?;

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
            let block = Block {
                k: 5,
                data_size: 1,
                kmer: bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                data: vec![1, 2, 3],
                minimizer_offset: 1,
                offset: 0,
            };

            let mut writable = Vec::new();

            block.write_minimizer(&mut writable, 3, 100)?;

            assert_eq!(writable, vec![3, 1, 0b00111101, 1, 2, 3]);

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

        write_nb_kmer(&mut writable, u8::MAX as u64, u8::MAX as u64)?;
        assert_eq!(writable, vec![255]);
        writable.clear();

        write_nb_kmer(&mut writable, (u8::MAX as u64) + 1, (u8::MAX as u64) + 1)?;
        assert_eq!(writable, vec![1, 0]);
        writable.clear();

        write_nb_kmer(&mut writable, u16::MAX as u64, u16::MAX as u64)?;
        assert_eq!(writable, vec![255, 255]);
        writable.clear();

        write_nb_kmer(&mut writable, (u16::MAX as u64) + 1, &(u16::MAX as u64) + 1)?;
        assert_eq!(writable, vec![0, 1, 0, 0]);
        writable.clear();

        write_nb_kmer(&mut writable, u32::MAX as u64, u32::MAX as u64)?;
        assert_eq!(writable, vec![255, 255, 255, 255]);
        writable.clear();

        write_nb_kmer(&mut writable, (u32::MAX as u64) + 1, (u32::MAX as u64) + 1)?;
        assert_eq!(writable, vec![0, 0, 0, 1, 0, 0, 0, 0]);
        writable.clear();

        write_nb_kmer(&mut writable, u64::MAX as u64, u64::MAX as u64)?;
        assert_eq!(writable, vec![255, 255, 255, 255, 255, 255, 255, 255]);
        writable.clear();

        Ok(())
    }
}
