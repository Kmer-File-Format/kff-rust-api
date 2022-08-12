//! Representation of a KFF kmer block

/* std use */

/* crate use */

/* project use */
use crate::error;

/// Represent a kmer
pub type Kmer = bitvec::boxed::BitBox<u8, bitvec::order::Msb0>;

/// Represent data associate to a kmer
pub type Data = Vec<u8>;

/// Trait share by two type of block
pub trait Block: std::iter::Iterator<Item = (Kmer, Data)> {
    /// Consume inner to populate block
    fn read<R>(&mut self, inner: &mut R) -> error::Result<()>
    where
        R: std::io::Read + crate::KffRead;

    /// Write content of block in outer
    fn write<W>(&self, outer: &mut W) -> error::Result<()>
    where
        W: std::io::Write + crate::KffWrite;

    /// Get the next kmer of the block
    fn next_kmer(&mut self) -> std::option::Option<(Kmer, Data)> {
        if self.offset() == self.nb_kmer() as usize {
            None
        } else {
            let k_range = self.offset() * 2..(self.offset() + self.k() as usize) * 2;
            let d_range = self.offset() * self.data_size() as usize
                ..(self.offset() + 1) * self.data_size() as usize;

            self.offset_inc();
            Some((
                bitvec::boxed::BitBox::from_bitslice(&self.kmer()[k_range]),
                self.data()[d_range].to_vec(),
            ))
        }
    }

    /// K
    fn k(&self) -> u64;

    /// Offset
    fn offset(&self) -> usize;

    /// Increment offset
    fn offset_inc(&mut self);

    /// Number of_kmer
    fn nb_kmer(&self) -> usize;

    /// Data size
    fn data_size(&self) -> usize;

    /// Kmer
    fn kmer(&self) -> &Kmer;

    /// Data
    fn data(&self) -> &Data;
}

/// Struct to represente a KFF Raw block
#[derive(getset::Getters, std::fmt::Debug)]
#[getset(get = "pub")]
pub struct Raw {
    /// Size of kmer
    pub(crate) k: u64,

    /// Maximum number in block
    pub(crate) max: u64,

    /// Size of data associate (in bytes) to each kmer
    pub(crate) data_size: usize,

    /// Number of kmer in block
    pub(crate) nb_kmer: usize,

    /// Bit field store all kmer of this block
    pub(crate) kmer: Kmer,

    /// Array store data associate to kmer of this block
    pub(crate) data: Data,

    /// Actual position of next kmer
    #[getset(skip)]
    pub(crate) offset: usize,
}

impl Raw {
    /// Create a new block
    pub fn new(k: u64, max: u64, data_size: usize) -> Self {
        Self {
            k,
            max,
            data_size,
            nb_kmer: 0,
            kmer: Kmer::default(),
            data: Data::default(),
            offset: 0,
        }
    }
}

impl Block for Raw {
    fn read<R>(&mut self, inner: &mut R) -> error::Result<()>
    where
        R: std::io::Read + crate::KffRead,
    {
        self.nb_kmer = if self.max <= 1 {
            1
        } else {
            read_nb_kmer(inner, self.max)? as usize
        };

        self.kmer = inner
            .read_2bits(self.nb_kmer + self.k as usize - 1)?
            .into_boxed_bitslice();

        self.data = inner.read_n_bytes_dyn((self.nb_kmer * self.data_size) as usize)?;

        Ok(())
    }

    fn write<W>(&self, outer: &mut W) -> error::Result<()>
    where
        W: std::io::Write + crate::KffWrite,
    {
        if self.max > 1 {
            write_nb_kmer(outer, self.max, self.nb_kmer as u64)?;
        }
        outer.write_bytes(self.kmer.as_raw_slice())?;
        outer.write_bytes(self.data.as_slice())?;

        Ok(())
    }

    fn k(&self) -> u64 {
        self.k
    }

    fn offset(&self) -> usize {
        self.offset
    }

    /// Increment offset
    fn offset_inc(&mut self) {
        self.offset += 1
    }

    /// Number of_kmer
    fn nb_kmer(&self) -> usize {
        self.nb_kmer
    }

    /// Data size
    fn data_size(&self) -> usize {
        self.data_size
    }

    /// Kmer
    fn kmer(&self) -> &Kmer {
        &self.kmer
    }

    /// Data
    fn data(&self) -> &Data {
        &self.data
    }
}

impl std::iter::Iterator for Raw {
    type Item = (Kmer, Data);

    fn next(&mut self) -> Option<Self::Item> {
        self.next_kmer()
    }
}

/// Struct to represente a KFF Raw block
#[derive(getset::Getters, std::fmt::Debug)]
#[getset(get = "pub")]
pub struct Minimizer {
    /// Size of kmer
    pub(crate) k: u64,

    /// Size of minimizer
    pub(crate) m: u64,

    /// Maximum number in block
    pub(crate) max: u64,

    /// Size of data associate (in bytes) to each kmer
    pub(crate) data_size: usize,

    /// Number of kmer in block
    pub(crate) nb_kmer: usize,

    /// Bit field store all kmer of this block
    pub(crate) kmer: Kmer,

    /// Array store data associate to kmer of this block
    pub(crate) data: Data,

    /// Minimizer sequence
    pub(crate) minimizer: Kmer,

    /// Minimizer offset
    pub(crate) minimizer_offset: usize,

    /// Actual position of next kmer
    #[getset(skip)]
    pub(crate) offset: usize,
}

impl Minimizer {
    /// Create a new Minimizer block
    pub fn new(k: u64, m: u64, max: u64, data_size: usize, minimizer: Kmer) -> Self {
        Self {
            k,
            m,
            max,
            data_size,
            nb_kmer: 0,
            kmer: Kmer::default(),
            data: Data::default(),
            minimizer,
            minimizer_offset: 0,
            offset: 0,
        }
    }
}

impl Block for Minimizer {
    fn read<R>(&mut self, inner: &mut R) -> error::Result<()>
    where
        R: std::io::Read + crate::KffRead,
    {
        self.nb_kmer = if self.max <= 1 {
            1
        } else {
            read_nb_kmer(inner, self.max)? as usize
        };
        println!("nb_kmer {}", self.nb_kmer);

        self.minimizer_offset =
            read_nb_kmer(inner, std::cmp::min(self.k + self.max - 1, u64::MAX))? as usize;
        println!(
            "{}, {}",
            self.k + self.max - 1,
            std::cmp::min(self.k + self.max - 1, u64::MAX)
        );
        println!("{}", self.minimizer_offset);

        let kmer_without_minimizer =
            inner.read_2bits(self.nb_kmer + self.k as usize - 1 - self.m as usize)?;
        println!("{}", kmer_without_minimizer);

        let mut kmer = bitvec::vec::BitVec::from_bitslice(
            &kmer_without_minimizer[..(self.minimizer_offset as usize * 2)],
        );
        println!("{:?}", kmer);
        kmer.extend_from_bitslice(&self.minimizer);
        println!("{:?}", kmer);
        kmer.extend_from_bitslice(&kmer_without_minimizer[(self.minimizer_offset as usize * 2)..]);
        println!("{:?}", kmer);
        self.kmer = kmer.into_boxed_bitslice();

        self.data = inner.read_n_bytes_dyn((self.nb_kmer * self.data_size) as usize)?;

        Ok(())
    }

    fn write<W>(&self, outer: &mut W) -> error::Result<()>
    where
        W: std::io::Write + crate::KffWrite,
    {
        if self.max > 1 {
            write_nb_kmer(outer, self.max, self.nb_kmer as u64)?;
        }
        write_nb_kmer(
            outer,
            std::cmp::min(self.k + self.max - 1, u64::MAX),
            self.minimizer_offset as u64,
        )?;

        let mut kmer =
            bitvec::vec::BitVec::from_bitslice(&self.kmer[..(self.minimizer_offset as usize * 2)]);
        kmer.extend_from_bitslice(
            &self.kmer[((self.minimizer_offset * 2 + self.minimizer.len()) as usize)..],
        );
        kmer.resize(
            (self.minimizer_offset * 2 + self.minimizer.len()) as usize,
            false,
        );

        outer.write_bytes(kmer.as_raw_slice())?;
        outer.write_bytes(self.data.as_slice())?;

        Ok(())
    }

    fn k(&self) -> u64 {
        self.k
    }

    fn offset(&self) -> usize {
        self.offset
    }

    fn offset_inc(&mut self) {
        self.offset += 1
    }

    fn nb_kmer(&self) -> usize {
        self.nb_kmer
    }

    fn data_size(&self) -> usize {
        self.data_size
    }

    fn kmer(&self) -> &Kmer {
        &self.kmer
    }

    fn data(&self) -> &Data {
        &self.data
    }
}

impl std::iter::Iterator for Minimizer {
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

            let mut block = Raw::new(5, 255, 1);
            block.read(&mut readable)?;

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

            let mut block = Raw::new(5, 255, 0);
            block.read(&mut readable)?;

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

            let mut block = Raw::new(5, 1, 1);
            block.read(&mut readable)?;

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
            let block = Raw {
                k: 5,
                max: 255,
                data_size: 1,
                nb_kmer: 3,
                kmer: bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0],
                data: vec![1, 2, 3],
                offset: 0,
            };

            let mut writable = Vec::new();

            block.write(&mut writable)?;

            assert_eq!(writable, vec![3, 0b00011011, 0b11110100, 1, 2, 3]);

            Ok(())
        }
    }

    mod minimizer {
        use super::*;

        #[test]
        fn full() -> error::Result<()> {
            let mut readable: &[u8] = &[3, 1, 0b00111101, 1, 2, 3];

            let mut block = Minimizer::new(
                5,
                3,
                200,
                1,
                bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1],
            );
            block.read(&mut readable)?;

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

            let mut block = Minimizer::new(
                5,
                3,
                100,
                0,
                bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1],
            );
            block.read(&mut readable)?;

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

            let mut block = Minimizer::new(
                5,
                3,
                1,
                1,
                bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1],
            );
            block.read(&mut readable)?;

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
            let block = Raw {
                k: 5,
                max: 255,
                data_size: 1,
                nb_kmer: 3,
                kmer: bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0],
                data: vec![1, 2, 3],
                offset: 0,
            };

            let mut writable = Vec::new();

            block.write(&mut writable)?;

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
