//! Read and Write Minimizer section

/* std use */

/* crate use */

/* project use */
use crate::error;
use crate::section;
use crate::{Kmer, Seq2Bit};

/// Struct to Read and Write Raw section
#[derive(getset::Getters, getset::Setters, getset::MutGetters, std::default::Default)]
#[getset(get = "pub", set = "pub", get_mut = "pub")]
pub struct Minimizer {
    /// Size of kmer
    k: u64,

    /// Size of minimizer
    m: u64,

    /// Kmer are sort by lexicographic order
    ordered: bool,

    /// Max number of kmer per block
    max: u64,

    /// Size in bytes of data associate to each kmer
    data_size: u64,
}

impl Minimizer {
    /// Intialize a minimizer section with Values
    pub fn new(values: &section::Values) -> error::Result<Self> {
        Ok(Self {
            k: values
                .get("k")
                .cloned()
                .ok_or_else(|| error::Kff::FieldIsMissing("k".to_string()))?,
            m: values
                .get("m")
                .cloned()
                .ok_or_else(|| error::Kff::FieldIsMissing("m".to_string()))?,
            ordered: true,
            // ordered: values
            //     .get("ordered")
            //     .cloned()
            //     .ok_or_else(|| error::Kff::FieldIsMissing("ordered".to_string()))
            //     .map(|x| x != 0)?,
            max: values
                .get("max")
                .cloned()
                .ok_or_else(|| error::Kff::FieldIsMissing("max".to_string()))?,
            data_size: values
                .get("data_size")
                .cloned()
                .ok_or_else(|| error::Kff::FieldIsMissing("data_size".to_string()))?,
        })
    }

    /// Read a Minimizer section, section flag must be already read
    pub fn read<R>(&self, inner: &mut R) -> error::Result<Vec<Kmer>>
    where
        R: std::io::Read + crate::KffRead,
    {
        let mut output = Vec::new();

        let minimizer = inner.read_2bits(self.m as usize)?.into_boxed_bitslice();

        let nb_block = inner.read_u64()?;

        for _ in 0..nb_block {
            let block = section::Block::read_minimizer(
                inner,
                self.k,
                self.m,
                self.data_size as usize,
                self.max,
                &minimizer,
            )?;

            output.extend(block);
        }

        Ok(output)
    }

    /// Write a Raw section, section flag isn't read
    pub fn write<W>(
        &self,
        outer: &mut W,
        minimizer: Seq2Bit,
        blocks: &[section::block::Block],
    ) -> error::Result<()>
    where
        W: std::io::Write + crate::KffWrite,
    {
        outer.write_bytes(minimizer.as_raw_slice())?;
        outer.write_u64(&(blocks.len() as u64))?;

        for block in blocks {
            block.write_minimizer(outer, self.m as usize, self.max)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use section::values::AbcValues;

    #[test]
    fn creation() -> error::Result<()> {
        let mut values = section::Values::with_capacity(4);
        assert!(Minimizer::new(&values).is_err());

        values.insert("k".to_string(), 5);
        assert!(Minimizer::new(&values).is_err());

        values.insert("m".to_string(), 2);
        assert!(Minimizer::new(&values).is_err());

        values.insert("ordered".to_string(), false as u64);
        assert!(Minimizer::new(&values).is_err());

        values.insert("max".to_string(), 255);
        assert!(Minimizer::new(&values).is_err());

        values.insert("data_size".to_string(), 1);
        assert!(Minimizer::new(&values).is_ok());

        Ok(())
    }

    #[test]
    fn read() -> error::Result<()> {
        let mut values = section::Values::with_capacity(5);

        values.insert("k".to_string(), 5);
        values.insert("m".to_string(), 3);
        values.insert("ordered".to_string(), false as u64);
        values.insert("max".to_string(), 100);
        values.insert("data_size".to_string(), 1);

        let minimizer = Minimizer::new(&values)?;

        let mut data: &[u8] = &[
            0b01101100, // minimizer sequence
            0, 0, 0, 0, 0, 0, 0, 3, // number of block
            3, 1, 0b00111101, 1, 2, 3, // one block with 3 kmer and 1 bytes data
            2, 1, 0b00111111, 1, 2, // one block with 2 kmer and 1 bytes data
            1, 1, 0b00110000, 1, // one block with 1 kmer and 1 bytes data
        ];

        let kmers = minimizer.read(&mut data)?;

        assert_eq!(
            kmers,
            vec![
                Kmer::new(
                    bitvec::bitbox![u8, bitvec::order::Msb0;  0, 0, 0, 1, 1, 0, 1, 1, 1, 1],
                    vec![1]
                ),
                Kmer::new(
                    bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1, 1, 1, 1, 1],
                    vec![2]
                ),
                Kmer::new(
                    bitvec::bitbox![u8, bitvec::order::Msb0; 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                    vec![3]
                ),
                Kmer::new(
                    bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1],
                    vec![1]
                ),
                Kmer::new(
                    bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1, 1, 1, 1, 1],
                    vec![2]
                ),
                Kmer::new(
                    bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1],
                    vec![1]
                )
            ]
        );

        Ok(())
    }

    #[test]
    fn write() -> error::Result<()> {
        let mut values = section::Values::with_capacity(4);

        values.insert("k".to_string(), 5);
        values.insert("m".to_string(), 3);
        values.insert("ordered".to_string(), false as u64);
        values.insert("max".to_string(), 100);
        values.insert("data_size".to_string(), 1);

        let minimizer = Minimizer::new(&values)?;

        let mut writable = Vec::new();

        minimizer.write(
            &mut writable,
	    bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1],
            &[
                section::block::Block{
                    k: 5,
                    data_size: 1,
                    kmer: Kmer::new(bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                     vec![1, 2, 3]),
		    minimizer_offset: 1,
		    offset: 0,
                },
                section::block::Block {
		    k: 5,
		    data_size: 1,
		    kmer: Kmer::new(bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1],
                    vec![1, 2]),
		    minimizer_offset: 1,
		    offset: 0,
                },
                section::block::Block {
		    k: 5,
                    data_size: 1,
		    kmer: Kmer::new(bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1],
		    vec![1]),
		    minimizer_offset: 1,
		    offset: 0,
                }
            ],
        )?;

        assert_eq!(
            writable,
            vec![
                0b01101100, // minimizer sequence
                0, 0, 0, 0, 0, 0, 0, 3, // number of block
                3, 1, 0b00111101, 1, 2, 3, // one block with 3 kmer and 1 bytes data
                2, 1, 0b00111100, 1, 2, // one block with 2 kmer and 1 bytes data
                1, 1, 0b00110000, 1, // one block with 1 kmer and 1 bytes data
            ]
        );

        Ok(())
    }

    #[test]
    fn write_single_block() -> error::Result<()> {
        let mut values = section::Values::with_capacity(4);

        values.insert("k".to_string(), 5);
        values.insert("m".to_string(), 4);
        values.insert("ordered".to_string(), false as u64);
        values.insert("max".to_string(), 100);
        values.insert("data_size".to_string(), 1);

        let minimizer = Minimizer::new(&values)?;

        let mut writable = Vec::new();

        let minimizer_val = bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1, 0, 1];
        let block = section::block::Block {
            k: 5,
            data_size: 1,
            kmer: Kmer::new(
                bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1],
                vec![1, 2, 3],
            ),
            minimizer_offset: 4,
            offset: 0,
        };
        minimizer.write(&mut writable, minimizer_val, &[block])?;

        assert_eq!(
            writable,
            vec![
                0b01101101, // minimizer sequence
                0, 0, 0, 0, 0, 0, 0, 1, // number of block
                5, 4, 0b00101101, 0b11000000, 1, 2, // kmer without minimizer
                3, // one block with 3 kmer and 1 bytes data
            ]
        );

        Ok(())
    }
}
