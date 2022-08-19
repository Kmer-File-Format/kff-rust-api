//! Read and Write Raw section

/* std use */

/* crate use */

/* project use */
use crate::error;
use crate::section;
use crate::Kmer;

/// Struct to Read and Write Raw section
#[derive(getset::Getters, getset::Setters, getset::MutGetters, std::default::Default)]
#[getset(get = "pub", set = "pub", get_mut = "pub")]
pub struct Raw {
    /// Size of kmer
    k: u64,

    /// Kmer are sort by lexicographic order
    ordered: bool,

    /// Max number of kmer per block
    max: u64,

    /// Size in bytes of data associate to each kmer
    data_size: u64,
}

impl Raw {
    /// Intialize a raw section with Values
    pub fn new(values: &section::Values) -> error::Result<Self> {
        Ok(Self {
            k: values
                .get("k")
                .cloned()
                .ok_or_else(|| error::Kff::FieldIsMissing("k".to_string()))?,
            ordered: values
                .get("ordered")
                .cloned()
                .ok_or_else(|| error::Kff::FieldIsMissing("ordered".to_string()))
                .map(|x| x != 0)?,
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

    /// Read a Raw section, section flag must be already read
    pub fn read<R>(&self, inner: &mut R) -> error::Result<Vec<Kmer>>
    where
        R: std::io::Read + crate::KffRead,
    {
        let mut output = Vec::new();

        let nb_block = inner.read_u64()?;

        for _ in 0..nb_block {
            let block =
                section::block::Block::read_raw(inner, self.k, self.data_size as usize, self.max)?;

            output.extend(block);
        }

        Ok(output)
    }

    /// Write a Raw section, section flag isn't write
    pub fn write<W>(&self, outer: &mut W, blocks: Vec<section::block::Block>) -> error::Result<()>
    where
        W: std::io::Write + crate::KffWrite,
    {
        outer.write_u64(&(blocks.len() as u64))?;

        for block in blocks {
            block.write_raw(outer, self.max)?;
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
        assert!(Raw::new(&values).is_err());

        values.insert("k".to_string(), 5);
        assert!(Raw::new(&values).is_err());

        values.insert("ordered".to_string(), false as u64);
        assert!(Raw::new(&values).is_err());

        values.insert("max".to_string(), 255);
        assert!(Raw::new(&values).is_err());

        values.insert("data_size".to_string(), 1);
        assert!(Raw::new(&values).is_ok());

        Ok(())
    }

    #[test]
    fn read() -> error::Result<()> {
        let mut values = section::Values::with_capacity(4);

        values.insert("k".to_string(), 5);
        values.insert("ordered".to_string(), false as u64);
        values.insert("max".to_string(), 255);
        values.insert("data_size".to_string(), 1);

        let raw = Raw::new(&values)?;

        let mut data: &[u8] = &[
            0, 0, 0, 0, 0, 0, 0, 3, // number of block
            3, 0b00011011, 0b11110100, 1, 2, 3, // one block with 3 kmer and 1 bytes data
            2, 0b00011011, 0b11110000, 1, 2, // one block with 2 kmer and 1 bytes data
            1, 0b00011011, 0b11000000, 1, // one block with 1 kmer and 1 bytes data
        ];

        let kmers = raw.read(&mut data)?;

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
        values.insert("ordered".to_string(), false as u64);
        values.insert("max".to_string(), 255);
        values.insert("data_size".to_string(), 1);

        let raw = Raw::new(&values)?;

        let mut writable = Vec::new();

        raw.write(
            &mut writable,
            vec![
                section::block::Block {
                    k: 5,
                    data_size: 1,
                    kmer: Kmer::new(bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                    vec![1, 2, 3]),
		    minimizer_offset: 0,
		    offset: 0,
                },
                section::block::Block{
                    k: 5,
                    data_size: 1,
                    kmer: Kmer::new(bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1],
                    vec![1, 2]),
		    minimizer_offset: 0,
		    offset: 0,
                },
                section::block::Block {
                    k: 5,
		    data_size: 1,
                    kmer: Kmer::new(bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1],
                    vec![1]),
		    minimizer_offset: 0,
		    offset: 0,
                },
            ],
        )?;

        assert_eq!(
            writable,
            vec![
                0, 0, 0, 0, 0, 0, 0, 3, // number of block
                3, 0b00011011, 0b11110100, 1, 2, 3, // one block with 3 kmer and 1 bytes data
                2, 0b00011011, 0b11110000, 1, 2, // one block with 2 kmer and 1 bytes data
                1, 0b00011011, 0b11000000, 1, // one block with 1 kmer and 1 bytes data
            ]
        );

        Ok(())
    }
}
