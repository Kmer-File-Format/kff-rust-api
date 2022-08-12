//! Read and Write Minimizer section

/* std use */

/* crate use */

/* project use */
use crate::error;
use crate::section;
use crate::section::block::Block;

/// Struct to Read and Write Raw section
#[derive(getset::Getters, getset::Setters, getset::MutGetters)]
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

    /// Read a Minimizer section
    pub fn read<R>(
        &self,
        inner: &mut R,
    ) -> error::Result<Vec<(section::block::Kmer, section::block::Data)>>
    where
        R: std::io::Read + crate::KffRead,
    {
        let mut output = Vec::new();

        let minimizer = inner.read_2bits(self.m as usize)?;
        println!("{:?}", minimizer);

        let nb_block = inner.read_u64()?;
        println!("{:?}", nb_block);

        for _ in 0..nb_block {
            let mut block = section::block::Minimizer::new(
                self.k,
                self.m,
                self.max,
                self.data_size as usize,
                minimizer.clone().into_boxed_bitslice(),
            );
            block.read(inner)?;
            println!("{:?}", block);

            for (kmer, data) in block {
                output.push((kmer, data))
            }
        }

        Ok(output)
    }

    /// Write a Raw section
    pub fn write<W>(
        &self,
        outer: &mut W,
        minimizer: section::block::Kmer,
        blocks: Vec<section::block::Minimizer>,
    ) -> error::Result<()>
    where
        W: std::io::Write + crate::KffWrite,
    {
        outer.write_bytes(minimizer.as_raw_slice())?;
        outer.write_u64(&(blocks.len() as u64))?;

        for block in blocks {
            block.write(outer)?;
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
                (
                    bitvec::bitbox![u8, bitvec::order::Msb0;  0, 0, 0, 1, 1, 0, 1, 1, 1, 1],
                    vec![1]
                ),
                (
                    bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1, 1, 1, 1, 1],
                    vec![2]
                ),
                (
                    bitvec::bitbox![u8, bitvec::order::Msb0; 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                    vec![3]
                ),
                (
                    bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1],
                    vec![1]
                ),
                (
                    bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1, 1, 1, 1, 1],
                    vec![2]
                ),
                (
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
            vec![
                section::block::Minimizer{
                    k: 5,
		    m: 3,
		    max: 100,
                    data_size: 1,
                    nb_kmer: 3,
                    kmer: bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0],
                    data: vec![1, 2, 3],
		    minimizer: bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1],
		    minimizer_offset: 1,
		    offset: 0,
                },
                section::block::Minimizer {
		    k: 5,
		    m: 3,
		    max: 100,                    data_size: 1,
                    nb_kmer: 2,
		    kmer: bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
                    data: vec![1, 2],
		    minimizer: bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1],
		    minimizer_offset: 1,
		    offset: 0,
                },
                section::block::Minimizer {
		    k: 5,
		    m: 3,
		    max: 100,
                    data_size: 1,
                    nb_kmer: 1,
		    kmer: bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
		    data: vec![1],
		    minimizer: bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1],
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
}