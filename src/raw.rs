//! Declaration of Raw section Reader and Writer

/* crate use */
use anyhow::Result;
use bitvec::prelude::*;
use byteorder::*;

/* local use */
use crate::data::Reader as DataReader;
use crate::data::Writer as DataWriter;
use crate::kff::Reader as KffReader;
use crate::kmer::Kmer;
use crate::utils::{BitBox, BitOrd, BitSlice};
use crate::variables::Variables;
use crate::variables::Variables1;
use crate::*;

/// A Raw section reader implement [data::Reader]
pub struct Reader<'input, R>
where
    R: std::io::Read,
{
    // Global
    k: u64,
    max: u64,
    data_size: u64,

    // Other
    reader: &'input mut KffReader<R>,

    // Section
    remaining_block: u32,

    // Block
    block_n: u64,
    block_seq: BitBox,
    block_data: Vec<u8>,
}

impl<'input, R> Reader<'input, R>
where
    R: std::io::Read,
{
    /// Create a new reader with a reference of kff::Reader
    pub fn new(reader: &'input mut KffReader<R>) -> Result<Self> {
        let k = reader.variables().k()?;
        let max = reader.variables().max()?;
        let data_size = reader.variables().data_size()?;

        let remaining_block = reader.input().read_u32::<utils::Order>()?;

        Ok(Self {
            k,
            max,
            data_size,
            reader,
            remaining_block,
            block_n: 0,
            block_seq: bitbox![BitOrd, u8; 0; 0],
            block_data: vec![],
        })
    }
}

impl<'input, R> Iterator for Reader<'input, R>
where
    R: std::io::Read,
{
    type Item = Result<Kmer>;

    fn next(&mut self) -> Option<Self::Item> {
        let tmp = self.next_kmer();

        match tmp {
            Err(ref _e) => Some(tmp),
            Ok(ref o) => {
                if o.is_empty() {
                    None
                } else {
                    Some(tmp)
                }
            }
        }
    }
}

impl<'input, R> DataReader<'input, R> for Reader<'input, R>
where
    R: std::io::Read,
{
    fn k(&self) -> u64 {
        self.k
    }

    fn input(&mut self) -> &mut R {
        self.reader.input()
    }

    fn max_kmer(&self) -> u64 {
        self.max
    }

    fn block_n(&self) -> u64 {
        self.block_n
    }

    fn data_size(&self) -> u64 {
        self.data_size
    }

    fn block_seq(&self) -> &BitSlice {
        &self.block_seq
    }

    fn block_data(&self) -> &[u8] {
        &self.block_data
    }

    fn rev_encoding(&self) -> u8 {
        self.reader.rev_encoding()
    }

    fn decrease_n(&mut self) {
        self.block_n -= 1;
    }

    fn read_block(&mut self) -> Result<usize> {
        if self.remaining_block == 0 {
            return Ok(0);
        }

        self.block_n = self.read_n()?;

        // Read sequence
        self.block_seq = self.read_seq(self.block_n + self.k - 1)?;

        // Read data
        self.block_data = self.read_data()?;

        self.remaining_block -= 1;

        if self.block_n == 1 {
            Ok((utils::ceil_to_8(((self.block_n + self.k - 1) * 2) / 8)
                + self.block_n * self.data_size) as usize)
        } else {
            Ok((1
                + utils::ceil_to_8(((self.block_n + self.k - 1) * 2) / 8)
                + (self.block_n * self.data_size)) as usize)
        }
    }
}

/// A Raw section writer implement [data::Writer]
pub struct Writer<'output, W>
where
    W: std::io::Write + std::io::Seek + 'output,
{
    // Global
    k: u64,
    max: u64,
    data_size: u64,

    nb_block_offset: u64,
    nb_block: u32,

    // Other
    output: &'output mut W,
    is_close: bool,
    encoding: u8,
}

impl<'output, W> Writer<'output, W>
where
    W: std::io::Write + std::io::Seek + 'output,
{
    /// Create a new Raw section writer
    pub fn new(variables: &Variables, encoding: u8, output: &'output mut W) -> Result<Self> {
        let k = variables.k()?;
        let max = variables.max()?;
        let data_size = variables.data_size()?;

        let nb_block_offset = output.seek(std::io::SeekFrom::Current(0))?;

        output.write_u32::<utils::Order>(0)?;

        Ok(Self {
            k,
            max,
            data_size,
            nb_block_offset,
            nb_block: 0,
            output,
            is_close: false,
            encoding,
        })
    }

    /// Write a raw block
    pub fn write_block(&mut self, seq: &BitSlice, data: &[u8]) -> Result<usize> {
        self.increment_nb_block()?;

        let nb_kmer = self.check_block(seq.len(), data.len())? as u64;

        let mut bytes_write = 0;

        if self.max != 1 {
            utils::write_dynamic_size_field(self.output, nb_kmer, self.max)?;
            bytes_write += utils::bytes_to_store_n(self.max) as usize;
        }

        let mut write_seq = if seq.len() % 8 == 0 {
            bitvec![Msb0, u8; 0; 0]
        } else {
            bitvec![Msb0, u8; 0; 8 - (seq.len() % 8)]
        };

        write_seq.extend(seq);
        self.output.write_all(write_seq.as_raw_slice())?;
        bytes_write += seq.as_slice().len();

        self.output.write_all(data)?;
        bytes_write += data.len();

        Ok(bytes_write)
    }

    /// Write a raw block, where sequence is encode in ASCII
    pub fn write_seq_block(&mut self, seq: &[u8], data: &[u8]) -> Result<usize> {
        self.write_block(&utils::seq2bits(seq, self.encoding)[..], data)
    }
}

impl<'output, W> DataWriter<'output, W> for Writer<'output, W>
where
    W: std::io::Write + std::io::Seek + 'output,
{
    // Getter
    fn data_size(&self) -> u64 {
        self.data_size
    }

    fn max(&self) -> u64 {
        self.max
    }

    fn nb_block_offset(&self) -> u64 {
        self.nb_block_offset
    }

    fn nb_block(&self) -> u32 {
        self.nb_block
    }

    fn is_close(&self) -> bool {
        self.is_close
    }

    fn output(&mut self) -> &mut W {
        self.output
    }

    // Setter
    fn set_nb_block(&mut self, value: u32) {
        self.nb_block = value;
    }

    fn set_close(&mut self, value: bool) {
        self.is_close = value;
    }

    // Computation
    fn nb_kmer(&self, seq_len: usize) -> usize {
        (seq_len / 2) - self.k as usize + 1
    }
}

impl<'output, W> Drop for Writer<'output, W>
where
    W: std::io::Write + std::io::Seek,
{
    fn drop(&mut self) {
        self.close().unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq2bits::Bits2Nuc;
    use crate::variables::Variables;

    #[test]
    fn init() {
        {
            let mut input: &[u8] = &[
                1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
                0, 0, 0, 1,
            ];

            let mut reader = KffReader::new(&mut input).unwrap();

            assert!(super::Reader::new(&mut reader).is_err());
        }

        {
            let mut input: &[u8] = &[
                1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
                0, 0, 0, 1,
            ];

            let mut reader = KffReader::new(&mut input).unwrap();

            reader.variables().insert("k".to_string(), 15);
            assert!(super::Reader::new(&mut reader).is_err());
        }

        {
            let mut input: &[u8] = &[
                1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
                0, 0, 0, 1,
            ];

            let mut reader = KffReader::new(&mut input).unwrap();

            reader.variables().insert("k".to_string(), 15);
            reader.variables().insert("max".to_string(), 256);
            assert!(super::Reader::new(&mut reader).is_err());
        }

        {
            let mut input: &[u8] = &[
                1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
                0, 0, 0, 1,
            ];

            let mut reader = KffReader::new(&mut input).unwrap();

            reader.variables().insert("k".to_string(), 15);
            reader.variables().insert("max".to_string(), 256);
            reader.variables().insert("data_size".to_string(), 0);
            assert!(super::Reader::new(&mut reader).is_ok());
        }
    }

    #[test]
    fn next_kmer() {
        let mut block: &[u8] = &[
            1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
            1, 0, 0, 0, 5, 0b00000011, 0b11011110, 0b00111101, 1, 2, 3, 4,
            5, // One block, 5 kmer in block,
        ];

        let mut reader = KffReader::new(&mut block).unwrap();
        reader.variables().insert("k".to_string(), 5);
        reader.variables().insert("max".to_string(), 5);
        reader.variables().insert("data_size".to_string(), 1);

        let mut raw = super::Reader::new(&mut reader).unwrap();

        assert_eq!(
            raw.next_kmer().unwrap(),
            Kmer::new(
                bitbox![Msb0, u8; 1, 1, 1, 1, 0, 1, 1, 1, 1, 0],
                vec![1].into_boxed_slice()
            )
        );

        assert_eq!(
            raw.next_kmer().unwrap(),
            Kmer::new(
                bitbox![Msb0, u8; 1, 1, 0, 1, 1, 1, 1, 0, 0, 0],
                vec![2].into_boxed_slice()
            )
        );

        assert_eq!(
            raw.next_kmer().unwrap(),
            Kmer::new(
                bitbox![Msb0, u8; 0, 1, 1, 1, 1, 0, 0, 0, 1, 1],
                vec![3].into_boxed_slice()
            )
        );

        assert_eq!(
            raw.next_kmer().unwrap(),
            Kmer::new(
                bitbox![Msb0, u8; 1, 1, 1, 0, 0, 0, 1, 1, 1, 1],
                vec![4].into_boxed_slice()
            )
        );

        assert_eq!(
            raw.next_kmer().unwrap(),
            Kmer::new(
                bitbox![Msb0, u8; 1, 0, 0, 0, 1, 1, 1, 1, 0, 1],
                vec![5].into_boxed_slice()
            )
        );

        assert_eq!(
            raw.next_kmer().unwrap(),
            Kmer::new(bitbox![Msb0, u8; 0; 0], vec![].into_boxed_slice())
        );

        assert_eq!(
            raw.next_kmer().unwrap(),
            Kmer::new(bitbox![Msb0, u8; 0; 0], vec![].into_boxed_slice())
        );
    }

    #[test]
    fn bin() {
        let mut block: &[u8] = &[
            1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
            1, 0, 0, 0, 5, 0b00000011, 0b11011110, 0b00111101, 1, 2, 3, 4,
            5, // One block, 5 kmer in block,
        ];

        let mut reader = KffReader::new(&mut block).unwrap();
        reader.variables().insert("k".to_string(), 5);
        reader.variables().insert("max".to_string(), 5);
        reader.variables().insert("data_size".to_string(), 1);

        let raw = super::Reader::new(&mut reader).unwrap();

        let mut it = raw.into_iter();

        assert_eq!(
            it.next().unwrap().unwrap(),
            Kmer::new(
                bitbox![Msb0, u8; 1, 1, 1, 1, 0, 1, 1, 1, 1, 0],
                vec![1].into_boxed_slice()
            )
        );

        assert_eq!(
            it.next().unwrap().unwrap(),
            Kmer::new(
                bitbox![Msb0, u8; 1, 1, 0, 1, 1, 1, 1, 0, 0, 0],
                vec![2].into_boxed_slice()
            )
        );

        assert_eq!(
            it.next().unwrap().unwrap(),
            Kmer::new(
                bitbox![Msb0, u8; 0, 1, 1, 1, 1, 0, 0, 0, 1, 1],
                vec![3].into_boxed_slice()
            )
        );

        assert_eq!(
            it.next().unwrap().unwrap(),
            Kmer::new(
                bitbox![Msb0, u8; 1, 1, 1, 0, 0, 0, 1, 1, 1, 1],
                vec![4].into_boxed_slice()
            )
        );

        assert_eq!(
            it.next().unwrap().unwrap(),
            Kmer::new(
                bitbox![Msb0, u8; 1, 0, 0, 0, 1, 1, 1, 1, 0, 1],
                vec![5].into_boxed_slice()
            )
        );

        assert!(it.next().is_none());

        assert!(it.next().is_none());
    }

    #[test]
    fn seq() {
        let mut block: &[u8] = &[
            1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
            1, 0, 0, 0, 5, 0b00000011, 0b11011110, 0b00111101, 1, 2, 3, 4,
            5, // One block, 5 kmer in block,
        ];

        let mut reader = KffReader::new(&mut block).unwrap();
        reader.variables().insert("k".to_string(), 5);
        reader.variables().insert("max".to_string(), 5);
        reader.variables().insert("data_size".to_string(), 1);

        let raw = super::Reader::new(&mut reader).unwrap();

        let mut it = raw.into_iter();

        let mut value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq().into_nuc(utils::rev_encoding(0b00011011)),
            vec![b'G', b'G', b'C', b'G', b'T'].into_boxed_slice(),
        );
        assert_eq!(value.data(), &[1]);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq().into_nuc(utils::rev_encoding(0b00011011)),
            vec![b'G', b'C', b'G', b'T', b'A'].into_boxed_slice(),
        );
        assert_eq!(value.data(), &[2]);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq().into_nuc(utils::rev_encoding(0b00011011)),
            vec![b'C', b'G', b'T', b'A', b'G'].into_boxed_slice(),
        );
        assert_eq!(value.data(), &[3]);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq().into_nuc(utils::rev_encoding(0b00011011)),
            vec![b'G', b'T', b'A', b'G', b'G'].into_boxed_slice(),
        );
        assert_eq!(value.data(), &[4]);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq().into_nuc(utils::rev_encoding(0b00011011)),
            vec![b'T', b'A', b'G', b'G', b'C'].into_boxed_slice(),
        );
        assert_eq!(value.data(), &[5]);

        assert!(it.next().is_none());

        assert!(it.next().is_none());
    }

    #[test]
    fn write() {
        let mut variables = Variables::new();
        variables.insert("k".to_string(), 5);
        variables.insert("max".to_string(), 5);
        variables.insert("data_size".to_string(), 1);

        let mut out = vec![0u8; 0];
        let mut output = std::io::Cursor::new(&mut out);
        {
            let mut writer = Writer::new(&variables, 0b00011011, &mut output).unwrap();

            writer.write_seq_block(b"GGCGTAG", &[10, 8, 9]).unwrap();
            writer.write_seq_block(b"GCGAT", &[1]).unwrap();

            writer.close().unwrap();
        }

        assert_eq!(
            out,
            [2, 0, 0, 0, 3, 0b00111101, 0b11100011, 10, 8, 9, 1, 0b00000011, 0b01110010, 1]
        );
    }

    #[test]
    fn write_data_size_0() {
        let mut variables = Variables::new();
        variables.insert("k".to_string(), 5);
        variables.insert("max".to_string(), 5);
        variables.insert("data_size".to_string(), 0);

        let mut out = vec![0u8; 0];
        let mut output = std::io::Cursor::new(&mut out);
        {
            let mut writer = Writer::new(&variables, 0b00011011, &mut output).unwrap();

            writer.write_seq_block(b"GGCGTAG", &[]).unwrap();
            writer.write_seq_block(b"GCGAT", &[]).unwrap();

            writer.close().unwrap();
        }

        assert_eq!(
            out,
            [2, 0, 0, 0, 3, 0b00111101, 0b11100011, 1, 0b00000011, 0b01110010]
        );
    }

    #[test]
    fn write_n1() {
        let mut variables = Variables::new();
        variables.insert("k".to_string(), 5);
        variables.insert("max".to_string(), 1);
        variables.insert("data_size".to_string(), 1);

        let mut output = std::io::Cursor::new(vec![0u8; 0]);
        {
            let mut writer = Writer::new(&variables, 0b00011011, &mut output).unwrap();

            writer.write_seq_block(b"GAGTT", &[10]).unwrap();
            writer.write_seq_block(b"GCGAT", &[1]).unwrap();

            writer.close().unwrap();
        }

        assert_eq!(
            output.into_inner(),
            [2, 0, 0, 0, 0b00000011, 0b00111010, 10, 0b00000011, 0b01110010, 1]
        );
    }

    #[test]
    fn write_read() {
        let mut variables = Variables::new();
        variables.insert("k".to_string(), 5);
        variables.insert("max".to_string(), 255);
        variables.insert("data_size".to_string(), 1);

        let mut buffer = std::io::Cursor::new(vec![1u8, 0, 30, 0, 0, 0, 0]);
        buffer.set_position(7);
        {
            let mut writer = Writer::new(&variables, 0b00011011, &mut buffer).unwrap();

            writer.write_seq_block(b"GAGTTAC", &[10, 8, 9]).unwrap();
            writer.write_seq_block(b"GCGAT", &[1]).unwrap();
        }

        let inp = buffer.into_inner();
        assert_eq!(
            inp,
            [
                1, 0, 30, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0b00110011, 0b10100001, 10, 8, 9, 1,
                0b00000011, 0b01110010, 1
            ]
        );

        let mut input = std::io::Cursor::new(inp);

        let mut reader = KffReader::new(&mut input).unwrap();
        reader.variables().insert("k".to_string(), 5);
        reader.variables().insert("max".to_string(), 5);
        reader.variables().insert("data_size".to_string(), 1);

        let raw = super::Reader::new(&mut reader).unwrap();
        let mut it = raw.into_iter();

        let mut value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq().into_nuc(utils::rev_encoding(0b00011011)),
            vec![b'G', b'A', b'G', b'T', b'T'].into_boxed_slice(),
        );
        assert_eq!(value.data(), &[10]);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq().into_nuc(utils::rev_encoding(0b00011011)),
            vec![b'A', b'G', b'T', b'T', b'A'].into_boxed_slice(),
        );

        assert_eq!(value.data(), &[8]);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq().into_nuc(utils::rev_encoding(0b00011011)),
            vec![b'G', b'T', b'T', b'A', b'C'].into_boxed_slice(),
        );

        assert_eq!(value.data(), &[9]);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq().into_nuc(utils::rev_encoding(0b00011011)),
            vec![b'G', b'C', b'G', b'A', b'T'].into_boxed_slice(),
        );

        assert_eq!(value.data(), &[1]);

        assert!(it.next().is_none());
    }
}
