/* crate use */
use bitvec::prelude::*;
use byteorder::*;

/* local use */
use crate::data::Reader as DataReader;
use crate::data::Writer as DataWriter;
use crate::error;
use crate::kff::Reader as KffReader;
use crate::kmer::Kmer;
use crate::utils::{BitBox, BitSlice, BitVec};
use crate::variables::BaseVariables;
use crate::variables::Variables;
use crate::*;

pub struct Reader<'input, R>
where
    R: std::io::Read,
{
    // Global
    k: u64,
    m: u64,
    max: u64,
    data_size: u64,

    // Other
    reader: &'input mut KffReader<R>,

    // Section
    minimizer: BitBox,
    remaining_block: u32,

    // Block
    block_n: u64,
    block_idx: u64,
    block_seq: BitBox,
    block_data: Vec<u8>,
}

impl<'input, R> Reader<'input, R>
where
    R: std::io::Read,
{
    pub fn new(reader: &'input mut KffReader<R>) -> crate::Result<Self> {
        let k = reader.variables().k()?;
        let m = reader.variables().m()?;
        let max = reader.variables().max()?;
        let data_size = reader.variables().data_size()?;

        let mut buffer = vec![0u8; utils::ceil_to_8(m * 2) as usize / 8];
        reader.input().read_exact(&mut buffer)?;
        let mut tmp = BitVec::from_vec(buffer);
        tmp.resize((m * 2) as usize, false);
        let minimizer = tmp.into_boxed_bitslice();

        let remaining_block = reader.input().read_u32::<LittleEndian>()?;

        Ok(Self {
            k,
            m,
            max,
            data_size,
            reader,
            minimizer,
            remaining_block,
            block_n: 0,
            block_idx: 0,
            block_seq: bitbox![Lsb0, u8; 0; 0],
            block_data: vec![],
        })
    }
}

impl<'input, R> Iterator for Reader<'input, R>
where
    R: std::io::Read,
{
    type Item = crate::Result<Kmer>;

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

    fn read_block(&mut self) -> crate::Result<usize> {
        if self.remaining_block == 0 {
            return Ok(0);
        }

        self.block_n = self.read_n()?;

        self.block_idx = utils::read_dynamic_size_field(
            self.reader.input(),
            std::cmp::min(self.k + self.m - 1, u64::MAX),
        )?;

        let seq_without_mini = self.read_seq(self.block_n + self.k - 1 - self.m)?;

        let mut seq = bitvec![Lsb0, u8; 0; 0];
        seq.extend(&seq_without_mini[0..(self.block_idx as usize * 2)]);
        seq.extend(self.minimizer.iter());
        seq.extend(&seq_without_mini[(self.block_idx as usize * 2)..]);

        self.block_seq = seq.into_boxed_bitslice();

        self.block_data = self.read_data()?;

        self.remaining_block -= 1;

        if self.block_n == 1 {
            Ok(
                (utils::ceil_to_8(((self.block_n + self.k - 1 - self.m) * 2) / 8)
                    + self.block_n * self.data_size) as usize,
            )
        } else {
            Ok((1
                + utils::ceil_to_8(((self.block_n + self.k - 1 - self.m) * 2) / 8)
                + (self.block_n * self.data_size)) as usize)
        }
    }
}

pub struct Writer<'output, W>
where
    W: std::io::Write + std::io::Seek + 'output,
{
    // Global
    k: u64,
    m: u64,
    max: u64,
    data_size: u64,

    minimizer_len: usize,
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
    pub fn new(
        variables: &Variables,
        minimizer: &[u8],
        encoding: u8,
        output: &'output mut W,
    ) -> crate::Result<Self> {
        let k = *variables.get("k").ok_or(error::Data::KMissing)?;
        let m = *variables.get("m").ok_or(error::Data::MMissing)?;
        let max = *variables.get("max").ok_or(error::Data::MaxMissing)?;
        let data_size = *variables
            .get("data_size")
            .ok_or(error::Data::DataSizeMissing)?;

        if m != minimizer.len() as u64 {
            return Err(Box::new(error::Minimizer::MinimizerSizeMDiff));
        }

        output.write_all(utils::seq2bits(minimizer, encoding).as_slice())?;
        let nb_block_offset = output.seek(std::io::SeekFrom::Current(0))?;
        output.write_u32::<LittleEndian>(0)?;

        Ok(Self {
            k,
            m,
            max,
            data_size,
            minimizer_len: minimizer.len(),
            nb_block_offset,
            nb_block: 0,
            output,
            is_close: false,
            encoding,
        })
    }

    pub fn write_block(
        &mut self,
        minimizer_idx: u64,
        seq: &BitSlice,
        data: &[u8],
    ) -> crate::Result<usize> {
        self.increment_nb_block()?;

        let nb_kmer = self.check_block(seq.len(), data.len())? as u64;

        let mut bytes_write = 0;

        if self.max != 1 {
            utils::write_dynamic_size_field(self.output, nb_kmer, self.max)?;
            bytes_write += utils::bytes_to_store_n(self.max) as usize;
        }

        utils::write_dynamic_size_field(
            self.output,
            minimizer_idx,
            std::cmp::min(self.k + self.m - 1, u64::MAX),
        )?;
        bytes_write +=
            utils::bytes_to_store_n(std::cmp::min(self.k + self.m - 1, u64::MAX)) as usize;

        self.output.write_all(seq.as_slice())?;
        bytes_write += seq.as_slice().len();

        self.output.write_all(data)?;
        bytes_write += data.len();

        Ok(bytes_write)
    }

    pub fn write_seq_block(
        &mut self,
        minimizer_idx: u64,
        seq: &[u8],
        data: &[u8],
    ) -> crate::Result<usize> {
        self.write_block(
            minimizer_idx,
            &utils::seq2bits(seq, self.encoding)[..],
            data,
        )
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
        (seq_len) / 2 + self.minimizer_len - self.k as usize + 1
    }
}

impl<'output, W> Drop for Writer<'output, W>
where
    W: std::io::Write + std::io::Seek + 'output,
{
    fn drop(&mut self) {
        self.close().unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::Reader;

    #[test]
    fn init() {
        {
            let mut input: &[u8] = &[
                1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
                0b10100101, 0, 0, 0, 4, // minimizer -> CCTT, 4 block
            ];

            let mut reader = KffReader::new(&mut input).unwrap();

            assert!(super::Reader::new(&mut reader).is_err());
        }

        {
            let mut input: &[u8] = &[
                1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
                0b10100101, 0, 0, 0, 4, // minimizer -> CCTT, 4 block
            ];

            let mut reader = KffReader::new(&mut input).unwrap();
            reader.variables().insert("k".to_string(), 15);

            assert!(super::Reader::new(&mut reader).is_err());
        }

        {
            let mut input: &[u8] = &[
                1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
                0b10100101, 0, 0, 0, 4, // minimizer -> CCTT, 4 block
            ];

            let mut reader = KffReader::new(&mut input).unwrap();
            reader.variables().insert("k".to_string(), 15);
            reader.variables().insert("m".to_string(), 4);

            assert!(super::Reader::new(&mut reader).is_err());
        }

        {
            let mut input: &[u8] = &[
                1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
                0b10100101, 0, 0, 0, 4, // minimizer -> CCTT, 4 block
            ];

            let mut reader = KffReader::new(&mut input).unwrap();
            reader.variables().insert("k".to_string(), 15);
            reader.variables().insert("m".to_string(), 4);
            reader.variables().insert("max".to_string(), 255);

            assert!(super::Reader::new(&mut reader).is_err());
        }

        {
            let mut input: &[u8] = &[
                1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
                0b10100101, 0, 0, 0, 4, // minimizer -> CCTT, 4 block
            ];

            let mut reader = KffReader::new(&mut input).unwrap();
            reader.variables().insert("k".to_string(), 15);
            reader.variables().insert("m".to_string(), 4);
            reader.variables().insert("max".to_string(), 255);

            assert!(super::Reader::new(&mut reader).is_err());
        }

        {
            let mut input: &[u8] = &[
                1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
                0b10100101, 0, 0, 0, 4, // minimizer -> CCTT, 4 block
            ];

            let mut reader = KffReader::new(&mut input).unwrap();
            reader.variables().insert("k".to_string(), 15);
            reader.variables().insert("m".to_string(), 4);
            reader.variables().insert("max".to_string(), 255);
            reader.variables().insert("data_size".to_string(), 0);

            assert!(super::Reader::new(&mut reader).is_ok());
        }
    }

    #[test]
    fn next_kmer() {
        let mut block: &[u8] = &[
            1,
            0,
            30,
            0,
            0,
            0,
            0, // version 1.0 encoding 0b00011011
            0b10100101,
            1,
            0,
            0,
            0,
            5,
            2,
            0b10011100,
            0b0000000011,
            1,
            2,
            3,
            4,
            5, // minimizer -> CCTT, 1 block, k -> 5, minimizer index -> 2
        ];

        let mut reader = KffReader::new(&mut block).unwrap();
        reader.variables().insert("k".to_string(), 5);
        reader.variables().insert("m".to_string(), 4);
        reader.variables().insert("max".to_string(), 5);
        reader.variables().insert("data_size".to_string(), 1);

        let mut minimizer = super::Reader::new(&mut reader).unwrap();

        assert_eq!(
            minimizer.next_kmer().unwrap(),
            Kmer::new(
                bitbox![Lsb0, u8; 0, 0, 1, 1, 1, 0, 1, 0, 0, 1],
                vec![1].into_boxed_slice()
            )
        );

        assert_eq!(
            minimizer.next_kmer().unwrap(),
            Kmer::new(
                bitbox![Lsb0, u8; 1, 1, 1, 0, 1, 0, 0, 1, 0, 1],
                vec![2].into_boxed_slice()
            )
        );
        assert_eq!(
            minimizer.next_kmer().unwrap(),
            Kmer::new(
                bitbox![Lsb0, u8; 1, 0, 1, 0, 0, 1, 0, 1, 1, 0],
                vec![3].into_boxed_slice()
            )
        );
        assert_eq!(
            minimizer.next_kmer().unwrap(),
            Kmer::new(
                bitbox![Lsb0, u8; 1, 0, 0, 1, 0, 1, 1, 0, 0, 1],
                vec![4].into_boxed_slice()
            )
        );
        assert_eq!(
            minimizer.next_kmer().unwrap(),
            Kmer::new(
                bitbox![Lsb0, u8; 0, 1, 0, 1, 1, 0, 0, 1, 1, 1],
                vec![5].into_boxed_slice()
            )
        );

        assert_eq!(
            minimizer.next_kmer().unwrap(),
            Kmer::new(bitbox![Lsb0, u8; 0; 0], vec![].into_boxed_slice())
        );

        assert_eq!(
            minimizer.next_kmer().unwrap(),
            Kmer::new(bitbox![Lsb0, u8; 0; 0], vec![].into_boxed_slice())
        );
    }

    #[test]
    fn bin() {
        let mut block: &[u8] = &[
            1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
            0b10100101, 1, 0, 0, 0, 5, 2, 0b10011100, 0b00000011, 1, 2, 3, 4,
            5, // minimizer -> CCTT, 1 block, k -> 5, minimizer index -> 2
        ];

        let mut reader = KffReader::new(&mut block).unwrap();
        reader.variables().insert("k".to_string(), 5);
        reader.variables().insert("m".to_string(), 4);
        reader.variables().insert("max".to_string(), 5);
        reader.variables().insert("data_size".to_string(), 1);

        let minimizer = super::Reader::new(&mut reader).unwrap();

        let mut it = minimizer.into_iter();

        assert_eq!(
            it.next().unwrap().unwrap(),
            Kmer::new(
                bitbox![Lsb0, u8; 0, 0, 1, 1, 1, 0, 1, 0, 0, 1],
                vec![1].into_boxed_slice()
            )
        );
        assert_eq!(
            it.next().unwrap().unwrap(),
            Kmer::new(
                bitbox![Lsb0, u8; 1, 1, 1, 0, 1, 0, 0, 1, 0, 1],
                vec![2].into_boxed_slice()
            )
        );
        assert_eq!(
            it.next().unwrap().unwrap(),
            Kmer::new(
                bitbox![Lsb0, u8; 1, 0, 1, 0, 0, 1, 0, 1, 1, 0],
                vec![3].into_boxed_slice()
            )
        );
        assert_eq!(
            it.next().unwrap().unwrap(),
            Kmer::new(
                bitbox![Lsb0, u8; 1, 0, 0, 1, 0, 1, 1, 0, 0, 1],
                vec![4].into_boxed_slice()
            )
        );
        assert_eq!(
            it.next().unwrap().unwrap(),
            Kmer::new(
                bitbox![Lsb0, u8; 0, 1, 0, 1, 1, 0, 0, 1, 1, 1],
                vec![5].into_boxed_slice()
            )
        );

        assert!(it.next().is_none());

        assert!(it.next().is_none());
    }

    #[test]
    fn seq() {
        let mut block: &[u8] = &[
            1,
            0,
            30,
            0,
            0,
            0,
            0, // version 1.0 encoding 0b00011011
            0b10100101,
            1,
            0,
            0,
            0,
            5,
            2,
            0b10011100,
            0b0000000011,
            1,
            2,
            3,
            4,
            5, // minimizer -> CCTT, 1 block, k -> 5, minimizer index -> 2
        ];

        let mut reader = KffReader::new(&mut block).unwrap();
        reader.variables().insert("k".to_string(), 5);
        reader.variables().insert("m".to_string(), 4);
        reader.variables().insert("max".to_string(), 5);
        reader.variables().insert("data_size".to_string(), 1);

        let minimizer = super::Reader::new(&mut reader).unwrap();

        let mut it = Box::new(minimizer.into_iter());

        let mut value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq(utils::rev_encoding(0b00011011)),
            vec![b'A', b'G', b'C', b'C', b'T'].into_boxed_slice()
        );

        assert_eq!(value.data(), &[1u8],);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq(utils::rev_encoding(0b00011011)),
            vec![b'G', b'C', b'C', b'T', b'T'].into_boxed_slice(),
        );
        assert_eq!(value.data(), &[2],);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq(utils::rev_encoding(0b00011011)),
            vec![b'C', b'C', b'T', b'T', b'C'].into_boxed_slice()
        );
        assert_eq!(value.data(), &[3],);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq(utils::rev_encoding(0b00011011)),
            vec![b'C', b'T', b'T', b'C', b'T'].into_boxed_slice(),
        );
        assert_eq!(value.data(), &[4],);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq(utils::rev_encoding(0b00011011)),
            vec![b'T', b'T', b'C', b'T', b'G'].into_boxed_slice(),
        );
        assert_eq!(value.data(), &[5],);

        assert!(it.next().is_none());

        assert!(it.next().is_none());
    }

    #[test]
    fn write() {
        let mut variables = Variables::new();
        variables.insert("k".to_string(), 5);
        variables.insert("m".to_string(), 2);
        variables.insert("max".to_string(), 255);
        variables.insert("data_size".to_string(), 1);

        // minimizer sequence AC
        let minimizer = b"AC";

        let mut buffer = std::io::Cursor::new(vec![0u8; 0]);
        {
            let mut writer = Writer::new(&variables, minimizer, 0b00011011, &mut buffer).unwrap();

            writer.write_seq_block(4, b"GAGTT", &[10, 8, 9]).unwrap();
            writer.write_seq_block(0, b"GAT", &[1]).unwrap();
        }

        assert_eq!(
            buffer.into_inner(),
            [4, 2, 0, 0, 0, 3, 4, 179, 2, 10, 8, 9, 1, 0, 35, 1]
        );
    }

    #[test]
    fn write_n1() {
        let mut variables = Variables::new();
        variables.insert("k".to_string(), 5);
        variables.insert("m".to_string(), 2);
        variables.insert("max".to_string(), 1);
        variables.insert("data_size".to_string(), 1);

        // minimizer sequence AC
        let minimizer = b"AC";

        let mut buffer = std::io::Cursor::new(vec![0u8; 0]);
        {
            let mut writer = Writer::new(&variables, minimizer, 0b00011011, &mut buffer).unwrap();

            writer.write_seq_block(1, b"GAG", &[10]).unwrap();
            writer.write_seq_block(0, b"GAT", &[1]).unwrap();
        }

        assert_eq!(buffer.into_inner(), [4, 2, 0, 0, 0, 1, 51, 10, 0, 35, 1]);
    }

    #[test]
    fn write_read() {
        let mut variables = Variables::new();
        variables.insert("k".to_string(), 5);
        variables.insert("m".to_string(), 2);
        variables.insert("max".to_string(), 255);
        variables.insert("data_size".to_string(), 1);

        // minimizer sequence AC
        let minimizer = b"AC";

        let mut buffer = std::io::Cursor::new(vec![1u8, 0, 30, 0]);
        buffer.set_position(7);
        {
            let mut writer = Writer::new(&variables, minimizer, 0b00011011, &mut buffer).unwrap();

            writer.write_seq_block(4, b"GAGTT", &[10, 8, 9]).unwrap();
            writer.write_seq_block(0, b"GAT", &[1]).unwrap();
        }

        let inp = buffer.into_inner();
        assert_eq!(
            inp,
            [1, 0, 30, 0, 0, 0, 0, 4, 2, 0, 0, 0, 3, 4, 179, 2, 10, 8, 9, 1, 0, 35, 1]
        );

        let mut input = std::io::Cursor::new(inp);

        let mut reader = KffReader::new(&mut input).unwrap();
        reader.variables().insert("k".to_string(), 5);
        reader.variables().insert("m".to_string(), 2);
        reader.variables().insert("max".to_string(), 255);
        reader.variables().insert("data_size".to_string(), 1);

        let reader = super::Reader::new(&mut reader).unwrap();
        let mut it = reader.into_iter();

        let mut value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq(utils::rev_encoding(0b00011011)),
            vec![b'G', b'A', b'G', b'T', b'A'].into_boxed_slice(),
        );
        assert_eq!(value.data(), &[10]);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq(utils::rev_encoding(0b00011011)),
            vec![b'A', b'G', b'T', b'A', b'C'].into_boxed_slice()
        );
        assert_eq!(value.data(), &[8]);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq(utils::rev_encoding(0b00011011)),
            vec![b'G', b'T', b'A', b'C', b'T'].into_boxed_slice(),
        );
        assert_eq!(value.data(), &[9]);

        value = it.next().unwrap().unwrap();
        assert_eq!(
            value.seq(utils::rev_encoding(0b00011011)),
            vec![b'A', b'C', b'G', b'A', b'T'].into_boxed_slice(),
        );
        assert_eq!(value.data(), &[1]);

        assert!(it.next().is_none());
    }
}
