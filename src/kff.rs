//! Declaration of Kff::Reader and Kff::Writer

/* crate use */
use anyhow::Result;
use byteorder::*;

/* local use */
use crate::data;
use crate::data::Writer as DataWriter;
use crate::error;
use crate::minimizer::Reader as MinimizerReader;
use crate::minimizer::Writer as MinimizerWriter;
use crate::raw::Reader as RawReader;
use crate::raw::Writer as RawWriter;
use crate::seq2bits::Seq2Bits;
use crate::utils;
use crate::utils::switch_56_n_78;
use crate::variables::Reader as VariablesReader;
use crate::variables::Variables;
use crate::variables::Writer as VariablesWriter;

/// A Kff Reader
pub struct Reader<R>
where
    R: std::io::Read,
{
    // header
    major: u8,
    minor: u8,
    encoding: u8,
    rev_encoding: u8,
    comment: Box<[u8]>,

    //
    input: R,

    // global variables
    variables: Variables,
}

impl<R> Reader<R>
where
    R: std::io::Read,
{
    /// Create a new Kff reader from a reading stream
    pub fn new(mut input: R) -> Result<Self> {
        let major = input.read_u8()?;
        let minor = input.read_u8()?;
        let encoding = utils::valid_encoding(utils::switch_56_n_78(input.read_u8()?))?;
        let rev_encoding = utils::rev_encoding(encoding);

        let comment_len = input.read_u32::<utils::Order>()? as usize;

        let mut comment = vec![0; comment_len].into_boxed_slice();

        input.read_exact(&mut comment)?;
        let comment = comment;

        if major < 1 {
            return Err(error::Error::Kff(error::Kff::NotSupportVersionNumber).into());
        }

        let variables = Variables::default();

        Ok(Self {
            major,
            minor,
            encoding,
            rev_encoding,
            comment,
            input,
            variables,
        })
    }

    // Getter
    /// Get the major verion number
    pub fn major(&self) -> u8 {
        self.major
    }

    /// Get the minor verion number
    pub fn minor(&self) -> u8 {
        self.minor
    }

    /// Get encoding used
    pub fn encoding(&self) -> u8 {
        self.encoding
    }

    /// Get reverse encoding used
    pub fn rev_encoding(&self) -> u8 {
        self.rev_encoding
    }

    /// Get comment
    pub fn comment(&self) -> &[u8] {
        &self.comment
    }

    /// Get a mutable reference to input stream
    pub fn input(&mut self) -> &'_ mut R {
        &mut self.input
    }

    /// Get a mutable reference of global variable
    pub fn variables(&mut self) -> &mut Variables {
        &mut self.variables
    }

    /// Get the next kmer section we have to parse
    pub fn next_section(&mut self) -> Result<Box<dyn data::Reader<R> + '_>> {
        match self.input.read_u8()? {
            b'r' => Ok(Box::new(RawReader::new(self)?)),
            b'm' => Ok(Box::new(MinimizerReader::new(self)?)),
            b'v' => {
                self.variables.deserialize(&mut self.input)?;
                self.next_section()
            }
            _ => Err(error::Error::Kff(error::Kff::UnknowSectionType).into()),
        }
    }
}

/// A Kff Writer
pub struct Writer<W>
where
    W: std::io::Write,
{
    output: W,
    encoding: u8,
    variables_buffer: Variables,
    variables: Variables,
}

impl<W> Writer<W>
where
    W: std::io::Write + std::io::Seek,
{
    /// Create a Kff Writer from a seekable output, encoding and comment
    pub fn new(mut output: W, encoding: u8, comment: &[u8]) -> Result<Self> {
        // write header
        output.write_all(&[1u8, 0, switch_56_n_78(encoding)])?;
        output.write_u32::<BigEndian>(comment.len() as u32)?;
        output.write_all(comment)?;

        Ok(Self {
            output,
            encoding,
            variables_buffer: Variables::default(),
            variables: Variables::default(),
        })
    }

    /// Get encoding used
    pub fn encoding(&self) -> u8 {
        self.encoding
    }

    /// Get a mutable reference of global variable
    pub fn variables(&mut self) -> &mut Variables {
        &mut self.variables_buffer
    }

    /// Write variable section
    pub fn write_variables(&mut self) -> Result<()> {
        self.output.write_u8(b'v')?;

        self.variables_buffer.serialize(&mut self.output)?;

        self.variables.extend(self.variables_buffer.drain());

        Ok(())
    }

    /// Write a raw section, with sequence encode in 2 bits
    pub fn write_raw_section<A: AsRef<[u8]>>(
        &mut self,
        seqs: &[Seq2Bits],
        datas: &[A],
    ) -> Result<usize> {
        self.output.write_u8(b'r')?;

        let mut raw = RawWriter::new(&self.variables, self.encoding, &mut self.output)?;

        let mut nb_bytes = 0;

        for (seq, data) in seqs.iter().zip(datas) {
            nb_bytes += raw.write_block(&seq, &data.as_ref())?
        }

        raw.close()?;

        Ok(nb_bytes)
    }

    /// Write a raw section with sequence encode in ASCII
    pub fn write_raw_seq_section<A: AsRef<[u8]>, B: AsRef<[u8]>>(
        &mut self,
        seqs: &[A],
        datas: &[B],
    ) -> Result<usize> {
        // Todo it's ugly
        let tmp: Vec<Seq2Bits> = seqs
            .iter()
            .map(|x| utils::seq2bits(x.as_ref(), self.encoding))
            .collect();

        let tmp2: Vec<&[u8]> = datas.iter().map(|x| x.as_ref()).collect();

        self.write_raw_section(&tmp[..], &tmp2[..])
    }

    /// Write a minimizer section, with sequence encode in 2 bits
    pub fn write_minimizer_section<A: AsRef<[u8]>>(
        &mut self,
        minimizer: &[u8],
        mini_index: &[u64],
        seqs: &[Seq2Bits],
        datas: &[A],
    ) -> Result<usize> {
        self.output.write_u8(b'm')?;

        let mut minimizer =
            MinimizerWriter::new(&self.variables, minimizer, self.encoding, &mut self.output)?;

        let mut nb_bytes = 0;

        for (index, (seq, data)) in mini_index.iter().zip(seqs.iter().zip(datas)) {
            nb_bytes += minimizer.write_block(*index, &seq, &data.as_ref())?
        }

        minimizer.close()?;

        Ok(nb_bytes)
    }

    /// Write a raw section with sequence not encode in ASCII
    pub fn write_minimizer_seq_section<A: AsRef<[u8]>, B: AsRef<[u8]>>(
        &mut self,
        minimizer: &[u8],
        mini_index: &[u64],
        seqs: &[A],
        datas: &[B],
    ) -> Result<usize> {
        // Todo it's ugly

        let tmp: Vec<Seq2Bits> = seqs
            .iter()
            .map(|x| utils::seq2bits(x.as_ref(), self.encoding))
            .collect();

        let tmp2: Vec<&[u8]> = datas.iter().map(|x| x.as_ref()).collect();

        self.write_minimizer_section(minimizer, mini_index, &tmp[..], &tmp2[..])
    }
}

#[cfg(test)]
mod tests {
    use crate::seq2bits::Bits2Nuc;

    #[test]
    fn read() {
        let mut input: &[u8] = &[
            1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
            b'v', 3, 0, 0, 0, 0, 0, 0, 0, // 3 variables
            109, 97, 120, 0, 5, 0, 0, 0, 0, 0, 0, 0, // max -> 5
            100, 97, 116, 97, 95, 115, 105, 122, 101, 0, 1, 0, 0, 0, 0, 0, 0,
            0, // data size -> 1
            107, 0, 5, 0, 0, 0, 0, 0, 0, 0, // variable k -> 5
            b'r', 1, 0, 0, 0, 5, 0b00000011, 0b11011110, 0b00111101, 1, 2, 3, 4, 5,
            // Raw block, 5 kmer in block,
            b'v', 1, 0, 0, 0, 0, 0, 0, 0, // 1 variable
            109, 0, 4, 0, 0, 0, 0, 0, 0, 0, // m -> 4
            b'm', 0b10100101, 1, 0, 0, 0, 5, 2, 0b00000010, 0b01110011, 1, 2, 3, 4,
            5,
            // Minimizer block, minimizer -> TTCC, 1 block, 5 kmer, minimizer index -> 2
        ];

        let mut reader = super::Reader::new(&mut input).unwrap();

        let rev_encoding = reader.rev_encoding();

        {
            let section = reader.next_section().unwrap();
            let mut it = section.into_iter();

            let mut value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'G', b'G', b'C', b'G', b'T'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[1]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'G', b'C', b'G', b'T', b'A'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[2]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'C', b'G', b'T', b'A', b'G'].into_boxed_slice()
            );
            assert_eq!(value.data(), &[3]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'G', b'T', b'A', b'G', b'G'].into_boxed_slice()
            );
            assert_eq!(value.data(), &[4]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'T', b'A', b'G', b'G', b'C'].into_boxed_slice()
            );
            assert_eq!(value.data(), &[5]);
            assert!(it.next().is_none());
        }

        {
            let section = reader.next_section().unwrap();
            let mut it = section.into_iter();

            let mut value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'T', b'C', b'T', b'T', b'C'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[1]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'C', b'T', b'T', b'C', b'C'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[2]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'T', b'T', b'C', b'C', b'G'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[3]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'T', b'C', b'C', b'G', b'A'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[4]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'C', b'C', b'G', b'A', b'G'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[5]);

            assert!(it.next().is_none());
        }
        assert!(reader.next_section().is_err());
    }

    #[test]
    fn realistic_bin_read() {
        let mut input: &[u8] = &[
            1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011
            b'v', 3, 0, 0, 0, 0, 0, 0, 0, // 3 variables
            109, 97, 120, 0, 5, 0, 0, 0, 0, 0, 0, 0, // max -> 5
            100, 97, 116, 97, 95, 115, 105, 122, 101, 0, 1, 0, 0, 0, 0, 0, 0,
            0, // data size -> 1
            107, 0, 5, 0, 0, 0, 0, 0, 0, 0, // variable k -> 5
            b'r', 1, 0, 0, 0, 5, 0b00000011, 0b11011110, 0b00111101, 1, 2, 3, 4, 5,
            // Raw block, 5 kmer in block,
            b'v', 1, 0, 0, 0, 0, 0, 0, 0, // 1 variable
            109, 0, 4, 0, 0, 0, 0, 0, 0, 0, // m -> 4
            b'm', 0b10100101, 1, 0, 0, 0, 5, 2, 0b00000010, 0b01110011, 1, 2, 3, 4,
            5,
            // Minimizer block, minimizer -> TTCC, 1 block, 5 kmer, minimizer index -> 2
        ];

        let mut reader = super::Reader::new(&mut input).unwrap();

        let mut nb_kmer = 0;

        while let Ok(section) = reader.next_section() {
            let mut it = section.into_iter();
            while let Some(Ok(_)) = it.next() {
                nb_kmer += 1;
            }
        }

        assert_eq!(nb_kmer, 10);
    }

    #[test]
    fn write_read() {
        let mut output = vec![0u8; 0];
        let buffer = std::io::Cursor::new(&mut output);

        let mut writer = super::Writer::new(buffer, 0b00011011, b"").unwrap();

        writer.variables().insert("k".to_string(), 5);
        writer.variables().insert("max".to_string(), 5);
        writer.variables().insert("data_size".to_string(), 1);

        writer.write_variables().unwrap();

        writer
            .write_raw_seq_section(&[b"GCGGGGATC"], &[vec![1u8, 2, 3, 4, 5]])
            .unwrap();

        writer.variables().insert("m".to_string(), 4);

        writer.write_variables().unwrap();

        writer
            .write_minimizer_seq_section(b"CCTT", &[2], &[b"AGCTG"], &[&[6, 7, 8, 9, 10]])
            .unwrap();

        let mut input = std::io::Cursor::new(output);
        let mut reader = super::Reader::new(&mut input).unwrap();

        let rev_encoding = reader.rev_encoding();

        {
            let section = reader.next_section().unwrap();
            let mut it = section.into_iter();

            let mut value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'G', b'C', b'G', b'G', b'G'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[1]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'C', b'G', b'G', b'G', b'G'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[2]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'G', b'G', b'G', b'G', b'A'].into_boxed_slice()
            );
            assert_eq!(value.data(), &[3]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'G', b'G', b'G', b'A', b'T'].into_boxed_slice()
            );
            assert_eq!(value.data(), &[4]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'G', b'G', b'A', b'T', b'C'].into_boxed_slice()
            );
            assert_eq!(value.data(), &[5]);
            assert!(it.next().is_none());
        }

        {
            let section = reader.next_section().unwrap();
            let mut it = section.into_iter();

            let mut value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'A', b'G', b'C', b'C', b'T'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[6]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'G', b'C', b'C', b'T', b'T'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[7]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'C', b'C', b'T', b'T', b'C'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[8]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'C', b'T', b'T', b'C', b'T'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[9]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq().into_nuc(rev_encoding),
                vec![b'T', b'T', b'C', b'T', b'G'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[10]);

            assert!(it.next().is_none());
        }

        assert!(reader.next_section().is_err());
    }
}
