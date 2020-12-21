/* std use */

/* crate use */
use byteorder::*;

/* local use */
use crate::data;
use crate::data::Writer as DataWriter;
use crate::error;
use crate::header::Header;
use crate::metadata::Metadata;
use crate::minimizer;
use crate::raw;
use crate::utils;
use crate::utils::BitBox;
use crate::variables::Variables;

pub struct Reader<R>
where
    R: std::io::Read,
{
    input: R,
    header: Header,
    variables: Variables,
}

impl<R> Reader<R>
where
    R: std::io::Read,
{
    pub fn new(mut input: R) -> crate::Result<Self> {
        let mut header = Header::default();

        header.deserialize(&mut input)?;

        if header.major() < 1 {
            return Err(Box::new(error::Kff::NotSupportVersionNumber));
        }

        let variables = Variables::default();

        Ok(Self {
            input,
            header,
            variables,
        })
    }

    pub fn input(&mut self) -> &mut R {
        &mut self.input
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn variables(&mut self) -> &mut Variables {
        &mut self.variables
    }

    pub fn encoding(&self) -> u8 {
        self.header.encoding()
    }

    pub fn rev_encoding(&self) -> u8 {
        utils::rev_encoding(self.header.encoding())
    }

    pub fn next_section(&mut self) -> crate::Result<Box<dyn data::Reader<R> + '_>> {
        match self.input.read_u8()? {
            b'r' => Ok(Box::new(raw::Reader::new(self)?)),
            b'm' => Ok(Box::new(minimizer::Reader::new(self)?)),
            b'v' => {
                self.variables.deserialize(&mut self.input)?;
                self.next_section()
            }
            _ => Err(Box::new(error::Kff::UnknowSectionType)),
        }
    }
}

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
    pub fn new(mut output: W, encoding: u8, comment: &[u8]) -> crate::Result<Self> {
        let header = Header::new(1, 0, encoding, Box::from(comment));
        header.serialize(&mut output)?;

        Ok(Self {
            output,
            encoding,
            variables_buffer: Variables::default(),
            variables: Variables::default(),
        })
    }

    pub fn encoding(&self) -> u8 {
        self.encoding
    }

    pub fn variables(&mut self) -> &mut Variables {
        &mut self.variables_buffer
    }

    pub fn write_variables(&mut self) -> crate::Result<()> {
        self.output.write_u8(b'v')?;

        self.variables_buffer.serialize(&mut self.output)?;

        self.variables.extend(self.variables_buffer.drain());

        Ok(())
    }

    pub fn write_raw_block(&mut self, seqs: &[BitBox], datas: &[&[u8]]) -> crate::Result<usize> {
        self.output.write_u8(b'r')?;

        let mut raw = raw::Writer::new(&self.variables, self.encoding, &mut self.output)?;

        let mut nb_bytes = 0;

        for (seq, data) in seqs.iter().zip(datas) {
            nb_bytes += raw.write_block(&seq, &data)?
        }

        raw.close()?;

        Ok(nb_bytes)
    }

    pub fn write_raw_seq_block(&mut self, seqs: &[&[u8]], datas: &[&[u8]]) -> crate::Result<usize> {
        let tmp: Vec<BitBox> = seqs
            .iter()
            .map(|x| utils::seq2bits(x, self.encoding))
            .collect();
        self.write_raw_block(&tmp[..], datas)
    }

    pub fn write_minimizer_block(
        &mut self,
        minimizer: &[u8],
        mini_index: &[u64],
        seqs: &[BitBox],
        datas: &[&[u8]],
    ) -> crate::Result<usize> {
        self.output.write_u8(b'm')?;

        let mut minimizer =
            minimizer::Writer::new(&self.variables, minimizer, self.encoding, &mut self.output)?;

        let mut nb_bytes = 0;

        for (index, (seq, data)) in mini_index.iter().zip(seqs.iter().zip(datas)) {
            nb_bytes += minimizer.write_block(*index, &seq, &data)?
        }

        minimizer.close()?;

        Ok(nb_bytes)
    }

    pub fn write_minimizer_seq_block(
        &mut self,
        minimizer: &[u8],
        mini_index: &[u64],
        seqs: &[&[u8]],
        datas: &[&[u8]],
    ) -> crate::Result<usize> {
        let tmp: Vec<BitBox> = seqs
            .iter()
            .map(|x| utils::seq2bits(x, self.encoding))
            .collect();
        self.write_minimizer_block(minimizer, mini_index, &tmp[..], datas)
    }
}

#[cfg(test)]
mod tests {
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
            // varibale m -> 4
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
                value.seq(rev_encoding),
                vec![b'G', b'G', b'C', b'G', b'T'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[1]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'G', b'C', b'G', b'T', b'A'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[2]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'C', b'G', b'T', b'A', b'G'].into_boxed_slice()
            );
            assert_eq!(value.data(), &[3]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'G', b'T', b'A', b'G', b'G'].into_boxed_slice()
            );
            assert_eq!(value.data(), &[4]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
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
                value.seq(rev_encoding),
                vec![b'T', b'C', b'T', b'T', b'C'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[1]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'C', b'T', b'T', b'C', b'C'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[2]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'T', b'T', b'C', b'C', b'G'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[3]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'T', b'C', b'C', b'G', b'A'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[4]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
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
            1, 0, 30, 0, 0, 0, 0, // version 1.0 encoding 0b00011011 no comment
            b'v', 3, 0, 0, 0, 0, 0, 0, 0, 109, 97, 120, 0, 5, 0, 0, 0, 0, 0, 0, 0, 100, 97, 116,
            97, 95, 115, 105, 122, 101, 0, 1, 0, 0, 0, 0, 0, 0, 0, 107, 0, 5, 0, 0, 0, 0, 0, 0, 0,
            // variable k -> 5, max -> 5, data_size -> 1
            b'r', 1, 0, 0, 0, 5, 0b11110111, 0b10001111, 0b0000001, 1, 2, 3, 4, 5,
            // Raw block, 5 kmer in block,
            b'v', 1, 0, 0, 0, 0, 0, 0, 0, 109, 0, 4, 0, 0, 0, 0, 0, 0, 0,
            // varibale m -> 4
            b'm', 0b10100101, 1, 0, 0, 0, 5, 2, 0b10011100, 0b00000011, 1, 2, 3, 4,
            5,
            // Minimizer block, minimizer -> CCTT, 1 block, k -> 5, minimizer index -> 2
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
            .write_raw_seq_block(&[b"GCGGGGATC"], &[&[1u8, 2, 3, 4, 5]])
            .unwrap();

        writer.variables().insert("m".to_string(), 4);

        writer.write_variables().unwrap();

        writer
            .write_minimizer_seq_block(b"CCTT", &[2], &[b"AGCTG"], &[&[6, 7, 8, 9, 10]])
            .unwrap();

        let mut input = std::io::Cursor::new(output);
        let mut reader = super::Reader::new(&mut input).unwrap();

        let rev_encoding = reader.rev_encoding();

        {
            let section = reader.next_section().unwrap();
            let mut it = section.into_iter();

            let mut value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'G', b'C', b'G', b'G', b'G'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[1]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'C', b'G', b'G', b'G', b'G'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[2]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'G', b'G', b'G', b'G', b'A'].into_boxed_slice()
            );
            assert_eq!(value.data(), &[3]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'G', b'G', b'G', b'A', b'T'].into_boxed_slice()
            );
            assert_eq!(value.data(), &[4]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
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
                value.seq(rev_encoding),
                vec![b'A', b'G', b'C', b'C', b'T'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[6]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'G', b'C', b'C', b'T', b'T'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[7]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'C', b'C', b'T', b'T', b'C'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[8]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'C', b'T', b'T', b'C', b'T'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[9]);

            value = it.next().unwrap().unwrap();
            assert_eq!(
                value.seq(rev_encoding),
                vec![b'T', b'T', b'C', b'T', b'G'].into_boxed_slice(),
            );
            assert_eq!(value.data(), &[10]);

            assert!(it.next().is_none());
        }

        assert!(reader.next_section().is_err());
    }
}
