/* crate use */
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

/* local use */
use crate::error;
use crate::metadata::Metadata;
use crate::*;

pub struct Header {
    major: u8,
    minor: u8,
    encoding: u8,
    rev_encoding: u8,
    comment: Box<[u8]>,
}

impl Header {
    pub fn default() -> Self {
        Header {
            major: 1,
            minor: 0,
            encoding: 0b00011011,
            rev_encoding: 0b00011011,
            comment: b"".to_vec().into_boxed_slice(),
        }
    }

    pub fn new(major: u8, minor: u8, encoding: u8, comment: Box<[u8]>) -> Self {
        Header {
            major,
            minor,
            encoding,
            rev_encoding: utils::rev_encoding(encoding),
            comment,
        }
    }

    pub fn major(&self) -> u8 {
        self.major
    }

    pub fn minor(&self) -> u8 {
        self.minor
    }

    pub fn encoding(&self) -> u8 {
        self.encoding
    }

    pub fn rev_encoding(&self) -> u8 {
        self.rev_encoding
    }

    pub fn comment(&self) -> &[u8] {
        &self.comment
    }
}

impl Metadata for Header {
    fn deserialize<R>(&mut self, input: &mut R) -> crate::Result<usize>
    where
        R: std::io::Read,
    {
        self.major = input.read_u8()?;
        self.minor = input.read_u8()?;
        self.encoding = valid_encoding(switch_56_n_78(input.read_u8()?))?;
        self.rev_encoding = utils::rev_encoding(self.encoding);

        let mut comment = vec![0; input.read_u32::<LittleEndian>()? as usize].into_boxed_slice();
        input.read_exact(&mut comment)?;
        self.comment = comment;

        Ok(4 + self.comment.len())
    }

    fn serialize<W: std::io::Write>(&self, output: &mut W) -> crate::Result<usize> {
        output.write_all(&[self.major, self.minor])?;
        output.write_u8(switch_56_n_78(self.encoding))?;
        output.write_u32::<LittleEndian>(self.comment.len() as u32)?;
        output.write_all(&self.comment)?;

        Ok(7 + self.comment.len())
    }
}

fn switch_56_n_78(input: u8) -> u8 {
    // Internal encoding order is ACTG kff order is ACGT
    (input & 0b11110000) ^ ((input & 0b00000011) << 2) ^ ((input & 0b00001100) >> 2)
}

fn valid_encoding(encoding: u8) -> crate::Result<u8> {
    let a = encoding >> 6;
    let c = (encoding >> 4) & 0b11;
    let t = (encoding >> 2) & 0b11;
    let g = encoding & 0b11;

    if a != c && a != t && a != g && c != t && t != g {
        Ok(encoding)
    } else {
        Err(Box::new(error::Header::BadEncoding))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn getter() {
        let header = Header {
            major: 12,
            minor: 2,
            encoding: 0b00011011,
            rev_encoding: 0b00011011,
            comment: b"test 1 2 1 2".to_vec().into_boxed_slice(),
        };

        assert_eq!(header.major(), 12);
        assert_eq!(header.minor(), 2);
        assert_eq!(header.encoding(), 0b00011011);
        assert_eq!(header.comment(), b"test 1 2 1 2");
    }

    #[test]
    fn write() {
        let mut output = Vec::new();

        let header = Header::new(
            78,
            46,
            0b00011011,
            b"test 1 2 1 2".to_vec().into_boxed_slice(),
        );

        header.serialize(&mut output).unwrap();

        assert_eq!(
            output,
            [78, 46, 30, 12, 0, 0, 0, 116, 101, 115, 116, 32, 49, 32, 50, 32, 49, 32, 50]
        );
    }

    #[test]
    fn read() {
        let mut input: &[u8] = &[
            78, 46, 30, 12, 0, 0, 0, 116, 101, 115, 116, 32, 49, 32, 50, 32, 49, 32, 50,
        ];

        let mut header = Header::default();

        header.deserialize(&mut input).unwrap();

        let mut output = Vec::new();

        header.serialize(&mut output).unwrap();

        assert_eq!(
            [78, 46, 30, 12, 0, 0, 0, 116, 101, 115, 116, 32, 49, 32, 50, 32, 49, 32, 50,].to_vec(),
            output
        );
    }

    #[test]
    fn kff2internal_order() {
        assert_eq!(switch_56_n_78(0b01011100), 0b01010011);
        assert_eq!(switch_56_n_78(0b00000011), 0b00001100);
    }

    #[test]
    fn validate_encoding() {
        assert!(valid_encoding(0b00011011).is_ok());
        assert!(valid_encoding(0b00011111).is_err());
    }
}
