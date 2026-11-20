//! Nucleotide encoding and 2-bit sequences
//!
//! A KFF file stores the nucleotides with 2 bits each. The header of the file gives the 2-bit
//! value of each nucleotide, in the order A, C, G, T (the first two bits of the encoding byte are
//! the value of A).
//!
//! A sequence of `n` nucleotides is stored big endian in `ceil(n / 4)` bytes, right aligned: the
//! last nucleotide is in the 2 lowest bits of the last byte and the padding bits (if any) are the
//! highest bits of the first byte. This crate writes padding bits set to 0.

/* project use */
use crate::error;

/// Encoding of the nucleotides: the 2-bit value of A, C, G and T.
///
/// The default encoding is the one of the C++ API: A=0, C=1, G=3, T=2.
#[derive(Clone, Copy, std::fmt::Debug, PartialEq, Eq, Hash)]
pub struct Encoding {
    /// 2-bit values of A, C, G, T
    values: [u8; 4],
    /// Nucleotide of each 2-bit value
    letters: [u8; 4],
}

impl Encoding {
    /// Build an encoding from the 2-bit values of A, C, G and T.
    ///
    /// As in the C++ API, only the 2 lowest bits of each value are used. The 4 values must be
    /// different.
    pub fn new(a: u8, c: u8, g: u8, t: u8) -> error::Result<Self> {
        let values = [a & 0b11, c & 0b11, g & 0b11, t & 0b11];
        let byte = (values[0] << 6) | (values[1] << 4) | (values[2] << 2) | values[3];

        let mut letters = [0u8; 4];
        let mut seen = [false; 4];
        for (value, letter) in values.iter().zip(*b"ACGT") {
            if seen[*value as usize] {
                return Err(error::Error::InvalidEncoding(byte));
            }
            seen[*value as usize] = true;
            letters[*value as usize] = letter;
        }

        Ok(Self { values, letters })
    }

    /// Build an encoding from the byte of a KFF header
    pub fn from_byte(byte: u8) -> error::Result<Self> {
        Self::new(byte >> 6, byte >> 4, byte >> 2, byte)
    }

    /// Byte of the encoding in a KFF header
    pub fn as_byte(&self) -> u8 {
        (self.values[0] << 6) | (self.values[1] << 4) | (self.values[2] << 2) | self.values[3]
    }

    /// The 2-bit values of A, C, G and T (in this order), as `get_encoding()` in the C++ API
    pub fn values(&self) -> [u8; 4] {
        self.values
    }

    /// 2-bit value of a nucleotide (A, C, G or T, upper or lower case)
    pub fn value(&self, nucleotide: u8) -> Option<u8> {
        match nucleotide {
            b'A' | b'a' => Some(self.values[0]),
            b'C' | b'c' => Some(self.values[1]),
            b'G' | b'g' => Some(self.values[2]),
            b'T' | b't' => Some(self.values[3]),
            _ => None,
        }
    }

    /// Nucleotide (upper case) of a 2-bit value. Only the 2 lowest bits are used.
    pub fn nucleotide(&self, value: u8) -> u8 {
        self.letters[(value & 0b11) as usize]
    }

    /// Encode a sequence of nucleotides with 2 bits per nucleotide, right aligned in
    /// `ceil(len / 4)` bytes with null padding bits.
    pub fn encode(&self, sequence: &[u8]) -> error::Result<Vec<u8>> {
        let mut bytes = vec![0; crate::nucleotides::bytes_for_nucleotides(sequence.len() as u64)];
        let padding = (4 - sequence.len() % 4) % 4;

        for (i, nucleotide) in sequence.iter().enumerate() {
            let value = self
                .value(*nucleotide)
                .ok_or(error::Error::InvalidNucleotide(*nucleotide as char))?;
            let position = padding + i;
            bytes[position / 4] |= value << (2 * (3 - position % 4));
        }

        Ok(bytes)
    }

    /// Decode the `nb_nucleotides` nucleotides of a right aligned 2-bit sequence. The padding bits
    /// are ignored.
    ///
    /// # Panics
    /// If `bytes` is shorter than `ceil(nb_nucleotides / 4)`.
    pub fn decode(&self, bytes: &[u8], nb_nucleotides: usize) -> Vec<u8> {
        let padding = (4 - nb_nucleotides % 4) % 4;

        (padding..padding + nb_nucleotides)
            .map(|position| self.nucleotide(bytes[position / 4] >> (2 * (3 - position % 4))))
            .collect()
    }
}

impl std::default::Default for Encoding {
    fn default() -> Self {
        Self {
            values: [0, 1, 3, 2],
            letters: *b"ACTG",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_the_cpp_one() {
        assert_eq!(Encoding::default(), Encoding::new(0, 1, 3, 2).unwrap());
        assert_eq!(Encoding::default().as_byte(), 0b00011110);
    }

    #[test]
    fn byte_order_is_acgt() {
        // README of the format: 0b00101101 means A=0, C=2, G=3, T=1
        let encoding = Encoding::from_byte(0b00101101).unwrap();
        assert_eq!(encoding.values(), [0, 2, 3, 1]);
        assert_eq!(encoding.value(b'G'), Some(3));
        assert_eq!(encoding.value(b't'), Some(1));
        assert_eq!(encoding.nucleotide(3), b'G');
        assert_eq!(encoding.nucleotide(1), b'T');
        assert_eq!(encoding.as_byte(), 0b00101101);
    }

    #[test]
    fn all_permutations() {
        let mut nb_valid = 0;
        for byte in 0..=255u8 {
            if let Ok(encoding) = Encoding::from_byte(byte) {
                nb_valid += 1;
                assert_eq!(encoding.as_byte(), byte);
                for nucleotide in *b"ACGT" {
                    assert_eq!(
                        encoding.nucleotide(encoding.value(nucleotide).unwrap()),
                        nucleotide
                    );
                }
            }
        }
        assert_eq!(nb_valid, 24);
    }

    #[test]
    fn invalid_encodings() {
        assert!(Encoding::new(0, 0, 1, 2).is_err());
        assert!(Encoding::new(0, 1, 0, 2).is_err());
        assert!(Encoding::new(0, 1, 2, 0).is_err());
        assert!(Encoding::new(0, 1, 1, 2).is_err());
        assert!(Encoding::new(0, 1, 2, 1).is_err());
        assert!(Encoding::new(0, 1, 2, 2).is_err());
        assert!(Encoding::new(3, 2, 1, 0).is_ok());
        // Only the 2 lowest bits are used: 0xff 0xfe 0xfc 0xfd = 3 2 0 1
        assert_eq!(
            Encoding::new(0xff, 0xfe, 0xfc, 0xfd).unwrap().values(),
            [3, 2, 0, 1]
        );
    }

    #[test]
    fn encode_readme_example() {
        // README of the format: A=0, C=2, G=3, T=1
        let encoding = Encoding::new(0, 2, 3, 1).unwrap();
        assert_eq!(
            encoding.encode(b"ACTAAACTGATT").unwrap(),
            vec![0x24, 0x09, 0xc5]
        );
        assert_eq!(
            encoding.encode(b"AAACTGATCG").unwrap(),
            vec![0x00, 0x27, 0x1b]
        );
        assert_eq!(
            encoding.encode(b"CTAAACTGATT").unwrap(),
            vec![0x24, 0x09, 0xc5]
        );
        assert_eq!(encoding.encode(b"").unwrap(), Vec::<u8>::new());
        assert!(encoding.encode(b"ACNT").is_err());
    }

    #[test]
    fn decode_ignores_padding() {
        let encoding = Encoding::new(0, 2, 3, 1).unwrap();
        assert_eq!(encoding.decode(&[0xf0, 0x27, 0x1b], 10), b"AAACTGATCG");
        assert_eq!(encoding.decode(&[0x24, 0x09, 0xc5], 12), b"ACTAAACTGATT");
        assert_eq!(encoding.decode(&[0xe4], 3), b"CTA");
    }

    #[test]
    fn round_trip_lower_case() {
        let encoding = Encoding::default();
        let bytes = encoding.encode(b"acgtACGTa").unwrap();
        assert_eq!(encoding.decode(&bytes, 9), b"ACGTACGTA");
    }
}
