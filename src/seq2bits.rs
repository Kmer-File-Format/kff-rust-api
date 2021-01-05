/* standard use */
use std::ops::Index;

/* local use */
use crate::utils::{bits2seq, seq2bits, BitBox, BitSlice, BitVec};

pub trait Nuc2Bits {
    fn from_nuc(input: &[u8], encoding: u8) -> Self;

    fn from_bitslice(input: &BitSlice) -> Self;

    fn from_bits(input: &[u8], nb_nuc: usize) -> Self;
}

pub trait Bits2Nuc {
    fn into_nuc(&self, rev_encoding: u8) -> Box<[u8]>;
}

pub type Seq2Slice = BitSlice;
pub type Seq2Bits = BitBox;

pub trait RangeNuc {
    fn range_nuc(&self, range: std::ops::Range<usize>) -> &Seq2Slice;
}

impl Nuc2Bits for Seq2Bits {
    fn from_nuc(input: &[u8], encoding: u8) -> Self {
        // Todo add check
        seq2bits(&input, encoding)
    }

    fn from_bitslice(input: &BitSlice) -> Self {
        // Todo add check
        BitBox::from_bitslice(input)
    }

    fn from_bits(input: &[u8], nb_nuc: usize) -> Self {
        // Todo add check
        let vec = BitVec::from_vec(input.to_vec());

        BitBox::from_bitslice(&vec[..nb_nuc * 2])
    }
}

impl Bits2Nuc for Seq2Bits {
    fn into_nuc(&self, rev_encoding: u8) -> Box<[u8]> {
        // Todo add check
        bits2seq(&self, rev_encoding)
    }
}

impl RangeNuc for Seq2Bits {
    fn range_nuc(&self, range: std::ops::Range<usize>) -> &Seq2Slice {
        // Todo add check
        &self.as_bitslice().index((range.start * 2)..(range.end * 2))
    }
}

impl Bits2Nuc for Seq2Slice {
    fn into_nuc(&self, rev_encoding: u8) -> Box<[u8]> {
        // Todo add check
        bits2seq(&self, rev_encoding)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitvec::prelude::*;

    #[test]
    fn from_nuc() {
        assert_eq!(
            bitbox![Msb0, u8; 0, 1, 0, 0],
            Seq2Bits::from_nuc(b"CA", 0b00011011)
        );
        assert_eq!(
            bitbox![Msb0, u8; 0, 1, 0, 0],
            Seq2Bits::from_nuc(b"AC", 0b01001011)
        );
    }

    #[test]
    fn from_bitslicec() {
        assert_eq!(
            bitbox![Msb0, u8; 0, 1, 0, 0],
            Seq2Bits::from_bitslice(&bitbox![Msb0, u8; 0, 1, 0, 0][..])
        );
        assert_eq!(
            bitbox![Msb0, u8; 0, 1, 1, 1],
            Seq2Bits::from_bitslice(&bitbox![Msb0, u8; 0, 1, 1, 1][..])
        );
    }

    #[test]
    fn from_bits() {
        assert_eq!(
            bitbox![Msb0, u8; 0, 1, 0, 0],
            Seq2Bits::from_bits(&[0b01000000], 2)
        );
        assert_eq!(
            bitbox![Msb0, u8; 0, 1, 0, 0],
            Seq2Bits::from_bits(&[0b01001111], 2)
        );
    }

    #[test]
    fn into_nuc() {
        assert_eq!(
            &Seq2Bits::from_nuc(b"CA", 0b00011011)
                .into_nuc(0b00011011)
                .into_vec(),
            b"CA"
        );
        assert_eq!(
            &Seq2Bits::from_nuc(b"CATG", 0b00011011)
                .into_nuc(0b00011011)
                .into_vec(),
            b"CATG"
        );
        assert_eq!(
            &Seq2Bits::from_nuc(b"CATGA", 0b00011011)
                .into_nuc(0b00011011)
                .into_vec(),
            b"CATGA"
        );
    }
}
