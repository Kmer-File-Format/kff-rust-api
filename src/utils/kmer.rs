//! Representation of a KFF kmer block

/* std use */

/* crate use */

/* project use */

/// Represent a sequence in 2 bit
pub type Seq2Bit = bitvec::boxed::BitBox<u8, bitvec::order::Msb0>;

/// Represent data associate to a kmer
pub type Data = Vec<u8>;

/// Represent a Kmer with data
#[derive(
    getset::Getters,
    getset::Setters,
    getset::MutGetters,
    std::fmt::Debug,
    std::cmp::PartialEq,
    std::cmp::Eq,
    std::default::Default,
)]
#[getset(get = "pub", set = "pub", get_mut = "pub")]
pub struct Kmer {
    /// 2 bit field for sequenc
    seq2bit: Seq2Bit,
    /// Data associate to kmer
    data: Data,
}

impl Kmer {
    /// Create a new kmer
    pub fn new(seq2bit: Seq2Bit, data: Data) -> Self {
        Self { seq2bit, data }
    }

    /// Create a Kmer for ascii sequence
    pub fn from_ascii(seq: &[u8], data: Data, encoding: u8) -> Self {
        Self {
            seq2bit: seq2bits(seq, encoding),
            data,
        }
    }

    /// Get seq in ascii
    pub fn seq(&self, encoding: u8) -> Vec<u8> {
        bits2seq(&self.seq2bit, encoding)
    }
}

/// Convert a nucleotide in internal encoding
#[inline]
fn nuc2internal(nuc: u8) -> u8 {
    (nuc >> 1) & 0b11
}

const INTERNAL2NUC: [u8; 4] = [b'A', b'C', b'T', b'G'];

/// Convert internal encoding to nucleotide
#[inline]
fn internal2nuc(internal: u8) -> u8 {
    INTERNAL2NUC[internal as usize]
}

/// Convert nucleotide in 2bit encoding
#[inline]
fn nuc2encoding(nuc: u8, encoding: u8) -> u8 {
    let index = nuc2internal(nuc) * 2;

    (encoding << index) & 0b11000000
}

/// Convert 2bit encoding in nucleotide
#[inline]
fn encoding2nuc(bits: u8, rev_encoding: u8) -> u8 {
    internal2nuc((rev_encoding >> (6 - ((bits >> 6) * 2))) & 0b11)
}

/// Convert a nucleotide in 2bit encoding
#[inline]
fn nuc2bits(nuc: u8, encoding: u8) -> Seq2Bit {
    let mut tmp = bitvec::vec::BitVec::from_vec(vec![nuc2encoding(nuc, encoding)]);

    tmp.resize(2, false);

    tmp.into_boxed_bitslice()
}

/// Convert a sequence of nucleotide in Seq2Bit
pub fn seq2bits(seq: &[u8], encoding: u8) -> Seq2Bit {
    let mut bits = bitvec::vec::BitVec::with_capacity(seq.len() * 2);

    for nuc in seq {
        bits.extend_from_bitslice(&nuc2bits(*nuc, encoding));
    }

    bits.into_boxed_bitslice()
}

/// Convert a Seq2Bit in sequence of nucleotide
pub fn bits2seq(bits: &Seq2Bit, encoding: u8) -> Vec<u8> {
    let rev_encoding = rev_encoding(encoding);

    let mut ret = Vec::with_capacity(bits.len());

    for bit in bits.chunks(2) {
        ret.push(encoding2nuc(
            ((bit[0] as u8) << 7) ^ (bit[1] as u8) << 6,
            rev_encoding,
        ))
    }

    ret
}

/// Convert an encoding in reverse version
#[inline]
fn rev_encoding(encoding: u8) -> u8 {
    let mut rev = 0;

    rev ^= 0b00 << (6 - ((encoding >> 6) * 2));
    rev ^= 0b01 << (6 - (((encoding >> 4) & 0b11) * 2));
    rev ^= 0b10 << (6 - (((encoding >> 2) & 0b11) * 2));
    rev ^= 0b11 << (6 - ((encoding & 0b11) * 2));

    rev
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::error;

    #[test]
    fn create() -> error::Result<()> {
        let mut kmer = Kmer::new(
            bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1],
            vec![1],
        );

        assert_eq!(
            kmer.seq2bit(),
            &bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1]
        );
        kmer.set_seq2bit(bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0]);
        assert_eq!(
            kmer.seq2bit(),
            &bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0]
        );
        kmer.seq2bit_mut().as_raw_mut_slice()[0] = 40;
        assert_eq!(
            kmer.seq2bit(),
            &bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 1, 0, 1, 0]
        );

        Ok(())
    }

    #[test]
    fn seq2bit() -> error::Result<()> {
        let encoding = 0b00011011;

        let kmer = Kmer::from_ascii(b"ACTG", vec![1], encoding);

        assert_eq!(kmer.seq(encoding), b"ACTG");

        Ok(())
    }

    #[test]
    fn internal_encoding() {
        assert_eq!(nuc2internal(b'A'), 0);
        assert_eq!(nuc2internal(b'C'), 1);
        assert_eq!(nuc2internal(b'T'), 2);
        assert_eq!(nuc2internal(b'G'), 3);
    }

    #[test]
    fn internal_decoding() {
        assert_eq!(internal2nuc(0), b'A');
        assert_eq!(internal2nuc(1), b'C');
        assert_eq!(internal2nuc(2), b'T');
        assert_eq!(internal2nuc(3), b'G');
    }

    #[test]
    fn encoding() {
        let encoding = 0b11100100;

        assert_eq!(nuc2encoding(b'A', encoding), 0b11000000);
        assert_eq!(nuc2encoding(b'C', encoding), 0b10000000);
        assert_eq!(nuc2encoding(b'T', encoding), 0b01000000);
        assert_eq!(nuc2encoding(b'G', encoding), 0b00000000);
    }

    #[test]
    fn decoding() {
        let mut rencoding = rev_encoding(0b00011011);

        assert_eq!(encoding2nuc(0b00000000, rencoding), b'A');
        assert_eq!(encoding2nuc(0b01000000, rencoding), b'C');
        assert_eq!(encoding2nuc(0b10000000, rencoding), b'T');
        assert_eq!(encoding2nuc(0b11000000, rencoding), b'G');

        rencoding = rev_encoding(0b11100100);

        assert_eq!(encoding2nuc(0b11000000, rencoding), b'A');
        assert_eq!(encoding2nuc(0b10000000, rencoding), b'C');
        assert_eq!(encoding2nuc(0b01000000, rencoding), b'T');
        assert_eq!(encoding2nuc(0b00000000, rencoding), b'G');

        rencoding = rev_encoding(0b01110010);

        assert_eq!(encoding2nuc(0b01000000, rencoding), b'A');
        assert_eq!(encoding2nuc(0b11000000, rencoding), b'C');
        assert_eq!(encoding2nuc(0b00000000, rencoding), b'T');
        assert_eq!(encoding2nuc(0b10000000, rencoding), b'G');
    }

    #[test]
    fn nuc2bits_() {
        let mut encoding = 0b00011011;

        assert_eq!(
            nuc2bits(b'A', encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0]
        );
        assert_eq!(
            nuc2bits(b'C', encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1]
        );
        assert_eq!(
            nuc2bits(b'T', encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0; 1, 0]
        );
        assert_eq!(
            nuc2bits(b'G', encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0; 1, 1]
        );

        encoding = 0b01110010;

        assert_eq!(
            nuc2bits(b'A', encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1]
        );
        assert_eq!(
            nuc2bits(b'C', encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0; 1, 1]
        );
        assert_eq!(
            nuc2bits(b'T', encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0]
        );
        assert_eq!(
            nuc2bits(b'G', encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0; 1, 0]
        );
    }

    #[test]
    fn seq2bits_() {
        let encoding = 0b00011011;

        assert_eq!(
            seq2bits(b"AC", encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1]
        );
        assert_eq!(
            seq2bits(b"ACG", encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 1]
        );
        assert_eq!(
            seq2bits(b"ACGTA", encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 1, 1, 0, 0, 0]
        );
    }

    #[test]
    fn bits2seq_() {
        let encoding = 0b00011011;

        assert_eq!(
            seq2bits(b"AC", encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0;0, 0, 0, 1]
        );

        assert_eq!(
            seq2bits(b"ACGT", encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0;0, 0, 0, 1, 1, 1, 1, 0]
        );

        assert_eq!(
            seq2bits(b"ACGTG", encoding),
            bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 1, 1, 0, 1, 1]
        );

        assert_eq!(
            bits2seq(&seq2bits(b"AC", encoding), rev_encoding(encoding)),
            vec![b'A', b'C']
        );
        assert_eq!(
            bits2seq(&seq2bits(b"ACGT", encoding), rev_encoding(encoding)),
            vec![b'A', b'C', b'G', b'T']
        );
        assert_eq!(
            bits2seq(&seq2bits(b"ACGTG", encoding), rev_encoding(encoding)),
            vec![b'A', b'C', b'G', b'T', b'G']
        );
    }

    #[test]
    fn rev_encoding_() {
        assert_eq!(rev_encoding(0b00011011), 0b00011011);

        assert_eq!(rev_encoding(0b11100100), 0b11100100);

        assert_eq!(rev_encoding(0b01110010), 0b10001101);
    }
}
