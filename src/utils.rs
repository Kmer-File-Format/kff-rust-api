/* crate use */
use anyhow::Result;
use byteorder::*;

/* local use */
use crate::error;

/// Read order of bytes in file
pub type Order = byteorder::LittleEndian;

/// Order of bit for bitvec
pub type BitOrd = bitvec::order::Msb0;

/// Syntaxic sugar around bitvec::BitBox
pub type BitBox = bitvec::boxed::BitBox<BitOrd, u8>;

/// Syntaxic sugar around bitvec::BitVec
pub type BitVec = bitvec::vec::BitVec<bitvec::order::Msb0, u8>;

/// Syntaxic sugar around bitvec::BitSlice
pub type BitSlice = bitvec::slice::BitSlice<BitOrd, u8>;

pub(crate) fn read_dynamic_size_field<R>(input: &mut R, max_value: u64) -> Result<u64>
where
    R: std::io::Read + ?Sized,
{
    let mut buffer = vec![0u8; bytes_to_store_n(max_value) as usize];

    input.read_exact(&mut buffer)?;

    buffer.resize((9 - bytes_to_store_n(max_value)) as usize, 0);

    Ok(Order::read_u64(&buffer))
}

pub(crate) fn write_dynamic_size_field<W>(
    output: &mut W,
    value: u64,
    max_value: u64,
) -> Result<usize>
where
    W: std::io::Write,
{
    let mut buffer = vec![0u8; 0];
    buffer.write_u64::<Order>(value)?;
    output.write_all(&buffer[..bytes_to_store_n(max_value) as usize])?;

    Ok(bytes_to_store_n(max_value) as usize)
}

#[inline]
pub(crate) fn nuc2internal(nuc: u8) -> u8 {
    (nuc as u8 >> 1) & 0b11
}

const INTERNAL2NUC: [u8; 4] = [b'A', b'C', b'T', b'G'];

#[inline]
pub(crate) fn internal2nuc(internal: u8) -> u8 {
    INTERNAL2NUC[internal as usize]
}

#[inline]
pub(crate) fn nuc2encoding(nuc: u8, encoding: u8) -> u8 {
    let index = nuc2internal(nuc) * 2;

    (encoding << index) & 0b11000000
}

#[inline]
pub(crate) fn encoding2nuc(bits: u8, rev_encoding: u8) -> u8 {
    internal2nuc((rev_encoding >> (6 - ((bits >> 6) * 2))) & 0b11)
}

#[inline]
pub(crate) fn nuc2bits(nuc: u8, encoding: u8) -> BitBox {
    let mut tmp = BitVec::from_vec(vec![nuc2encoding(nuc, encoding)]);

    tmp.resize(2, false);

    tmp.into_boxed_bitslice()
}

pub(crate) fn seq2bits(seq: &[u8], encoding: u8) -> BitBox {
    let mut bits = BitVec::with_capacity(seq.len() * 2);

    for nuc in seq {
        bits.extend_from_bitslice(&nuc2bits(*nuc, encoding));
    }

    bits.into_boxed_bitslice()
}

pub(crate) fn bits2seq(bits: &BitSlice, rev_encoding: u8) -> Box<[u8]> {
    let mut ret = Vec::with_capacity(bits.len());

    for bit in bits.chunks(2) {
        ret.push(encoding2nuc(
            ((bit[0] as u8) << 7) ^ (bit[1] as u8) << 6,
            rev_encoding,
        ))
    }

    ret.into_boxed_slice()
}

#[inline]
pub(crate) fn rev_encoding(encoding: u8) -> u8 {
    let mut rev = 0;

    rev ^= 0b00 << (6 - ((encoding >> 6) * 2));
    rev ^= 0b01 << (6 - (((encoding >> 4) & 0b11) * 2));
    rev ^= 0b10 << (6 - (((encoding >> 2) & 0b11) * 2));
    rev ^= 0b11 << (6 - ((encoding & 0b11) * 2));

    rev
}

#[inline]
pub(crate) fn bytes_to_store_n(n: u64) -> u64 {
    let nb_bytes_needs = ((n + 1) as f64).log2().ceil() as u64;
    ceil_to_8(nb_bytes_needs) / 8
}

#[inline]
pub(crate) fn ceil_to_8(n: u64) -> u64 {
    (n + 7) & !(7)
}

pub(crate) fn switch_56_n_78(input: u8) -> u8 {
    // Internal encoding order is ACTG kff order is ACGT
    (input & 0b11110000) ^ ((input & 0b00000011) << 2) ^ ((input & 0b00001100) >> 2)
}

pub(crate) fn valid_encoding(encoding: u8) -> Result<u8> {
    let a = encoding >> 6;
    let c = (encoding >> 4) & 0b11;
    let t = (encoding >> 2) & 0b11;
    let g = encoding & 0b11;

    if a != c && a != t && a != g && c != t && t != g {
        Ok(encoding)
    } else {
        Err(error::Error::Header(error::Header::BadEncoding).into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitvec::prelude::*;

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

        assert_eq!(nuc2bits(b'A', encoding), bitbox![0, 0]);
        assert_eq!(nuc2bits(b'C', encoding), bitbox![0, 1]);
        assert_eq!(nuc2bits(b'T', encoding), bitbox![1, 0]);
        assert_eq!(nuc2bits(b'G', encoding), bitbox![1, 1]);

        encoding = 0b01110010;

        assert_eq!(nuc2bits(b'A', encoding), bitbox![0, 1]);
        assert_eq!(nuc2bits(b'C', encoding), bitbox![1, 1]);
        assert_eq!(nuc2bits(b'T', encoding), bitbox![0, 0]);
        assert_eq!(nuc2bits(b'G', encoding), bitbox![1, 0]);
    }

    #[test]
    fn seq2bits_() {
        let encoding = 0b00011011;

        assert_eq!(seq2bits(b"AC", encoding), bitbox![0, 0, 0, 1]);
        assert_eq!(seq2bits(b"ACG", encoding), bitbox![0, 0, 0, 1, 1, 1]);
        assert_eq!(
            seq2bits(b"ACGTA", encoding),
            bitbox![0, 0, 0, 1, 1, 1, 1, 0, 0, 0]
        );
    }

    #[test]
    fn bits2seq_() {
        let encoding = 0b00011011;

        assert_eq!(seq2bits(b"AC", encoding), bitbox![0, 0, 0, 1]);

        assert_eq!(seq2bits(b"ACGT", encoding), bitbox![0, 0, 0, 1, 1, 1, 1, 0]);

        assert_eq!(
            seq2bits(b"ACGTG", encoding),
            bitbox![0, 0, 0, 1, 1, 1, 1, 0, 1, 1]
        );

        assert_eq!(
            bits2seq(&seq2bits(b"AC", encoding), rev_encoding(encoding)),
            vec![b'A', b'C'].into_boxed_slice()
        );
        assert_eq!(
            bits2seq(&seq2bits(b"ACGT", encoding), rev_encoding(encoding)),
            vec![b'A', b'C', b'G', b'T'].into_boxed_slice()
        );
        assert_eq!(
            bits2seq(&seq2bits(b"ACGTG", encoding), rev_encoding(encoding)),
            vec![b'A', b'C', b'G', b'T', b'G'].into_boxed_slice()
        );
    }

    #[test]
    fn rev_encoding_() {
        assert_eq!(rev_encoding(0b00011011), 0b00011011);

        assert_eq!(rev_encoding(0b11100100), 0b11100100);

        assert_eq!(rev_encoding(0b01110010), 0b10001101);
    }

    #[test]
    fn padding_computation() {
        assert_eq!(bytes_to_store_n(2), 1);
        assert_eq!(bytes_to_store_n(255), 1);
        assert_eq!(bytes_to_store_n(256), 2);
    }

    #[test]
    fn up_to_next_8() {
        assert_eq!(ceil_to_8(1), 8);
        assert_eq!(ceil_to_8(2), 8);
        assert_eq!(ceil_to_8(15), 16);
        assert_eq!(ceil_to_8(16), 16);
        assert_eq!(ceil_to_8(25), 32);
        assert_eq!(ceil_to_8(46), 48);
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
