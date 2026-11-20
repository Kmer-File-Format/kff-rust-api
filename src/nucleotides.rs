//! Operations on right aligned 2-bit sequences
//!
//! The nucleotide `i` of a byte array is stored in the bits `7 - 2(i % 4)` and `6 - 2(i % 4)` of
//! its byte `i / 4`. A sequence of `n` nucleotides is right aligned: its first nucleotide is at
//! position `(4 - n % 4) % 4` of the array.

/// Number of bytes used by a sequence of `nb` nucleotides
#[inline]
pub fn bytes_for_nucleotides(nb: u64) -> usize {
    ((nb + 3) / 4) as usize
}

/// Position of the first nucleotide of a right aligned sequence of `nb` nucleotides
#[inline]
pub(crate) fn padding(nb: u64) -> usize {
    ((4 - nb % 4) % 4) as usize
}

/// Copy of the subsequence of `len` nucleotides starting at nucleotide `start` of a right aligned
/// sequence of `seq_size` nucleotides. The result is right aligned with null padding bits.
///
/// # Panics
/// If the subsequence is not in the sequence or if `seq` is shorter than `ceil(seq_size / 4)`.
pub fn subsequence(seq: &[u8], seq_size: u64, start: u64, len: u64) -> Vec<u8> {
    assert!(start + len <= seq_size, "subsequence out of the sequence");
    let mut out = vec![0; bytes_for_nucleotides(len)];
    copy_nucleotides(
        &mut out,
        padding(len),
        seq,
        padding(seq_size) + start as usize,
        len as usize,
    );
    out
}

/// Value of the `nb` last nucleotides (`nb <= 32`) of a right aligned sequence. The padding bits are
/// ignored.
#[inline]
pub(crate) fn load_nucleotides(seq: &[u8], nb: u64) -> u64 {
    let mut value = 0u64;
    for byte in &seq[..bytes_for_nucleotides(nb)] {
        value = (value << 8) | *byte as u64;
    }
    if nb >= 32 {
        value
    } else {
        value & ((1u64 << (2 * nb)) - 1)
    }
}

/// Write a value of `nb` nucleotides (`nb <= 32`) right aligned in `ceil(nb / 4)` bytes
#[inline]
pub(crate) fn store_nucleotides(seq: &mut [u8], nb: u64, mut value: u64) {
    for byte in seq[..bytes_for_nucleotides(nb)].iter_mut().rev() {
        *byte = value as u8;
        value >>= 8;
    }
}

/// Left shift of `nb` nucleotides, 0 when all the bits are shifted out
#[inline]
fn shift_nucleotides(value: u64, nb: u64) -> u64 {
    if nb >= 32 {
        0
    } else {
        value << (2 * nb)
    }
}

/// Copy `nb` nucleotides of `src`, starting at nucleotide `src_pos`, into `dst` at nucleotide
/// `dst_pos`. The bits of the copied area must be 0 in `dst` (bitwise or).
pub(crate) fn copy_nucleotides(
    dst: &mut [u8],
    dst_pos: usize,
    src: &[u8],
    src_pos: usize,
    nb: usize,
) {
    if nb == 0 {
        return;
    }
    let src = &src[src_pos / 4..];
    let dst = &mut dst[dst_pos / 4..];
    let src_bit = 2 * (src_pos % 4);
    let dst_bit = 2 * (dst_pos % 4);
    let bits = 2 * nb;
    let src_bytes = (src_bit + bits + 7) / 8;
    let dst_bytes = (dst_bit + bits + 7) / 8;
    // Mask of the bits of the last destination byte that are part of the copy
    let last_bits = dst_bit + bits - 8 * (dst_bytes - 1);
    let last_mask = (0xFFu16 << (8 - last_bits)) as u8;

    if src_bit == dst_bit {
        // Same alignment: whole bytes are copied
        if dst_bytes == 1 {
            dst[0] |= src[0] & (0xFF >> dst_bit) & last_mask;
            return;
        }
        dst[0] |= src[0] & (0xFF >> dst_bit);
        dst[1..dst_bytes - 1].copy_from_slice(&src[1..dst_bytes - 1]);
        dst[dst_bytes - 1] |= src[dst_bytes - 1] & last_mask;
        return;
    }

    // 8 bits of the source from the bit p (only the bytes of the source are read)
    let get8 = |p: usize| -> u8 {
        let byte = p / 8;
        let offset = p % 8;
        let mut word = (src[byte] as u16) << 8;
        if offset != 0 && byte + 1 < src_bytes {
            word |= src[byte + 1] as u16;
        }
        (word >> (8 - offset)) as u8
    };

    // First destination byte: bits dst_bit..7 receive the source bits from src_bit
    let first = if dst_bit > src_bit {
        get8(0) >> (dst_bit - src_bit)
    } else {
        get8(src_bit - dst_bit)
    };
    let mut first_mask = 0xFFu8 >> dst_bit;
    if dst_bytes == 1 {
        first_mask &= last_mask;
    }
    dst[0] |= first & first_mask;
    if dst_bytes == 1 {
        return;
    }
    // Following bytes: the destination bit 8j comes from the source bit src_bit + 8j - dst_bit
    for (j, byte) in dst.iter_mut().enumerate().take(dst_bytes - 1).skip(1) {
        *byte = get8(src_bit + 8 * j - dst_bit);
    }
    dst[dst_bytes - 1] |= get8(src_bit + 8 * (dst_bytes - 1) - dst_bit) & last_mask;
}

/// Write in `out` the right aligned sequence `seq` of `seq_size` nucleotides without its `m`
/// nucleotides starting at `pos`. `out` must be zeroed and contain at least
/// `ceil((seq_size - m) / 4)` bytes.
pub(crate) fn remove_nucleotides(out: &mut [u8], seq: &[u8], seq_size: u64, pos: u64, m: u64) {
    let out_size = seq_size - m;
    if seq_size <= 32 {
        // Short sequences: 64 bits arithmetic
        let suffix_size = seq_size - pos - m;
        let full = load_nucleotides(seq, seq_size);
        let suffix = full & (shift_nucleotides(1, suffix_size).wrapping_sub(1));
        let prefix = if pos == 0 {
            0
        } else {
            full >> (2 * (suffix_size + m))
        };
        store_nucleotides(
            out,
            out_size,
            shift_nucleotides(prefix, suffix_size) | suffix,
        );
    } else {
        // Prefix then suffix of the removed part
        let seq_offset = padding(seq_size);
        let out_offset = padding(out_size);
        copy_nucleotides(out, out_offset, seq, seq_offset, pos as usize);
        copy_nucleotides(
            out,
            out_offset + pos as usize,
            seq,
            seq_offset + (pos + m) as usize,
            (seq_size - pos - m) as usize,
        );
    }
}

/// Write in `out` the right aligned sequence `seq` of `seq_size` nucleotides where the `m`
/// nucleotides of `inserted` are inserted at position `pos`. `out` must be zeroed and contain at
/// least `ceil((seq_size + m) / 4)` bytes.
pub(crate) fn insert_nucleotides(
    out: &mut [u8],
    seq: &[u8],
    seq_size: u64,
    pos: u64,
    inserted: &[u8],
    m: u64,
) {
    let out_size = seq_size + m;
    if out_size <= 32 {
        // Short sequences: 64 bits arithmetic
        let suffix_size = seq_size - pos;
        let no_mini = load_nucleotides(seq, seq_size);
        let suffix = no_mini & (shift_nucleotides(1, suffix_size).wrapping_sub(1));
        let prefix = if pos == 0 {
            0
        } else {
            no_mini >> (2 * suffix_size)
        };
        let mini = load_nucleotides(inserted, m);
        let full = shift_nucleotides(prefix, suffix_size + m)
            | shift_nucleotides(mini, suffix_size)
            | suffix;
        store_nucleotides(out, out_size, full);
    } else {
        let out_offset = padding(out_size);
        let seq_offset = padding(seq_size);
        let mini_offset = padding(m);
        copy_nucleotides(out, out_offset, seq, seq_offset, pos as usize);
        copy_nucleotides(
            out,
            out_offset + pos as usize,
            inserted,
            mini_offset,
            m as usize,
        );
        copy_nucleotides(
            out,
            out_offset + (pos + m) as usize,
            seq,
            seq_offset + pos as usize,
            (seq_size - pos) as usize,
        );
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Deterministic pseudo-random generator (xorshift64*) for the tests
    pub(crate) struct Random(u64);

    impl Random {
        pub(crate) fn new(seed: u64) -> Self {
            Self(seed.wrapping_mul(0x9E3779B97F4A7C15) | 1)
        }

        pub(crate) fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545F4914F6CDD1D)
        }

        /// Inclusive bounds
        pub(crate) fn between(&mut self, min: u64, max: u64) -> u64 {
            min + self.next() % (max - min + 1)
        }

        pub(crate) fn bytes(&mut self, len: usize) -> Vec<u8> {
            (0..len).map(|_| self.next() as u8).collect()
        }

        pub(crate) fn sequence(&mut self, len: usize) -> Vec<u8> {
            (0..len)
                .map(|_| b"ACGT"[self.between(0, 3) as usize])
                .collect()
        }
    }

    /// Nucleotide by nucleotide reference implementation
    fn get(seq: &[u8], pos: usize) -> u8 {
        (seq[pos / 4] >> (2 * (3 - pos % 4))) & 0b11
    }

    fn naive_copy(dst: &mut [u8], dst_pos: usize, src: &[u8], src_pos: usize, nb: usize) {
        for i in 0..nb {
            let pos = dst_pos + i;
            dst[pos / 4] |= get(src, src_pos + i) << (2 * (3 - pos % 4));
        }
    }

    #[test]
    fn bytes_and_padding() {
        assert_eq!(bytes_for_nucleotides(0), 0);
        assert_eq!(bytes_for_nucleotides(1), 1);
        assert_eq!(bytes_for_nucleotides(4), 1);
        assert_eq!(bytes_for_nucleotides(5), 2);
        assert_eq!(padding(0), 0);
        assert_eq!(padding(1), 3);
        assert_eq!(padding(3), 1);
        assert_eq!(padding(8), 0);
    }

    #[test]
    fn copy_matches_naive_copy() {
        let mut rnd = Random::new(1);
        for _ in 0..20000 {
            let src_len = rnd.between(1, 40) as usize;
            let src = rnd.bytes(src_len);
            let src_pos = rnd.between(0, (src.len() * 4 - 1) as u64) as usize;
            let nb = rnd.between(0, (src.len() * 4 - src_pos) as u64) as usize;
            let dst_len = rnd.between(1, 45) as usize;
            if nb > dst_len * 4 {
                continue;
            }
            let dst_pos = rnd.between(0, (dst_len * 4 - nb) as u64) as usize;

            let mut expected = vec![0; dst_len];
            naive_copy(&mut expected, dst_pos, &src, src_pos, nb);
            let mut dst = vec![0; dst_len];
            copy_nucleotides(&mut dst, dst_pos, &src, src_pos, nb);
            assert_eq!(dst, expected, "{src_pos} {dst_pos} {nb}");
        }
    }

    #[test]
    fn remove_and_insert_are_inverse() {
        let encoding = crate::Encoding::default();
        let mut rnd = Random::new(2);
        for _ in 0..20000 {
            let size = rnd.between(1, 100);
            let m = rnd.between(0, size);
            let pos = rnd.between(0, size - m);
            let seq = rnd.sequence(size as usize);
            let bytes = encoding.encode(&seq).unwrap();

            let mut without_seq = seq[..pos as usize].to_vec();
            without_seq.extend_from_slice(&seq[(pos + m) as usize..]);
            let without = encoding.encode(&without_seq).unwrap();
            let mini = encoding
                .encode(&seq[pos as usize..(pos + m) as usize])
                .unwrap();

            // Padding bits of the input are ignored
            let mut noisy = bytes.clone();
            if !noisy.is_empty() && padding(size) > 0 {
                noisy[0] |= 0xFF << (8 - 2 * padding(size));
            }
            let mut removed = vec![0; bytes_for_nucleotides(size - m)];
            remove_nucleotides(&mut removed, &noisy, size, pos, m);
            assert_eq!(removed, without, "remove {size} {m} {pos}");

            let mut inserted = vec![0; bytes_for_nucleotides(size)];
            insert_nucleotides(&mut inserted, &without, size - m, pos, &mini, m);
            assert_eq!(inserted, bytes, "insert {size} {m} {pos}");
        }
    }

    #[test]
    fn subsequences() {
        let encoding = crate::Encoding::default();
        let mut rnd = Random::new(3);
        for _ in 0..5000 {
            let size = rnd.between(1, 80);
            let seq = rnd.sequence(size as usize);
            let bytes = encoding.encode(&seq).unwrap();
            let start = rnd.between(0, size - 1);
            let len = rnd.between(0, size - start);
            assert_eq!(
                subsequence(&bytes, size, start, len),
                encoding
                    .encode(&seq[start as usize..(start + len) as usize])
                    .unwrap()
            );
        }
    }
}
