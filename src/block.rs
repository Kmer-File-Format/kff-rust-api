//! Layout of the blocks of the raw and minimizer sections

/* project use */
use crate::error;
use crate::nucleotides::{bytes_for_nucleotides, subsequence};
use crate::section::{bytes_for_max, Values};

/// Field widths and variables of the blocks of a raw or minimizer section
#[derive(Clone, std::fmt::Debug, PartialEq, Eq)]
pub(crate) struct BlockFormat {
    /// True for a minimizer section
    pub minimizer: bool,
    pub k: u64,
    /// Size of the minimizer, 0 for raw sections
    pub m: u64,
    pub max: u64,
    pub data_size: u64,
    /// Width of the number of kmers of a block (0 when max = 1)
    pub nb_kmers_bytes: usize,
    /// Width of the position of the minimizer
    pub mini_pos_bytes: usize,
    /// Size of the minimizer in bytes
    pub mini_bytes: usize,
}

fn variable(vars: &Values, name: &'static str) -> error::Result<u64> {
    vars.get(name)
        .copied()
        .ok_or(error::Error::MissingVariable(name))
}

fn invalid(name: &'static str, reason: &str) -> error::Error {
    error::Error::InvalidVariable {
        name,
        reason: reason.to_string(),
    }
}

impl BlockFormat {
    /// Format of a raw section with the variables of the current values section
    pub fn raw(vars: &Values) -> error::Result<Self> {
        Self::new(vars, false)
    }

    /// Format of a minimizer section with the variables of the current values section
    pub fn minimizer(vars: &Values) -> error::Result<Self> {
        Self::new(vars, true)
    }

    fn new(vars: &Values, minimizer: bool) -> error::Result<Self> {
        let k = variable(vars, "k")?;
        let m = if minimizer { variable(vars, "m")? } else { 0 };
        let max = variable(vars, "max")?;
        let data_size = variable(vars, "data_size")?;

        if k == 0 {
            return Err(invalid("k", "k must be at least 1"));
        }
        if max == 0 {
            return Err(invalid("max", "max must be at least 1"));
        }
        if m > k {
            return Err(invalid("m", "the minimizer can't be longer than k"));
        }
        if k.checked_add(max).is_none() || data_size.checked_mul(max).is_none() {
            return Err(invalid("max", "the blocks are too large"));
        }

        Ok(Self {
            minimizer,
            k,
            m,
            max,
            data_size,
            nb_kmers_bytes: bytes_for_max(max),
            mini_pos_bytes: if minimizer {
                bytes_for_max(k + max - 1)
            } else {
                0
            },
            mini_bytes: bytes_for_nucleotides(m),
        })
    }

    /// Number of nucleotides of a block of `nb_kmers` kmers
    #[inline]
    pub fn sequence_size(&self, nb_kmers: u64) -> u64 {
        nb_kmers + self.k - 1
    }

    /// Number of nucleotides stored in a block of `nb_kmers` kmers (without the minimizer)
    #[inline]
    pub fn stored_size(&self, nb_kmers: u64) -> u64 {
        nb_kmers + self.k - 1 - self.m
    }

    /// Number of bytes of the data of a block of `nb_kmers` kmers
    #[inline]
    pub fn data_bytes(&self, nb_kmers: u64) -> usize {
        (nb_kmers * self.data_size) as usize
    }

    /// Check the number of kmers of a block read or to write
    pub fn check_nb_kmers(&self, nb_kmers: u64) -> error::Result<()> {
        if nb_kmers == 0 || nb_kmers > self.max {
            return Err(error::Error::InvalidBlock(format!(
                "{} kmers in a block, the section accepts 1 to {} kmers (max)",
                nb_kmers, self.max
            )));
        }
        // Without field (max = 1), the number of kmers is implicit
        if self.nb_kmers_bytes > 0
            && self.nb_kmers_bytes < 8
            && nb_kmers >> (8 * self.nb_kmers_bytes) != 0
        {
            return Err(error::Error::InvalidBlock(format!(
                "{} kmers do not fit in the {} bytes of the number of kmers (ceil(log2(max)) bits, max = {})",
                nb_kmers, self.nb_kmers_bytes, self.max
            )));
        }
        Ok(())
    }

    /// Check the position of the minimizer of a block of `nb_kmers` kmers
    pub fn check_mini_pos(&self, nb_kmers: u64, mini_pos: u64) -> error::Result<()> {
        if mini_pos > self.stored_size(nb_kmers) {
            return Err(error::Error::InvalidBlock(format!(
                "minimizer at position {} in a sequence of {} nucleotides (m = {})",
                mini_pos,
                self.sequence_size(nb_kmers),
                self.m
            )));
        }
        if self.mini_pos_bytes < 8 && mini_pos >> (8 * self.mini_pos_bytes) != 0 {
            return Err(error::Error::InvalidBlock(format!(
                "minimizer position {} does not fit in {} bytes",
                mini_pos, self.mini_pos_bytes
            )));
        }
        Ok(())
    }
}

/// Check that a buffer contains at least `needed` bytes
#[inline]
pub(crate) fn check_buffer(needed: usize, available: usize) -> error::Result<()> {
    if available < needed {
        Err(error::Error::BufferTooSmall { needed, available })
    } else {
        Ok(())
    }
}

/// A block read from a raw or minimizer section: a sequence and the data of its kmers
#[derive(Clone, std::fmt::Debug, PartialEq, Eq, Default)]
pub struct Block {
    /// Size of the kmers
    pub k: u64,

    /// Number of kmers of the block
    pub nb_kmers: u64,

    /// Size of the data of one kmer in bytes
    pub data_size: u64,

    /// Sequence of the block (`nb_kmers + k - 1` nucleotides, minimizer included), right aligned
    /// 2-bit nucleotides with null padding bits
    pub sequence: Vec<u8>,

    /// Data of the kmers, `data_size` bytes per kmer
    pub data: Vec<u8>,
}

impl Block {
    /// Number of nucleotides of the sequence
    pub fn sequence_size(&self) -> u64 {
        if self.nb_kmers == 0 {
            0
        } else {
            self.nb_kmers + self.k - 1
        }
    }

    /// The `i`-th kmer of the block (right aligned 2-bit nucleotides)
    ///
    /// # Panics
    /// If `i >= nb_kmers`.
    pub fn kmer(&self, i: u64) -> Vec<u8> {
        assert!(i < self.nb_kmers, "kmer {} out of the block", i);
        subsequence(&self.sequence, self.sequence_size(), i, self.k)
    }

    /// Data of the `i`-th kmer of the block
    ///
    /// # Panics
    /// If `i >= nb_kmers`.
    pub fn kmer_data(&self, i: u64) -> &[u8] {
        assert!(i < self.nb_kmers, "kmer {} out of the block", i);
        let start = (i * self.data_size) as usize;
        &self.data[start..start + self.data_size as usize]
    }

    /// Kmers of the block with their data
    pub fn kmers(&self) -> impl Iterator<Item = (Vec<u8>, &[u8])> + '_ {
        (0..self.nb_kmers).map(move |i| (self.kmer(i), self.kmer_data(i)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vars(list: &[(&str, u64)]) -> Values {
        list.iter().map(|(n, v)| (n.to_string(), *v)).collect()
    }

    #[test]
    fn missing_variables() {
        assert!(matches!(
            BlockFormat::raw(&vars(&[("k", 5), ("max", 5)])),
            Err(error::Error::MissingVariable("data_size"))
        ));
        assert!(matches!(
            BlockFormat::raw(&vars(&[("max", 5), ("data_size", 5)])),
            Err(error::Error::MissingVariable("k"))
        ));
        assert!(matches!(
            BlockFormat::raw(&vars(&[("k", 5), ("data_size", 5)])),
            Err(error::Error::MissingVariable("max"))
        ));
        assert!(matches!(
            BlockFormat::minimizer(&vars(&[("k", 5), ("max", 5), ("data_size", 0)])),
            Err(error::Error::MissingVariable("m"))
        ));
        assert!(BlockFormat::raw(&vars(&[("k", 5), ("max", 5), ("data_size", 0)])).is_ok());
    }

    #[test]
    fn invalid_variables() {
        assert!(BlockFormat::raw(&vars(&[("k", 0), ("max", 5), ("data_size", 0)])).is_err());
        assert!(BlockFormat::raw(&vars(&[("k", 5), ("max", 0), ("data_size", 0)])).is_err());
        assert!(
            BlockFormat::minimizer(&vars(&[("k", 5), ("m", 6), ("max", 5), ("data_size", 0)]))
                .is_err()
        );
        assert!(BlockFormat::raw(&vars(&[("k", 5), ("max", u64::MAX), ("data_size", 0)])).is_err());
    }

    #[test]
    fn widths() {
        // ceil(log2(k + max - 1)) bits for the minimizer position (tests of the C++ API)
        let cases: &[(u64, u64, usize)] = &[
            (1, 1, 0),
            (2, 1, 1),
            (3, 255, 2),
            (10, 247, 1),
            (10, 248, 2),
            (31, 226, 1),
            (31, 227, 2),
            (31, 65506, 2),
            (31, 65507, 3),
        ];
        for (k, max, bytes) in cases {
            let format = BlockFormat::minimizer(&vars(&[
                ("k", *k),
                ("m", 1),
                ("max", *max),
                ("data_size", 0),
            ]))
            .unwrap();
            assert_eq!(format.mini_pos_bytes, *bytes, "k {k} max {max}");
        }
    }

    #[test]
    fn nb_kmers_limits() {
        let format = BlockFormat::raw(&vars(&[("k", 5), ("max", 256), ("data_size", 0)])).unwrap();
        assert_eq!(format.nb_kmers_bytes, 1);
        assert!(format.check_nb_kmers(255).is_ok());
        // 256 does not fit in ceil(log2(256)) = 8 bits
        assert!(format.check_nb_kmers(256).is_err());
        assert!(format.check_nb_kmers(0).is_err());
        let format = BlockFormat::raw(&vars(&[("k", 5), ("max", 1), ("data_size", 0)])).unwrap();
        assert_eq!(format.nb_kmers_bytes, 0);
        assert!(format.check_nb_kmers(1).is_ok());
        assert!(format.check_nb_kmers(2).is_err());
    }

    #[test]
    fn block_kmers() {
        let encoding = crate::Encoding::default();
        let block = Block {
            k: 3,
            nb_kmers: 3,
            data_size: 1,
            sequence: encoding.encode(b"ACGTA").unwrap(),
            data: vec![1, 2, 3],
        };
        let kmers: Vec<(Vec<u8>, Vec<u8>)> = block
            .kmers()
            .map(|(kmer, data)| (encoding.decode(&kmer, 3), data.to_vec()))
            .collect();
        assert_eq!(
            kmers,
            vec![
                (b"ACG".to_vec(), vec![1]),
                (b"CGT".to_vec(), vec![2]),
                (b"GTA".to_vec(), vec![3])
            ]
        );
    }
}
