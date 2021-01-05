//! Define a kmer type

/* local use */
use crate::seq2bits::Seq2Bits;

/// Kmer type
#[derive(Debug, PartialEq)]
pub struct Kmer {
    seq: Seq2Bits,
    data: Box<[u8]>,
}

impl Kmer {
    /// Create a new kmer
    pub fn new(seq: Seq2Bits, data: Box<[u8]>) -> Self {
        Self { seq, data }
    }

    // Getter
    /// Get sequence in 2 bits encoding
    pub fn seq(&self) -> &Seq2Bits {
        &self.seq
    }

    /// Get data
    pub fn data(&self) -> &[u8] {
        &self.data
    }

    /// Get length of kmer
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Return true if kmer length is 0
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}
