/* local use */
use crate::seq2bits::Seq2Bits;

#[derive(Debug, PartialEq)]
pub struct Kmer {
    seq: Seq2Bits,
    data: Box<[u8]>,
}

impl Kmer {
    pub fn new(seq: Seq2Bits, data: Box<[u8]>) -> Self {
        Self { seq, data }
    }

    // Getter
    pub fn seq(&self) -> &Seq2Bits {
        &self.seq
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}
