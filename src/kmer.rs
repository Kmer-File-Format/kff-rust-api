/* local use */
use crate::utils;

#[derive(Debug, PartialEq)]
pub struct Kmer {
    bits: utils::BitBox,
    data: Box<[u8]>,
}

impl Kmer {
    pub fn new(bits: utils::BitBox, data: Box<[u8]>) -> Self {
        Self { bits, data }
    }

    // Getter
    pub fn bits(&self) -> &utils::BitSlice {
        &self.bits
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }

    pub fn seq(&self, rev_encoding: u8) -> Box<[u8]> {
        utils::bits2seq(&self.bits, rev_encoding)
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}
