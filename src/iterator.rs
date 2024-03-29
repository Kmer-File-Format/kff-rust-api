//! An iterator over all kmer present in a Kff file

/* std use */

/* crate use */

/* project use */
use crate::error;
use crate::Kff;
use crate::Kmer;

/* mod declaration */

/// A Kmer Iterator that consume Kff file
pub struct KmerIterator<R>
where
    R: std::io::Read + std::io::BufRead + crate::KffRead,
{
    inner: Kff<R>,
    front: Option<Box<dyn Iterator<Item = Kmer>>>,
    back: Option<Box<dyn Iterator<Item = Kmer>>>,
}

impl<R> KmerIterator<R>
where
    R: std::io::Read + std::io::BufRead + crate::KffRead,
{
    /// Build a KmerIterator with a Kff object
    pub fn new(inner: Kff<R>) -> Self {
        KmerIterator {
            inner,
            front: None,
            back: None,
        }
    }
}

impl<R> Iterator for KmerIterator<R>
where
    R: std::io::Read + std::io::BufRead + crate::KffRead,
{
    type Item = error::Result<Kmer>;

    fn next(&mut self) -> std::option::Option<Self::Item> {
        loop {
            if let Some(ref mut inner) = self.front {
                match inner.next() {
                    None => self.front = None,
                    Some(elt) => return Some(Ok(elt)),
                }
            }
            match self.inner.next_kmer_section() {
                None => match self.back.as_mut()?.next() {
                    None => {
                        self.back = None;
                        return None;
                    }
                    Some(elt) => return Some(Ok(elt)),
                },
                Some(Ok(inner)) => self.front = Some(Box::new(inner.into_iter())),
                Some(Err(e)) => return Some(Err(e)),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::Seq2Bit;

    const KFF_FILE: &[u8] = &[
        75, 70, 70, 1, 0, 30, 0, 0, 0, 0, 0, 0, 118, 0, 0, 0, 0, 0, 0, 0, 4, 100, 97, 116, 97, 95,
        115, 105, 122, 101, 0, 0, 0, 0, 0, 0, 0, 0, 0, 107, 0, 0, 0, 0, 0, 0, 0, 0, 31, 109, 97,
        120, 0, 0, 0, 0, 0, 0, 0, 0, 255, 111, 114, 100, 101, 114, 101, 100, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 114, 0, 0, 0, 0, 0, 0, 0, 10, 18, 0, 40, 137, 120, 204, 137, 243, 74, 136, 103, 244,
        8, 13, 0, 156, 203, 186, 149, 112, 37, 152, 137, 165, 65, 10, 169, 165, 151, 173, 202, 148,
        201, 132, 210, 91, 3, 0, 171, 159, 211, 30, 240, 136, 83, 211, 12, 7, 208, 106, 177, 81,
        104, 31, 65, 170, 197, 69, 5, 34, 95, 27, 33, 148, 249, 76, 132, 147, 17, 49, 44, 116, 226,
        247, 227, 74, 233, 63, 164, 215, 254, 3, 2, 17, 128, 147, 131, 175, 77, 36, 92, 2, 132, 96,
        36, 224, 235, 211, 73, 21, 2, 62, 203, 109, 65, 161, 178, 163, 184, 105, 0, 0, 0, 0, 0, 0,
        0, 2, 118, 255, 255, 255, 255, 255, 255, 255, 38, 114, 255, 255, 255, 255, 255, 255, 255,
        103, 0, 0, 0, 0, 0, 0, 0, 0, 118, 0, 0, 0, 0, 0, 0, 0, 2, 102, 105, 114, 115, 116, 95, 105,
        110, 100, 101, 120, 0, 0, 0, 0, 0, 0, 0, 0, 195, 102, 111, 111, 116, 101, 114, 95, 115,
        105, 122, 101, 0, 0, 0, 0, 0, 0, 0, 0, 49, 75, 70, 70,
    ];

    #[test]
    fn read_kmer() -> error::Result<()> {
        let reader = Kff::<std::io::BufReader<&[u8]>>::read(std::io::BufReader::new(KFF_FILE))?;

        let kmers: Vec<Seq2Bit> = reader
            .kmers()
            .map(|x| x.unwrap().seq2bit().clone())
            .collect();

        assert_eq!(
            kmers,
            vec![
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1,
                    1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1,
                    0, 0, 1, 1, 0, 1, 0, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1,
                    1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0,
                    1, 1, 0, 1, 0, 0, 1, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0,
                    0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1,
                    0, 1, 0, 0, 1, 0, 1, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0,
                    1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1,
                    0, 0, 1, 0, 1, 0, 1, 0, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1,
                    0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0,
                    1, 0, 1, 0, 1, 0, 0, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0,
                    1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0,
                    1, 0, 1, 0, 0, 0, 1, 0, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1,
                    0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0,
                    1, 0, 0, 0, 1, 0, 0, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
                    1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0,
                    0, 0, 1, 0, 0, 0, 0, 1, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                    0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0,
                    1, 0, 0, 0, 0, 1, 1, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0,
                    1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0,
                    0, 0, 0, 1, 1, 0, 0, 1, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0,
                    0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0,
                    0, 1, 1, 0, 0, 1, 1, 1, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1,
                    1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1,
                    1, 0, 0, 1, 1, 1, 1, 1, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1,
                    1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0,
                    0, 1, 1, 1, 1, 1, 1, 1, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1,
                    0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1,
                    1, 1, 1, 1, 1, 1, 0, 1, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0,
                    1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1,
                    1, 1, 1, 1, 0, 1, 0, 0, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1,
                    0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1,
                    1, 1, 0, 1, 0, 0, 0, 0, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1,
                    0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1,
                    0, 1, 0, 0, 0, 0, 0, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0,
                    1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1,
                    0, 0, 0, 0, 0, 0, 1, 0, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0,
                    1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0,
                    0, 1, 0, 1, 1, 0, 0, 1, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1,
                    1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1,
                    0, 1, 1, 0, 0, 1, 1, 0, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0,
                    1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1,
                    1, 0, 0, 1, 1, 0, 0, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0,
                    1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0,
                    0, 1, 1, 0, 0, 0, 1, 0, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0,
                    0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1,
                    1, 0, 0, 0, 1, 0, 0, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1,
                    0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0,
                    0, 0, 1, 0, 0, 0, 1, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1,
                    0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0,
                    1, 0, 0, 0, 1, 0, 0, 1, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1,
                    0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0,
                    0, 0, 1, 0, 0, 1, 1, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1,
                    1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0,
                    1, 0, 0, 1, 1, 0, 1, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1,
                    0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0,
                    0, 1, 1, 0, 1, 0, 0, 1, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0,
                    0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 1, 0, 1, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0,
                    0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0,
                    1, 0, 0, 1, 0, 1, 0, 1, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0,
                    0, 1, 0, 1, 0, 1, 0, 0, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0,
                    1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0,
                    1, 0, 0, 1, 1, 0, 0, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0,
                    1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                    0, 1, 1, 0, 0, 0, 0, 1, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1,
                    0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1,
                    1, 0, 0, 0, 0, 1, 0, 0, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1,
                    1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0,
                    0, 0, 0, 1, 0, 0, 1, 1, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1,
                    0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0,
                    0, 1, 0, 0, 1, 1, 0, 1, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0,
                    1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1,
                    0, 0, 1, 1, 0, 1, 0, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0,
                    1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0,
                    1, 1, 0, 1, 0, 0, 1, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0,
                    1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1,
                    0, 1, 0, 0, 1, 0, 0, 1, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0,
                    0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1,
                    0, 0, 1, 0, 0, 1, 0, 1, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1,
                    0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0,
                    1, 0, 0, 1, 0, 1, 1, 0, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1,
                    0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0,
                    1, 0, 0, 0, 0, 1, 0, 1, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1,
                    0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0,
                    0, 0, 0, 1, 0, 1, 0, 0, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0,
                    1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
                    0, 1, 0, 1, 0, 0, 1, 1, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0,
                    1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1,
                    1, 1, 1, 1, 0, 1, 0, 0, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1,
                    0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1,
                    1, 1, 0, 1, 0, 0, 0, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0,
                    0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
                    0, 1, 0, 0, 0, 0, 0, 1, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1,
                    0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1,
                    0, 0, 0, 0, 0, 1, 1, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1,
                    0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0,
                    0, 0, 0, 1, 1, 0, 1, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1,
                    0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0,
                    0, 1, 1, 0, 1, 0, 1, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0,
                    0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1,
                    1, 0, 1, 0, 1, 0, 1, 0, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1,
                    0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0,
                    1, 0, 1, 0, 1, 0, 1, 1, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1,
                    1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0,
                    1, 0, 1, 0, 1, 1, 0, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0,
                    1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0,
                    1, 0, 1, 1, 0, 0, 0, 1, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0,
                    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0,
                    1, 1, 0, 0, 0, 1, 0, 1, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0,
                    0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1,
                    0, 0, 0, 1, 0, 1, 0, 1, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0,
                    1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0,
                    1, 1, 0, 0, 1, 0, 0, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0,
                    0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1,
                    0, 0, 1, 0, 0, 0, 0, 1, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0,
                    0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0,
                    1, 0, 0, 0, 0, 1, 0, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1,
                    1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                    0, 0, 0, 1, 0, 0, 1, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0,
                    0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0,
                    0, 1, 0, 0, 1, 0, 0, 1, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1,
                    1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0,
                    1, 0, 1, 0, 1, 1, 1, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0,
                    0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0,
                    1, 0, 1, 1, 1, 0, 1, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0,
                    1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0,
                    1, 1, 1, 0, 1, 0, 0, 1, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0,
                    1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1,
                    1, 0, 1, 0, 0, 1, 0, 0, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1,
                    1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0,
                    1, 0, 0, 1, 0, 0, 1, 1, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1,
                    0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0,
                    0, 1, 0, 0, 1, 1, 1, 1, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1,
                    1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1,
                    0, 0, 1, 1, 1, 1, 1, 1, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1,
                    1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0,
                    1, 1, 1, 1, 1, 1, 1, 0, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
                    1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1,
                    1, 1, 1, 1, 1, 0, 1, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0,
                    0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1,
                    1, 1, 1, 0, 1, 0, 0, 1, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0,
                    1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1,
                    1, 0, 1, 0, 0, 1, 0, 0, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1,
                    0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0,
                    1, 0, 0, 1, 0, 0, 1, 1, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1,
                    0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0,
                    0, 1, 0, 0, 1, 1, 0, 1, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0,
                    1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1,
                    0, 0, 1, 1, 0, 1, 0, 1, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0,
                    1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0,
                    1, 1, 0, 1, 0, 1, 1, 1, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0,
                    1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1,
                    0, 1, 0, 1, 1, 1, 1, 1, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1,
                    1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1,
                    0, 1, 1, 1, 1, 1, 1, 1, 1, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                    0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0,
                    1, 1, 0, 1, 0, 0, 1, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1,
                    0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1,
                    0, 1, 0, 0, 1, 0, 0, 1, 0, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0,
                    1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1,
                    0, 0, 1, 0, 0, 1, 0, 0, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1,
                    1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0,
                    1, 0, 0, 1, 0, 0, 0, 1, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0,
                    0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0,
                    0, 1, 0, 0, 0, 1, 0, 1, 0, 1
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1,
                    0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0,
                    0, 0, 1, 1, 1, 0, 1, 1, 1, 0
                ],
                bitvec::bitbox![u8, bitvec::order::Msb0;
                    1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0,
                    0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0,
                    1, 1, 1, 0, 1, 1, 1, 0, 0, 0
                ],
            ]
        );

        Ok(())
    }
}
