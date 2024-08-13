//! Build a global index of a kff file

/* std use */

/* crate use */

/* project use */
use crate::error;
use crate::section;
use crate::KffRead;

/// Struct that manage and build a global index of a Kff file
#[derive(getset::Getters, getset::Setters, getset::MutGetters, std::default::Default)]
#[getset(get = "pub")]
pub struct GlobalIndex {
    /// Pair of section type and position from begin of file
    pair: Vec<(u8, u64)>,
}

impl GlobalIndex {
    /// Create a GlobalIndex by scan all index in file
    ///
    /// At the end of scan position in file is the start of the first_index
    pub fn new<R>(inner: &mut R, first_index: u64) -> error::Result<Self>
    where
        R: std::io::Read + std::io::Seek + KffRead,
    {
        let mut pair = Vec::new();

        inner.seek(std::io::SeekFrom::Start(first_index))?;

        let section_type = inner.read_u8()?;
        if section_type != b'i' {
            return Err(error::Kff::NotAnIndex.into());
        }

        loop {
            let local_index = section::Index::read(inner)?;

            let relative_to = inner.stream_position()?;

            pair.extend(
                local_index
                    .pair()
                    .iter()
                    .map(|(t, pos)| (*t, (relative_to as i64 + *pos) as u64)),
            );

            if local_index.next_index() == &0 {
                break;
            } else {
                inner.seek(std::io::SeekFrom::Current(
                    (*local_index.next_index()) as i64,
                ))?;
            }
        }

        pair.sort_unstable_by_key(|k| k.1);

        Ok(Self { pair })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create() -> error::Result<()> {
        let mut file = std::io::Cursor::new(vec![
            b'i', // Index start
            0, 0, 0, 0, 0, 0, 0, 3, // Number of value indexed
            b'r', 0, 0, 0, 0, 0, 0, 0, 1, // Raw section
            b'i', 0, 0, 0, 0, 0, 0, 0, 0, // Index refere to them self
            b'm', 0, 0, 0, 0, 0, 0, 0, 3, // Minimizer section
            0, 0, 0, 0, 0, 0, 0, 1,    // Next index section
            b'i', // Index start
            0, 0, 0, 0, 0, 0, 0, 3, // Number of value indexed
            b't', 255, 255, 255, 255, 255, 255, 255, 253, // a T section with value in past
            b'r', 0, 0, 0, 0, 0, 0, 0, 1, // Raw section
            b'm', 0, 0, 0, 0, 0, 0, 0, 3, // Minimizer section
            0, 0, 0, 0, 0, 0, 0, 9, // Next index section
            0, 0, 0, 0, 0, 0, 0, 0,    // empty part
            b'i', // Index start
            0, 0, 0, 0, 0, 0, 0, 3, // Number of value indexed
            b'm', 0, 0, 0, 0, 0, 0, 0, 3, // Minimizer section
            b'r', 0, 0, 0, 0, 0, 0, 0, 1, // Raw section
            b't', 255, 255, 255, 255, 255, 255, 255, 253, // a T section with value in past
            0, 0, 0, 0, 0, 0, 0, 0, // Next index section
        ]);

        let index = GlobalIndex::new(&mut file, 0)?;

        assert_eq!(
            index.pair(),
            &vec![
                (105, 44),
                (114, 45),
                (109, 47),
                (116, 85),
                (114, 89),
                (109, 91),
                (116, 137),
                (114, 141),
                (109, 143)
            ]
        );

        Ok(())
    }
}
