//! Representation of a KFF index

/* std use */

/* crate use */

/* project use */
use crate::error;

/// Struct to Read and Write Index section
#[derive(getset::Getters, getset::Setters, getset::MutGetters)]
#[getset(get = "pub", set = "pub", get_mut = "pub")]
pub struct Index {
    /// Vector of block type and relative position
    pair: Vec<(u8, i64)>,
    /// Position of the next index
    next_index: u64,
}

impl Index {
    /// Read an Index section
    pub fn read<R>(inner: &mut R) -> error::Result<Self>
    where
        R: std::io::Read + crate::KffRead,
    {
        let mut output = Vec::new();

        let nb_block = inner.read_u64()?;

        for _ in 0..nb_block {
            let section_type = inner.read_u8()?;
            let delta = inner.read_i64()?;

            output.push((section_type, delta));
        }

        Ok(Self {
            pair: output,
            next_index: inner.read_u64()?,
        })
    }

    /// Write an Index section
    pub fn write<W>(&self, outer: &mut W) -> error::Result<()>
    where
        W: std::io::Write + crate::KffWrite,
    {
        outer.write_u64(&(self.pair.len() as u64))?;

        for (section_type, delta) in &self.pair {
            outer.write_u8(section_type)?;
            outer.write_i64(delta)?;
        }

        outer.write_u64(&self.next_index)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read() -> error::Result<()> {
        let mut data: &[u8] = &[
            0, 0, 0, 0, 0, 0, 0, 3, // number of pair
            b'r', 0, 0, 0, 0, 0, 0, 55, 255, // Raw section
            b't', 255, 0, 0, 0, 0, 0, 0, 255, // a T section with value in past
            b'm', 0, 0, 0, 0, 0, 255, 0, 255, // Minimizer section
            0, 0, 0, 0, 0, 45, 33, 0, // Next index section
        ];

        let index = Index::read(&mut data)?;

        assert_eq!(
            index.pair(),
            &[(b'r', 14335), (b't', -72057594037927681), (b'm', 16711935)]
        );
        assert_eq!(index.next_index(), &2957568);

        Ok(())
    }

    #[test]
    fn write() -> error::Result<()> {
        let mut data: &[u8] = &[
            0, 0, 0, 0, 0, 0, 0, 3, // number of pair
            b'r', 0, 0, 0, 0, 0, 0, 55, 255, // Raw section
            b't', 255, 0, 0, 0, 0, 0, 0, 255, // a T section with value in past
            b'm', 0, 0, 0, 0, 0, 255, 0, 255, // Minimizer section
            0, 0, 0, 0, 0, 45, 33, 0, // Next index section
        ];

        let index = Index::read(&mut data)?;

        let mut output = Vec::new();
        index.write(&mut output)?;

        assert_eq!(
            output,
            &[
                0, 0, 0, 0, 0, 0, 0, 3, // number of pair
                b'r', 0, 0, 0, 0, 0, 0, 55, 255, // Raw section
                b't', 255, 0, 0, 0, 0, 0, 0, 255, // a T section with value in past
                b'm', 0, 0, 0, 0, 0, 255, 0, 255, // Minimizer section
                0, 0, 0, 0, 0, 45, 33, 0, // Next index section
            ]
        );

        Ok(())
    }
}
