//! Parse, manage and write Values section

/* std use */

/* crate use */

/* project use */
use crate::error;

/// Struct to parse, manage and write Values section
pub type Values = rustc_hash::FxHashMap<String, u64>;

trait AbcValues: Sized {
    /// Build an empty Values
    fn new() -> Self;

    /// Build an empty Values but with a minimal capacity
    fn with_capacity(capacity: usize) -> Self;

    /// Build a Values from a readable key are overwrite
    fn read<R>(inner: &mut R) -> error::Result<Self>
    where
        R: std::io::Read + crate::KffRead;

    /// Write contents of Values in writables
    fn write<W>(&self, outer: &mut W) -> error::Result<()>
    where
        W: std::io::Write;
}

impl AbcValues for Values {
    fn new() -> Self {
        Values::default()
    }

    fn with_capacity(capacity: usize) -> Self {
        Values::with_capacity_and_hasher(
            capacity,
            core::hash::BuildHasherDefault::<rustc_hash::FxHasher>::default(),
        )
    }

    fn read<R>(inner: &mut R) -> error::Result<Self>
    where
        R: std::io::Read + crate::KffRead,
    {
        let nb_variable = inner.read_u64()?;
        let mut obj = Self::with_capacity(nb_variable as usize);

        for _ in 0..nb_variable {
            let key = String::from_utf8(inner.read_ascii()?)?;
            let value = inner.read_u64()?;

            obj.insert(key, value);
        }

        Ok(obj)
    }

    fn write<W>(&self, outer: &mut W) -> error::Result<()>
    where
        W: std::io::Write + crate::KffWrite,
    {
        outer.write_u64(&(self.len() as u64))?;
        for (key, value) in self.iter() {
            outer.write_ascii(key.as_bytes())?;
            outer.write_u64(value)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creation() -> error::Result<()> {
        let default = Values::new();
        assert_eq!(default.capacity(), 0);

        let capacity = Values::with_capacity(20);
        assert!(capacity.capacity() > 20);

        Ok(())
    }
}
