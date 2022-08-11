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

    #[test]
    fn read() -> error::Result<()> {
        let mut input_file: std::io::BufReader<&[u8]> = std::io::BufReader::new(&[
            0, 0, 0, 0, 0, 0, 0, 4, 109, 97, 120, 0, 0, 0, 0, 0, 0, 0, 0, 255, 100, 97, 116, 97,
            95, 115, 105, 122, 101, 0, 0, 0, 0, 0, 0, 0, 0, 1, 111, 114, 100, 101, 114, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 107, 0, 0, 0, 0, 0, 0, 0, 0, 15,
        ]);

        let values = Values::read(&mut input_file)?;

        let mut keys: Vec<String> = values.keys().cloned().collect();
        let mut values: Vec<u64> = values.values().cloned().collect();

        keys.sort();
        values.sort();

        assert_eq!(
            keys,
            [
                "data_size".to_string(),
                "k".to_string(),
                "max".to_string(),
                "order".to_string(),
            ]
        );

        assert_eq!(values, [0, 1, 15, 255]);

        Ok(())
    }

    #[test]
    fn write() -> error::Result<()> {
        let mut values = Values::with_capacity(20);

        values.insert("k".to_string(), 15);
        values.insert("order".to_string(), 0);
        values.insert("max".to_string(), 255);
        values.insert("data_size".to_string(), 1);

        let mut outer = Vec::new();
        values.write(&mut outer)?;

        assert_eq!(
            outer,
            &[
                0, 0, 0, 0, 0, 0, 0, 4, 109, 97, 120, 0, 0, 0, 0, 0, 0, 0, 0, 255, 100, 97, 116,
                97, 95, 115, 105, 122, 101, 0, 0, 0, 0, 0, 0, 0, 0, 1, 111, 114, 100, 101, 114, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 107, 0, 0, 0, 0, 0, 0, 0, 0, 15
            ]
        );

        Ok(())
    }
}