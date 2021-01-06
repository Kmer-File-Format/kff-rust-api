//! Management of global variable

/* crate use */
use anyhow::Result;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

/* local use */
use crate::error;

/// Variables is a specialisation of HashMap
pub type Variables = std::collections::HashMap<String, u64>;

/// Trait to read global variable section
pub trait Reader {
    /// Read variable from input
    fn deserialize<R>(&mut self, input: &mut R) -> Result<usize>
    where
        R: std::io::Read;
}

/// Trait to write global variable section
pub trait Writer {
    /// Write variable in output
    fn serialize<W>(&self, output: &mut W) -> Result<usize>
    where
        W: std::io::Write;
}

/// Trait of variable needed by KFF 1.0
pub trait Variables1 {
    /// Get value of k
    fn k(&self) -> Result<u64>;

    /// Get value of m
    fn m(&self) -> Result<u64>;

    /// Get value of max
    fn max(&self) -> Result<u64>;

    /// Get value of data_size
    fn data_size(&self) -> Result<u64>;
}

impl Variables1 for Variables {
    fn k(&self) -> Result<u64> {
        self.get("k")
            .copied()
            .ok_or_else(|| error::Error::Variables(error::Variables::KMissing).into())
    }

    fn m(&self) -> Result<u64> {
        self.get("m")
            .copied()
            .ok_or_else(|| error::Error::Variables(error::Variables::MMissing).into())
    }

    fn max(&self) -> Result<u64> {
        self.get("max")
            .copied()
            .ok_or_else(|| error::Error::Variables(error::Variables::MaxMissing).into())
    }

    fn data_size(&self) -> Result<u64> {
        self.get("data_size")
            .copied()
            .ok_or_else(|| error::Error::Variables(error::Variables::DataSizeMissing).into())
    }
}

impl Reader for Variables {
    fn deserialize<R>(&mut self, input: &mut R) -> Result<usize>
    where
        R: std::io::Read,
    {
        let nb_variables = input.read_u64::<LittleEndian>()?;
        let mut name = Vec::new();
        let mut char;

        let mut nb_bytes = 1;

        for _ in 0..nb_variables {
            char = input.read_u8()?;
            nb_bytes += 1;
            while char != 0 {
                name.push(char);
                char = input.read_u8()?;
                nb_bytes += 1;

                if char == 0 {
                    break;
                }
            }

            self.insert(
                String::from_utf8(name.clone())?,
                input.read_u64::<LittleEndian>()?,
            );
            nb_bytes += 8;
            name.clear();
        }

        Ok(nb_bytes)
    }
}

impl Writer for Variables {
    fn serialize<W>(&self, output: &mut W) -> Result<usize>
    where
        W: std::io::Write,
    {
        output.write_u64::<LittleEndian>(self.len() as u64)?;

        let mut nb_bytes = 8;
        for (key, value) in self.iter() {
            output.write_all(key.as_bytes())?;
            output.write_u8(0)?;
            output.write_u64::<LittleEndian>(*value)?;

            nb_bytes += key.len() + 1 + 8;
        }

        Ok(nb_bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_write_read_access() {
        let mut variables = Variables::new();

        variables.insert("k".to_string(), 15);
        variables.insert("max".to_string(), 255);
        variables.insert("data_size".to_string(), 0);

        let mut output = Vec::new();

        variables.serialize(&mut output).unwrap();

        let mut var2 = Variables::new();

        let mut input: &[u8] = &output;
        var2.deserialize(&mut input).unwrap();

        assert_eq!(var2.get("k"), Some(&15u64));
        assert_eq!(var2.get("k"), variables.get("k"));

        assert_eq!(var2.get("max"), Some(&255u64));
        assert_eq!(var2.get("max"), variables.get("max"));

        assert_eq!(var2.get("data_size"), Some(&0u64));
        assert_eq!(var2.get("data_size"), variables.get("data_size"));
    }

    #[test]
    fn over_rigth() {
        let mut variables = Variables::new();

        variables.insert("k".to_string(), 15);

        assert_eq!(variables.get("k"), Some(&15u64));

        variables.insert("k".to_string(), 12);

        assert_eq!(variables.get("k"), Some(&12u64));
    }
}
