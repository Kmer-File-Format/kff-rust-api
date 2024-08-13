//! Utils function to write KFF

/* std use */

/* crate use */

/* project use */
use crate::error;

/// Define trait containts utils function to write kff
pub trait KffWrite {
    /// Function that write all bytes
    fn write_bytes(&mut self, bytes: &[u8]) -> error::Result<()>;

    /// Function that write bytes plus a '\0' at end
    fn write_ascii(&mut self, ascii: &[u8]) -> error::Result<()> {
        self.write_bytes(ascii)?;
        self.write_bytes(b"\0")
    }

    /// Function that write one bit and convert it as bool
    fn write_bool(&mut self, value: &bool) -> error::Result<()> {
        self.write_bytes(&((*value as u8).to_be_bytes()))
    }

    /// Function that write u8
    fn write_u8(&mut self, value: &u8) -> error::Result<()> {
        self.write_bytes(&value.to_be_bytes())
    }

    /// Function that write u16
    fn write_u16(&mut self, value: &u16) -> error::Result<()> {
        self.write_bytes(&value.to_be_bytes())
    }

    /// Function that write u32
    fn write_u32(&mut self, value: &u32) -> error::Result<()> {
        self.write_bytes(&value.to_be_bytes())
    }

    /// Function that write u64
    fn write_u64(&mut self, value: &u64) -> error::Result<()> {
        self.write_bytes(&value.to_be_bytes())
    }

    /// Function that write i64
    fn write_i64(&mut self, value: &i64) -> error::Result<()> {
        self.write_bytes(&value.to_be_bytes())
    }
}

impl<T> KffWrite for T
where
    T: std::io::Write,
{
    fn write_bytes(&mut self, bytes: &[u8]) -> error::Result<()> {
        self.write_all(bytes)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn write_bytes() -> error::Result<()> {
        let mut writer = Vec::new();

        writer.write_bytes(b"Lorem")?;

        assert_eq!(writer, b"Lorem");

        Ok(())
    }

    #[test]
    fn write_ascii() -> error::Result<()> {
        let mut writer = Vec::new();

        writer.write_ascii(b"Lorem")?;

        assert_eq!(writer, &[b'L', b'o', b'r', b'e', b'm', 0]);

        Ok(())
    }

    #[test]
    fn write_bool() -> error::Result<()> {
        let mut writer = Vec::new();

        writer.write_bool(&true)?;
        writer.write_bool(&false)?;

        assert_eq!(writer, &[1, 0]);

        Ok(())
    }

    #[test]
    fn write_u8() -> error::Result<()> {
        let mut writer = Vec::new();

        writer.write_u8(&1)?;
        writer.write_u8(&128)?;
        writer.write_u8(&130)?;

        assert_eq!(writer, &[1, 128, 130]);

        Ok(())
    }

    #[test]
    fn write_u16() -> error::Result<()> {
        let mut writer = Vec::new();

        writer.write_u16(&250)?;
        writer.write_u16(&500)?;

        assert_eq!(writer, &[0, 250, 1, 244]);

        Ok(())
    }

    #[test]
    fn write_u32() -> error::Result<()> {
        let mut writer = Vec::new();

        writer.write_u32(&3511)?;
        writer.write_u32(&511110)?;

        assert_eq!(writer, &[0, 0, 13, 183, 0, 7, 204, 134]);

        Ok(())
    }

    #[test]
    fn write_u64() -> error::Result<()> {
        let mut writer = Vec::new();

        writer.write_u64(&35191823831)?;
        writer.write_u64(&444799335191823831)?;

        assert_eq!(
            writer,
            &[0, 0, 0, 8, 49, 152, 157, 215, 6, 44, 62, 163, 130, 112, 69, 215]
        );

        Ok(())
    }
}
