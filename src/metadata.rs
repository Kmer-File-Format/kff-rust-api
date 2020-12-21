/* local use */
use crate::error::Result;

pub trait Metadata {
    fn deserialize<R>(&mut self, input: &mut R) -> Result<usize>
    where
        Self: Sized,
        R: std::io::Read;

    fn serialize<W: std::io::Write>(&self, output: &mut W) -> Result<usize>;
}

