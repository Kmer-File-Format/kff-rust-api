//! Some utils function to read and write in stream

/* std use */

/* crate use */

/* project use */
use crate::*;

pub fn read_fix<R, const S: usize>(input: &mut R) -> error::Result<[u8; S]>
where
    R: std::io::Read,
{
    let mut buffer = [0u8; S];

    input.read_exact(&mut buffer)?;

    Ok(buffer)
}

pub fn read_dynamic<R>(input: &mut R, size: usize) -> error::Result<Vec<u8>>
where
    R: std::io::Read,
{
    let mut buffer = vec![0u8; size];

    input.read_exact(&mut buffer)?;

    Ok(buffer)
}
