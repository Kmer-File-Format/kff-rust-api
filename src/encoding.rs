//! Module to manage encoding and decoding

/* std use */

/* crate use */

/* project use */
use crate::*;

/// Struct to store validate and use encoding
pub struct Encoding {
    value: u8,
}

impl Encoding {
    pub fn new(encoding: u8) -> error::Result<Self> {
        Ok(Self {
            value: valid_encoding(encoding)?,
        })
    }
}

pub(crate) fn valid_encoding(encoding: u8) -> error::Result<u8> {
    let a = encoding >> 6;
    let c = (encoding >> 4) & 0b11;
    let t = (encoding >> 2) & 0b11;
    let g = encoding & 0b11;

    if a != c && a != t && a != g && c != t && t != g {
        Ok(encoding)
    } else {
        Err(error::Error::EncodingNotValid)
    }
}
