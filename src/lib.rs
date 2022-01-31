//! KFF rust crate to handle, read and write [KFF format](https://github.com/Kmer-File-Format/kff-reference)

/* std use */

/* crate use */

/* project use */

/* mod declaration */
mod error;
mod header;
mod io;

pub type Endianess = byteorder::BigEndian;
