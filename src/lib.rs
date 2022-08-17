//! Kmer File Format Rust parser and writer

#![warn(missing_docs)]

/* std use */

/* crate use */

/* project use */

/* mod declaration */

pub mod error;
pub mod iterator;
pub mod kff;
pub mod section;
pub mod utils;

pub use self::kff::Kff;
pub use iterator::KmerIterator;
pub use utils::*;
