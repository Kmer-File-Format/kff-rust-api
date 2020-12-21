pub mod error;
pub mod header;
pub mod minimizer;
pub mod raw;
pub mod utils;
pub mod variables;
pub mod data;
pub mod kff;
pub mod kmer;
pub mod metadata;

pub use crate::error::Result;

pub use crate::kff::*;
