//! # kff-rust-api
//!
//! A rust API to read and write Kmer File Format file. [KFF specification](https://github.com/Kmer-File-Format/kff-reference)
//!
//! ## Read
//!
//! ```rust
//! use kff::error::Error;
//! use kff::Writer;
//! use kff::Reader;
//! use kff::error::LocalResult;
//! use crate::kff::seq2bits::Bits2Nuc;
//!
//! # fn main() -> Result<(), Error> {
//!
//! let mut output = vec![0u8; 0];
//! let buffer = std::io::Cursor::new(&mut output);
//!
//! // Define 2 bits encoding and writer
//! let encoding = 0b00011011; // 00 -> A, 01 -> C, 10 -> T, 11 -> G
//! let mut writer = Writer::new(buffer, 0b00011011, b"").unwrap();
//!
//! // Define required globale kff variables
//! writer.variables().insert("k".to_string(), 5);
//! writer.variables().insert("max".to_string(), 5);
//! writer.variables().insert("data_size".to_string(), 1);
//!
//! // Write these variables, before use variables you must call write_variables
//! writer.write_variables()?;
//!
//! // Write a Raw section
//! writer.write_raw_seq_section(&[b"GCGGGGATC"], &[&[1u8, 2, 3, 4, 5]])?;
//!
//! // Add a new variables and write a new Variable Section
//! writer.variables().insert("m".to_string(), 4);
//! writer.write_variables().unwrap();
//!
//! // Write a Minimizer section, minimizer -> CTT
//! writer.write_minimizer_seq_section(b"CCTT", &[2], &[b"AGCTG"], &[&[6, 7, 8, 9, 10]])?;
//!
//! // Create Reader
//! let mut input = std::io::Cursor::new(output);
//! let mut reader = Reader::new(&mut input).unwrap();
//!
//! // Get a copy of information need to perform reverse encoding
//! let rev_encoding = reader.rev_encoding();
//!
//! while let Ok(section) = reader.next_section() {
//!     let mut it = section.into_iter();
//!     while let Some(Ok(kmer)) = it.next() {
//!        println!("{} {:?}", std::str::from_utf8(&kmer.seq().into_nuc(rev_encoding)).map_local()?, kmer.data());
//!     }
//! }
//!
//! # Ok(())
//! # }
//! ```

pub mod data;
pub mod error;
pub mod kff;
pub mod kmer;
pub mod minimizer;
pub mod raw;
pub mod seq2bits;
pub mod utils;
pub mod variables;

pub use crate::error::Result;

pub use crate::kff::*;
