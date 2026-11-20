//! Reader and writer of KFF files (Kmer File Format).
//!
//! The format is described in the [KFF reference](https://github.com/Kmer-File-Format/kff-reference).
//! This crate offers the features of the [C++ API](https://github.com/Kmer-File-Format/kff-cpp-api)
//! and reads and writes the same bytes:
//!
//! - High level reading: [KmerReader] enumerates the kmers (or the blocks) of a file with their
//!   data, without dealing with the sections (`Kff_reader` in the C++ API).
//! - Low level reading: [Reader] reads the header, the footer, the index and each section
//!   (`Kff_file` in "r" mode and the `Section_*` classes).
//! - Low level writing: [Writer] writes the header and each section, and by default an index of the
//!   sections with a footer (`Kff_file` in "w" mode and the `Section_*` classes).
//!
//! The sequences are given as right aligned 2-bit nucleotides (see [encoding]): [Encoding]
//! converts them from and to text.
//!
//! ```no_run
//! # fn main() -> kff::error::Result<()> {
//! // Write a file
//! let mut writer = kff::Writer::create("file.kff")?;
//! let encoding = writer.encoding();
//! let mut values = writer.values_section()?;
//! values.write_var("k", 5)?;
//! values.write_var("max", 10)?;
//! values.write_var("data_size", 1)?;
//! values.close()?;
//! let mut raw = writer.raw_section()?;
//! raw.write_compacted_sequence(&encoding.encode(b"ACGTACG")?, 7, &[1, 2, 3])?;
//! raw.close()?;
//! writer.close()?;
//!
//! // Read its kmers
//! for kmer in kff::KmerReader::open("file.kff")?.kmers() {
//!     let kmer = kmer?;
//!     println!("{} {:?}", kmer, kmer.data);
//! }
//! # Ok(())
//! # }
//! ```

#![warn(missing_docs)]

/* mod declaration */
mod block;
pub mod encoding;
pub mod error;
mod header;
mod kmer_reader;
pub mod nucleotides;
pub mod reader;
pub mod section;
pub mod writer;

/* pub use */
pub use block::Block;
pub use encoding::Encoding;
pub use header::Header;
pub use kmer_reader::{Kmer, KmerReader, KmerView, Kmers};
pub use reader::{
    BlockReader, BlockSectionReader, MinimizerSectionReader, RawSectionReader, Reader, Section,
};
pub use section::{IndexSection, SectionType, Values};
pub use writer::{
    IndexSectionWriter, MinimizerSectionWriter, RawSectionWriter, ValuesSectionWriter, Writer,
};

/// Major version of the KFF format written in the files. The files of a higher version are not
/// read.
pub const KFF_FORMAT_VERSION_MAJOR: u8 = 1;

/// Minor version of the KFF format written in the files
pub const KFF_FORMAT_VERSION_MINOR: u8 = 0;
