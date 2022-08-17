//! Kmer File Format Rust parser

#![warn(missing_docs)]

/* std use */

/* crate use */

use std::io::Read;

use clap::Parser as _;

/* project use */
use kff::error;

/// Example: Kmer File Format Rust parser
#[derive(clap::Parser, std::fmt::Debug)]
#[clap(
    name = "kff",
    version = "0.1",
    author = "Pierre Marijon <pierre@marijon.fr>"
)]
pub struct Command {
    /// Silence all output
    #[clap(short = 'q', long = "quiet")]
    pub quiet: bool,

    /// Verbose mode (-v, -vv, -vvv, etc)
    #[clap(short = 'v', long = "verbosity", parse(from_occurrences))]
    pub verbosity: usize,

    /// Timestamp (sec, ms, ns, none)
    #[clap(short = 'T', long = "timestamp")]
    pub ts: Option<stderrlog::Timestamp>,

    /// Kff input file
    #[clap(short = 'i', long = "input")]
    pub input: std::path::PathBuf,
}

fn main() -> error::Result<()> {
    // parse cli
    let params = Command::parse();

    // Setup logger
    stderrlog::new()
        .quiet(params.quiet)
        .verbosity(params.verbosity)
        .timestamp(params.ts.unwrap_or(stderrlog::Timestamp::Off))
        .init()
        .unwrap();

    log::trace!("Open file");
    let file = kff::read::Kff::<std::io::BufReader<std::fs::File>>::open(params.input)?;

    let mut iter = file.kmers();
    while let Some(Ok((kmer, _data))) = iter.next() {
        println!("{}", kmer);
    }

    Ok(())
}
