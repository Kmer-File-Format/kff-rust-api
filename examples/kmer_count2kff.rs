//! Convert a CSV file of kmer counts ("kmer,count" lines) into a KFF file
//!
//! Each kmer is written in its own block of a raw section, with its count on one byte.

/* std use */
use std::io::BufRead as _;

/* crate use */
use clap::Parser as _;

/* project use */
use kff::error;

/// Convert a CSV file of kmer counts into a KFF file
#[derive(clap::Parser, std::fmt::Debug)]
#[clap(name = "kmer_count2kff")]
pub struct Command {
    /// CSV input file, one "kmer,count" per line (count lower than 256)
    #[clap(short = 'i', long = "input")]
    pub input: std::path::PathBuf,

    /// KFF output file
    #[clap(short = 'o', long = "output")]
    pub output: std::path::PathBuf,
}

fn main() -> error::Result<()> {
    let params = Command::parse();

    let input = std::fs::File::open(params.input).map(std::io::BufReader::new)?;
    let mut kmers = Vec::new();
    for line in input.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let (kmer, count) = line.split_once(',').expect("lines must be \"kmer,count\"");
        let count: u8 = count.trim().parse().expect("counts must be lower than 256");
        kmers.push((kmer.as_bytes().to_vec(), count));
    }
    let k = kmers
        .first()
        .map(|(kmer, _)| kmer.len() as u64)
        .unwrap_or(1);

    let mut writer = kff::Writer::create(params.output)?;
    writer.set_uniqueness(true)?;
    writer.write_metadata(b"producer: kmer_count2kff")?;
    let encoding = writer.encoding();

    let mut values = writer.values_section()?;
    values.write_var("k", k)?;
    values.write_var("max", 1)?;
    values.write_var("data_size", 1)?;
    values.write_var("ordered", 0)?;
    values.close()?;

    let mut raw = writer.raw_section()?;
    for (kmer, count) in kmers {
        raw.write_compacted_sequence(&encoding.encode(&kmer)?, kmer.len() as u64, &[count])?;
    }
    raw.close()?;

    writer.close()?;
    Ok(())
}
