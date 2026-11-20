//! Print the kmers of a KFF file, one per line, with their data in hexadecimal

/* crate use */
use clap::Parser as _;

/* project use */
use kff::error;

/// Print the kmers of a KFF file, one per line
#[derive(clap::Parser, std::fmt::Debug)]
#[clap(name = "kff2kmer")]
pub struct Command {
    /// KFF input file
    #[clap(short = 'i', long = "input")]
    pub input: std::path::PathBuf,

    /// Print the data of the kmers in hexadecimal, after a tabulation
    #[clap(short = 'd', long = "data")]
    pub data: bool,
}

fn main() -> error::Result<()> {
    let params = Command::parse();

    let mut reader = kff::KmerReader::open(params.input)?;
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    while let Some(kmer) = reader.next_kmer()? {
        use std::io::Write as _;
        out.write_all(&kmer.to_ascii())?;
        if params.data {
            out.write_all(b"\t")?;
            for byte in kmer.data {
                write!(out, "{:02x}", byte)?;
            }
        }
        out.write_all(b"\n")?;
    }

    Ok(())
}
