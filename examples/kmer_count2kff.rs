//! Read a kff file and write in stdout kmer

/* std use */
use std::io::BufRead as _;

/* crate use */
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
    #[clap(short = 'v', long = "verbosity", action = clap::ArgAction::Count)]
    pub verbosity: u8,

    /// Timestamp (sec, ms, ns, none)
    #[clap(short = 'T', long = "timestamp")]
    pub ts: Option<stderrlog::Timestamp>,

    /// Kff input file
    #[clap(short = 'i', long = "input")]
    pub input: std::path::PathBuf,

    /// Output
    #[clap(short = 'o', long = "output")]
    pub output: std::path::PathBuf,
}

fn main() -> error::Result<()> {
    // parse cli
    let params = Command::parse();

    // Setup logger
    stderrlog::new()
        .quiet(params.quiet)
        .verbosity(params.verbosity as usize)
        .timestamp(params.ts.unwrap_or(stderrlog::Timestamp::Off))
        .init()
        .unwrap();

    log::info!("Open input file");
    let mut kmer_size = 0;
    let input = std::fs::File::open(params.input).map(std::io::BufReader::new)?;
    let mut kmers = vec![];
    for l in input.lines() {
        let line = l?;
        if line.is_empty() {
            break;
        }

        let values: Vec<&str> = line.split(',').collect();
        let mut seq = bitvec::vec::BitVec::<u8, bitvec::order::Msb0>::new();

        kmer_size = values[0].len();

        for nuc in values[0].bytes() {
            let bits = (nuc >> 1) & 0b11;
            seq.push((bits >> 1) != 0);
            seq.push((bits & 0b1) != 0);
        }

        kmers.push(kff::section::Block::new(
            kmer_size as u64,
            1,
            kff::Kmer::new(
                seq.into_boxed_bitslice(),
                vec![values[1].parse::<u8>().unwrap()],
            ),
            0,
        ));
    }
    log::info!("End of input file reading");

    log::info!("Start create kff file");
    let header = kff::section::Header::new(
        1,
        0,
        0b00011110,
        true,
        true,
        b"producer: kmer_count2kff".to_vec(),
    )?;

    let mut kff = kff::Kff::create(params.output, header)?;
    let mut values = kff::section::Values::default();
    values.insert("k".to_string(), kmer_size as u64);
    values.insert("ordered".to_string(), true as u64);
    values.insert("max".to_string(), 200);
    values.insert("data_size".to_string(), 1);

    kff.write_values(values.clone())?;
    kff.write_raw(kff::section::Raw::new(&values)?, &kmers)?;
    kff.finalize()?;
    log::info!("End create kff file");

    Ok(())
}
