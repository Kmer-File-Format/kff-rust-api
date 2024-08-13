//! Read a kff file show index or show content of a section

/* std use */

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
    #[clap(short = 'v', long = "verbosity",  action = clap::ArgAction::Count)]
    pub verbosity: usize,

    /// Timestamp (sec, ms, ns, none)
    #[clap(short = 'T', long = "timestamp")]
    pub ts: Option<stderrlog::Timestamp>,

    /// Kff input file
    #[clap(short = 'i', long = "input")]
    pub input: std::path::PathBuf,

    /// Flag to show index content
    #[clap(short = 's', long = "show-index")]
    pub show: bool,

    /// Get target section
    #[clap(short = 't', long = "target")]
    pub target: Option<usize>,
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

    let mut kff = kff::Kff::with_index(params.input)?;

    if params.show {
        if let Some(index) = kff.index() {
            println!("index, type, offset");
            for (i, (t, p)) in index.pair().iter().enumerate() {
                println!("{},{},{}", i, *t as char, p);
            }
        } else {
            log::error!("This kff file didn't containts index");
            return Ok(());
        }
    } else if let Some(target) = params.target {
        let section = kff.kmer_of_section(target)?;
        for kmer in section {
            println!("{}", String::from_utf8(kmer.seq(*kff.header().encoding()))?);
        }
    } else {
        log::error!("You must set option `show` or `target`")
    }

    Ok(())
}
