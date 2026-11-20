//! Print the index of a KFF file, or the kmers of one of its indexed sections

/* crate use */
use clap::Parser as _;

/* project use */
use kff::{error, BlockReader, SectionType};

/// Print the index of a KFF file or the kmers of one of its sections
#[derive(clap::Parser, std::fmt::Debug)]
#[clap(name = "index_read")]
pub struct Command {
    /// KFF input file
    #[clap(short = 'i', long = "input")]
    pub input: std::path::PathBuf,

    /// Print the sections registered in the index: number, type and position in the file
    #[clap(short = 's', long = "show-index")]
    pub show: bool,

    /// Print the kmers of the section of this number
    #[clap(short = 't', long = "target")]
    pub target: Option<usize>,
}

fn main() -> error::Result<()> {
    let params = Command::parse();

    let mut reader = kff::Reader::open(params.input)?;
    if !reader.is_indexed() {
        eprintln!("This KFF file has no index");
        std::process::exit(1);
    }
    let sections: Vec<(u64, SectionType)> = reader.indexed_sections().into_iter().collect();

    if params.show {
        println!("index, type, offset");
        for (i, (position, section_type)) in sections.iter().enumerate() {
            println!("{},{},{}", i, section_type.as_byte() as char, position);
        }
    } else if let Some(target) = params.target {
        let (position, section_type) = match sections.get(target) {
            Some(section) => *section,
            None => {
                eprintln!("The index contains {} sections", sections.len());
                std::process::exit(1);
            }
        };
        if !section_type.has_kmers() {
            eprintln!("Section {} is not a raw or minimizer section", target);
            std::process::exit(1);
        }

        // The variables of the section are in the last values section before it
        let values = sections[..target]
            .iter()
            .rev()
            .find(|(_, section_type)| *section_type == SectionType::Values);
        if let Some((values_position, _)) = values {
            reader.jump_to(*values_position)?;
            reader.read_values()?;
        }

        reader.jump_to(position)?;
        let encoding = reader.encoding();
        let mut section = reader.open_block_section()?.unwrap();
        while let Some(block) = section.read_block()? {
            for (kmer, _) in block.kmers() {
                println!(
                    "{}",
                    String::from_utf8_lossy(&encoding.decode(&kmer, block.k as usize))
                );
            }
        }
    } else {
        eprintln!("Use the option --show-index or --target");
        std::process::exit(1);
    }

    Ok(())
}
