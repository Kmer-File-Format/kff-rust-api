//! Compatibility with the C++ API (kff-cpp-api)
//!
//! The files of tests/data/cpp are written by the C++ API (tests/cpp/generate_fixtures.cpp) with
//! the list of their kmers read by its high level reader. This crate must read the same kmers, and
//! write the same bytes when the sections are written again with the same calls.

/* std use */
use std::io::Cursor;
use std::path::PathBuf;

/* project use */
use kff::{BlockReader, KmerReader, Reader, SectionType, Writer};

const FIXTURES: &[&str] = &[
    "empty",
    "raw",
    "raw_sizes",
    "minimizer",
    "mixed",
    "chained_index",
];

fn path(name: &str, extension: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cpp")
        .join(format!("{}.{}", name, extension))
}

fn expected_kmers(name: &str) -> Vec<String> {
    std::fs::read_to_string(path(name, "txt"))
        .unwrap()
        .lines()
        .map(|line| line.to_string())
        .collect()
}

fn hex(data: &[u8]) -> String {
    data.iter().map(|byte| format!("{:02x}", byte)).collect()
}

#[test]
fn next_kmer_reads_the_kmers_of_the_cpp_reader() -> kff::error::Result<()> {
    for name in FIXTURES {
        let mut reader = KmerReader::open(path(name, "kff"))?;
        let encoding = reader.encoding();
        let mut kmers = Vec::new();
        while let Some(kmer) = reader.next_kmer()? {
            assert_eq!(kmer.sequence.len(), (kmer.k as usize + 3) / 4);
            kmers.push(format!(
                "{}\t{}",
                String::from_utf8(encoding.decode(kmer.sequence, kmer.k as usize)).unwrap(),
                hex(kmer.data)
            ));
            assert_eq!(kmer.data.len() as u64, reader.data_size());
        }
        assert!(!reader.has_next()?);
        assert_eq!(kmers, expected_kmers(name), "file {}", name);
    }
    Ok(())
}

#[test]
fn kmers_iterator_and_blocks_read_the_same_kmers() -> kff::error::Result<()> {
    for name in FIXTURES {
        let kmers = KmerReader::open(path(name, "kff"))?
            .kmers()
            .map(|kmer| kmer.map(|kmer| format!("{}\t{}", kmer, hex(&kmer.data))))
            .collect::<kff::error::Result<Vec<String>>>()?;
        assert_eq!(kmers, expected_kmers(name), "file {}", name);

        let mut reader = KmerReader::open(path(name, "kff"))?;
        let encoding = reader.encoding();
        let mut kmers = Vec::new();
        while let Some(block) = reader.next_block()? {
            for (kmer, data) in block.kmers() {
                kmers.push(format!(
                    "{}\t{}",
                    String::from_utf8(encoding.decode(&kmer, block.k as usize)).unwrap(),
                    hex(data)
                ));
            }
        }
        assert_eq!(kmers, expected_kmers(name), "file {}", name);
    }
    Ok(())
}

/// Write again every section of a file with this crate. The index and the footer written by the
/// C++ API are not copied: the writer writes them on closing.
///
/// With `full_minimizer`, the blocks of the minimizer sections are read with their minimizer and
/// written with `write_compacted_sequence` (the writer removes the minimizer), otherwise they are
/// copied without their minimizer.
fn rewrite(name: &str, full_minimizer: bool) -> kff::error::Result<Vec<u8>> {
    let mut reader = Reader::open(path(name, "kff"))?;
    let mut writer = Writer::new(Cursor::new(Vec::new()))?;

    let auto_index = reader.footer().is_some();
    writer.set_indexation(auto_index);
    writer.set_encoding(reader.encoding())?;
    writer.set_uniqueness(reader.uniqueness())?;
    writer.set_canonicity(reader.canonicity())?;
    writer.write_metadata(reader.metadata())?;
    let footer_position = reader.footer_position();

    while let Some(section_type) = reader.read_section_type()? {
        match section_type {
            SectionType::Values => {
                let position = reader.tell();
                let values = reader.read_values()?;
                if Some(position) != footer_position {
                    writer.write_values(&values)?;
                }
            }
            SectionType::Index => {
                let index = reader.read_index()?;
                if !auto_index {
                    let mut section = writer.index_section()?;
                    for (position, section_type) in &index.entries {
                        section.register_section(*section_type, *position);
                    }
                    section.set_next_index(index.next_index);
                    section.close()?;
                }
            }
            SectionType::Raw => reader.open_raw()?.copy_to(&mut writer)?,
            SectionType::Minimizer if !full_minimizer => {
                reader.open_minimizer()?.copy_to(&mut writer)?
            }
            SectionType::Minimizer => {
                let mut section = reader.open_minimizer()?;
                let mut out = writer.minimizer_section()?;
                out.write_minimizer(section.minimizer())?;
                let k = section.k();
                let mut seq = vec![0; section.max_sequence_bytes()];
                let mut data = vec![0; section.max_data_bytes()];
                while section.remaining_blocks() > 0 {
                    let (nb_kmers, mini_pos) =
                        section.read_compacted_sequence_without_mini(&mut seq, &mut data)?;
                    section.add_minimizer(nb_kmers, &mut seq, mini_pos)?;
                    out.write_compacted_sequence(&seq, nb_kmers + k - 1, mini_pos, &data)?;
                }
                out.close()?;
                section.close()?;
            }
            SectionType::Other(byte) => panic!("unknown section {}", byte),
        }
    }

    Ok(writer.close()?.into_inner())
}

#[test]
fn rewritten_files_are_the_files_of_the_cpp_api() -> kff::error::Result<()> {
    for name in FIXTURES {
        let expected = std::fs::read(path(name, "kff"))?;
        assert!(rewrite(name, false)? == expected, "file {}", name);
        assert!(
            rewrite(name, true)? == expected,
            "file {} (full minimizers)",
            name
        );
    }
    Ok(())
}

#[test]
fn header_footer_and_index_of_the_fixtures() -> kff::error::Result<()> {
    let reader = Reader::open(path("raw", "kff"))?;
    assert_eq!(reader.metadata(), b"kff-rust-api compatibility test");
    assert_eq!(reader.encoding(), kff::Encoding::default());
    assert!(!reader.uniqueness() && !reader.canonicity());
    let footer = reader.footer().unwrap();
    assert_eq!(footer.get("footer_size"), Some(&49));
    assert_eq!(reader.index().len(), 1);
    assert_eq!(
        reader.index()[0].beginning,
        *footer.get("first_index").unwrap()
    );
    let sections: Vec<SectionType> = reader.indexed_sections().values().copied().collect();
    assert_eq!(sections, vec![SectionType::Values, SectionType::Raw]);

    let reader = Reader::open(path("raw_sizes", "kff"))?;
    assert_eq!(reader.encoding().values(), [2, 0, 1, 3]);
    assert!(reader.uniqueness() && reader.canonicity());
    assert!(reader.footer().is_none());
    assert!(!reader.is_indexed());

    let reader = Reader::open(path("chained_index", "kff"))?;
    assert!(reader.footer().is_none());
    assert_eq!(reader.index().len(), 2);
    assert_eq!(reader.index()[0].beginning, 15);
    assert_eq!(reader.index()[0].next_index, 9);
    assert_eq!(reader.index()[1].beginning, 50);
    let positions: Vec<(u64, SectionType)> = reader.indexed_sections().into_iter().collect();
    assert_eq!(
        positions,
        vec![(41, SectionType::Values), (85, SectionType::Values)]
    );

    Ok(())
}
