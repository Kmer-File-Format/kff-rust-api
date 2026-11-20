//! The examples give the expected results
//!
//! The expected kmers of tests/data/*.kmers.txt are the output of `kff-tools outstr`, which reads
//! the files with the C++ API.

/* std use */
use std::path::PathBuf;
use std::process::Command;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// Run an example and return its standard output
fn run_example(name: &str, args: &[&str]) -> String {
    let output = Command::new(env!("CARGO"))
        .args(["run", "--quiet", "--example", name, "--"])
        .args(args)
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .expect("could not run cargo");
    assert!(
        output.status.success(),
        "example {} failed\nstdout: {}\nstderr: {}",
        name,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn kff2kmer() {
    for name in ["test", "index"] {
        let kff = data(&format!("{}.kff", name));
        let output = run_example("kff2kmer", &["-i", kff.to_str().unwrap()]);
        let expected = std::fs::read_to_string(data(&format!("{}.kmers.txt", name))).unwrap();
        assert_eq!(output, expected, "file {}.kff", name);
    }
}

#[test]
fn index_read_show_index() {
    let kff = data("index.kff");
    let output = run_example("index_read", &["-i", kff.to_str().unwrap(), "-s"]);
    let expected = std::fs::read_to_string(data("index.sections.txt")).unwrap();
    assert_eq!(output, expected);
}

#[test]
fn index_read_kmers_of_a_section() {
    let kff = data("index.kff");
    let output = run_example("index_read", &["-i", kff.to_str().unwrap(), "-t", "5"]);
    assert_eq!(
        output,
        "TATCCGGACTGATAC\nATCCGGACTGATACT\nTCCGGACTGATACTC\nCCGGACTGATACTCC\nCGGACTGATACTCCA\n"
    );
    // The kmers of the section are kmers of the file
    let all = std::fs::read_to_string(data("index.kmers.txt")).unwrap();
    let all: Vec<&str> = all.lines().collect();
    assert!(output.lines().all(|kmer| all.contains(&kmer)));
}

#[test]
fn kmer_count2kff() -> kff::error::Result<()> {
    let dir = tempfile::tempdir()?;
    let output = dir.path().join("counts.kff");
    let csv = data("kmer_count.csv");
    run_example(
        "kmer_count2kff",
        &["-i", csv.to_str().unwrap(), "-o", output.to_str().unwrap()],
    );

    let expected: Vec<(String, u8)> = std::fs::read_to_string(csv)?
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let (kmer, count) = line.split_once(',').unwrap();
            (kmer.to_string(), count.parse().unwrap())
        })
        .collect();

    let reader = kff::KmerReader::open(&output)?;
    assert!(reader.reader().uniqueness());
    assert_eq!(reader.reader().metadata(), b"producer: kmer_count2kff");
    let kmers = reader
        .kmers()
        .map(|kmer| kmer.map(|kmer| (kmer.to_string(), kmer.data[0])))
        .collect::<kff::error::Result<Vec<(String, u8)>>>()?;
    assert_eq!(kmers, expected);

    Ok(())
}
