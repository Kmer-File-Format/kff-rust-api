//! Check example produce the good output

/* std use */
use std::io::Read;
use std::process::{Command, Stdio};

/* crate use */
use tempfile;

#[test]
fn kff2kmers() -> kff::error::Result<()> {
    let tmp_file = tempfile::NamedTempFile::new().unwrap();

    let args = vec![
        "run",
        "--example",
        "kmer_count2kff",
        "--",
        "-i",
        "tests/data/kmer_count.csv",
        "-o",
        tmp_file.path().to_str().unwrap(),
    ];

    let mut child = Command::new("cargo")
        .args(&args)
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Couldn't create cargo example child");

    if !child.wait().expect("Error durring yacrd run").success() {
        let mut stdout = String::new();
        let mut stderr = String::new();

        child.stdout.unwrap().read_to_string(&mut stdout)?;
        child.stderr.unwrap().read_to_string(&mut stderr)?;

        panic!("\nstdout: {}\nstderr: {}", stdout, stderr);
    }

    let mut file = std::fs::File::open(tmp_file.path())?;
    let mut data = Vec::new();

    file.read_to_end(&mut data)?;

    assert_eq!(
        data,
        vec![
            b'K', b'F', b'F', // KFF
            1, 0,  // Version number
            30, // Encoding
            1, 1, // Uniq, Canonical
            0, 0, 0, 24, // Free space size length
            112, 114, 111, 100, 117, 99, 101, 114, 58, 32, 107, 109, 101, 114, 95, 99, 111, 117,
            110, 116, 50, 107, 102, 102, 0, b'v', 0, 0, 0, 0, 0, 0, 0, 4, // Four value
            b'o', b'r', b'd', b'e', b'r', b'e', b'd', 0, 0, 0, 0, 0, 0, 0, 0, 1, //
            b'd', b'a', b't', b'a', b'_', b's', b'i', b'z', b'e', 0, 0, 0, 0, 0, 0, 0, 0,
            1, //
            b'k', 0, 0, 0, 0, 0, 0, 0, 0, 31, //
            b'm', b'a', b'x', 0, 0, 0, 0, 0, 0, 0, 0, 200, //
            b'r', 0, 0, 0, 0, 0, 0, 0, 10, // Ten block
            1, 0, 60, 205, 108, 136, 205, 162, 76, 20, // One kmer in block data 20
            1, 0, 243, 53, 178, 35, 54, 137, 60, 40, // One kmer in block data 40
            1, 3, 204, 214, 200, 140, 218, 36, 252, 60, // One kmer in block data 60
            1, 15, 51, 91, 34, 51, 104, 147, 240, 1, // One kmer in block data 1
            1, 60, 205, 108, 136, 205, 162, 79, 204, 2, // One kmer in block data 2
            1, 243, 53, 178, 35, 54, 137, 63, 48, 3, // One kmer in block data 3
            1, 204, 214, 200, 140, 218, 36, 252, 196, 4, // One kmer in block data 4
            1, 51, 91, 34, 51, 104, 147, 243, 28, 5, // One kmer in block data 5
            1, 205, 108, 136, 205, 162, 79, 204, 116, 6, // One kmer in block data 6
            1, 53, 178, 35, 54, 137, 63, 49, 216, 7, // One kmer in block data 7
            b'K', b'F', b'F', // KFF
        ]
    );

    Ok(())
}
