//! Check example produce the good output

use std::io::Read;
use std::process::{Command, Stdio};

#[test]
#[cfg(not(tarpaulin))]
fn kff2kmers() -> kff::error::Result<()> {
    let args = vec![
        "run",
        "--example",
        "kff2kmer",
        "--",
        "-i",
        "tests/data/test.kff",
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

    let mut value_str = String::new();
    child.stdout.unwrap().read_to_string(&mut value_str)?;
    value_str.pop().unwrap();

    let value: Vec<Vec<u8>> = value_str
        .split('\n')
        .map(|x| x.bytes().collect::<Vec<u8>>())
        .collect();

    assert_eq!(
        value,
        vec![
            b"AAAAAGGAGAGCCTGATATAGAGCTTATCAG",
            b"AAAAGGAGAGCCTGATATAGAGCTTATCAGG",
            b"AAAGGAGAGCCTGATATAGAGCTTATCAGGG",
            b"AAGGAGAGCCTGATATAGAGCTTATCAGGGA",
            b"AGGAGAGCCTGATATAGAGCTTATCAGGGAG",
            b"GGAGAGCCTGATATAGAGCTTATCAGGGAGA",
            b"GAGAGCCTGATATAGAGCTTATCAGGGAGAC",
            b"AGAGCCTGATATAGAGCTTATCAGGGAGACG",
            b"GAGCCTGATATAGAGCTTATCAGGGAGACGC",
            b"AGCCTGATATAGAGCTTATCAGGGAGACGCT",
            b"GCCTGATATAGAGCTTATCAGGGAGACGCTT",
            b"CCTGATATAGAGCTTATCAGGGAGACGCTTT",
            b"CTGATATAGAGCTTATCAGGGAGACGCTTTC",
            b"TGATATAGAGCTTATCAGGGAGACGCTTTCA",
            b"GATATAGAGCTTATCAGGGAGACGCTTTCAA",
            b"ATATAGAGCTTATCAGGGAGACGCTTTCAAA",
            b"TATAGAGCTTATCAGGGAGACGCTTTCAAAG",
            b"ATAGAGCTTATCAGGGAGACGCTTTCAAAGA",
            b"AAAAGCTATAGTGTGGGCCCCTAAAGCCGCG",
            b"AAAGCTATAGTGTGGGCCCCTAAAGCCGCGA",
            b"AAGCTATAGTGTGGGCCCCTAAAGCCGCGAG",
            b"AGCTATAGTGTGGGCCCCTAAAGCCGCGAGA",
            b"GCTATAGTGTGGGCCCCTAAAGCCGCGAGAG",
            b"CTATAGTGTGGGCCCCTAAAGCCGCGAGAGC",
            b"TATAGTGTGGGCCCCTAAAGCCGCGAGAGCG",
            b"ATAGTGTGGGCCCCTAAAGCCGCGAGAGCGG",
            b"TAGTGTGGGCCCCTAAAGCCGCGAGAGCGGC",
            b"AGTGTGGGCCCCTAAAGCCGCGAGAGCGGCC",
            b"GTGTGGGCCCCTAAAGCCGCGAGAGCGGCCC",
            b"TGTGGGCCCCTAAAGCCGCGAGAGCGGCCCA",
            b"GTGGGCCCCTAAAGCCGCGAGAGCGGCCCAA",
            b"GGGCGGCCGCCTGGTCTAGGGCCATAGCGAC",
            b"GGCGGCCGCCTGGTCTAGGGCCATAGCGACA",
            b"GCGGCCGCCTGGTCTAGGGCCATAGCGACAT",
            b"CGGCCGCCTGGTCTAGGGCCATAGCGACATC",
            b"GGCCGCCTGGTCTAGGGCCATAGCGACATCA",
            b"GCCGCCTGGTCTAGGGCCATAGCGACATCAG",
            b"CCGCCTGGTCTAGGGCCATAGCGACATCAGC",
            b"CGCCTGGTCTAGGGCCATAGCGACATCAGCC",
            b"GCCTGGTCTAGGGCCATAGCGACATCAGCCG",
            b"CCTGGTCTAGGGCCATAGCGACATCAGCCGT",
            b"AAAAGGGTGCTTTCATACTGTTAAGAGACCA",
            b"AAAGGGTGCTTTCATACTGTTAAGAGACCAT",
            b"AAGGGTGCTTTCATACTGTTAAGAGACCATT",
            b"AACTTCAACGGGGTACCCACCGGAACTTCAA",
            b"ACTTCAACGGGGTACCCACCGGAACTTCAAC",
            b"CTTCAACGGGGTACCCACCGGAACTTCAACG",
            b"TTCAACGGGGTACCCACCGGAACTTCAACGG",
            b"TCAACGGGGTACCCACCGGAACTTCAACGGG",
            b"CAACGGGGTACCCACCGGAACTTCAACGGGG",
            b"AACGGGGTACCCACCGGAACTTCAACGGGGT",
            b"ACGGGGTACCCACCGGAACTTCAACGGGGTA",
            b"CGGGGTACCCACCGGAACTTCAACGGGGTAC",
            b"GGGGTACCCACCGGAACTTCAACGGGGTACC",
            b"GGGTACCCACCGGAACTTCAACGGGGTACCC",
            b"GGTACCCACCGGAACTTCAACGGGGTACCCA",
            b"AGAGCCTTACGTAGACGCCATTGCCATAGAC",
            b"GAGCCTTACGTAGACGCCATTGCCATAGACA",
            b"AGCCTTACGTAGACGCCATTGCCATAGACAG",
            b"GCCTTACGTAGACGCCATTGCCATAGACAGC",
            b"CCTTACGTAGACGCCATTGCCATAGACAGCA",
            b"ATACAGTACTCATGAGTTCTTGATCAGGTGG",
            b"TACAGTACTCATGAGTTCTTGATCAGGTGGC",
            b"ACAGTACTCATGAGTTCTTGATCAGGTGGCA",
            b"CAGTACTCATGAGTTCTTGATCAGGTGGCAT",
            b"AGTACTCATGAGTTCTTGATCAGGTGGCATT",
            b"GTACTCATGAGTTCTTGATCAGGTGGCATTT",
            b"TACTCATGAGTTCTTGATCAGGTGGCATTTG",
            b"ACTCATGAGTTCTTGATCAGGTGGCATTTGG",
            b"CTCATGAGTTCTTGATCAGGTGGCATTTGGC",
            b"TCATGAGTTCTTGATCAGGTGGCATTTGGCA",
            b"CATGAGTTCTTGATCAGGTGGCATTTGGCAT",
            b"ATGAGTTCTTGATCAGGTGGCATTTGGCATC",
            b"TGAGTTCTTGATCAGGTGGCATTTGGCATCC",
            b"GAGTTCTTGATCAGGTGGCATTTGGCATCCT",
            b"AGTTCTTGATCAGGTGGCATTTGGCATCCTT",
            b"GTTCTTGATCAGGTGGCATTTGGCATCCTTT",
            b"TTCTTGATCAGGTGGCATTTGGCATCCTTTT",
            b"AAAGACACGAAAGCATGAATGGTTCATCAGC",
            b"AAGACACGAAAGCATGAATGGTTCATCAGCA",
            b"AGACACGAAAGCATGAATGGTTCATCAGCAC",
            b"GACACGAAAGCATGAATGGTTCATCAGCACC",
            b"ACACGAAAGCATGAATGGTTCATCAGCACCC",
            b"ATTGTAGTCGTCCAACGGACGTAGGGATGTG",
            b"TTGTAGTCGTCCAACGGACGTAGGGATGTGA",
        ]
    );

    Ok(())
}
