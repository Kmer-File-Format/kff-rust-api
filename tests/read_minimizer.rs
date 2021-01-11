use kff::seq2bits::Bits2Nuc;

const TRUTH: &[&[u8]] = &[
    b"AAAAATCAACT",
    b"AAAACGATTGC",
    b"AAAATCAACTG",
    b"AAACACATACA",
    b"AAACACTGACA",
    b"AAACACGGCTA",
    b"AAACCGTCAAG",
    b"AAACCGGCTTT",
    b"AAACTTGCCAT",
    b"AAACGATTGCA",
    b"AAATAACCAGC",
    b"AAATAGCTAGC",
    b"AAATCCTTCCT",
    b"AAAGAAACCAG",
    b"AAAGACGCGAA",
    b"AAAGATGTAAA",
    b"AAAGCAAGCCC",
    b"AAAGCCCTAAA",
    b"AAAGCCGAATC",
    b"AAAGTAGGCCC",
    b"AAAGGATTGAG",
    b"AAAGGCTTGGT",
    b"AACAATATGAT",
    b"AACACATAAAC",
    b"AACACATACAG",
    b"AACACCAGGGC",
    b"AACACTGACAG",
    b"AACACGGTTCT",
    b"AACATTCAACC",
    b"AACAGCGCGAT",
    b"AACAGCGTGTG",
    b"AACAGGGGTTC",
    b"AACCAAGCCTT",
    b"AACCACAAAAT",
    b"AACCAGCGGAC",
    b"AACCCCTGTTT",
    b"AACCTAATGAA",
    b"AACCTTACGCA",
    b"AACCGTCAAGA",
    b"AACCGTGTTTG",
    b"AACCGGCTTTG",
    b"AACTCCATCTT",
    b"AACTCTCAAGT",
    b"AACTTAGTCTG",
    b"AACTTGCCATG",
    b"AACTTGTGTAC",
    b"AACTTGGACCT",
    b"AACGAACAGGG",
    b"AACGAACGCAC",
    b"AACGATTTCGC",
];

#[test]
#[should_panic]
fn kmers2kff() {
    let mut truth = std::collections::HashSet::new();

    for k in TRUTH {
        truth.insert(k.to_vec().into_boxed_slice());
    }

    let mut input =
        std::io::BufReader::new(std::fs::File::open("tests/data/kmers2kff.kff").unwrap());
    let mut reader = kff::Reader::new(&mut input).unwrap();

    let rev_encoding = reader.rev_encoding();

    let mut reads = std::collections::HashSet::new();
    while let Ok(section) = reader.next_section() {
        let mut it = section.into_iter();
        while let Some(Ok(kmer)) = it.next() {
            reads.insert(kmer.seq().into_nuc(rev_encoding));
        }
    }

    assert_eq!(truth, reads);
}
