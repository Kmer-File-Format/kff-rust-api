use kff;
use kff::seq2bits::Bits2Nuc;

fn read_file(path: &str) -> (Vec<Vec<u8>>, Vec<u8>) {
    let mut input = std::io::BufReader::new(std::fs::File::open(path).unwrap());

    let mut reader = kff::Reader::new(&mut input).unwrap();
    let rev_encoding = reader.rev_encoding();

    let mut kmers: Vec<Vec<u8>> = Vec::new();
    let mut datas: Vec<u8> = Vec::new();

    while let Ok(section) = reader.next_section() {
        let mut it = section.into_iter();
        while let Some(Ok(kmer)) = it.next() {
            kmers.push(kmer.seq().into_nuc(rev_encoding).to_vec());
            datas.push(kmer.data()[0]);
        }
    }

    (kmers, datas)
}

#[test]
fn example() {
    let (seqs, datas) = read_file("tests/data/example.kff");

    assert_eq!(
        seqs,
        [
            b"ACTAAACTGA",
            b"CTAAACTGAT",
            b"TAAACTGATT",
            b"AAACTGATCG",
            b"CTAAACTGAT",
            b"TAAACTGATT",
            b"ACTAAACTGA",
            b"CTAAACTGAT",
            b"TAAACTGATT",
            b"AAACTGATCG",
            b"CTAAACTGAT",
            b"TAAACTGATT",
        ]
    );

    assert_eq!(datas, [32, 47, 1, 12, 1, 47, 32, 47, 1, 12, 1, 47]);
}

#[test]
fn example_tcga() {
    let (seqs, datas) = read_file("tests/data/example_TCGA.kff");

    assert_eq!(
        seqs,
        [
            b"ACTAAACTGA",
            b"CTAAACTGAT",
            b"TAAACTGATT",
            b"AAACTGATCG",
            b"CTAAACTGAT",
            b"TAAACTGATT",
            b"ACTAAACTGA",
            b"CTAAACTGAT",
            b"TAAACTGATT",
            b"AAACTGATCG",
            b"CTAAACTGAT",
            b"TAAACTGATT",
        ]
    );

    assert_eq!(datas, [32, 47, 1, 12, 1, 47, 32, 47, 1, 12, 1, 47]);
}

#[test]
fn m_1() {
    let (seqs, datas) = read_file("tests/data/m_1.kff");

    assert_eq!(
        seqs,
        [
            b"ACTAAACTGA",
            b"CTAAACTGAT",
            b"TAAACTGATT",
            b"AAACTGATCG",
            b"CTAAACTGAT",
            b"TAAACTGATT",
        ]
    );

    assert_eq!(datas, [32, 47, 1, 12, 1, 47]);
}

#[test]
fn r_0() {
    let (seqs, datas) = read_file("tests/data/r_0.kff");

    assert_eq!(
        seqs,
        [
            b"ACTAAACTGA",
            b"CTAAACTGAT",
            b"TAAACTGATT",
            b"AAACTGATCG",
            b"CTAAACTGAT",
            b"TAAACTGATT",
        ]
    );

    assert_eq!(datas, [32, 47, 1, 12, 1, 47]);
}
