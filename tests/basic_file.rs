use kff;

fn read_file(path: &str) -> (Vec<Vec<u8>>, Vec<u8>) {
    let mut input = std::io::BufReader::new(std::fs::File::open(path).unwrap());

    let mut reader = kff::Reader::new(&mut input).unwrap();
    let rev_encoding = reader.rev_encoding();

    println!();
    println!("rev_encoding {:08b}", rev_encoding);
    
    let mut kmers: Vec<Vec<u8>> = Vec::new();
    let mut datas: Vec<u8> = Vec::new();

    while let Ok(section) = reader.next_section() {
        let mut it = section.into_iter();
        while let Some(Ok(kmer)) = it.next() {
	    println!("{} {} {}", kmer.bits(), String::from_utf8(kmer.seq(rev_encoding).to_vec()).unwrap(), kmer.data()[0]);
	    
            kmers.push(kmer.seq(rev_encoding).to_vec());
            datas.push(kmer.data()[0]);
        }
    }

    (kmers, datas)
}

#[test]
fn example() {
    let (_seqs, _datas) = read_file("tests/data/example.kff");

    /*
    assert_eq!(
        kmers,
        [
            b"ACTAAACTGA".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
            b"AAACTGATCG".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
            b"ACTAAACTGA".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
            b"AAACTGATCG".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
        ]
    );

    assert_eq!(datas, [32, 47, 1, 12, 1, 47, 32, 47, 1, 12, 1, 47,]);
    */
}

#[test]
fn example_tcga() {
    let (_seqs, _datas) = read_file("tests/data/example_TCGA.kff");

    /*
    assert_eq!(
        kmers,
        [
            b"ACTAAACTGA".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
            b"AAACTGATCG".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
            b"ACTAAACTGA".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
            b"AAACTGATCG".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
        ]
    );

    assert_eq!(datas, [32, 47, 1, 12, 1, 47, 32, 47, 1, 12, 1, 47,]);
    */
}

#[test]
fn m_1() {
    let (_seqs, _datas) = read_file("tests/data/m_1.kff");

    /*
    assert_eq!(
        kmers,
        [
            b"ACTAAACTGA".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
            b"AAACTGATCG".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
            b"ACTAAACTGA".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
            b"AAACTGATCG".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
        ]
    );

    assert_eq!(datas, [32, 47, 1, 12, 1, 47, 32, 47, 1, 12, 1, 47,]);
    */
}

#[test]
fn r_0() {
    let (_seqs, _datas) = read_file("tests/data/r_0.kff");

    /*
    assert_eq!(
        kmers,
        [
            b"ACTAAACTGA".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
            b"AAACTGATCG".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
            b"ACTAAACTGA".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
            b"AAACTGATCG".to_vec(),
            b"CTAAACTGAT".to_vec(),
            b"TAAACTGATT".to_vec(),
        ]
    );

    assert_eq!(datas, [32, 47, 1, 12, 1, 47, 32, 47, 1, 12, 1, 47,]);
    */
}
