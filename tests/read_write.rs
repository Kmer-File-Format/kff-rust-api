use std::io::Read;
use std::io::Write;

use kff;
use kff::seq2bits::Bits2Nuc;
use kff::variables::Variables1;

fn read_and_write(i: &str, o: &str, k: u64, m: u64, data_size: u64) {
    let mut input = std::io::BufReader::new(std::fs::File::open(i).unwrap());
    let mut output = std::fs::File::create(o).unwrap();
    //let mut output = std::io::BufWriter::new(std::fs::File::create(output).unwrap());

    let mut reader = kff::Reader::new(&mut input).unwrap();
    let mut writer = kff::Writer::new(&mut output, reader.encoding(), reader.comment()).unwrap();

    writer.variables().insert("k".to_string(), k);
    writer.variables().insert("max".to_string(), m);
    writer
        .variables()
        .insert("data_size".to_string(), data_size);

    writer.write_variables().unwrap();

    while let Ok(mut section) = reader.next_section() {
        let mut seqs = Vec::new();
        let mut data = Vec::new();

        loop {
            if section.read_block().unwrap() == 0 {
                break;
            }

            seqs.push(kff::seq2bits::Seq2Bits::from_bitslice(section.block_seq()));
            data.push(Vec::from(section.block_data()));
        }
        writer.write_raw_section(&seqs[..], &data).unwrap();
    }
}

#[test]
fn r_0() {
    read_and_write("tests/data/r_0.kff", "tests/temp_r_0.kff", 10, 240, 1);

    let mut truth = Vec::new();
    let mut my = Vec::new();

    let mut input = std::io::BufReader::new(std::fs::File::open("tests/data/r_0.kff").unwrap());
    input.read_to_end(&mut truth).unwrap();

    input = std::io::BufReader::new(std::fs::File::open("tests/temp_r_0.kff").unwrap());
    input.read_to_end(&mut my).unwrap();

    std::fs::remove_file("tests/temp_r_0.kff");

    assert_eq!(truth[67..], my[67..]);
}

#[test]
fn kmers2kff() {
    let mut input =
        std::io::BufReader::new(std::fs::File::open("tests/data/kmers2kff.kff").unwrap());
    let mut reader = kff::Reader::new(&mut input).unwrap();

    let rev_encoding = reader.rev_encoding();

    while let Ok(section) = reader.next_section() {
        println!("section");
        let mut it = section.into_iter();
        while let Some(Ok(kmer)) = it.next() {
            println!(
                "{}",
                std::str::from_utf8(&kmer.seq().into_nuc(rev_encoding)).unwrap()
            );
        }
    }
}
