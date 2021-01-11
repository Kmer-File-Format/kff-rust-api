use std::io::Read;

use kff;
use kff::seq2bits::Bits2Nuc;

fn compare_file(truth_path: &str, test_path: &str, position: usize) {
    let mut truth = Vec::new();
    let mut my = Vec::new();

    let mut input = std::io::BufReader::new(std::fs::File::open(truth_path).unwrap());
    input.read_to_end(&mut truth).unwrap();

    input = std::io::BufReader::new(std::fs::File::open(test_path).unwrap());
    input.read_to_end(&mut my).unwrap();

    assert_eq!(truth[position..], my[position..]);
}

#[test]
fn r_0() {
    let truth_path = "tests/data/r_0.kff";
    let test_path = "tests/temp_r_0.kff";

    let mut input = std::io::BufReader::new(std::fs::File::open(truth_path).unwrap());
    let mut output = std::fs::File::create(test_path).unwrap();

    let mut reader = kff::Reader::new(&mut input).unwrap();
    let mut writer = kff::Writer::new(&mut output, reader.encoding(), reader.comment()).unwrap();

    writer.variables().insert("k".to_string(), 10);
    writer.variables().insert("max".to_string(), 240);
    writer.variables().insert("data_size".to_string(), 1);

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

    compare_file(truth_path, test_path, 67);

    std::fs::remove_file(test_path).unwrap();
}

#[test]
fn m_1() {
    let truth_path = "tests/data/m_1.kff";
    let test_path = "tests/temp_m_1.kff";

    let mut input = std::io::BufReader::new(std::fs::File::open(truth_path).unwrap());
    let mut output = std::fs::File::create(test_path).unwrap();

    let mut reader = kff::Reader::new(&mut input).unwrap();
    let mut writer = kff::Writer::new(&mut output, reader.encoding(), reader.comment()).unwrap();

    let rev_encode = reader.rev_encoding();

    writer.variables().insert("k".to_string(), 10);
    writer.variables().insert("m".to_string(), 8);
    writer.variables().insert("max".to_string(), 240);
    writer.variables().insert("data_size".to_string(), 1);

    writer.write_variables().unwrap();

    let mini = "AAACTGAT";

    while let Ok(mut section) = reader.next_section() {
        let mut poss = Vec::new();
        let mut seqs = Vec::new();
        let mut data = Vec::new();

        loop {
            if section.read_block().unwrap() == 0 {
                break;
            }

            let pos = String::from_utf8(section.block_seq().into_nuc(rev_encode).into_vec())
                .unwrap()
                .find(mini)
                .unwrap() as u64;

            let mut seq = bitvec::vec::BitVec::new();

            seq.extend(section.block_seq()[..(pos as usize * 2)].iter());
            seq.extend(section.block_seq()[((pos as usize + mini.len()) * 2)..].iter());

            poss.push(pos);
            seqs.push(seq.into_boxed_bitslice());
            data.push(Vec::from(section.block_data()));
        }
        writer
            .write_minimizer_section(mini.as_bytes(), &poss, &seqs[..], &data)
            .unwrap();
    }

    compare_file(truth_path, test_path, 77);

    std::fs::remove_file(test_path).unwrap();
}
