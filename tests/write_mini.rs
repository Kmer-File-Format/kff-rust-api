use kff;

#[test]
fn write_minimizer() {
    // seq GCTTCAAAGGATTGAG
    // min      AAAGGA

    let mut output = std::io::Cursor::new(Vec::new());
    let mut writer = kff::Writer::new(&mut output, 0b00011011, b"").unwrap();

    let sequence = b"GCTTCTTGAG";
    let mini = b"AAAGGA";

    writer.variables().insert("k".to_string(), 11);
    writer.variables().insert("m".to_string(), 6);
    writer.variables().insert("max".to_string(), 255);
    writer.variables().insert("data_size".to_string(), 0);

    writer.write_variables().unwrap();

    let data: Vec<u8> = Vec::new();

    writer
        .write_minimizer_seq_section(mini, &[5u64], &[sequence], &[data])
        .unwrap();
    let out = output.into_inner();

    assert_eq!(vec![109, 0, 60, 1, 0, 0, 0, 6, 5, 13, 166, 179], out[66..]);
}
