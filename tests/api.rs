//! Tests of the API, ported from the unit tests of the C++ API (kff-cpp-api/tests).
//!
//! The expected bytes are built from the specification of the format (kff-reference/README.md),
//! independently of the writer.

/* std use */
use std::io::Cursor;

/* project use */
use kff::error::Error;
use kff::{BlockReader, Encoding, KmerReader, Reader, Section, SectionType, Writer};

type Bytes = Vec<u8>;

// ----- Helpers -----

/// Deterministic pseudo-random generator (splitmix64)
struct Random(u64);

impl Random {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    /// Inclusive bounds
    fn between(&mut self, min: u64, max: u64) -> u64 {
        min + self.next() % (max - min + 1)
    }
    fn sequence(&mut self, size: u64) -> Vec<u8> {
        (0..size)
            .map(|_| b"ACGT"[self.between(0, 3) as usize])
            .collect()
    }
    fn bytes(&mut self, size: u64) -> Bytes {
        (0..size).map(|_| self.next() as u8).collect()
    }
}

#[derive(Clone, Debug)]
struct Block {
    /// Full sequence (with the minimizer for 'm' sections)
    sequence: Vec<u8>,
    data: Bytes,
    mini_pos: u64,
}

impl Block {
    fn without_mini(&self, m: u64) -> Vec<u8> {
        let mut without = self.sequence[..self.mini_pos as usize].to_vec();
        without.extend_from_slice(&self.sequence[(self.mini_pos + m) as usize..]);
        without
    }
}

fn random_block(
    rnd: &mut Random,
    k: u64,
    nb_kmers: u64,
    data_size: u64,
    minimizer: &[u8],
) -> Block {
    let mut sequence = rnd.sequence(k + nb_kmers - 1);
    let data = rnd.bytes(data_size * nb_kmers);
    let mut mini_pos = 0;
    if !minimizer.is_empty() {
        mini_pos = rnd.between(0, sequence.len() as u64 - minimizer.len() as u64);
        sequence[mini_pos as usize..mini_pos as usize + minimizer.len()].copy_from_slice(minimizer);
    }
    Block {
        sequence,
        data,
        mini_pos,
    }
}

fn random_blocks(rnd: &mut Random, k: u64, max: u64, data_size: u64, nb: usize) -> Vec<Block> {
    (0..nb)
        .map(|_| {
            let nb_kmers = rnd.between(1, max);
            random_block(rnd, k, nb_kmers, data_size, b"")
        })
        .collect()
}

fn encode(sequence: &[u8]) -> Bytes {
    Encoding::default().encode(sequence).unwrap()
}

fn decode(bytes: &[u8], size: usize) -> Vec<u8> {
    Encoding::default().decode(bytes, size)
}

fn big_endian(value: u64, nb_bytes: usize) -> Bytes {
    (0..nb_bytes)
        .rev()
        .map(|i| if i < 8 { (value >> (8 * i)) as u8 } else { 0 })
        .collect()
}

fn c_string(s: &str) -> Bytes {
    let mut bytes = s.as_bytes().to_vec();
    bytes.push(0);
    bytes
}

/// Number of bytes of a value that can go up to max (ceil(log2(max)) bits)
fn bytes_for_max(max: u64) -> usize {
    let mut bits = 0;
    while bits < 64 && (1u64 << bits) < max {
        bits += 1;
    }
    (bits + 7) / 8
}

fn header_bytes(encoding: [u8; 4], unique: bool, canonical: bool, metadata: &[u8]) -> Bytes {
    let mut b = vec![
        b'K',
        b'F',
        b'F',
        1,
        0,
        (encoding[0] << 6) | (encoding[1] << 4) | (encoding[2] << 2) | encoding[3],
        unique as u8,
        canonical as u8,
    ];
    b.extend(big_endian(metadata.len() as u64, 4));
    b.extend_from_slice(metadata);
    b
}

fn default_header() -> Bytes {
    header_bytes([0, 1, 3, 2], false, false, b"")
}

/// Variables are written sorted by name
fn variables_bytes(sorted_vars: &[(&str, u64)]) -> Bytes {
    let mut b = vec![b'v'];
    b.extend(big_endian(sorted_vars.len() as u64, 8));
    for (name, value) in sorted_vars {
        b.extend(c_string(name));
        b.extend(big_endian(*value, 8));
    }
    b
}

fn raw_section_bytes(blocks: &[Block], k: u64, max: u64, encoding: &Encoding) -> Bytes {
    let mut b = vec![b'r'];
    b.extend(big_endian(blocks.len() as u64, 8));
    for block in blocks {
        b.extend(big_endian(
            block.sequence.len() as u64 - k + 1,
            bytes_for_max(max),
        ));
        b.extend(encoding.encode(&block.sequence).unwrap());
        b.extend_from_slice(&block.data);
    }
    b
}

fn minimizer_section_bytes(
    minimizer: &[u8],
    blocks: &[Block],
    k: u64,
    max: u64,
    encoding: &Encoding,
) -> Bytes {
    let m = minimizer.len() as u64;
    let mut b = vec![b'm'];
    b.extend(encoding.encode(minimizer).unwrap());
    b.extend(big_endian(blocks.len() as u64, 8));
    for block in blocks {
        b.extend(big_endian(
            block.sequence.len() as u64 - k + 1,
            bytes_for_max(max),
        ));
        b.extend(big_endian(block.mini_pos, bytes_for_max(k + max - 1)));
        b.extend(encoding.encode(&block.without_mini(m)).unwrap());
        b.extend_from_slice(&block.data);
    }
    b
}

/// Size of the automatic footer: values section with first_index and footer_size
const FOOTER_SIZE: u64 = 1 + 8 + (12 + 8) + (12 + 8);

fn new_writer() -> Writer<Cursor<Bytes>> {
    Writer::new(Cursor::new(Vec::new())).unwrap()
}

fn bytes_of(writer: Writer<Cursor<Bytes>>) -> Bytes {
    writer.close().unwrap().into_inner()
}

fn reader_of(bytes: &[u8]) -> kff::error::Result<Reader<Cursor<Bytes>>> {
    Reader::new(Cursor::new(bytes.to_vec()))
}

/// Write a values section, returns its position
fn write_vars<W>(writer: &mut Writer<W>, vars: &[(&str, u64)]) -> u64
where
    W: std::io::Write + std::io::Seek,
{
    let mut section = writer.values_section().unwrap();
    let beginning = section.beginning();
    for (name, value) in vars {
        section.write_var(name, *value).unwrap();
    }
    section.close().unwrap();
    beginning
}

// ----- Header -----

#[test]
fn header_default_header_of_a_new_file() {
    let bytes = bytes_of(new_writer());
    let mut expected = default_header();
    // Index without section
    expected.push(b'i');
    expected.extend(big_endian(0, 8));
    expected.extend(big_endian(0, 8));
    // Footer
    expected.extend(variables_bytes(&[
        ("first_index", 12),
        ("footer_size", FOOTER_SIZE),
    ]));
    expected.extend_from_slice(b"KFF");
    assert_eq!(bytes, expected);
    assert_eq!(bytes.len(), 12 + 17 + 49 + 3);
}

#[test]
fn header_encodings_flags_and_metadata_round_trip() {
    let metadata = b"hello\0\xffworld";
    let mut permutation = 0;
    for byte in 0..=255u8 {
        let encoding = match Encoding::from_byte(byte) {
            Ok(encoding) => encoding,
            Err(_) => continue,
        };
        let unique = permutation % 2 == 0;
        let canonical = permutation % 3 == 0;
        let values = encoding.values();

        let mut writer = new_writer();
        writer
            .set_encoding(Encoding::new(values[0], values[1], values[2], values[3]).unwrap())
            .unwrap();
        writer.set_uniqueness(unique).unwrap();
        writer.set_canonicity(canonical).unwrap();
        writer.write_metadata(metadata).unwrap();
        let bytes = bytes_of(writer);
        assert_eq!(
            bytes[..12 + metadata.len()],
            header_bytes(values, unique, canonical, metadata)[..]
        );

        let reader = reader_of(&bytes).unwrap();
        assert_eq!(reader.encoding().values(), values);
        assert_eq!(reader.uniqueness(), unique);
        assert_eq!(reader.canonicity(), canonical);
        assert_eq!(reader.header().major_version, 1);
        assert_eq!(reader.header().minor_version, 0);
        assert_eq!(reader.metadata(), metadata);
        assert_eq!(reader.first_section_position(), 12 + metadata.len() as u64);
        permutation += 1;
    }
    assert_eq!(permutation, 24);
}

#[test]
fn header_flags_can_change_after_the_sections() {
    let mut writer = new_writer();
    write_vars(&mut writer, &[("k", 3)]);
    writer.set_canonicity(true).unwrap();
    writer
        .set_encoding(Encoding::new(3, 2, 1, 0).unwrap())
        .unwrap();
    assert!(matches!(
        writer.write_metadata(b"late"),
        Err(Error::MetadataAfterSection)
    ));
    let reader = reader_of(&bytes_of(writer)).unwrap();
    assert!(reader.canonicity());
    assert!(!reader.uniqueness());
    assert_eq!(reader.encoding().values(), [3, 2, 1, 0]);
}

#[test]
fn header_invalid_encoding_is_rejected_when_read() {
    let mut bytes = bytes_of(new_writer());
    bytes[5] = 0b00011100; // A=0 C=1 G=3 T=0
    assert!(matches!(reader_of(&bytes), Err(Error::InvalidEncoding(_))));
}

#[test]
fn header_missing_signatures_are_detected() {
    let bytes = bytes_of(new_writer());

    let mut start = bytes.clone();
    start[1] = b'X';
    assert!(matches!(
        reader_of(&start),
        Err(Error::MissingSignature("beginning"))
    ));

    let mut end = bytes.clone();
    *end.last_mut().unwrap() = b'X';
    assert!(matches!(
        reader_of(&end),
        Err(Error::MissingSignature("end"))
    ));

    assert!(reader_of(b"KF").is_err());
    assert!(reader_of(b"").is_err());
}

#[test]
fn header_files_from_a_newer_format_version_are_rejected() {
    let bytes = bytes_of(new_writer());

    let mut major = bytes.clone();
    major[3] = 2;
    assert!(matches!(
        reader_of(&major),
        Err(Error::UnsupportedVersion { major: 2, minor: 0 })
    ));

    let mut minor = bytes.clone();
    minor[4] = 1;
    assert!(matches!(
        reader_of(&minor),
        Err(Error::UnsupportedVersion { major: 1, minor: 1 })
    ));

    // Older versions are read
    let mut older = bytes;
    older[3] = 0;
    older[4] = 9;
    assert!(reader_of(&older).is_ok());
}

#[test]
fn header_a_missing_file_cannot_be_opened() {
    let dir = tempfile::tempdir().unwrap();
    assert!(matches!(
        Reader::open(dir.path().join("absent.kff")),
        Err(Error::Io(_))
    ));
}

#[test]
fn header_smallest_valid_file_without_section_nor_footer() {
    let mut bytes = default_header();
    bytes.extend_from_slice(b"KFF");
    assert_eq!(bytes.len(), 15);

    let mut reader = reader_of(&bytes).unwrap();
    assert!(reader.footer().is_none());
    assert_eq!(reader.tell(), 12);
    assert_eq!(reader.end_position(), 12);
    assert_eq!(reader.read_section_type().unwrap(), None);

    let mut reader = KmerReader::new(reader_of(&bytes).unwrap());
    assert!(!reader.has_next().unwrap());
}

// ----- Variables -----

#[test]
fn variables_names_are_written_in_sorted_order() {
    let mut writer = new_writer();
    writer.set_indexation(false);
    let mut section = writer.values_section().unwrap();
    assert_eq!(section.beginning(), 12);
    section.write_var("max", 255).unwrap();
    section.write_var("k", 10).unwrap();
    section.write_var("data_size", 1).unwrap();
    assert_eq!(section.vars().len(), 3);
    section.close().unwrap();

    let mut expected = default_header();
    expected.extend(variables_bytes(&[
        ("data_size", 1),
        ("k", 10),
        ("max", 255),
    ]));
    expected.extend_from_slice(b"KFF");
    assert_eq!(bytes_of(writer), expected);
}

#[test]
fn variables_64_bits_values_and_arbitrary_names_round_trip() {
    let values: &[(&str, u64)] = &[
        ("zero", 0),
        ("one", 1),
        ("2^32", 1 << 32),
        ("2^63", 1 << 63),
        ("max", u64::MAX),
        ("name with spaces", 42),
        ("x", 0x0123456789abcdef),
    ];
    let mut writer = new_writer();
    write_vars(&mut writer, values);
    let mut reader = reader_of(&bytes_of(writer)).unwrap();
    assert_eq!(
        reader.read_section_type().unwrap(),
        Some(SectionType::Values)
    );
    let vars = reader.read_values().unwrap();
    assert_eq!(vars.len(), values.len());
    for (name, value) in values {
        assert_eq!(vars.get(*name), Some(value));
        assert_eq!(reader.var(name), Some(*value));
    }
}

#[test]
fn variables_a_new_section_resets_the_scope() {
    let mut writer = new_writer();
    write_vars(&mut writer, &[("k", 10), ("m", 5)]);
    assert_eq!(writer.global_vars().len(), 2);
    let section = writer.values_section().unwrap();
    drop(section);
    assert_eq!(writer.global_vars().len(), 0);
    write_vars(&mut writer, &[("k", 12)]);
    assert_eq!(writer.var("k"), Some(12));

    let mut reader = reader_of(&bytes_of(writer)).unwrap();
    reader.read_values().unwrap();
    assert_eq!(reader.global_vars().len(), 2);
    reader.read_values().unwrap();
    assert_eq!(reader.global_vars().len(), 0);
    reader.read_values().unwrap();
    assert_eq!(reader.global_vars().len(), 1);
    assert_eq!(reader.var("k"), Some(12));
    assert_eq!(reader.var("m"), None);
}

#[test]
fn variables_empty_section() {
    let mut writer = new_writer();
    writer.set_indexation(false);
    write_vars(&mut writer, &[]);
    let bytes = bytes_of(writer);
    let mut expected = default_header();
    expected.extend(variables_bytes(&[]));
    expected.extend_from_slice(b"KFF");
    assert_eq!(bytes, expected);

    let mut reader = reader_of(&bytes).unwrap();
    assert!(reader.read_values().unwrap().is_empty());
}

#[test]
fn variables_sections_written_in_any_order_are_read() {
    let mut bytes = default_header();
    bytes.push(b'v');
    bytes.extend(big_endian(3, 8));
    for (name, value) in [("k", 10), ("max", 255), ("data_size", 1)] {
        bytes.extend(c_string(name));
        bytes.extend(big_endian(value, 8));
    }
    bytes.extend_from_slice(b"KFF");

    let mut reader = reader_of(&bytes).unwrap();
    let vars = reader.read_values().unwrap();
    assert_eq!(vars.get("k"), Some(&10));
    assert_eq!(vars.get("max"), Some(&255));
    assert_eq!(vars.get("data_size"), Some(&1));
    assert_eq!(reader.tell(), bytes.len() as u64 - 3);
}

#[test]
fn variables_the_footer_is_discovered_when_the_file_is_opened() {
    let mut writer = new_writer();
    write_vars(&mut writer, &[("k", 21)]);
    let mut reader = reader_of(&bytes_of(writer)).unwrap();
    let footer = reader.footer().unwrap();
    assert_eq!(footer.get("footer_size"), Some(&49));
    assert_eq!(footer.get("first_index"), Some(&(12 + 19)));
    // The footer is not in the scope of the first section
    assert!(reader.global_vars().is_empty());
    assert_eq!(
        reader.read_section_type().unwrap(),
        Some(SectionType::Values)
    );
}

#[test]
fn variables_a_variable_written_twice_is_declared_once() {
    let mut writer = new_writer();
    write_vars(&mut writer, &[("k", 10), ("k", 12), ("max", 1)]);
    let mut reader = reader_of(&bytes_of(writer)).unwrap();
    let vars = reader.read_values().unwrap();
    assert_eq!(vars.len(), 2);
    assert_eq!(vars.get("k"), Some(&12));
    assert_eq!(
        reader.read_section_type().unwrap(),
        Some(SectionType::Index)
    );
}

#[test]
fn variables_names_cannot_contain_a_null_byte() {
    let mut writer = new_writer();
    let mut section = writer.values_section().unwrap();
    assert!(matches!(
        section.write_var("a\0b", 1),
        Err(Error::InvalidVariableName(_))
    ));
}

#[test]
fn variables_copy_to_another_file() {
    let mut source = new_writer();
    source.set_indexation(false);
    write_vars(&mut source, &[]);
    write_vars(&mut source, &[("k", 31), ("data_size", 2)]);
    let mut reader = reader_of(&bytes_of(source)).unwrap();

    let mut copy = new_writer();
    copy.set_indexation(false);
    // As with the C++ API, the copy of an empty section is not written: copy only non empty sections
    for _ in 0..2 {
        let vars = reader.read_values().unwrap();
        if !vars.is_empty() {
            copy.write_values(&vars).unwrap();
        }
    }
    let mut expected = default_header();
    expected.extend(variables_bytes(&[("data_size", 2), ("k", 31)]));
    expected.extend_from_slice(b"KFF");
    assert_eq!(bytes_of(copy), expected);
}

// ----- Raw sections -----

/// Write a file made of one 'v' section and one 'r' section (no index). Returns the bytes and the
/// position of the raw section.
fn write_raw_file(k: u64, max: u64, data_size: u64, blocks: &[Block]) -> (Bytes, u64) {
    let mut writer = new_writer();
    writer.set_indexation(false);
    write_vars(
        &mut writer,
        &[("k", k), ("max", max), ("data_size", data_size)],
    );
    let mut section = writer.raw_section().unwrap();
    let beginning = section.beginning();
    for block in blocks {
        section
            .write_compacted_sequence(
                &encode(&block.sequence),
                block.sequence.len() as u64,
                &block.data,
            )
            .unwrap();
    }
    assert_eq!(section.nb_blocks(), blocks.len() as u64);
    section.close().unwrap();
    (bytes_of(writer), beginning)
}

#[test]
fn raw_section_bytes_and_round_trip_for_many_k_max_and_data_sizes() {
    let mut rnd = Random(1);
    let encoding = Encoding::default();
    for k in [
        1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127,
    ] {
        // Powers of 256 are avoided: n == max does not fit in the counter field
        for max in [1, 2, 3, 5, 100, 255, 300] {
            for data_size in [0, 1, 2, 3, 8] {
                let mut blocks = random_blocks(&mut rnd, k, max, data_size, 20);
                blocks.push(random_block(&mut rnd, k, max, data_size, b"")); // a full block
                let (bytes, beginning) = write_raw_file(k, max, data_size, &blocks);

                // Exact bytes
                assert_eq!(
                    bytes[beginning as usize..bytes.len() - 3],
                    raw_section_bytes(&blocks, k, max, &encoding)[..]
                );

                // Read back
                let mut reader = reader_of(&bytes).unwrap();
                reader.read_values().unwrap();
                let mut section = reader.open_raw().unwrap();
                assert_eq!(section.nb_blocks(), blocks.len() as u64);
                assert_eq!(section.k(), k);
                assert_eq!(section.max(), max);
                assert_eq!(section.data_size(), data_size);
                let mut seq = vec![0; section.max_sequence_bytes()];
                let mut data = vec![0; section.max_data_bytes()];
                for block in &blocks {
                    let nb_kmers = section
                        .read_compacted_sequence(&mut seq, &mut data)
                        .unwrap();
                    assert_eq!(nb_kmers, block.sequence.len() as u64 - k + 1);
                    assert_eq!(decode(&seq, block.sequence.len()), block.sequence);
                    assert_eq!(data[..block.data.len()], block.data[..]);
                }
                assert_eq!(section.remaining_blocks(), 0);
                section.close().unwrap();
            }
        }
    }
}

#[test]
fn raw_width_of_the_kmer_counter_depends_on_max() {
    let widths: &[(u64, usize)] = &[
        (1, 0),
        (2, 1),
        (3, 1),
        (255, 1),
        (256, 1),
        (257, 2),
        (65535, 2),
        (65536, 2),
        (65537, 3),
        (1 << 24, 3),
        ((1 << 24) + 1, 4),
        (1 << 32, 4),
        ((1 << 32) + 1, 5),
    ];
    let mut rnd = Random(2);
    for (max, width) in widths {
        let mut blocks = vec![
            random_block(&mut rnd, 5, 1, 1, b""),
            random_block(&mut rnd, 5, 3, 1, b""),
        ];
        if *max < 3 {
            blocks.pop();
        }
        let (bytes, beginning) = write_raw_file(5, *max, 1, &blocks);
        let expected_size: usize = 9 + blocks
            .iter()
            .map(|block| width + 2 + block.data.len())
            .sum::<usize>();
        assert_eq!(bytes.len() - 3 - beginning as usize, expected_size);
        assert_eq!(bytes_for_max(*max), *width);
        // Read back
        let mut reader = reader_of(&bytes).unwrap();
        reader.read_values().unwrap();
        let mut section = reader.open_raw().unwrap();
        for block in &blocks {
            let read = section.read_block().unwrap().unwrap();
            assert_eq!(decode(&read.sequence, block.sequence.len()), block.sequence);
            assert_eq!(read.data, block.data);
        }
        assert!(section.read_block().unwrap().is_none());
    }
}

#[test]
fn raw_read_compacted_sequence_joined_puts_the_data_after_the_sequence() {
    let mut rnd = Random(3);
    for k in [3, 4, 5, 6, 21] {
        let blocks = random_blocks(&mut rnd, k, 50, 2, 30);
        let (bytes, _) = write_raw_file(k, 50, 2, &blocks);
        let mut reader = reader_of(&bytes).unwrap();
        reader.read_values().unwrap();
        let mut section = reader.open_raw().unwrap();
        let mut buffer = vec![0; section.max_sequence_bytes() + section.max_data_bytes()];
        for block in &blocks {
            let nb_kmers = section.read_compacted_sequence_joined(&mut buffer).unwrap();
            let seq_bytes = (block.sequence.len() + 3) / 4;
            assert_eq!(nb_kmers, block.sequence.len() as u64 - k + 1);
            assert_eq!(decode(&buffer, block.sequence.len()), block.sequence);
            assert_eq!(
                buffer[seq_bytes..seq_bytes + block.data.len()],
                block.data[..]
            );
        }
    }
}

#[test]
fn raw_jump_over_blocks_and_sections() {
    let mut rnd = Random(4);
    let blocks = random_blocks(&mut rnd, 11, 20, 1, 10);
    let (bytes, _) = write_raw_file(11, 20, 1, &blocks);

    // jump_sequence skips exactly one block, close jumps over the remaining blocks
    {
        let mut reader = reader_of(&bytes).unwrap();
        reader.read_values().unwrap();
        let mut section = reader.open_raw().unwrap();
        section.jump_sequence().unwrap();
        section.jump_sequence().unwrap();
        let block = section.read_block().unwrap().unwrap();
        assert_eq!(
            decode(&block.sequence, blocks[2].sequence.len()),
            blocks[2].sequence
        );
        assert_eq!(section.remaining_blocks(), 7);
        section.close().unwrap();
        assert_eq!(reader.tell(), reader.end_position());
    }
    // A dropped section is jumped too
    {
        let mut reader = reader_of(&bytes).unwrap();
        reader.read_values().unwrap();
        let mut section = reader.open_raw().unwrap();
        section.jump_sequence().unwrap();
        drop(section);
        assert_eq!(reader.tell(), reader.end_position());
    }
    // jump_section
    {
        let mut reader = reader_of(&bytes).unwrap();
        reader.read_values().unwrap();
        let mut section = reader.open_raw().unwrap();
        section.jump_section().unwrap();
        assert!(section.jump_sequence().is_err());
        drop(section);
        assert_eq!(reader.tell(), reader.end_position());
    }
    // jump_next_section
    {
        let mut reader = reader_of(&bytes).unwrap();
        assert!(!reader.jump_next_section().unwrap()); // a 'v' section is not skipped
        reader.read_values().unwrap();
        assert!(reader.jump_next_section().unwrap());
        assert_eq!(reader.tell(), reader.end_position());
        assert!(!reader.jump_next_section().unwrap());
    }
    // skip_section goes over any section
    {
        let mut reader = reader_of(&bytes).unwrap();
        assert_eq!(reader.skip_section().unwrap(), Some(SectionType::Values));
        assert_eq!(reader.var("k"), Some(11));
        assert_eq!(reader.skip_section().unwrap(), Some(SectionType::Raw));
        assert_eq!(reader.skip_section().unwrap(), None);
    }
}

#[test]
fn raw_section_without_block() {
    let (bytes, beginning) = write_raw_file(7, 10, 0, &[]);
    assert_eq!(
        bytes[beginning as usize..bytes.len() - 3],
        raw_section_bytes(&[], 7, 10, &Encoding::default())[..]
    );
    let mut reader = reader_of(&bytes).unwrap();
    reader.read_values().unwrap();
    let section = reader.open_raw().unwrap();
    assert_eq!(section.nb_blocks(), 0);
    section.close().unwrap();
    assert_eq!(reader.tell(), reader.end_position());
}

#[test]
fn raw_missing_variables_are_detected() {
    let mut writer = new_writer();
    write_vars(&mut writer, &[("k", 5), ("max", 5)]);
    assert!(matches!(
        writer.raw_section(),
        Err(Error::MissingVariable("data_size"))
    ));
    write_vars(&mut writer, &[("max", 5), ("data_size", 5)]);
    assert!(matches!(
        writer.raw_section(),
        Err(Error::MissingVariable("k"))
    ));
    write_vars(&mut writer, &[("k", 5), ("data_size", 5)]);
    assert!(matches!(
        writer.raw_section(),
        Err(Error::MissingVariable("max"))
    ));

    // Same when reading a file where a variable is missing
    let mut bytes = default_header();
    bytes.extend(variables_bytes(&[("data_size", 1), ("k", 5)]));
    bytes.extend_from_slice(&[b'r', 0, 0, 0, 0, 0, 0, 0, 0]);
    bytes.extend_from_slice(b"KFF");
    let mut reader = reader_of(&bytes).unwrap();
    reader.read_values().unwrap();
    assert!(matches!(
        reader.open_raw(),
        Err(Error::MissingVariable("max"))
    ));
}

#[test]
fn raw_invalid_blocks_are_rejected() {
    let mut writer = new_writer();
    write_vars(&mut writer, &[("k", 5), ("max", 3), ("data_size", 1)]);
    let mut section = writer.raw_section().unwrap();
    let seq = encode(b"ACGTACGTAC");
    // Too short for a kmer
    assert!(section.write_compacted_sequence(&seq, 4, &[1]).is_err());
    // More kmers than max
    assert!(section.write_compacted_sequence(&seq, 8, &[1; 4]).is_err());
    // Missing data
    assert!(matches!(
        section.write_compacted_sequence(&seq, 7, &[1, 2]),
        Err(Error::BufferTooSmall {
            needed: 3,
            available: 2
        })
    ));
    // Missing nucleotides
    assert!(section
        .write_compacted_sequence(&seq[..1], 7, &[1; 3])
        .is_err());
    section.write_compacted_sequence(&seq, 7, &[1; 3]).unwrap();
    assert_eq!(section.nb_blocks(), 1);

    // n == max does not fit when max is a power of 256 (ceil(log2(max)) bits)
    let mut writer = new_writer();
    write_vars(&mut writer, &[("k", 1), ("max", 256), ("data_size", 0)]);
    let mut section = writer.raw_section().unwrap();
    let seq = vec![0; 64];
    assert!(section.write_compacted_sequence(&seq, 256, &[]).is_err());
    section.write_compacted_sequence(&seq, 255, &[]).unwrap();
}

#[test]
fn raw_a_too_small_buffer_does_not_consume_the_block() {
    let mut rnd = Random(6);
    let blocks = random_blocks(&mut rnd, 9, 30, 2, 5);
    let (bytes, _) = write_raw_file(9, 30, 2, &blocks);
    let mut reader = reader_of(&bytes).unwrap();
    reader.read_values().unwrap();
    let mut section = reader.open_raw().unwrap();
    let mut small = vec![0; 1];
    let mut data = vec![0; section.max_data_bytes()];
    assert!(matches!(
        section.read_compacted_sequence(&mut small, &mut data),
        Err(Error::BufferTooSmall { .. })
    ));
    assert_eq!(section.remaining_blocks(), 5);
    let block = section.read_block().unwrap().unwrap();
    assert_eq!(
        decode(&block.sequence, blocks[0].sequence.len()),
        blocks[0].sequence
    );
}

#[test]
fn raw_copy_to_another_file() {
    let mut rnd = Random(5);
    let blocks = random_blocks(&mut rnd, 9, 40, 3, 50);
    let (source, _) = write_raw_file(9, 40, 3, &blocks);

    let mut reader = reader_of(&source).unwrap();
    let mut copy = new_writer();
    copy.set_indexation(false);
    let vars = reader.read_values().unwrap();
    copy.write_values(&vars).unwrap();
    reader.open_raw().unwrap().copy_to(&mut copy).unwrap();
    assert_eq!(bytes_of(copy), source);
}

#[test]
fn raw_the_readme_example() {
    // k=10, max=255, data_size=1, encoding A=0 C=2 G=3 T=1
    let encoding = Encoding::new(0, 2, 3, 1).unwrap();
    let mut writer = new_writer();
    writer.set_indexation(false);
    writer.set_encoding(encoding).unwrap();
    write_vars(&mut writer, &[("k", 10), ("max", 255), ("data_size", 1)]);
    let mut section = writer.raw_section().unwrap();
    let beginning = section.beginning() as usize;
    for (sequence, data) in [
        (&b"ACTAAACTGATT"[..], &[32, 47, 1][..]),
        (b"AAACTGATCG", &[12]),
        (b"CTAAACTGATT", &[1, 47]),
    ] {
        section
            .write_compacted_sequence(
                &encoding.encode(sequence).unwrap(),
                sequence.len() as u64,
                data,
            )
            .unwrap();
    }
    section.close().unwrap();
    let bytes = bytes_of(writer);
    assert_eq!(bytes[5], 0x2d);
    // Bytes given by the README
    let expected: Bytes = vec![
        b'r', 0, 0, 0, 0, 0, 0, 0, 3, //
        3, 0x24, 0x09, 0xc5, 0x20, 0x2f, 0x01, //
        1, 0x00, 0x27, 0x1b, 0x0c, //
        2, 0x24, 0x09, 0xc5, 0x01, 0x2f,
    ];
    assert_eq!(bytes[beginning..bytes.len() - 3], expected[..]);
}

// ----- Minimizer sections -----

struct MiniFile {
    bytes: Bytes,
    beginning: u64,
    minimizer: Vec<u8>,
    blocks: Vec<Block>,
}

/// Write a file with one 'v' and one 'm' section. Odd blocks are written with
/// write_compacted_sequence (minimizer removed by the API), even blocks with
/// write_compacted_sequence_without_mini.
fn write_mini_file(
    rnd: &mut Random,
    k: u64,
    m: u64,
    max: u64,
    data_size: u64,
    nb_blocks: usize,
) -> MiniFile {
    let minimizer = rnd.sequence(m);
    let mut writer = new_writer();
    writer.set_indexation(false);
    write_vars(
        &mut writer,
        &[("k", k), ("m", m), ("max", max), ("data_size", data_size)],
    );
    let mut section = writer.minimizer_section().unwrap();
    let beginning = section.beginning();
    section.write_minimizer(&encode(&minimizer)).unwrap();
    let mut blocks = Vec::new();
    for i in 0..nb_blocks {
        let nb_kmers = if i == 0 { max } else { rnd.between(1, max) };
        let block = random_block(rnd, k, nb_kmers, data_size, &minimizer);
        if i % 2 == 1 {
            section
                .write_compacted_sequence(
                    &encode(&block.sequence),
                    block.sequence.len() as u64,
                    block.mini_pos,
                    &block.data,
                )
                .unwrap();
        } else {
            let without = block.without_mini(m);
            section
                .write_compacted_sequence_without_mini(
                    &encode(&without),
                    without.len() as u64,
                    block.mini_pos,
                    &block.data,
                )
                .unwrap();
        }
        blocks.push(block);
    }
    assert_eq!(section.nb_blocks(), nb_blocks as u64);
    section.close().unwrap();
    MiniFile {
        bytes: bytes_of(writer),
        beginning,
        minimizer,
        blocks,
    }
}

#[test]
fn minimizer_section_bytes_and_round_trip_for_many_k_m_max_and_data_sizes() {
    let mut rnd = Random(10);
    let encoding = Encoding::default();
    for m in [1, 2, 3, 4, 5, 7, 8, 9, 13, 16, 17, 31] {
        for delta in [0, 1, 2, 3, 4, 5, 9, 20] {
            for max in [1, 3, 7, 200] {
                for data_size in [0, 1, 3] {
                    let k = m + delta;
                    let mf = write_mini_file(&mut rnd, k, m, max, data_size, 12);
                    // Exact bytes
                    assert_eq!(
                        mf.bytes[mf.beginning as usize..mf.bytes.len() - 3],
                        minimizer_section_bytes(&mf.minimizer, &mf.blocks, k, max, &encoding)[..]
                    );

                    // 1 - Read without the minimizer
                    let mut reader = reader_of(&mf.bytes).unwrap();
                    reader.read_values().unwrap();
                    let mut section = reader.open_minimizer().unwrap();
                    assert_eq!(section.k(), k);
                    assert_eq!(section.m(), m);
                    assert_eq!(decode(section.minimizer(), m as usize), mf.minimizer);
                    let mut seq = vec![0; section.max_sequence_bytes()];
                    let mut data = vec![0; section.max_data_bytes()];
                    for block in &mf.blocks {
                        let (nb_kmers, mini_pos) = section
                            .read_compacted_sequence_without_mini(&mut seq, &mut data)
                            .unwrap();
                        let without = block.without_mini(m);
                        assert_eq!(nb_kmers, block.sequence.len() as u64 - k + 1);
                        assert_eq!(mini_pos, block.mini_pos);
                        assert_eq!(seq[..(without.len() + 3) / 4], encode(&without)[..]);
                        assert_eq!(data[..block.data.len()], block.data[..]);
                        // The minimizer is inserted back
                        section.add_minimizer(nb_kmers, &mut seq, mini_pos).unwrap();
                        let seq_bytes = (block.sequence.len() + 3) / 4;
                        assert_eq!(seq[..seq_bytes], encode(&block.sequence)[..]);
                    }
                    section.close().unwrap();

                    // 2 - Read with the minimizer inserted back, padding bits set to 0
                    let mut reader = reader_of(&mf.bytes).unwrap();
                    reader.read_values().unwrap();
                    let mut section = reader.open_minimizer().unwrap();
                    for block in &mf.blocks {
                        let nb_kmers = section
                            .read_compacted_sequence(&mut seq, &mut data)
                            .unwrap();
                        let seq_bytes = (block.sequence.len() + 3) / 4;
                        assert_eq!(nb_kmers, block.sequence.len() as u64 - k + 1);
                        assert_eq!(seq[..seq_bytes], encode(&block.sequence)[..]);
                        assert_eq!(data[..block.data.len()], block.data[..]);
                    }
                    section.close().unwrap();

                    // 3 - Sequence and data in the same buffer
                    let mut reader = reader_of(&mf.bytes).unwrap();
                    reader.read_values().unwrap();
                    let mut section = reader.open_minimizer().unwrap();
                    let mut buffer = vec![0; seq.len() + data.len()];
                    for block in &mf.blocks {
                        let nb_kmers = section.read_compacted_sequence_joined(&mut buffer).unwrap();
                        let seq_bytes = (block.sequence.len() + 3) / 4;
                        assert_eq!(nb_kmers, block.sequence.len() as u64 - k + 1);
                        assert_eq!(buffer[..seq_bytes], encode(&block.sequence)[..]);
                        assert_eq!(
                            buffer[seq_bytes..seq_bytes + block.data.len()],
                            block.data[..]
                        );
                    }
                    section.close().unwrap();

                    // 4 - Owned blocks
                    let mut reader = reader_of(&mf.bytes).unwrap();
                    reader.read_values().unwrap();
                    let mut section = reader.open_minimizer().unwrap();
                    for block in &mf.blocks {
                        let read = section.read_block().unwrap().unwrap();
                        assert_eq!(read.sequence, encode(&block.sequence));
                        assert_eq!(read.data, block.data);
                    }
                    assert!(section.read_block().unwrap().is_none());
                }
            }
        }
    }
}

#[test]
fn minimizer_width_of_the_minimizer_position() {
    let cases: &[(u64, u64)] = &[
        (1, 1),
        (2, 1),
        (3, 255),
        (10, 247),
        (10, 248),
        (31, 226),
        (31, 227),
        (31, 65506),
        (31, 65507),
    ];
    let mut rnd = Random(11);
    let encoding = Encoding::default();
    for (k, max) in cases {
        let mf = write_mini_file(&mut rnd, *k, 1, *max, 0, 1);
        assert_eq!(
            mf.bytes[mf.beginning as usize..mf.bytes.len() - 3],
            minimizer_section_bytes(&mf.minimizer, &mf.blocks, *k, *max, &encoding)[..]
        );
    }
}

#[test]
fn minimizer_jump_over_blocks() {
    let mut rnd = Random(12);
    let mf = write_mini_file(&mut rnd, 21, 7, 30, 2, 10);
    let mut reader = reader_of(&mf.bytes).unwrap();
    reader.read_values().unwrap();
    let mut section = reader.open_minimizer().unwrap();
    section.jump_sequence().unwrap();
    let block = section.read_block().unwrap().unwrap();
    assert_eq!(
        decode(&block.sequence, mf.blocks[1].sequence.len()),
        mf.blocks[1].sequence
    );
    section.close().unwrap();
    assert_eq!(reader.tell(), reader.end_position());
}

#[test]
fn minimizer_copy_to_another_file() {
    let mut rnd = Random(13);
    let mf = write_mini_file(&mut rnd, 25, 9, 40, 1, 30);
    let mut reader = reader_of(&mf.bytes).unwrap();
    let mut copy = new_writer();
    copy.set_indexation(false);
    let vars = reader.read_values().unwrap();
    copy.write_values(&vars).unwrap();
    reader.open_minimizer().unwrap().copy_to(&mut copy).unwrap();
    assert_eq!(bytes_of(copy), mf.bytes);
}

#[test]
fn minimizer_missing_m_is_detected() {
    let mut writer = new_writer();
    write_vars(&mut writer, &[("k", 5), ("max", 5), ("data_size", 0)]);
    assert!(matches!(
        writer.minimizer_section(),
        Err(Error::MissingVariable("m"))
    ));
}

#[test]
fn minimizer_invalid_minimizer_positions_are_rejected() {
    let mut writer = new_writer();
    write_vars(
        &mut writer,
        &[("k", 5), ("m", 3), ("max", 10), ("data_size", 0)],
    );
    let mut section = writer.minimizer_section().unwrap();
    let seq = encode(b"ACGTACG");
    assert!(section.write_compacted_sequence(&seq, 7, 5, &[]).is_err());
    assert!(section
        .write_compacted_sequence_without_mini(&seq, 4, 5, &[])
        .is_err());
    section.write_compacted_sequence(&seq, 7, 4, &[]).unwrap();
    section
        .write_compacted_sequence_without_mini(&seq, 4, 4, &[])
        .unwrap();
}

#[test]
fn minimizer_the_readme_example() {
    // k=10, m=8, max=255, data_size=1, encoding A=0 C=2 G=3 T=1
    let encoding = Encoding::new(0, 2, 3, 1).unwrap();
    let mut writer = new_writer();
    writer.set_indexation(false);
    writer.set_encoding(encoding).unwrap();
    write_vars(
        &mut writer,
        &[("k", 10), ("m", 8), ("max", 255), ("data_size", 1)],
    );
    let mut section = writer.minimizer_section().unwrap();
    let beginning = section.beginning() as usize;
    section
        .write_minimizer(&encoding.encode(b"AAACTGAT").unwrap())
        .unwrap();
    for (sequence, mini_pos, data) in [
        (&b"ACTAAACTGATT"[..], 3, &[32, 47, 1][..]),
        (b"AAACTGATCG", 0, &[12]),
        (b"CTAAACTGATT", 2, &[1, 47]),
    ] {
        section
            .write_compacted_sequence(
                &encoding.encode(sequence).unwrap(),
                sequence.len() as u64,
                mini_pos,
                data,
            )
            .unwrap();
    }
    section.close().unwrap();
    let bytes = bytes_of(writer);
    // Same blocks as the README, except that the minimizer position uses
    // ceil(log2(k + max - 1)) = 9 bits, i.e. 2 bytes (the README draws 1 byte).
    let expected: Bytes = vec![
        b'm', 0x02, 0x71, 0, 0, 0, 0, 0, 0, 0, 3, //
        3, 0, 3, 0x25, 0x20, 0x2f, 0x01, //
        1, 0, 0, 0x0b, 0x0c, //
        2, 0, 2, 0x25, 0x01, 0x2f,
    ];
    assert_eq!(bytes[beginning..bytes.len() - 3], expected[..]);

    // The full sequences are rebuilt when reading
    let mut reader = reader_of(&bytes).unwrap();
    reader.read_values().unwrap();
    let mut section = reader.open_minimizer().unwrap();
    for sequence in [&b"ACTAAACTGATT"[..], b"AAACTGATCG", b"CTAAACTGATT"] {
        let block = section.read_block().unwrap().unwrap();
        assert_eq!(encoding.decode(&block.sequence, sequence.len()), sequence);
    }
}

// ----- Index -----

/// Writes v, r, v, m, r sections and returns the bytes and the beginning of each section
fn write_indexed_file(indexed: bool) -> (Bytes, Vec<(u64, SectionType)>) {
    let mut sections = Vec::new();
    let mut rnd = Random(20);
    let mut writer = new_writer();
    writer.set_indexation(indexed);

    let position = write_vars(&mut writer, &[("k", 15), ("max", 10), ("data_size", 1)]);
    sections.push((position, SectionType::Values));
    let mut raw = writer.raw_section().unwrap();
    sections.push((raw.beginning(), SectionType::Raw));
    for i in 0..5 {
        let block = random_block(&mut rnd, 15, 1 + i, 1, b"");
        raw.write_compacted_sequence(
            &encode(&block.sequence),
            block.sequence.len() as u64,
            &block.data,
        )
        .unwrap();
    }
    raw.close().unwrap();

    let position = write_vars(
        &mut writer,
        &[("k", 15), ("m", 6), ("max", 10), ("data_size", 1)],
    );
    sections.push((position, SectionType::Values));
    let mut mini = writer.minimizer_section().unwrap();
    sections.push((mini.beginning(), SectionType::Minimizer));
    mini.write_minimizer(&encode(b"ACGTAC")).unwrap();
    for i in 0..5 {
        let block = random_block(&mut rnd, 15, 1 + i, 1, b"ACGTAC");
        mini.write_compacted_sequence(
            &encode(&block.sequence),
            block.sequence.len() as u64,
            block.mini_pos,
            &block.data,
        )
        .unwrap();
    }
    mini.close().unwrap();

    let raw = writer.raw_section().unwrap();
    sections.push((raw.beginning(), SectionType::Raw));
    raw.close().unwrap();

    if indexed {
        assert_eq!(
            writer
                .section_positions()
                .iter()
                .map(|(p, t)| (*p, *t))
                .collect::<Vec<_>>(),
            sections
        );
    }
    (bytes_of(writer), sections)
}

#[test]
fn index_the_footer_index_lists_every_section() {
    let (bytes, sections) = write_indexed_file(true);
    let reader = reader_of(&bytes).unwrap();
    assert!(reader.footer().is_some());
    assert!(reader.is_indexed());
    assert_eq!(reader.index().len(), 1);
    let index = &reader.index()[0];
    assert_eq!(index.next_index, 0);
    assert_eq!(index.entries.len(), sections.len());

    // Positions are relative to the end of the index section
    let index_position = *reader.footer().unwrap().get("first_index").unwrap();
    assert_eq!(index.beginning, index_position);
    let index_end = index_position + 17 + 9 * sections.len() as u64;
    assert_eq!(index.end(), index_end);
    assert_eq!(bytes[index_position as usize], b'i');
    assert!(index.entries.keys().all(|position| *position < 0));
    assert_eq!(
        reader.indexed_sections().into_iter().collect::<Vec<_>>(),
        sections
    );
    for (position, section_type) in &sections {
        assert_eq!(bytes[*position as usize], section_type.as_byte());
    }

    // Exact bytes of the index section (entries sorted by relative position)
    let mut expected = vec![b'i'];
    expected.extend(big_endian(sections.len() as u64, 8));
    for (position, section_type) in &sections {
        expected.push(section_type.as_byte());
        expected.extend(big_endian((*position as i64 - index_end as i64) as u64, 8));
    }
    expected.extend(big_endian(0, 8));
    assert_eq!(
        bytes[index_position as usize..index_end as usize],
        expected[..]
    );

    // Then the footer and the signature
    let mut footer = variables_bytes(&[("first_index", index_position), ("footer_size", 49)]);
    footer.extend_from_slice(b"KFF");
    assert_eq!(bytes[index_end as usize..], footer[..]);
}

#[test]
fn index_no_footer_when_the_indexation_is_disabled() {
    let (bytes, sections) = write_indexed_file(false);
    // The last section (an empty raw section) is directly followed by the signature
    assert_eq!(bytes.len() as u64, sections.last().unwrap().0 + 9 + 3);
    let mut reader = reader_of(&bytes).unwrap();
    assert!(reader.footer().is_none());
    assert!(reader.index().is_empty());
    assert_eq!(
        reader.read_section_type().unwrap(),
        Some(SectionType::Values)
    );
    assert_eq!(reader.tell(), 12);
}

#[test]
fn index_chained_index_sections_are_all_loaded() {
    let mut writer = new_writer();
    writer.set_indexation(false);
    // First index at the beginning of the file, pointing to a second one
    let mut first = writer.index_section().unwrap();
    first.register_section(SectionType::Values, 0);
    // first index: 1 + 8 + 9 + 8 = 26 bytes, then a 'v' section of 9 bytes
    first.set_next_index(9);
    first.close().unwrap();
    write_vars(&mut writer, &[]);
    // Second index at 47, its end is at 47 + 26 = 73 and the 'v' section at 38
    let mut second = writer.index_section().unwrap();
    assert_eq!(second.beginning(), 47);
    second.register_section(SectionType::Values, 38 - 73);
    second.close().unwrap();

    let mut reader = reader_of(&bytes_of(writer)).unwrap();
    assert_eq!(reader.index().len(), 2);
    assert_eq!(reader.index()[0].next_index, 9);
    assert_eq!(
        reader.index()[0].entries.get(&0),
        Some(&SectionType::Values)
    );
    assert_eq!(reader.index()[1].next_index, 0);
    assert_eq!(
        reader.index()[1].entries.get(&(38 - 73)),
        Some(&SectionType::Values)
    );
    // Reading still starts at the first section
    assert_eq!(
        reader.read_section_type().unwrap(),
        Some(SectionType::Index)
    );
}

#[test]
fn index_a_loop_in_the_chain_is_detected() {
    let mut writer = new_writer();
    writer.set_indexation(false);
    let mut first = writer.index_section().unwrap();
    // Points to itself: 26 bytes before its end
    first.set_next_index(-26);
    first.close().unwrap();
    assert!(matches!(
        reader_of(&bytes_of(writer)),
        Err(Error::InvalidFooter(_))
    ));
}

#[test]
fn index_index_collisions_are_detected() {
    let mut bytes = default_header();
    bytes.push(b'i');
    bytes.extend(big_endian(2, 8));
    for _ in 0..2 {
        bytes.push(b'v');
        bytes.extend(big_endian(0, 8));
    }
    bytes.extend(big_endian(0, 8));
    bytes.extend_from_slice(b"KFF");
    assert!(matches!(reader_of(&bytes), Err(Error::IndexCollision)));
}

#[test]
fn index_sections_can_be_built_from_their_type() {
    let (bytes, _) = write_indexed_file(true);
    let mut reader = reader_of(&bytes).unwrap();
    let mut types = String::new();
    while let Some(section) = reader.read_section().unwrap() {
        types.push(match section {
            Section::Raw(section) => {
                section.close().unwrap();
                'r'
            }
            Section::Minimizer(section) => {
                section.close().unwrap();
                'm'
            }
            Section::Values(_) => 'v',
            Section::Index(_) => 'i',
        });
    }
    assert_eq!(types, "vrvmriv");
}

// ----- Files -----

/// Writes a raw section large enough to go through the write buffer several times
fn write_large_file<W>(writer: &mut Writer<W>, nb_blocks: usize) -> Vec<Block>
where
    W: std::io::Write + std::io::Seek,
{
    let mut rnd = Random(40);
    let mut blocks = Vec::new();
    write_vars(writer, &[("k", 31), ("max", 100), ("data_size", 1)]);
    let mut raw = writer.raw_section().unwrap();
    for _ in 0..nb_blocks {
        let nb_kmers = rnd.between(1, 100);
        let block = random_block(&mut rnd, 31, nb_kmers, 1, b"");
        raw.write_compacted_sequence(
            &encode(&block.sequence),
            block.sequence.len() as u64,
            &block.data,
        )
        .unwrap();
        blocks.push(block);
    }
    raw.close().unwrap();
    blocks
}

#[test]
fn file_large_files_are_written_through_the_buffer_and_read_back() {
    let dir = tempfile::tempdir().unwrap();
    let name = dir.path().join("large.kff");
    let mut writer = Writer::create(&name).unwrap();
    let blocks = write_large_file(&mut writer, 50000);
    writer.close().unwrap();

    let bytes = std::fs::read(&name).unwrap();
    assert!(bytes.len() > 3 * (1 << 20));
    let mut expected = default_header();
    expected.extend(variables_bytes(&[
        ("data_size", 1),
        ("k", 31),
        ("max", 100),
    ]));
    let raw_position = expected.len();
    expected.extend(raw_section_bytes(&blocks, 31, 100, &Encoding::default()));
    assert_eq!(bytes[..expected.len()], expected[..]);

    let mut reader = Reader::open(&name).unwrap();
    assert_eq!(
        reader.footer().unwrap().get("first_index"),
        Some(&(expected.len() as u64))
    );
    reader.read_values().unwrap();
    let mut section = reader.open_raw().unwrap();
    assert_eq!(section.beginning(), raw_position as u64);
    for block in &blocks {
        let read = section.read_block().unwrap().unwrap();
        assert_eq!(decode(&read.sequence, block.sequence.len()), block.sequence);
        assert_eq!(read.data, block.data);
    }
}

#[test]
fn file_releasing_the_file_does_not_change_it() {
    let dir = tempfile::tempdir().unwrap();
    let write = |name: &str, release: bool| {
        let path = dir.path().join(name);
        let mut writer = Writer::create(&path).unwrap();
        for section in 0..45 {
            if release {
                writer.release().unwrap();
                assert!(writer.is_released());
            }
            if section == 0 {
                writer.set_uniqueness(true).unwrap();
            }
            write_large_file(&mut writer, 1000);
        }
        writer.close().unwrap();
        std::fs::read(path).unwrap()
    };
    let normal = write("normal.kff", false);
    let released = write("released.kff", true);
    assert!(normal.len() > 3 * (1 << 20));
    assert!(normal == released);

    let mut stream = new_writer();
    assert!(matches!(stream.release(), Err(Error::CannotRelease)));
}

#[test]
fn file_write_at_on_the_stream_the_buffer_and_across_both() {
    let content = Random(41).bytes(3 * (1 << 20) + 12345);
    let mut expected = content.clone();
    let mut writer = new_writer();
    writer.set_indexation(false);
    for chunk in content.chunks(1000) {
        writer.write_bytes(chunk).unwrap();
    }
    let end = writer.tell();
    assert_eq!(end, 8 + content.len() as u64);
    // Patches of 16 bytes every 10 bytes: some are in the part already written in the stream,
    // some in the buffer and some across both, wherever the buffer was flushed
    let patch: Bytes = (1..=16).collect();
    for position in (8..end - 16).step_by(10) {
        writer.write_at(&patch, position).unwrap();
        expected[(position - 8) as usize..(position + 8) as usize].copy_from_slice(&patch);
    }
    // Writing after the end extends the file, writing out of the file fails
    assert!(writer.write_at(&patch, end + 1).is_err());
    writer.write_at(&patch, end - 8).unwrap();
    expected.truncate(expected.len() - 8);
    expected.extend_from_slice(&patch);
    assert_eq!(writer.tell(), end + 8);

    let bytes = bytes_of(writer);
    // No section was opened: the size of the metadata is written at the closing
    assert!(bytes[8..8 + expected.len()] == expected[..]);
    assert_eq!(bytes.len(), 8 + expected.len() + 4 + 3);
}

#[test]
fn file_a_discarded_file_is_deleted() {
    let dir = tempfile::tempdir().unwrap();
    let name = dir.path().join("discarded.kff");
    let mut writer = Writer::create(&name).unwrap();
    write_large_file(&mut writer, 100000);
    assert!(name.exists());
    writer.discard().unwrap();
    assert!(!name.exists());
}

#[test]
fn file_a_dropped_writer_is_closed() {
    let dir = tempfile::tempdir().unwrap();
    let name = dir.path().join("dropped.kff");
    {
        let mut writer = Writer::create(&name).unwrap();
        let mut section = writer.values_section().unwrap();
        section.write_var("k", 21).unwrap();
    }
    let mut reader = Reader::open(&name).unwrap();
    assert!(reader.footer().is_some());
    assert_eq!(reader.read_values().unwrap().get("k"), Some(&21));
}

#[test]
fn file_positions_jumps_and_reading_past_the_end() {
    let mut writer = new_writer();
    write_vars(&mut writer, &[("k", 21)]);
    let bytes = bytes_of(writer);
    let mut reader = reader_of(&bytes).unwrap();
    assert_eq!(reader.end_position(), bytes.len() as u64 - 3);
    assert_eq!(reader.file_size(), bytes.len() as u64);
    reader.jump_to(bytes.len() as u64 - 3).unwrap();
    let mut signature = [0; 3];
    reader.read_bytes(&mut signature).unwrap();
    assert_eq!(&signature, b"KFF");
    reader.jump(-10).unwrap();
    assert_eq!(reader.tell(), bytes.len() as u64 - 10);
    let mut buffer = [0; 16];
    assert!(reader.read_bytes(&mut buffer).is_err());
    assert!(reader.jump_to(bytes.len() as u64 + 1).is_err());
    assert!(reader.jump(-(bytes.len() as i64)).is_err());
}

// ----- High level reader -----

struct Kmer {
    sequence: Vec<u8>,
    data: Bytes,
    k: u64,
}

/// File with raw and minimizer sections that change k, max and data_size
fn write_mixed_file(encoding: Encoding) -> (Bytes, Vec<Kmer>) {
    let mut kmers = Vec::new();
    let mut rnd = Random(30);
    let mut writer = new_writer();
    writer.set_encoding(encoding).unwrap();
    let add_kmers = |kmers: &mut Vec<Kmer>, block: &Block, k: u64, data_size: u64| {
        for i in 0..=(block.sequence.len() - k as usize) {
            kmers.push(Kmer {
                sequence: block.sequence[i..i + k as usize].to_vec(),
                data: block.data[i * data_size as usize..(i + 1) * data_size as usize].to_vec(),
                k,
            });
        }
    };
    // k, m, max, data_size, type
    let params: &[(u64, u64, u64, u64, char)] = &[
        (5, 0, 4, 1, 'r'),
        (7, 0, 3, 0, 'r'),
        (6, 3, 5, 2, 'm'),
        (33, 0, 1, 3, 'r'),
        (6, 3, 5, 2, 'r'),
        (31, 11, 40, 1, 'm'),
        (64, 0, 9, 1, 'r'),
    ];
    for (k, m, max, data_size, section_type) in params {
        let mut vars = vec![("k", *k), ("max", *max), ("data_size", *data_size)];
        if *section_type == 'm' {
            vars.push(("m", *m));
        }
        write_vars(&mut writer, &vars);
        if *section_type == 'r' {
            let mut raw = writer.raw_section().unwrap();
            for _ in 0..25 {
                let nb_kmers = rnd.between(1, *max);
                let block = random_block(&mut rnd, *k, nb_kmers, *data_size, b"");
                raw.write_compacted_sequence(
                    &encoding.encode(&block.sequence).unwrap(),
                    block.sequence.len() as u64,
                    &block.data,
                )
                .unwrap();
                add_kmers(&mut kmers, &block, *k, *data_size);
            }
            raw.close().unwrap();
            writer.raw_section().unwrap().close().unwrap();
        } else {
            let minimizer = rnd.sequence(*m);
            let mut section = writer.minimizer_section().unwrap();
            section
                .write_minimizer(&encoding.encode(&minimizer).unwrap())
                .unwrap();
            for _ in 0..25 {
                let nb_kmers = rnd.between(1, *max);
                let block = random_block(&mut rnd, *k, nb_kmers, *data_size, &minimizer);
                section
                    .write_compacted_sequence(
                        &encoding.encode(&block.sequence).unwrap(),
                        block.sequence.len() as u64,
                        block.mini_pos,
                        &block.data,
                    )
                    .unwrap();
                add_kmers(&mut kmers, &block, *k, *data_size);
            }
            section.close().unwrap();
        }
    }
    (bytes_of(writer), kmers)
}

#[test]
fn reader_next_kmer_enumerates_every_kmer_and_its_data() {
    let encoding = Encoding::new(2, 0, 3, 1).unwrap();
    for indexed in [true, false] {
        let (bytes, expected) = write_mixed_file(encoding);
        let bytes = if indexed {
            bytes
        } else {
            // Same file without index nor footer: the last kmers are still read
            let mut reader = reader_of(&bytes).unwrap();
            let mut writer = new_writer();
            writer.set_indexation(false);
            writer.set_encoding(encoding).unwrap();
            let footer = reader.footer_position();
            while let Some(section_type) = reader.read_section_type().unwrap() {
                let position = reader.tell();
                match section_type {
                    SectionType::Values => {
                        let vars = reader.read_values().unwrap();
                        if Some(position) != footer {
                            writer.write_values(&vars).unwrap();
                        }
                    }
                    SectionType::Raw | SectionType::Minimizer => {
                        let mut section = reader.open_block_section().unwrap().unwrap();
                        section.copy_to(&mut writer).unwrap();
                    }
                    _ => {
                        reader.skip_section().unwrap();
                    }
                }
            }
            bytes_of(writer)
        };

        let mut reader = KmerReader::new(reader_of(&bytes).unwrap());
        assert_eq!(reader.encoding(), encoding);
        let mut idx = 0;
        while reader.has_next().unwrap() {
            let kmer = reader.next_kmer().unwrap().unwrap();
            assert!(idx < expected.len());
            assert_eq!(kmer.k, expected[idx].k);
            assert_eq!(kmer.to_ascii(), expected[idx].sequence);
            assert_eq!(kmer.data, &expected[idx].data[..]);
            assert_eq!(reader.k(), expected[idx].k);
            idx += 1;
        }
        assert_eq!(idx, expected.len());
        assert!(reader.next_kmer().unwrap().is_none());
    }
}

#[test]
fn reader_next_block_returns_whole_blocks() {
    let mut rnd = Random(31);
    let blocks = random_blocks(&mut rnd, 13, 50, 2, 40);
    let mut writer = new_writer();
    write_vars(&mut writer, &[("k", 13), ("max", 50), ("data_size", 2)]);
    let mut raw = writer.raw_section().unwrap();
    for block in &blocks {
        raw.write_compacted_sequence(
            &encode(&block.sequence),
            block.sequence.len() as u64,
            &block.data,
        )
        .unwrap();
    }
    raw.close().unwrap();
    let bytes = bytes_of(writer);

    let mut reader = KmerReader::new(reader_of(&bytes).unwrap());
    assert!(reader.has_next().unwrap());
    assert_eq!(reader.var("max"), Some(50));
    assert_eq!(reader.var("data_size"), Some(2));
    for block in &blocks {
        let read = reader.next_block().unwrap().unwrap();
        assert_eq!(read.nb_kmers, block.sequence.len() as u64 - 12);
        assert_eq!(decode(&read.sequence, block.sequence.len()), block.sequence);
        assert_eq!(read.data, block.data);
    }
    assert!(!reader.has_next().unwrap());
    assert!(reader.next_block().unwrap().is_none());

    // A block partially read with next_kmer is skipped by next_block
    let mut reader = KmerReader::new(reader_of(&bytes).unwrap());
    reader.next_kmer().unwrap().unwrap();
    let read = reader.next_block().unwrap().unwrap();
    assert_eq!(
        decode(&read.sequence, blocks[1].sequence.len()),
        blocks[1].sequence
    );
}

#[test]
fn reader_file_without_kmer() {
    let mut writer = new_writer();
    write_vars(&mut writer, &[("k", 13)]);
    let mut reader = KmerReader::new(reader_of(&bytes_of(writer)).unwrap());
    assert!(!reader.has_next().unwrap());
    assert!(reader.next_kmer().unwrap().is_none());
    assert_eq!(reader.reader().footer().unwrap().len(), 2);
}
