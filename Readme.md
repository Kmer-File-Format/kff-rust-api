# Kmer File Format Rust parser 🧬 💻

![tests](https://github.com/Kmer-File-Format/kff-rust-api/workflows/tests/badge.svg)
![lints](https://github.com/Kmer-File-Format/kff-rust-api/workflows/lints/badge.svg)
![msrv](https://github.com/Kmer-File-Format/kff-rust-api/workflows/msrv/badge.svg)
[![website](https://github.com/Kmer-File-Format/kff-rust-api/workflows/website/badge.svg)](https://kmer-file-format.github.io/kff-rust-api)
[![dev-doc](https://img.shields.io/badge/dev-doc-blue)](https://kmer-file-format.github.io/kff-rust-api/doc/kff)
[![license](https://img.shields.io/badge/license-AGPL-purple)](https://github.com/Kmer-File-Format/kff-rust-api//blob/main/LICENSE)
[![copier](https://img.shields.io/badge/copier-template-yellow)](https://github.com/natir/copier-rust)

Welcome in the rust library description for Reading or Writing kff files.

The crate offers the features of the [C++ API](https://github.com/Kmer-File-Format/kff-cpp-api)
and reads and writes the same bytes (the tests compare it with files written by the C++ API):
- high level reading: enumerate the kmers or the blocks of a file with their data (`KmerReader`);
- low level reading: header, footer, index and each section, with jumps and copies (`Reader`);
- low level writing: header, values, raw, minimizer and index sections, with the index and the
  footer written automatically (`Writer`).

## Usage

Add this in you Cargo.toml
```toml
kff = { version = "0.10" }
```

```rust
// Write a file
let mut writer = kff::Writer::create("file.kff")?;
let encoding = writer.encoding();
let mut values = writer.values_section()?;
values.write_var("k", 5)?;
values.write_var("max", 10)?;
values.write_var("data_size", 1)?;
values.close()?;
let mut raw = writer.raw_section()?;
raw.write_compacted_sequence(&encoding.encode(b"ACGTACG")?, 7, &[1, 2, 3])?;
raw.close()?;
writer.close()?;

// Read its kmers
let mut reader = kff::KmerReader::open("file.kff")?;
while let Some(kmer) = reader.next_kmer()? {
    println!("{} {:?}", String::from_utf8_lossy(&kmer.to_ascii()), kmer.data);
}
```

- [How to read a KFF file](https://kmer-file-format.github.io/kff-rust-api/read.html)
- [How to write a KFF file](https://kmer-file-format.github.io/kff-rust-api/write.html)
- You could also check [example](https://github.com/Kmer-File-Format/kff-rust-api/tree/main/examples)

## Documentation

An overview of features is [available](https://kmer-file-format.github.io/kff-rust-api/index.html), as well as detailed documentation of the [release version](https://docs.rs/kff/latest/kff/) and the version under [development](https://kmer-file-format.github.io/kff-rust-api/index.html).

## Contributors

- Pierre Marijon
- Lucas Robidou

## Minimum supported Rust version

Currently the minimum supported Rust version is 1.62.
