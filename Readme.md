# Kmer File Format Rust parser 🧬 💻

![tests](https://github.com/Kmer-File-Format/kff-rust-api/workflows/tests/badge.svg)
![lints](https://github.com/Kmer-File-Format/kff-rust-api/workflows/lints/badge.svg)
![msrv](https://github.com/Kmer-File-Format/kff-rust-api/workflows/msrv/badge.svg)
[![website](https://github.com/Kmer-File-Format/kff-rust-api/workflows/website/badge.svg)](https://kmer-file-format.github.io/kff-rust-api)
[![dev-doc](https://img.shields.io/badge/dev-doc-blue)](https://kmer-file-format.github.io/kff-rust-api/doc/kff)
[![license](https://img.shields.io/badge/license-AGPL-purple)](https://github.com/Kmer-File-Format/kff-rust-api//blob/main/LICENSE)
[![copier](https://img.shields.io/badge/copier-template-yellow)](https://github.com/natir/copier-rust)

Welcome in the rust library description for Reading or Writing kff files.
This code is not guaranty to be 100% bug-free, so please submit issues if encounter one.

## Usage

```
kff = { version = "0.9" }
```



## Minimum supported Rust version

Currently the minimum supported Rust version is 1.62.

## Read a file

### Open a file

`kff::Kff` is the main object in the library.
This is the object needed to manipulate a binary kff file.

```rust
let file = kff::Kff::<std::io::BufReader<std::fs::File>>::open(kff_path).expect("could not open kff file");
```

### Read header and encoding

When you open a file to read it, the header (including the encoding) is automatically read.
It is accessible via the file property *header*. You can access the encoding through the *encoding* header.

```rust
let header = file.header();
let encoding = file.header().encoding();
```

### Enumerating kmers from a file

This high-level reader API is made to be very easy to use.
This part of the API allows you to enumerate each pair of kmer/data through the whole file, hiding all the kff data structures.

```rust
let mut file = kff::Kff::<std::io::BufReader<std::fs::File>>::open(kff_path).expect("could not open kff file");
let encoding = *(file.header().encoding());
for kmer in file.kmers() {
    let kmer = kmer.expect("error reading the kmer");
    println!("{}", String::from_utf8(kmer.seq(encoding)).expect("could not parse utf 8"));
}
```

### Enumerating blocks

TODO


### How to know the properties of my kmers?

The values stored in the file (e.g. k or data_size) are accessible through the *values* method.

```rust
let vars: &kff::section::Values = file.values();
let m: Option<&u64> = vars.get("m");
```


## Write a file

### Open/close a file
Creating a kff file requires creating a header.
The header first contains the kff file version (minor and major).
You then have to fill in if the kmers are unique and/or canonical.
Then you have to write the encoding you are using.
Optionally, you can add some data in a vector of bytes.

```rust
const ENCODING: u8 = 0b00011011; // A C T G, 2 bits per letter
[...]
let major_version: u8 = 1;
let minor_version: u8 = 0;
let uniq_kmer: bool = ...;
let canonical_kmer: bool = ...;
let free_block: Vec<u8> = ...;
let header = kff::section::Header::new(major_version, minor_version, ENCODING, uniq_kmer, canonical_kmer, free_block).expect("invalid header");
let file = Kff::create(output_file, header).expect("unable to initiate kff");
```

### String sequences to binary

This crate provides two functions to transform a DNA string into a byte array coding for the binary sequence:
- `from_ascii` for a sequence with associated data (e.g. kmers)
- `seq2bits` for a sequence without any associated data (e.g. minimizers).
You can, of course, use another one of your own.

### Write values

As defined in the standard, writing values is necessary for writing some other sections. Please refer to the standard for more information on which value to write.

```rust
let mut file = Kff::create(/*output_file*/, /*header*/).expect("unable to initiate kff");

let mut values = kff::section::Values::default();
values.insert("k".to_string(), /* k */);
values.insert("m".to_string(), /* m */);
values.insert("ordered".to_string(), /* ordered */);
values.insert("max".to_string(), /* max_nb_kmers */);
values.insert("data_size".to_string(), /* data_size */);

file.write_values(values.clone()).expect("unable to write values");
```
### Raw sequences section

TODO


### Minimizer sequences section

Writing a minimizer section consists of creating one or multiple blocks and then writing them.

Creating a block:
```rust
const ENCODING: u8 = 0b00011011; // A C T G, 2 bits per letter
[...]
let k: u64 = ...;
let data_size: usize = ...;  // size of the data associated with a kmer (in byte)
let data: Vec<u8> = ...;  // data for all kmers, size should be `nb_kmers * data_size` 
let sequence: String = ...;
// encode the sequence
let kmer_seq = kff::Kmer::from_ascii(sequence.as_bytes(), data, ENCODING);
// build the block
let block = kff::section::Block::new(
    k,
    data_size,
    kmer_seq,
    *minimizer_start_pos,
);
```
Writing a vector of blocks:
```rust
const ENCODING: u8 = 0b00011011; // A C T G, 2 bits per letter
[...]
let values = ...;
let minimizer: String = ...;

let section = kff::section::Minimizer::new(values).expect("unable to create a minimizer section from the values");

// encode minimizer
let minimizer: BitBox<u8, Msb0> = seq2bits(minimizer.as_bytes(), ENCODING);
// create an array of blocks
let blocks: Vec<kff::section::Block> = ...;
// write these blocks
file.write_minimizer(section, bitbox, blocks).expect("error writing the minimizer section");
```

