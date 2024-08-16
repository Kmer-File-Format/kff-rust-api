# Write

## String sequences to binary

This crate provides two functions to transform a DNA string into a byte array coding for the binary sequence:
- `from_ascii` for a sequence with associated data (e.g. kmers)
- `seq2bits` for a sequence without any associated data (e.g. minimizers).
You can, of course, use another one of your own.

## Open a file
Creating a kff file requires creating a header.
The header first contains the kff file version (minor and major).
You then have to fill in if the kmers are unique and/or canonical.
Then you have to write the encoding you are using.
Optionally, you can add some data in a vector of bytes.

```rust
const ENCODING: u8 = 0b00011011; // A C T G, 2 bits per letter

// Build header of your kff file
let header = kff:section::Header::new(
	1,                       // major header version
	0,                       // minor header version
	ENCODING,              // nucleotide 2bit encoding A -> 00, C -> 01, T -> 11, G -> 10
	true,                    // are stored kmers unique ?
	true,                    // are stored kmers canonical ?
	b"producer: kff_example" // free header block
)?;

// Create file
let mut kff = kff::Kff::create("raw_kmer.kff", header)?;
```

## Close a file
Before dropping the file, be sure to call the `finalize` method.
```rust
kff.finalize()?;
```

## Write values

As defined in the standard, writing values is necessary for writing some other sections. Please refer to the standard for more information on which value to write.

```rust
// First Values section required
let mut values = kff::section::Values::default();
values.insert("k".to_string(), kmer_size);
values.insert("ordered".to_string(), true);
values.insert("max".to_string(), 200); // max number of kmer per block
values.insert("data_size".to_string(), 1); // number of bytes to store data associated to kmer

kff.write_values(values.clone())?;
```

## Raw kmer

```rust
let blocks = Vec::new();
loop {
	let seq = bitvec::vec::BitVec::<u8, bitvec::order::Msb0>::new();

	... set seq ...

	let data = Vec::<u8>::new();

	... set data ...

	blocks.push(kff::section::Block::new(
		kmer_size,
		1, // number of kmer in block
		kff::Kmer::new(
			seq,
			data,
		),
		0 // starting position of the minimizer in the kmer
	));
}

kff.write_raw(kff::section::Raw::new(&values), &blocks)?;

kff.finalize()?;  // be sure to call the finalize method
```

## Minimizer sequences section

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
let values = ...;   // see section above to create values
let minimizer: String = ...;  // let's suppose you have a minimizer as a string

// create a minimizer section
let section = kff::section::Minimizer::new(values)?;

// encode minimizer
let minimizer: BitBox<u8, Msb0> = seq2bits(minimizer.as_bytes(), ENCODING);
// create an array of blocks
let blocks: Vec<kff::section::Block> = ...;
// write these blocks in the section
kff.write_minimizer(section, bitbox, blocks)?;
[...]
kff.finalize()?;  // be sure to call the finalize method
```