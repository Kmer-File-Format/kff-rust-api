# Write

## Raw kmer

```rust
// Build header of your kff file
let header = kff:section::Header::new(
1,                       // major header version
0,                       // minor header version
0b00011110,              // nucleotide 2bit encoding A -> 00, C -> 01, T -> 11, G -> 10
true,                    // kmer store are uniq ?
true,                    // kmer store are canonical ?
b"producer: kff_example" // comment header block
)?;

// Create file
let mut kff = kff::Kff::create("raw_kmer.kff", header)?;

// First Values section required
let mut values = kff::section::Values::default();
values.insert("k".to_string(), kmer_size);
values.insert("ordered".to_string(), true);
values.insert("max".to_string(), 200); // max number of kmer per block
values.insert("data_size".to_string(), 1); // number of bytes to store data associated to kmer

kff.write_values(values.clone())?;

let blocks = Vec::new();
loop {
	let seq = bitvec::vec::BitVec::<u8, bitvec::order::Msb0>::new();

	... set seq ...

	let values = Vec::<u8>::new();

	... set values ...

	blocks.push(kff::section::Block::new(
		kmer_size,
		1, // number of kmer in block
		kff::Kmer::new(
			seq,
			values,
		),
		0 // minimizer offset for raw kmer set to 0
	));
}

kff.write_raw(kff::section::Raw::new(&values), &blocks)?;

kff.finalize()?;
```

## Minimizer kmer

```rust
// TODO
```
