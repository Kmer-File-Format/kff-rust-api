# Read

## Raw kmer

```rust
let file: kff::Kff::<std::io::BufReader<std::fs::File>> = kff::Kff::open("raw_kmer.kff")?;
let encoding = *(file.header().encoding()); // extract encoding

let mut iter = file.kmers();
while let Some(Ok(kmer)) = iter.next() {
	let kmer_seq: Vec<u8> = kmer.seq(encoding);
	let data = kmer.get_data();

	... Do what you want ...
}
```

## Minimizer kmer

```rust
// TODO
```
