# Read

## Open a file

`kff::Kff` is the main object in the library.
This is the object needed to manipulate a binary kff file.

```rust
let file = kff::Kff::<std::io::BufReader<std::fs::File>>::open(kff_path).expect("could not open kff file");
```

## How to know the properties of my kmers?

The values stored in the file (e.g. k or data_size) are accessible through the *values* method.

```rust
let vars: &kff::section::Values = file.values();
let m: Option<&u64> = vars.get("m");
```

## Read header and encoding

When you open a file to read it, the header (including the encoding) is automatically read.
It is accessible via the file property *header*. You can access the encoding through the *encoding* header.

```rust
let header = file.header();
let encoding = file.header().encoding();
```

## Enumerating kmers from a file

This high-level reader API is made to be very easy to use.
This part of the API allows you to enumerate each pair of kmer/data through the whole file, hiding all the kff data structures.

```rust
let file = kff::Kff::<std::io::BufReader<std::fs::File>>::open("test")
    .expect("could not open kff file");
let encoding = *(file.header().encoding());
for kmer in file.kmers() {
    let kmer = kmer.expect("error reading the kmer");
    let data = kmer.data();

    // use kmer and data, e.g.:
    let kmer = String::from_utf8(kmer.seq(encoding)).expect("could not parse utf 8");
    println!("{}", kmer);
}
```

## Enumerating kmers sections

```rust
let mut file = kff::Kff::<std::io::BufReader<std::fs::File>>::open(args.input_kff).expect("could not open kff file");
let encoding = *(file.header().encoding());

while let Some(kmer_section) = file.next_kmer_section() {
    let kmer_section = kmer_section.expect("could not read the kmer section");
    for kmer in kmer_section {
        let data = kmer.data();
		let kmer = kmer.seq(encoding);
		// use kmer and data
    }
}
```

## Enumerating blocks

```
// TODO
```