# Read

The crate offers two readers, as the C++ API:
- `kff::KmerReader`, the high level reader, enumerates the kmers or the blocks of a file without
  dealing with its sections (`Kff_reader` in the C++ API);
- `kff::Reader`, the low level reader, reads the header, the footer, the index and each section
  (`Kff_file` opened in "r" mode and the `Section_*` classes).

The sequences are right aligned 2-bit nucleotides: a sequence of `n` nucleotides uses `ceil(n / 4)`
bytes, its last nucleotide is in the 2 lowest bits of the last byte. `kff::Encoding` converts them
from and to text with the encoding of the file.

## Enumerating kmers

```rust
let mut reader = kff::KmerReader::open("file.kff")?;
while let Some(kmer) = reader.next_kmer()? {
    // kmer.sequence: 2-bit nucleotides, kmer.data: data of the kmer, kmer.k: size of the kmer
    println!("{} {:?}", String::from_utf8_lossy(&kmer.to_ascii()), kmer.data);
}
```

The kmer and its data are borrowed from the reader until the next call. To get owned kmers, use
the iterator:

```rust
for kmer in kff::KmerReader::open("file.kff")?.kmers() {
    let kmer = kmer?;
    println!("{} {:?}", kmer, kmer.data);
}
```

## Enumerating blocks

```rust
let mut reader = kff::KmerReader::open("file.kff")?;
while let Some(block) = reader.next_block()? {
    // block.sequence: nb_kmers + k - 1 nucleotides, block.data: data_size bytes per kmer
    for (kmer, data) in block.kmers() {
        // ...
    }
}
```

## Properties of the kmers

```rust
let mut reader = kff::KmerReader::open("file.kff")?;
reader.has_next()?;
// Variables of the current values section (get_var in the C++ API)
let m: Option<u64> = reader.var("m");
// Variables used by all the kmers
let (k, max, data_size) = (reader.k(), reader.max(), reader.data_size());
// Encoding: 2-bit values of A, C, G and T
let encoding: [u8; 4] = reader.encoding().values();
```

## Low level reader

### Open a file

```rust
let reader = kff::Reader::open("file.kff")?;
// Any std::io::Read + std::io::Seek can be read
let reader = kff::Reader::new(std::io::Cursor::new(bytes))?;
```

At the opening, the header is read and checked (signatures, version of the format, encoding), the
footer and the chain of index sections are loaded, and the reader is placed on the first section.

```rust
let header = reader.header();
let encoding = reader.encoding();
let unique = reader.uniqueness();
let canonical = reader.canonicity();
let metadata: &[u8] = reader.metadata();

if let Some(footer) = reader.footer() {
    for (name, value) in footer {
        println!("{}: {}", name, value);
    }
}
for index in reader.index() {
    println!("index section at {}", index.beginning);
    for (relative_position, section_type) in &index.entries {
        println!("{:?} at {} from the end of the index section", section_type, relative_position);
    }
}
// Absolute positions of all the indexed sections
let sections = reader.indexed_sections();
```

### Sections

`read_section_type` gives the type of the next section without moving (`None` at the end of the
file). Each type of section is then read with its function:

```rust
use kff::{BlockReader, SectionType};

while let Some(section_type) = reader.read_section_type()? {
    match section_type {
        SectionType::Values => {
            // The variables become the current scope
            let values = reader.read_values()?;
        }
        SectionType::Index => {
            let index = reader.read_index()?;
        }
        SectionType::Raw => {
            let mut section = reader.open_raw()?;
            let mut seq = vec![0; section.max_sequence_bytes()];
            let mut data = vec![0; section.max_data_bytes()];
            while section.remaining_blocks() > 0 {
                let nb_kmers = section.read_compacted_sequence(&mut seq, &mut data)?;
                // ...
            }
        }
        SectionType::Minimizer => {
            let mut section = reader.open_minimizer()?;
            let minimizer = section.minimizer().to_vec();
            let mut seq = vec![0; section.max_sequence_bytes()];
            let mut data = vec![0; section.max_data_bytes()];
            while section.remaining_blocks() > 0 {
                // Sequence without the minimizer, as stored in the file
                let (nb_kmers, mini_pos) =
                    section.read_compacted_sequence_without_mini(&mut seq, &mut data)?;
                // Insert the minimizer back if needed
                section.add_minimizer(nb_kmers, &mut seq, mini_pos)?;
            }
        }
        SectionType::Other(byte) => panic!("unknown section {}", byte),
    }
}
```

The raw and minimizer sections share the `kff::BlockReader` trait (`Block_section_reader` in the
C++ API): `read_compacted_sequence` (full sequence, the minimizer is inserted back),
`read_compacted_sequence_joined` (sequence and data in the same buffer), `read_block` (owned
block), `jump_sequence` and `jump_section`. `reader.open_block_section()` opens the next section if
it is a raw or minimizer section, and `reader.read_section()` reads or opens any section
(`SectionBuilder::build`).

When a section reader is closed or dropped, its remaining blocks are jumped. A whole section can
also be jumped with `reader.jump_next_section()` (raw and minimizer sections only, as in the C++
API) or `reader.skip_section()` (any section, the values sections are read).

### Copy

The sections can be copied to a `kff::Writer` (`copy` in the C++ API):

```rust
let values = reader.read_values()?;
writer.write_values(&values)?;
reader.open_raw()?.copy_to(&mut writer)?;
```

### Positions

`tell`, `jump`, `jump_to` and `read_bytes` give a direct access to the file.

## Correspondence with the C++ API

| C++ API | Rust |
| --- | --- |
| `Kff_reader(name)` | `KmerReader::open(name)` |
| `has_next`, `next_kmer`, `next_block` | `has_next`, `next_kmer`, `next_block` |
| `get_var`, `get_encoding`, `k`, `max`, `data_size` | `var`, `encoding().values()`, `k()`, `max()`, `data_size()` |
| `Kff_file(name, "r")` | `Reader::open(name)` |
| `encoding`, `uniqueness`, `canonicity`, `read_metadata` | `encoding()`, `uniqueness()`, `canonicity()`, `metadata()` |
| `footer`, `index`, `global_vars` | `footer()`, `index()`, `global_vars()` |
| `read_section_type`, `jump_next_section` | `read_section_type`, `jump_next_section` |
| `Section_GV(&file)` | `read_values()` |
| `Section_Index(&file)` | `read_index()` |
| `Section_Raw(&file)`, `Section_Minimizer(&file)` | `open_raw()`, `open_minimizer()` |
| `Block_section_reader::construct_section` | `open_block_section()` |
| `SectionBuilder::build` | `read_section()` |
| `read_compacted_sequence(seq, data)` | `read_compacted_sequence(&mut seq, &mut data)` |
| `read_compacted_sequence(seq_data)` | `read_compacted_sequence_joined(&mut seq_data)` |
| `read_compacted_sequence_without_mini` | `read_compacted_sequence_without_mini` (returns the minimizer position) |
| `copy(&file)` | `copy_to(&mut writer)` |
| `tellp`, `jump`, `jump_to`, `read` | `tell`, `jump`, `jump_to`, `read_bytes` |
