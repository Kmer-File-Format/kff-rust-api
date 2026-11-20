# Write

`kff::Writer` writes KFF files (`Kff_file` opened in "w" mode and the `Section_*` classes of the C++
API). With the same calls, it writes the same bytes as the C++ API.

## Sequences

The sequences are given as right aligned 2-bit nucleotides: a sequence of `n` nucleotides uses
`ceil(n / 4)` bytes, its last nucleotide is in the 2 lowest bits of the last byte. `kff::Encoding`
encodes text sequences with the encoding of the file:

```rust
let encoding = writer.encoding();
let bytes: Vec<u8> = encoding.encode(b"GATTACA")?;
```

## Open a file

```rust
let mut writer = kff::Writer::create("file.kff")?;
// Any std::io::Write + std::io::Seek can be written, from its current position
let mut writer = kff::Writer::new(std::io::Cursor::new(Vec::new()))?;
```

## Header

The signature and the version of the format are written at the creation. The default encoding is
A=0, C=1, G=3, T=2 and the flags are false. They can be changed at any time:

```rust
//                                     A  C  G  T
writer.set_encoding(kff::Encoding::new(0, 1, 3, 2)?)?;
writer.set_uniqueness(true)?;
writer.set_canonicity(false)?;
```

The metadata must be written before the first section. Without call, an empty metadata is written.

```rust
writer.write_metadata(b"<3 KFF <3")?;
```

## Values section

A values section starts a new scope: the variables of the previous one are forgotten. The variables
are written when the section is closed, sorted by name.

```rust
let mut values = writer.values_section()?;
values.write_var("k", 7)?;
values.write_var("m", 3)?;
values.write_var("max", 8)?;
values.write_var("data_size", 1)?;
values.close()?;

// Or from a map
let mut vars = kff::Values::new();
vars.insert("k".to_string(), 7);
writer.write_values(&vars)?;
```

## Raw section

The variables k, max and data_size must be defined in the current values section. Each block is a
sequence and the data of its kmers. The number of blocks is written when the section is closed.

```rust
let encoding = writer.encoding();
let mut raw = writer.raw_section()?;
// One 7-mer
raw.write_compacted_sequence(&encoding.encode(b"GATTACA")?, 7, &[3])?;
// 8 overlapping 7-mers: TGGTACA, GGTACAA, GTACAAG, ...
raw.write_compacted_sequence(
    &encoding.encode(b"TGGTACAAGTTACC")?,
    14,
    &[1, 221, 7, 1, 5, 6, 12, 198],
)?;
raw.close()?;
```

## Minimizer section

The variables k, m, max and data_size must be defined in the current values section. The minimizer
is written once for the section, and each block gives the position of the minimizer in its
sequence.

```rust
let mut section = writer.minimizer_section()?;
section.write_minimizer(&encoding.encode(b"ACA")?)?;
// Full sequence: the minimizer (position 4) is removed by the writer
section.write_compacted_sequence(&encoding.encode(b"GATTACA")?, 7, 4, &[3])?;
// Sequence GGTACAAGTTACCT already without its minimizer (position 3)
section.write_compacted_sequence_without_mini(
    &encoding.encode(b"GGTAGTTACCT")?,
    11,
    3,
    &[1, 221, 7, 1, 5, 6, 12, 198],
)?;
section.close()?;
```

## Index and footer

By default, as with the C++ API, the writer registers the position of each section and writes, on
closing, an index section of all the sections and a footer (a values section with "first_index"
and "footer_size"). This can be disabled:

```rust
writer.set_indexation(false);
```

Index sections can also be written by hand, with positions relative to the end of the index
section:

```rust
let mut index = writer.index_section()?;
index.register_section(kff::SectionType::Raw, -1463);
index.set_next_index(0);
index.close()?;
```

## Close the file

```rust
writer.close()?;
```

`close` writes the index and the footer (if enabled) and the final signature. A writer that is
dropped without being closed is closed, ignoring the errors; a section writer that is dropped is
closed too. `writer.discard()` stops the writing and deletes the file (`close(false)` in the C++
API).

When a lot of files are written at the same time, `writer.release()` writes the buffer and closes
the file, which is reopened at the next write (`tmp_close` in the C++ API).

## Correspondence with the C++ API

| C++ API | Rust |
| --- | --- |
| `Kff_file(name, "w")` | `Writer::create(name)` |
| `write_encoding(a, c, g, t)` | `set_encoding(Encoding::new(a, c, g, t)?)` |
| `set_uniqueness`, `set_canonicity`, `write_metadata` | `set_uniqueness`, `set_canonicity`, `write_metadata` |
| `set_indexation`, `register_position` | `set_indexation`, `register_position` |
| `Section_GV(&file)`, `write_var`, `close` | `values_section()`, `write_var`, `close` |
| `Section_Index(&file)`, `register_section`, `set_next_index` | `index_section()`, `register_section`, `set_next_index` |
| `Section_Raw(&file)`, `write_compacted_sequence` | `raw_section()`, `write_compacted_sequence` |
| `Section_Minimizer(&file)`, `write_minimizer` | `minimizer_section()`, `write_minimizer` |
| `write_compacted_sequence(seq, size, mini_pos, data)` | same |
| `write_compacted_sequence_without_mini` | same |
| `close()`, `close(false)` | `close()`, `discard()` |
| `tmp_close`, `reopen` | `release()` (reopened automatically) |
| `write`, `write_at`, `tellp` | `write_bytes`, `write_at`, `tell` |
