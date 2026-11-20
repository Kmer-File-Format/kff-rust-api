# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- New implementation with the features of the C++ API (kff-cpp-api 1.1), reading and writing the
  same bytes. The previous API (`Kff`, `KmerIterator`, `section::*`, `utils::*`) is replaced by
  `KmerReader` (high level reader), `Reader` (low level reader) and `Writer` (low level writer).
- Sequences are right aligned 2-bit nucleotides as in the format reference (they were left
  aligned).
- The encoding byte is read in the order A, C, G, T (G and T were swapped).
- The number of kmers of a block and the position of the minimizer use `ceil(log2(max))` and
  `ceil(log2(k + max - 1))` bits rounded to bytes, as the C++ API (they used 1, 2, 4 or 8 bytes).
- The variables of a values section are written sorted by name, as the C++ API.

### Added
- Index and footer written automatically, footer and index chain discovered when a file is opened.
- Minimizer sections: reading with or without the minimizer, writing with or without the
  minimizer.
- Jumps over blocks and sections, copy of sections from a reader to a writer.
- Release of the files of the writers (`tmp_close` in the C++ API), discard of a file being written.
- Compatibility tests with files written by the C++ API (tests/cpp).
