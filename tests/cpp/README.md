# Reference files of the C++ API

`generate_fixtures.cpp` writes the files of `tests/data/cpp` with the
[C++ API](https://github.com/Kmer-File-Format/kff-cpp-api): for each file, `NAME.kff` and
`NAME.txt`, the kmers read by the C++ high level reader (one line per kmer: nucleotides, a
tabulation and the data in hexadecimal).

`tests/cpp_compat.rs` checks that this crate reads the same kmers, and that it writes the same bytes
when the sections are written again.

The files were generated with kff-cpp-api 1.1. To generate them again:

```bash
KFF_CPP_API=path/to/kff-cpp-api
g++ -std=c++17 -O2 -I $KFF_CPP_API tests/cpp/generate_fixtures.cpp $KFF_CPP_API/kff_io.cpp -o generate
./generate tests/data/cpp
```

The generation is deterministic: the `.kff` files must not change.
