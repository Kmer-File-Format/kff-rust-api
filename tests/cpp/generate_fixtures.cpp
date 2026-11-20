/**
 * Writes the reference files of the compatibility tests with the C++ API (kff-cpp-api).
 *
 * For each file NAME, writes NAME.kff with the low level API of kff-cpp-api and NAME.txt, the
 * kmers of the file read with its high level reader (Kff_reader, block by block): one line per
 * kmer with the nucleotides and the data in hexadecimal, separated by a tabulation.
 *
 * Build and run (see tests/cpp/README.md):
 *   g++ -std=c++17 -O2 -I KFF_CPP_API tests/cpp/generate_fixtures.cpp KFF_CPP_API/kff_io.cpp -o generate
 *   ./generate tests/data/cpp
 */

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "kff_io.hpp"

using Bytes = std::vector<uint8_t>;

/** Deterministic pseudo-random generator (splitmix64) */
class Random {
	uint64_t state;
public:
	explicit Random(uint64_t seed) : state(seed) {}
	uint64_t next() {
		uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
		z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
		z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
		return z ^ (z >> 31);
	}
	uint64_t between(uint64_t min, uint64_t max) { return min + next() % (max - min + 1); }
	std::string sequence(size_t size) {
		std::string s(size, 'A');
		for (char & c : s)
			c = "ACGT"[between(0, 3)];
		return s;
	}
	Bytes bytes(size_t size) {
		Bytes b(size);
		for (uint8_t & byte : b)
			byte = (uint8_t)next();
		return b;
	}
};

static uint8_t nucleotide_value(char nucleotide, const uint8_t encoding[4]) {
	switch (nucleotide) {
		case 'A': return encoding[0];
		case 'C': return encoding[1];
		case 'G': return encoding[2];
		default: return encoding[3];
	}
}

/** 2 bits per nucleotide, right aligned, null padding bits */
static Bytes encode(const std::string & sequence, const uint8_t encoding[4]) {
	const size_t size = sequence.size();
	const size_t padding = (4 - size % 4) % 4;
	Bytes bytes((size + 3) / 4, 0);
	for (size_t i = 0; i < size; i++) {
		size_t position = padding + i;
		bytes[position / 4] |= nucleotide_value(sequence[i], encoding) << (2 * (3 - position % 4));
	}
	return bytes;
}

static std::string decode(const uint8_t * bytes, size_t size, const uint8_t encoding[4]) {
	char letters[4];
	letters[encoding[0]] = 'A';
	letters[encoding[1]] = 'C';
	letters[encoding[2]] = 'G';
	letters[encoding[3]] = 'T';
	const size_t padding = (4 - size % 4) % 4;
	std::string sequence(size, '?');
	for (size_t i = 0; i < size; i++) {
		size_t position = padding + i;
		sequence[i] = letters[(bytes[position / 4] >> (2 * (3 - position % 4))) & 0b11];
	}
	return sequence;
}

static const uint8_t DEFAULT_ENCODING[4] = {0, 1, 3, 2};

/** Raw section of nb_blocks random blocks */
static void raw_section(Kff_file & file, Random & rnd, uint64_t k, uint64_t max, uint64_t data_size,
                        size_t nb_blocks, const uint8_t encoding[4]) {
	Section_GV sgv(&file);
	sgv.write_var("k", k);
	sgv.write_var("max", max);
	sgv.write_var("data_size", data_size);
	sgv.close();

	Section_Raw sr(&file);
	for (size_t i = 0; i < nb_blocks; i++) {
		// The first block is full (except for the large max, to keep the files small)
		uint64_t nb_kmers = (i == 0 and max <= 1000) ? max : rnd.between(1, max < 40 ? max : 40);
		std::string sequence = rnd.sequence(k + nb_kmers - 1);
		Bytes seq = encode(sequence, encoding);
		Bytes data = rnd.bytes(nb_kmers * data_size + 1);
		sr.write_compacted_sequence(seq.data(), sequence.size(), data.data());
	}
	sr.close();
}

/** Minimizer section of nb_blocks random blocks. Odd blocks are written with the full sequence
 *  (the API removes the minimizer), even blocks without the minimizer. */
static void minimizer_section(Kff_file & file, Random & rnd, uint64_t k, uint64_t m, uint64_t max,
                              uint64_t data_size, size_t nb_blocks, const uint8_t encoding[4]) {
	Section_GV sgv(&file);
	sgv.write_var("k", k);
	sgv.write_var("m", m);
	sgv.write_var("max", max);
	sgv.write_var("data_size", data_size);
	sgv.write_var("ordered", 0);
	sgv.close();

	const std::string minimizer = rnd.sequence(m);
	Section_Minimizer sm(&file);
	Bytes mini = encode(minimizer, encoding);
	sm.write_minimizer(mini.data());
	for (size_t i = 0; i < nb_blocks; i++) {
		uint64_t nb_kmers = (i == 0 and max <= 1000) ? max : rnd.between(1, max < 40 ? max : 40);
		std::string sequence = rnd.sequence(k + nb_kmers - 1);
		uint64_t mini_pos = rnd.between(0, sequence.size() - m);
		sequence.replace(mini_pos, m, minimizer);
		Bytes data = rnd.bytes(nb_kmers * data_size + 1);
		if (i % 2 == 1) {
			Bytes seq = encode(sequence, encoding);
			sm.write_compacted_sequence(seq.data(), sequence.size(), mini_pos, data.data());
		} else {
			std::string without = sequence.substr(0, mini_pos) + sequence.substr(mini_pos + m);
			Bytes seq = encode(without, encoding);
			sm.write_compacted_sequence_without_mini(seq.data(), without.size(), mini_pos, data.data());
		}
	}
	sm.close();
}

/** Kmers of the file read with Kff_reader.
 *
 * The blocks are read with next_block: next_kmer of kff-cpp-api 1.1 loses the kmers of the last
 * block (except the first one) when the file ends without footer, because has_next() compares the
 * position in the file with its end while the last block is already read. */
static void dump(const std::string & kff, const std::string & txt) {
	Kff_reader reader(kff);
	std::ofstream out(txt);
	std::vector<uint8_t> seq_buffer(1 << 20), data_buffer(1 << 22);
	while (reader.has_next()) {
		uint8_t * seq = seq_buffer.data();
		uint8_t * data = data_buffer.data();
		const uint64_t nb_kmers = reader.next_block(seq, data);
		const uint64_t k = reader.k;
		const std::string sequence = decode(seq, nb_kmers + k - 1, reader.get_encoding());
		for (uint64_t i = 0; i < nb_kmers; i++) {
			out << sequence.substr(i, k) << '\t';
			for (uint64_t j = 0; j < reader.data_size; j++) {
				char hex[3];
				snprintf(hex, 3, "%02x", data[i * reader.data_size + j]);
				out << hex;
			}
			out << '\n';
		}
	}
}

int main(int argc, char ** argv) {
	if (argc != 2) {
		std::cerr << "usage: " << argv[0] << " OUTPUT_DIRECTORY" << std::endl;
		return 1;
	}
	const std::string dir = std::string(argv[1]) + "/";
	std::vector<std::string> names;

	// Empty file: header, empty index and footer
	{
		Kff_file file(dir + "empty.kff", "w");
		file.close();
		names.push_back("empty");
	}

	// One raw section, indexed, with metadata
	{
		Random rnd(1);
		Kff_file file(dir + "raw.kff", "w");
		const std::string meta = "kff-rust-api compatibility test";
		file.write_metadata(meta.size(), (const uint8_t *)meta.data());
		raw_section(file, rnd, 15, 50, 1, 30, DEFAULT_ENCODING);
		file.close();
		names.push_back("raw");
	}

	// Raw sections of many sizes, not indexed, other encoding and flags
	{
		Random rnd(2);
		const uint8_t encoding[4] = {2, 0, 1, 3};
		Kff_file file(dir + "raw_sizes.kff", "w");
		file.set_indexation(false);
		file.write_encoding(2, 0, 1, 3);
		file.set_uniqueness(true);
		file.set_canonicity(true);
		const uint64_t ks[] = {1, 2, 4, 5, 15, 16, 31, 32, 33, 63, 64, 65};
		const uint64_t maxs[] = {1, 2, 3, 255, 300, 70000};
		const uint64_t data_sizes[] = {0, 1, 3};
		size_t idx = 0;
		for (uint64_t k : ks)
			for (uint64_t max : maxs) {
				uint64_t data_size = data_sizes[idx++ % 3];
				raw_section(file, rnd, k, max, data_size, max > 1000 ? 2 : 6, encoding);
			}
		file.close();
		names.push_back("raw_sizes");
	}

	// Minimizer sections of many sizes (sequences shorter and longer than 32 nucleotides)
	{
		Random rnd(3);
		Kff_file file(dir + "minimizer.kff", "w");
		const uint64_t params[][4] = {
			// k, m, max, data_size
			{1, 1, 1, 0}, {5, 1, 3, 1}, {5, 5, 4, 2}, {10, 8, 255, 1}, {15, 7, 20, 0},
			{21, 9, 300, 2}, {31, 11, 40, 1}, {32, 13, 7, 1}, {33, 10, 1, 3}, {40, 31, 70, 0},
			{63, 15, 2, 1}, {64, 21, 30, 2}, {65, 32, 100, 1}, {100, 33, 10, 1}};
		for (const auto & p : params)
			minimizer_section(file, rnd, p[0], p[1], p[2], p[3], 7, DEFAULT_ENCODING);
		file.close();
		names.push_back("minimizer");
	}

	// Mixed sections, empty sections, extra variables, indexed
	{
		Random rnd(4);
		const uint8_t encoding[4] = {3, 2, 1, 0};
		Kff_file file(dir + "mixed.kff", "w");
		file.write_encoding(3, 2, 1, 0);
		file.set_canonicity(true);
		raw_section(file, rnd, 7, 10, 2, 5, encoding);
		{
			Section_Raw empty(&file);
			empty.close();
		}
		minimizer_section(file, rnd, 21, 6, 30, 1, 8, encoding);
		{
			Section_GV empty(&file);
			empty.close();
		}
		{
			Section_GV sgv(&file);
			sgv.write_var("k", 9);
			sgv.write_var("m", 4);
			sgv.write_var("max", 5);
			sgv.write_var("data_size", 0);
			sgv.write_var("nb_kmers", 1234);
			sgv.write_var("a variable name", 0xFFFFFFFFFFFFFFFFULL);
			sgv.close();
			Section_Minimizer empty(&file);
			Bytes mini = encode("ACGT", encoding);
			empty.write_minimizer(mini.data());
			empty.close();
			Section_Raw sr(&file);
			Bytes seq = encode("ACGTACGTACGTA", encoding);
			sr.write_compacted_sequence(seq.data(), 13, nullptr);
			sr.close();
		}
		raw_section(file, rnd, 31, 1, 4, 3, encoding);
		file.close();
		names.push_back("mixed");
	}

	// Chain of index sections written by hand at the beginning of the file, without footer
	{
		Random rnd(5);
		Kff_file file(dir + "chained_index.kff", "w");
		file.set_indexation(false);
		file.write_metadata(3, (const uint8_t *)"abc");
		// First index at 15 (26 bytes), an empty values section at 41, the second index at 50
		Section_Index first(&file);
		first.register_section('v', 0);
		first.set_next_index(9);
		first.close();
		Section_GV sgv(&file);
		sgv.close();
		// Second index at 50 (35 bytes): the values section at 41 - 85, the next values section at 85
		Section_Index second(&file);
		second.register_section('v', 41 - 85);
		second.register_section('v', 0);
		second.close();
		raw_section(file, rnd, 11, 20, 1, 4, DEFAULT_ENCODING);
		file.close();
		names.push_back("chained_index");
	}

	for (const std::string & name : names)
		dump(dir + name + ".kff", dir + name + ".txt");

	return 0;
}
