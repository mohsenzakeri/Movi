#ifndef UTILS_HPP
#define UTILS_HPP

#include <sys/stat.h> 

#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>
#include <cstring>
#include <random>

#include <zlib.h>
#include "kseq.h"

#include "commons.hpp"
#include "move_query.hpp"
#include "movi_options.hpp"
#include "version.hpp"

// Forward declarations
class MoveStructure;

// Magic number to identify BPF files (Base Profile Format)
const uint32_t BPF_MAGIC = 0x42504600;  // "BPF\0" in ASCII

// Magic number to identify MOVI index files
const uint32_t MOVI_MAGIC = 0x4D4F5649;  // "Movi" in ASCII

// Header for MOVI index files
struct MoviHeader {
    uint32_t magic;      // Magic number to identify MOVI files
    uint8_t version = MOVI_VERSION_MAJOR; // Version number for future compatibility
    uint8_t version_minor = MOVI_VERSION_MINOR; // Minor version number for future compatibility
    uint8_t version_patch = MOVI_VERSION_PATCH; // Patch version number for future compatibility
    char type;           // Index mode (LARGE, CONSTANT, REGULAR, etc.)
    uint8_t reserved;    // Reserved for future use

    uint64_t length;      // Length of the original string
    uint64_t r;           // Number of movi rows after splitting
    uint64_t original_r;  // Number of runs before splitting
    uint64_t end_bwt_idx; // End of BWT index

    // Initialize header with mode

    void init(char type_, uint64_t _length, uint64_t _r, uint64_t _original_r, uint64_t _end_bwt_idx) {
        magic = MOVI_MAGIC;
        type = type_;
        reserved = 0;
        length = _length;
        r = _r;
        original_r = _original_r;
        end_bwt_idx = _end_bwt_idx;
    }

    // Write header to file
    void write(std::ofstream& file) {
        file.write(reinterpret_cast<const char*>(this), sizeof(MoviHeader));
    }
};

// Header for base profile format files
struct BPFHeader {
    uint32_t magic;      // Magic number to identify BPF files
    uint8_t version = BPF_VERSION_MAJOR; // Major version number for future compatibility
    uint8_t version_minor = BPF_VERSION_MINOR; // Minor version number for future compatibility
    uint8_t version_patch = BPF_VERSION_PATCH; // Patch version number for future compatibility
    uint8_t entry_size;  // Size of each entry in bits (16, 32, or 64)
    uint16_t reserved;   // Reserved for future use

    // Initialize header with entry size
    void init(uint8_t size) {
        magic = BPF_MAGIC;
        entry_size = size;
        reserved = 0;
    }

    // Write header to file
    void write(std::ofstream& file) {
        file.write(reinterpret_cast<const char*>(this), sizeof(BPFHeader));
    }

    // Write BPF header to a file
    void write_bpf_header(std::ofstream& file, uint8_t entry_size);
};

// Enum for different data types that can be output
enum class DataType {
    match_length,
    color,
    sa_entry
};

// Struct to hold all output file streams used across the codebase
struct OutputFiles {
    std::ofstream mls_file;
    std::ofstream matches_file;
    std::ofstream costs_file;
    std::ofstream scans_file;
    std::ofstream fastforwards_file;
    std::ofstream colors_file;
    std::ofstream kmer_file;
    std::ofstream sa_entries_file;
    std::ofstream out_file;

    // Default constructor
    OutputFiles() = default;

    // Move constructor and assignment (ofstream is not copyable)
    OutputFiles(OutputFiles&&) = default;
    OutputFiles& operator=(OutputFiles&&) = default;

    // Disable copy constructor and assignment
    OutputFiles(const OutputFiles&) = delete;
    OutputFiles& operator=(const OutputFiles&) = delete;
};

#define my_prefetch_r(address) __builtin_prefetch((void *)address, 0, 3)
#define my_prefetch_w(address) __builtin_prefetch((void *)address, 1, 3)

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)
kseq_t* open_kseq(gzFile& fp, std::string file_address);
void close_kseq(kseq_t *seq, gzFile& fp);
#endif

extern uint32_t alphamap_3[4][4];

#define LARGE_INDEX MODE == 0
#define CONSTANT_INDEX MODE == 1
#define SPLIT_INDEX MODE == 4
#define REGULAR_INDEX MODE == 3
#define REGULAR_THRESHOLDS_INDEX MODE == 6
#define BLOCKED_INDEX MODE == 2
#define BLOCKED_THRESHOLDS_INDEX MODE == 8
#define TALLY_INDEX MODE == 5
#define TALLY_THRESHOLDS_INDEX MODE == 7

#define NO_SAMPPLED_ID MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4 or MODE == 8 or MODE == 3 or MODE == 6
#define USE_THRESHOLDS MODE == 0 or MODE == 1 or MODE == 4 or MODE == 6 or MODE == 7 or MODE == 8
#define NO_THRESHOLDS MODE == 2 or MODE == 3 or MODE == 5
#define THRESHOLDS_WITHOUT_NEXTS MODE == 0 or MODE == 4 or MODE == 8 or MODE == 7 or MODE == 6
#define USE_NEXT_POINTERS MODE == 1
#define SPLIT_THRESHOLDS_FALSE MODE == 0 or MODE == 1 or MODE == 4
#define SPLIT_THRESHOLDS_TRUE MODE == 8 or MODE == 7 or MODE == 6
#define SPLIT_MAX_RUN MODE == 3 or MODE == 6 or MODE == 2 or MODE == 8 or MODE == 5 or MODE == 7
#define SPLIT_ARRAY MODE == 1 or MODE == 4
#define NO_EXTRA_TABLE MODE == 0 or MODE == 1 or MODE == 4 or MODE == 3 or MODE == 6
#define REGULAR_MODES MODE == 3 or MODE == 6
#define BLOCKED_MODES MODE == 2 or MODE == 8
#define TALLY_MODES MODE == 5 or MODE == 7
#define MOVI1_STYLE MODE == 0 or MODE == 1 or MODE == 4

// Feature compatibility macros
#define SUPPORTS_SEPARATORS (MODE == 2 or MODE == 3 or MODE == 5 or MODE == 6 or MODE == 7 or MODE == 8)

#define END_CHARACTER 0
#define THRBYTES 5
#define MIN_MATCHING_LENGTH 3
#define NULL_READ_CHUNK 150
#define NUM_NULL_READS 800 // 150,000 = 150 bp * 1000 reads
#define NULL_READ_BOUND 1000

#define UNCLASSIFIED_THRESHOLD 0.4

// To be used for generating random numbers for each thread
struct ThreadRandom {
    std::mt19937 generator;
    std::uniform_int_distribution<int> dist;

    ThreadRandom()
        : generator(std::random_device{}()),
          dist(1, 100)
    {}

    int get_random() {
        return dist(generator);
    }
};

std::string program();

std::string query_type(MoviOptions& movi_options);

char complement(char c);

std::string reverse_complement(std::string& fw);

std::string reverse_complement_from_pos(MoveQuery& mq_fw, int32_t pos_on_r, uint64_t match_len);

std::string number_to_kmer(size_t j, size_t m, std::vector<unsigned char>& alphabet, std::vector<uint64_t>& alphamap);

uint64_t kmer_to_number(size_t k, std::string& r, int32_t pos, std::vector<uint64_t>& alphamap, bool rc = false) ;

uint8_t F_char(std::vector<uint64_t>& first_runs, uint64_t run);

void read_thresholds(std::string tmp_filename, std::vector<uint64_t>& thresholds);

void output_base_stats(DataType data_type, bool to_stdout, std::ofstream& output_file, MoveQuery& mq);

void output_counts(bool to_stdout, std::ofstream& count_file, size_t query_length, int32_t pos_on_r, uint64_t match_count, MoveQuery& mq);

void output_kmers(bool to_stdout, std::ofstream& kmer_file, size_t all_kmer_count, MoveQuery& mq);

void output_logs(std::ofstream& costs_file, std::ofstream& scans_file, std::ofstream& fastforwards_file, MoveQuery& mq);

void output_read(MoveQuery& mq);

void open_output_files(MoviOptions& movi_options, OutputFiles& output_files);

void close_output_files(MoviOptions& movi_options, OutputFiles& output_files);

void print_query_stats(MoviOptions& movi_options, uint64_t total_ff_count, MoveStructure& mv);

// Borrowed from spumoni written by Omar Ahmed: https://github.com/oma219/spumoni/tree/main
std::string parse_null_reads(const char* ref_file, const char* output_path);

#endif