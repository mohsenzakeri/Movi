#ifndef UTILS_HPP
#define UTILS_HPP

#include <sys/stat.h> 

#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>
#include <cstring>

#include <sdsl/int_vector.hpp>
#include <zlib.h>
#include "kseq.h"

#include "move_query.hpp"
#include "movi_options.hpp"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)
kseq_t* open_kseq(gzFile& fp, std::string file_address);
void close_kseq(kseq_t *seq, gzFile& fp);
#endif

extern uint32_t alphamap_3[4][4];

#define USE_THRESHOLDS MODE == 0 or MODE == 1 or MODE == 4 or MODE == 6 or MODE == 7 or MODE == 8
#define NO_THRESHOLDS MODE == 2 or MODE == 3 or MODE == 5
#define THRESHOLDS_WITHOUT_NEXTS MODE == 0 or MODE == 4 or MODE == 8 or MODE == 7 or MODE == 6
#define USE_NEXT_POINTERS MODE == 1
#define SPLIT_THRESHOLDS_FALSE MODE == 0 or MODE == 1 or MODE == 4
#define SPLIT_THRESHOLDS_TRUE MODE == 8 or MODE == 7 or MODE == 6
#define SPLIT_MAX_RUN MODE == 3 or MODE == 6 or MODE == 2 or MODE == 8 or MODE == 5 or MODE == 7
#define SPLIT_ARRAY MODE == 1 or MODE == 4
#define NO_SPLIT MODE == 0
#define NO_EXTRA_TABLE MODE == 0 or MODE == 1 or MODE == 4 or MODE == 3 or MODE == 6
#define CONSTANT_MODE MODE == 1
#define COMPACT_MODE MODE == 3 or MODE == 6 // to be changed to regular
#define BLOCKED_MODE MODE == 2 or MODE == 8
#define TALLY_MODE MODE == 5 or MODE == 7

#define END_CHARACTER 0
#define THRBYTES 5
#define MIN_MATCHING_LENGTH 3
#define NULL_READ_CHUNK 150
#define NUM_NULL_READS 800 // 150,000 = 150 bp * 1000 reads
#define NULL_READ_BOUND 1000

std::string program();

std::string query_type(MoviOptions& movi_options);

char complement(char c);

std::string reverse_complement(std::string& fw);

std::string reverse_complement_from_pos(MoveQuery& mq_fw, int32_t pos_on_r, uint64_t match_len);

std::string number_to_kmer(size_t j, size_t m, std::vector<unsigned char>& alphabet, std::vector<uint64_t>& alphamap);

uint64_t kmer_to_number(size_t k, std::string& r, int32_t pos, std::vector<uint64_t>& alphamap, bool rc = false) ;

uint8_t F_char(std::vector<uint64_t>& first_runs, uint64_t run);

void read_thresholds(std::string tmp_filename, sdsl::int_vector<>& thresholds);

void output_matching_lengths(bool to_stdout, std::ofstream& mls_file, std::string read_id, MoveQuery& mq);

void output_counts(bool to_stdout, std::ofstream& count_file, std::string read_id, size_t query_length, int32_t pos_on_r, uint64_t match_count);

void output_logs(std::ofstream& costs_file, std::ofstream& scans_file, std::ofstream& fastforwards_file, std::string read_id, MoveQuery& mq);

// Borrowed from spumoni written by Omar Ahmed: https://github.com/oma219/spumoni/tree/main
std::string parse_null_reads(const char* ref_file, const char* output_path);

#endif