#ifndef __MOVE_STRUCTURE__
#define __MOVE_STRUCTURE__

#include <fstream>
#include <cstdint>
#include <zlib.h>
#include <stdio.h>
#include <chrono>
#include <cstddef>
#include <unistd.h>
#include <sys/stat.h>
#include <vector>
#include <unordered_map>

#include "kseq.h"
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>

#include "move_row.hpp"
#include "move_query.hpp"

#define END_CHARACTER 0
#define THRBYTES 5 

class MoveStructure {
    public:
        MoveStructure() { }
        MoveStructure(bool verbose_ = false, bool logs_ = false, uint16_t splitting = false);
        MoveStructure(char* input_file_, bool bit1_ = false, bool verbose_ = false, bool logs_ = false, uint16_t splitting = false);

        void build(std::ifstream &bwt_file);
        uint64_t query_ms(MoveQuery& mq, bool random);
        // void all_lf_test(/*std::ifstream &bwt_file*/);
        // uint64_t random_lf_test();
        // std::string reconstruct();
        // std::string reconstruct_move();

        uint64_t LF(uint64_t row_number);
        // uint64_t LF_move(uint64_t& pointer, uint64_t& i);
        // uint64_t fast_forward(uint64_t pointer, uint64_t index);
        uint64_t fast_forward(uint64_t& offset, uint64_t index, uint64_t x);
        char compute_char(uint64_t idx);
        uint64_t compute_threshold(uint64_t r_idx, uint64_t pointer, char lookup_char);
        uint32_t compute_index(char row_char, char lookup_char);

        // uint64_t naive_lcp(uint64_t row1, uint64_t row2);
        // uint64_t naive_sa(uint64_t bwt_row);

        uint64_t jump_up(uint64_t idx, char c);
        uint64_t jump_down(uint64_t idx, char c);
        bool jump_thresholds(uint64_t& idx, uint64_t offset, char r_char);
        bool jump_randomly(uint64_t& idx, char r_char);
        // bool jump_naive_lcp(uint64_t& idx, uint64_t pointer, char r_char, uint64_t& lcp);
        void compute_nexts();

        void serialize(char* output_dir);
        void deserialize(char* index_dir);
        
        std::unordered_map<uint32_t, uint32_t> jumps;
        std::unordered_map<uint32_t, uint32_t> ff_counts;
        std::unordered_map<uint64_t, uint64_t> run_lengths;

        uint64_t get_n(uint64_t idx);
        uint64_t get_offset(uint64_t idx);
        uint64_t get_thresholds(uint64_t idx, uint32_t alphabet_index);

    private:
        bool bit1;
        uint16_t splitting;
        std::string bwt_string;
        std::string orig_string;
        bool reconstructed;
        uint64_t length;
        uint64_t r;
        uint64_t original_r;
        uint64_t end_bwt_idx;
        uint64_t end_bwt_idx_thresholds[4];
        bool verbose;
        bool logs;
	char* input_file;

        // Map from 2bit encoded character to the actual character
        // Example: alphabet[0] -> A, alphabet[1] -> C
        std::vector<unsigned char> alphabet;
        // Number of each character
        std::vector<uint64_t> counts;
        // Map from the character to the index of the character
        // Example: alphamap[A] -> 0, alphamap[C] -> 1
        std::vector<uint64_t> alphamap;

        std::vector<std::unique_ptr<sdsl::bit_vector> > occs;
        std::vector<std::unique_ptr<sdsl::rank_support_v<> > > occs_rank;
        sdsl::bit_vector bits;
        sdsl::rank_support_v<> rbits;
        sdsl::select_support_mcl<> sbits;
        sdsl::int_vector<> thresholds;

        std::vector<MoveRow> rlbwt;
        // std::vector<char> rlbwt_chars;
        uint64_t eof_row;
        uint64_t bit1_begin;
        uint64_t bit1_after_eof;

        // auxilary datastructures for the length, offset and thresholds overflow
        std::vector<uint64_t> n_overflow;
        std::vector<uint64_t> offset_overflow;
        std::vector<std::vector<uint64_t> > thresholds_overflow;
};

#endif
