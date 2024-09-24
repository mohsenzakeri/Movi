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
        MoveStructure(bool verbose_, bool logs_);
        MoveStructure(bool onebit_, bool verbose_, bool logs_, uint16_t splitting = 0, bool constant = false);
        MoveStructure(std::string input_file_, bool onebit_, bool verbose_, bool logs_, uint16_t splitting = 0, bool constant = false);

        bool check_mode();
        std::string index_type();
        void build(std::ifstream &bwt_file);
        void build_rlbwt(std::string input_file);
        uint64_t query_pml(MoveQuery& mq, bool random);
        uint64_t backward_search(std::string& R,  int32_t& pos_on_r);
        uint64_t exact_matches(MoveQuery& mq);
    
        void sequential_lf();
        void random_lf();
        std::string reconstruct_lf();

        // std::string reconstruct();
        // char compute_char(uint64_t idx);

        uint64_t LF(uint64_t row_number);
        uint16_t LF_move(uint64_t& pointer, uint64_t& i);
        uint64_t fast_forward(uint64_t& offset, uint64_t index, uint64_t x);

        uint64_t compute_threshold(uint64_t r_idx, uint64_t pointer, char lookup_char);
        uint32_t compute_index(char row_char, char lookup_char);
        void compute_nexts();

        // Finds SA entries of all rows in BWT.
        void find_all_SA();
        // Find which document a given SA entry belongs to.
        uint16_t find_document(uint64_t SA);
        void print_SA();
        // Builds document sets for each run in rlbwt.
        void build_doc_sets();
        // Builds document information for all rows.
        void build_doc_pats();
        // Writes frequencies of document sets to file.
        void write_doc_set_freqs(std::string fname);
        
        // uint64_t naive_lcp(uint64_t row1, uint64_t row2);
        // uint64_t naive_sa(uint64_t bwt_row);
        // bool jump_naive_lcp(uint64_t& idx, uint64_t pointer, char r_char, uint64_t& lcp);

        uint64_t jump_up(uint64_t idx, char c, uint64_t& scan_count);
        uint64_t jump_down(uint64_t idx, char c, uint64_t& scan_count);
        bool jump_thresholds(uint64_t& idx, uint64_t offset, char r_char, uint64_t& scan_count);
        bool jump_randomly(uint64_t& idx, char r_char, uint64_t& scan_count);
        
        void serialize_doc_pats(std::string output_dir);
        void deserialize_doc_pats(std::string index_dir);
        void serialize_doc_sets(std::string output_dir);
        void deserialize_doc_sets(std::string index_dir);
        void serialize(std::string output_dir);
        void deserialize(std::string index_dir);
        
        void print_stats();
        bool check_alphabet(char c);

        std::unordered_map<uint32_t, uint32_t> jumps;
        std::unordered_map<uint32_t, uint32_t> ff_counts;
        std::unordered_map<uint64_t, uint64_t> run_lengths;

        void set_use_doc_pats(bool val) { use_doc_pats = val; }
        int get_num_docs() { return num_docs; }
        uint64_t get_n(uint64_t idx);
        uint64_t get_n_ff(uint64_t idx);
        uint64_t get_offset(uint64_t idx);
        uint64_t get_thresholds(uint64_t idx, uint32_t alphabet_index);
        uint16_t get_rlbwt_thresholds(uint64_t idx, uint16_t i);
        void set_rlbwt_thresholds(uint64_t idx, uint16_t i, uint16_t value);
        void set_onebit();
        
        // Counts of genotype queries outputting each document.
        std::vector<uint32_t> genotype_cnts;
    
        friend class ReadProcessor;
    private:
        // Sorted vector of the start offsets of each document.  
        std::vector<uint64_t> doc_offsets;
        int num_docs;

        // Offset of run heads in the rlbwt. For experimental purposes.
        std::vector<uint64_t> run_offsets;
    
        // Vector of all SA entries. For experiment purposes.
        std::vector<uint64_t> SA_entries;

        // Document sets.
        std::vector<sdsl::bit_vector> unique_doc_sets;
        std::vector<uint32_t> doc_set_inds;

        // Mask for most frequent document sets.
        std::vector<bool> top_X_frequent;

        // Document patterns.
        std::vector<uint8_t> doc_pats;

        // Flag to determine which document method to use
        bool use_doc_pats;
    
        std::vector<uint64_t> first_runs;
        std::vector<uint64_t> first_offsets;
        std::vector<uint64_t> last_runs;
        std::vector<uint64_t> last_offsets;
        bool onebit;
        bool constant;
        uint16_t splitting;
        std::string bwt_string;
        std::string orig_string;
        bool reconstructed;
        uint64_t length;
        uint64_t r;
        uint64_t original_r;
        uint64_t end_bwt_idx;
        uint64_t end_bwt_idx_thresholds[4];
        uint64_t end_bwt_idx_next_up[4];
        uint64_t end_bwt_idx_next_down[4];
        bool verbose;
        bool logs;
	    std::string input_file;
    
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
        uint64_t eof_row;

        // auxilary datastructures for the length, offset and thresholds overflow
        std::vector<uint64_t> n_overflow;
        std::vector<uint64_t> offset_overflow;
        std::vector<std::vector<uint64_t> > thresholds_overflow;
};

class DocSet {
public:
    int size;
    sdsl::bit_vector bv;

    DocSet(int n) : size(n) { bv.resize(n); }
    bool operator==(const DocSet &o) const {
        for (int i = 0; i < size; i++) {
            if (bv[i] != o.bv[i]) {
                return false;
            }
        }
        return true;
    }
};

template<>
struct std::hash<DocSet> {
    const size_t MOD = 1000000007;
    std::size_t operator()(const DocSet &dc) const {
        size_t hash = 0;
        for (int i = 0; i < dc.size; i++) {
            hash = hash * 2 + dc.bv[i];
            if (hash >= MOD) {
                hash -= MOD;
            }
        }
        return hash;
    }
};

#endif
