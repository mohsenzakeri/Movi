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

#include "movi_options.hpp"
#include "move_row.hpp"
#include "move_query.hpp"

#define END_CHARACTER 0
#define THRBYTES 5 

struct MoveInterval {
    MoveInterval() {}

    MoveInterval& operator =(const MoveInterval& interval) {
        // Check self assignment
        if (this == &interval)
            return *this;
        run_start = interval.run_start;
        offset_start = interval.offset_start;
        run_end = interval.run_end;
        offset_end = interval.offset_end;
        return *this;
    }

    uint64_t run_start;
    uint64_t offset_start;
    uint64_t run_end;
    uint64_t offset_end;

    MoveInterval(uint64_t run_start_, uint64_t offset_start_, uint64_t run_end_, uint64_t offset_end_) {
        run_start = run_start_;
        offset_start = offset_start_;
        run_end = run_end_;
        offset_end = offset_end_;
    }

    void make_empty() {
        run_start = 1;
        offset_start = 0;
        run_end = 0;
        offset_end = 0;
    }

    bool is_empty() {
        return !((run_start < run_end) or (run_start == run_end and offset_start <= offset_end));
    }

    uint64_t count(std::vector<MoveRow>& rlbwt) {
        uint64_t row_count = 0;
        if (run_start == run_end) {
            row_count = offset_end - offset_start + 1;
        } else {
            row_count = (rlbwt[run_start].get_n() - offset_start) + (offset_end + 1);
            for (uint64_t k = run_start + 1; k < run_end; k ++) {
                row_count += rlbwt[k].get_n();
            }
        }
        return row_count;
    }

    friend std::ostream& operator<<(std::ostream& output, const MoveInterval& mi) {
        // output << "The matching statistics are:\n";
        output << mi.run_start << ":" << mi.offset_start << " --- " << mi.run_end << ":" << mi.offset_end;
        return output;
    }

    bool operator==(const MoveInterval &m) {
        if (run_start == m.run_start and offset_start == m.offset_start and
            run_end == m.run_end and offset_end == m.offset_end)
            return true;
        return false;
    }
};

class MoveStructure {
    public:
        MoveStructure(MoviOptions* movi_options_);
        MoveStructure(MoviOptions* movi_options_, bool onebit_, uint16_t splitting = 0, bool constant = false);

        bool check_mode();
        std::string index_type();
        void set_onebit();
        void build();
        void build_rlbwt();
        uint64_t query_pml(MoveQuery& mq, bool random);
        uint64_t query_backward_search(std::string& query_seq,  int32_t& pos_on_r);
        uint64_t query_zml(MoveQuery& mq);
        uint64_t backward_search(std::string& R,  int32_t& pos_on_r);
        uint64_t backward_search_step(std::string& R,  int32_t& pos_on_r, MoveInterval& interval);
        MoveInterval backward_search(std::string& R,  int32_t& pos_on_r, MoveInterval interval);
        MoveInterval initialize_backward_search(std::string& query_seq,  int32_t& pos_on_r);
        void update_interval(MoveInterval& interval, char next_char);

        void sequential_lf();
        void random_lf();
        std::string reconstruct_lf();

        uint64_t LF(uint64_t row_number);
        uint16_t LF_move(uint64_t& pointer, uint64_t& i);
        uint64_t fast_forward(uint64_t& offset, uint64_t index, uint64_t x);

        uint64_t compute_threshold(uint64_t r_idx, uint64_t pointer, char lookup_char);
        uint32_t compute_index(char row_char, char lookup_char);
        void compute_nexts();
        void compute_ftab();
        void write_ftab();
        void read_ftab();

        // The following are used during development only
        // std::string reconstruct();
        // char compute_char(uint64_t idx);
        // uint64_t naive_lcp(uint64_t row1, uint64_t row2);
        // uint64_t naive_sa(uint64_t bwt_row);
        // bool jump_naive_lcp(uint64_t& idx, uint64_t pointer, char r_char, uint64_t& lcp);

        uint64_t jump_up(uint64_t idx, char c, uint64_t& scan_count);
        uint64_t jump_down(uint64_t idx, char c, uint64_t& scan_count);
#if MODE == 0 or MODE == 1 or MODE == 2
        bool jump_thresholds(uint64_t& idx, uint64_t offset, char r_char, uint64_t& scan_count);
#endif
        bool jump_randomly(uint64_t& idx, char r_char, uint64_t& scan_count);

        void verify_lfs();
        void print_stats();
        void analyze_rows();
        bool check_alphabet(char c);

        void serialize();
        void deserialize();

        char get_char(uint64_t idx);
        uint64_t get_n(uint64_t idx);
        uint64_t get_offset(uint64_t idx);
#if MODE == 0 or MODE == 1 or MODE == 2
        uint64_t get_thresholds(uint64_t idx, uint32_t alphabet_index);
        uint16_t get_rlbwt_thresholds(uint64_t idx, uint16_t i);
        void set_rlbwt_thresholds(uint64_t idx, uint16_t i, uint16_t value);
#endif
        friend class ReadProcessor;
    private:
        MoviOptions* movi_options;
	    bool onebit;
        bool constant;
        uint16_t splitting;

        std::string bwt_string;
        uint64_t length;
        uint64_t r;
        uint64_t original_r;

        // The explicit values for the end bwt row
        uint64_t end_bwt_idx;
        uint64_t eof_row; // This pointer is the same as end_bwt_idx, should be removed.
        uint64_t end_bwt_idx_thresholds[4];
        uint64_t end_bwt_idx_next_up[4];
        uint64_t end_bwt_idx_next_down[4];

        // Map from 2bit encoded character to the actual character
        // Example: alphabet[0] -> A, alphabet[1] -> C
        std::vector<unsigned char> alphabet;
        // Number of each character
        std::vector<uint64_t> counts;
        // Map from the character to the index of the character
        // Example: alphamap[A] -> 0, alphamap[C] -> 1
        std::vector<uint64_t> alphamap;

        // The move structure rows
        std::vector<MoveRow> rlbwt;

        // auxilary datastructures for the length, offset and thresholds overflow
        std::vector<uint64_t> n_overflow;
        std::vector<uint64_t> offset_overflow;
        std::vector<std::vector<uint64_t> > thresholds_overflow;

        // Values used for starting the beckawrd search
        std::vector<uint64_t> first_runs;
        std::vector<uint64_t> first_offsets;
        std::vector<uint64_t> last_runs;
        std::vector<uint64_t> last_offsets;
        std::vector<MoveInterval> ftab;

        std::string orig_string;
        bool reconstructed;

        // Used for gathering statistics
        std::unordered_map<uint32_t, uint32_t> jumps;
        std::unordered_map<uint32_t, uint32_t> ff_counts;
        std::unordered_map<uint64_t, uint64_t> run_lengths;

        // The following are only used for construction, not stored in the index
        sdsl::int_vector<> thresholds;
        std::vector<uint64_t> all_p;
        std::vector<std::unique_ptr<sdsl::bit_vector> > occs;
        std::vector<std::unique_ptr<sdsl::rank_support_v<> > > occs_rank;
        sdsl::bit_vector bits;
        sdsl::rank_support_v<> rbits;
        // sdsl::select_support_mcl<> sbits;
};

#endif
