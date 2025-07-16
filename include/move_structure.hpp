#ifndef MOVE_STRUCTURE_HPP
#define MOVE_STRUCTURE_HPP

#include <fstream>
#include <cstdint>
#include <stdio.h>
#include <chrono>
#include <cstddef>
#include <unistd.h>
#include <sys/stat.h>
#include <vector>
#include <unordered_map>

#include <span>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>

#include "movi_options.hpp"
#include "move_row.hpp"
#include "move_query.hpp"
#include "sequitur.hpp"
#include "mem_finder.hpp"
#include "utils.hpp"

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

struct MoveBiInterval {
    MoveInterval fw_interval;
    MoveInterval rc_interval;
    uint64_t match_len;

    MoveBiInterval() {
        fw_interval = MoveInterval(1, 0, 0, 0);
        rc_interval = MoveInterval(1, 0, 0, 0);
        match_len = 0;
    }

    friend std::ostream& operator<<(std::ostream& output, const MoveBiInterval& mi_bi) {
        output << mi_bi.match_len << ": " << mi_bi.fw_interval << "\t" << mi_bi.rc_interval;
        return output;
    }
};

class MoveStructure {
    public:
        MoveStructure(MoviOptions* movi_options_);
        MoveStructure(MoviOptions* movi_options_, uint16_t splitting, bool constant);

        void build();
        void fill_bits_by_thresholds();
        void build_rlbwt();
        uint64_t query_pml(MoveQuery& mq);
        uint64_t query_backward_search(MoveQuery& mq, int32_t& pos_on_r);
        uint64_t query_zml(MoveQuery& mq);

        void query_all_kmers(MoveQuery& mq, bool kmer_counts = false);
        uint64_t query_kmers_from_bidirectional(MoveQuery& mq, int32_t& pos_on_r);
        uint64_t query_kmers_from(MoveQuery& mq, int32_t& pos_on_r, bool single = false);

        uint64_t query_mems(MoveQuery& mq);
        bool query_mem_bml(MoveQuery& mq, int32_t& pos_on_r,int32_t& min_mem_length, std::string& query_seq, size_t& ftab_k);
        uint64_t query_all_mems(MoveQuery& mq);

        MoveInterval try_ftab(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len, size_t ftab_k, bool rc = false);
        bool look_ahead_ftab(MoveQuery& mq, uint32_t pos_on_r, int32_t& step);
        bool look_ahead_backward_search(MoveQuery& mq, uint32_t pos_on_r, int32_t step);

        bool extend_bidirectional(char c_, MoveInterval& fw_interval, MoveInterval& rc_interval);
        bool extend_left(char c, MoveBiInterval& bi_interval);
        bool extend_right(char c, MoveBiInterval& bi_interval);
        MoveBiInterval backward_search_bidirectional(std::string& R, int32_t& pos_on_r, MoveBiInterval interval, int32_t max_length);
        MoveBiInterval initialize_bidirectional_search(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len);

        bool backward_search_step(char c, MoveInterval& interval);
        bool forward_search_step(char c, MoveInterval& rc_interval);
        uint64_t backward_search_step(std::string& R, int32_t& pos_on_r, MoveInterval& interval);
        MoveInterval backward_search(std::string& R, int32_t& pos_on_r, MoveInterval interval, int32_t max_length);
        MoveInterval initialize_backward_search(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len, bool rc = false);
        void update_interval(MoveInterval& interval, char next_char);

        void sequential_lf();
        void random_lf();
        std::string reconstruct_lf();

        uint64_t LF(uint64_t row_number, uint64_t alphabet_index);
        // The 3rd argument of LF_move is used in the latency_hiding_tally mode
        uint16_t LF_move(uint64_t& pointer, uint64_t& i, uint64_t id = std::numeric_limits<uint64_t>::max());
        uint64_t fast_forward(uint64_t& offset, uint64_t index, uint64_t x);

        uint64_t compute_threshold(uint64_t r_idx, uint64_t pointer, char lookup_char);
        uint32_t compute_index(char row_char, char lookup_char);

#if USE_NEXT_POINTERS
        void compute_nexts();
#endif

        void compute_ftab();
        void write_ftab();
        void compute_run_lcs();
        void read_ftab();

        // Finds SA entries of all rows in BWT.
        void find_sampled_SA_entries();
        uint64_t get_SA_entries(uint64_t idx, uint64_t offset);

        // The following are used during development only
        // std::string reconstruct();
        // char compute_char(uint64_t idx);
        // uint64_t naive_lcp(uint64_t row1, uint64_t row2);
        // uint64_t naive_sa(uint64_t bwt_row);
        // bool jump_naive_lcp(uint64_t& idx, uint64_t pointer, char r_char, uint64_t& lcp);

        uint64_t reposition_up(uint64_t idx, char c, uint64_t& scan_count);
        uint64_t reposition_down(uint64_t idx, char c, uint64_t& scan_count);
        bool reposition_randomly(uint64_t& idx, char r_char, uint64_t& scan_count);

        void verify_lfs();
        void print_stats();
        void print_ids();
        void analyze_rows();
        bool check_alphabet(char& c);

        void serialize();
        void deserialize();
        void serialize_sampled_SA();
        void deserialize_sampled_SA();

        char get_char(uint64_t idx);
        uint64_t get_n(uint64_t idx);
        uint64_t get_offset(uint64_t idx);
        uint64_t get_id(uint64_t idx);
        MoveRow& get_move_row(uint64_t idx);
#if USE_THRESHOLDS
        void compute_thresholds();
        bool reposition_thresholds(uint64_t& idx, uint64_t offset, char r_char, uint64_t& scan_count);
        uint64_t get_thresholds(uint64_t idx, uint32_t alphabet_index);
#endif
#if SPLIT_THRESHOLDS_FALSE
        uint16_t get_rlbwt_thresholds(uint64_t idx, uint16_t i);
        void set_rlbwt_thresholds(uint64_t idx, uint16_t i, uint16_t value);
#endif
        KmerStatistics kmer_stats;
        friend class ReadProcessor;
        std::vector<MoveRow> get_rlbwt();
    private:
        MoviOptions* movi_options;
	    bool onebit; // This is not used any more as the onebit modes is deprecated
        bool constant;
        uint16_t splitting;

        // Vector of sampled SA entries. For experiment purposes.
        std::vector<uint64_t> sampled_SA_entries;

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
        std::span<MoveRow> rlbwt_view;
#if TALLY_MODE
        uint32_t tally_checkpoints;
        std::vector<std::vector<MoveTally>> tally_ids;
#endif
        std::vector<std::vector<uint32_t>> id_blocks;

        // auxilary datastructures for the length, offset and thresholds overflow
        std::vector<uint64_t> n_overflow;
        std::vector<uint64_t> offset_overflow;
        std::vector<std::vector<uint64_t> > thresholds_overflow;

        // Values used for starting the beckawrd search
        std::vector<uint64_t> first_runs;
        std::vector<uint64_t> first_offsets;
        std::vector<uint64_t> last_runs;
        std::vector<uint64_t> last_offsets;
        std::vector<std::vector<MoveInterval> > ftabs;
        std::vector<MoveInterval> ftab;

        std::string orig_string;
        bool reconstructed;

        // Used for gathering statistics
        uint64_t no_ftab;
        uint64_t all_initializations;
        std::unordered_map<uint32_t, uint32_t> repositions;
        std::unordered_map<uint32_t, uint32_t> ff_counts;
        std::unordered_map<uint64_t, uint64_t> run_lengths;

        // The following are only used for construction, not stored in the index
        sdsl::int_vector<> thresholds;
        std::vector<uint64_t> all_p;
        std::vector<char> heads;
        std::vector<std::unique_ptr<sdsl::bit_vector> > occs;
        std::vector<std::unique_ptr<sdsl::rank_support_v<> > > occs_rank;
        sdsl::bit_vector bits;
        sdsl::rank_support_v<> rbits;
        // sdsl::select_support_mcl<> sbits;
};

#endif
