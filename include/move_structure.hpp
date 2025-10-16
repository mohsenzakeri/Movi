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
#include <map>
#include <sstream>
#include <filesystem>
#include <span>
#include <sys/stat.h>


#include <span>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "fastcluster.h"

#include "sdsl_wrapper.hpp"

#include "movi_options.hpp"
#include "move_row.hpp"
#include "move_row_colored.hpp"
#include "move_query.hpp"
#include "sequitur.hpp"
#include "utils.hpp"
#include "move_intervals.hpp"
#include "doc_set.hpp"

class Classifier;

struct ThresholdsRow {
    uint16_t values[4];
};

class MoveStructure {
    public:
        MoveStructure(MoviOptions* movi_options_);
        MoveStructure(MoviOptions* movi_options_, uint16_t nt_splitting_, bool constant_);

        std::vector<MoveRow> get_rlbwt();
        char get_char(uint64_t idx);
        uint64_t get_n(uint64_t idx);
        uint64_t get_offset(uint64_t idx);
        uint64_t get_id(uint64_t idx);
        bool use_separator();
        // TODO: The following is useful for mmaping, there is a slowdown though
        // MoveRow& get_move_row(uint64_t idx);

#if USE_THRESHOLDS
        uint64_t get_thresholds(uint64_t idx, uint32_t alphabet_index);
#endif

#if SPLIT_THRESHOLDS_FALSE
        uint16_t get_rlbwt_thresholds(uint64_t idx, uint16_t i);
        void set_rlbwt_thresholds(uint64_t idx, uint16_t i, uint16_t value);
#endif

        uint64_t get_SA_entries(uint64_t idx, uint64_t offset);

        uint64_t LF(uint64_t row_number, uint64_t alphabet_index);
        uint64_t LF_heads(uint64_t run_number, uint64_t alphabet_index);
        // The 3rd argument of LF_move is used in the latency_hiding_tally mode
        uint16_t LF_move(uint64_t& pointer, uint64_t& i, uint64_t id = std::numeric_limits<uint64_t>::max());
        uint64_t fast_forward(uint64_t& offset, uint64_t index, uint64_t x);

        bool check_alphabet(char& c);
        uint32_t compute_index(char row_char, char lookup_char);

        void analyze_rows();
        void print_stats();
        void print_basic_index_info();

/***************************************************************************/
/****  Beginning of functions implemented in move_structure_search.cpp  ****/
        MoveInterval try_ftab(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len, size_t ftab_k, bool rc = false);
        bool look_ahead_ftab(MoveQuery& mq, uint32_t pos_on_r, int32_t& step);
        bool look_ahead_backward_search(MoveQuery& mq, uint32_t pos_on_r, int32_t step);

        uint64_t query_backward_search(MoveQuery& mq, int32_t& pos_on_r);
        bool backward_search_step(char c, MoveInterval& interval);
        uint64_t backward_search_step(std::string& R, int32_t& pos_on_r, MoveInterval& interval);
        MoveInterval backward_search(std::string& R, int32_t& pos_on_r, MoveInterval interval, int32_t max_length);
        MoveInterval initialize_backward_search(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len, bool rc = false);
        void update_interval(MoveInterval& interval, char next_char);

        bool extend_bidirectional(char c_, MoveInterval& fw_interval, MoveInterval& rc_interval);
        bool extend_left(char c, MoveBiInterval& bi_interval);
        bool extend_right(char c, MoveBiInterval& bi_interval);
        MoveBiInterval backward_search_bidirectional(std::string& R, int32_t& pos_on_r, MoveBiInterval interval, int32_t max_length);
        MoveBiInterval initialize_bidirectional_search(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len);
/*******  End of functions implemented in move_structure_search.cpp  *******/
/***************************************************************************/

/***************************************************************************/
/****  Beginning of functions implemented in move_structure_build.cpp   ****/
        void build();
        uint64_t compute_length_from_bwt();
        void initialize_bits();
        void fill_bits_by_thresholds();
        void detect_move_row_boundaries();
        void add_detected_run(uint64_t scanned_bwt_length, uint64_t run_char, uint16_t& run_length);
        void find_run_heads_information();
        void build_move_rows();
        void find_base_interval_data();
        void build_alphabet(std::vector<uint64_t>& all_possible_chars);

#if USE_THRESHOLDS
        void compute_thresholds();
        // Helper for threshold calculation
        void debug_threshold_calculation(uint64_t i, uint64_t thr_i, uint64_t j, char rlbwt_c,
                                         const std::vector<uint64_t>& alphabet_thresholds);
        void set_threshold_for_one_character(uint64_t i, uint16_t threshold_index, uint64_t value,
                                             uint16_t value_split, std::vector<uint64_t>& current_thresholds);
#endif

#if USE_NEXT_POINTERS
        void compute_nexts();
#endif

#if BLOCKED_MODES
        void compute_blocked_ids(std::vector<uint64_t>& raw_ids);
#endif

        void build_ftab();

        // Finds SA entries of all rows in BWT.
        void find_sampled_SA_entries();
/*******  End of functions implemented in move_structure_build.cpp   *******/
/***************************************************************************/

/***************************************************************************/
/****  Beginning of functions implemented in move_structure_query.cpp   ****/
        void verify_lfs();
        void verify_lf_loop();
        void sequential_lf();
        void random_lf();
        void reconstruct_lf();

        uint64_t query_pml(MoveQuery& mq);
        uint64_t query_zml(MoveQuery& mq);

        uint64_t reposition_up(uint64_t idx, char c, uint64_t& scan_count);
        uint64_t reposition_down(uint64_t idx, char c, uint64_t& scan_count);
        bool reposition_randomly(uint64_t& idx, uint64_t& offset, char r_char, uint64_t& scan_count);

#if USE_THRESHOLDS
        bool reposition_thresholds(uint64_t& idx, uint64_t offset, char r_char, uint64_t& scan_count);
        // Helper for threshold calculation
        void handle_reposition_up(uint64_t& idx, uint64_t saved_idx, char r_char, uint16_t next_up, uint64_t& scan_count);
        void handle_reposition_down(uint64_t& idx, uint64_t saved_idx, char r_char, uint16_t next_down, uint64_t& scan_count);
#endif
/*******  End of functions implemented in move_structure_query.cpp  *******/
/***************************************************************************/

/***************************************************************************/
/**********  Beginning of functions implemented in sequitur.cpp  ***********/
        void query_all_kmers(MoveQuery& mq, bool kmer_counts = false);
        uint64_t query_kmers_from_bidirectional(MoveQuery& mq, int32_t& pos_on_r);
        uint64_t query_kmers_from(MoveQuery& mq, int32_t& pos_on_r, bool single = false);
/*******  End of functions implemented in sequitur.cpp  ***********/
/***************************************************************************/

/***************************************************************************/
/****  Beginning of functions implemented in move_structure_color.cpp  *****/
        int get_num_docs() { return num_docs; }
        int get_num_species() { return num_species; }

        // Fill the run offsets array (used for building colors among other things)
        void fill_run_offsets();
        // Builds document sets for each run in rlbwt.
        void build_doc_sets();
        uint32_t hash_collapse(std::unordered_map<DocSet, uint32_t> &keep_set, DocSet &bv);
        void build_tree_doc_sets();
        void build_doc_set_similarities();
        void compress_doc_sets();
        void compute_color_ids_from_flat();
        // Finds documents corresponding to rows in BWT.
        void build_doc_pats();
        // Initialize classify counts
        void initialize_classify_cnts();
        void set_classifier(Classifier *cl) { classifier = cl; }
        void set_output_files(OutputFiles *of) { output_files = of; }

        // Document tree functions
        bool is_ancestor(uint16_t x, uint16_t y);
        uint16_t LCA(uint16_t x, uint16_t y);
        void dfs_times(uint16_t cur, uint16_t &t);
        void compute_run_lcs();

        // The following is just used to test how storing colors in the rlbwt affects the performance
        void add_colors_to_rlbwt();
/******** End of functions implemented in move_structure_color.cpp *********/
/***************************************************************************/

/***************************************************************************/
/*****  Beginning of functions implemented in move_structure_io.cpp  *******/

        // Writes frequencies of document sets to file.
        void write_doc_set_freqs(std::string fname);

        void flat_and_serialize_colors_vectors();
        void deserialize_doc_sets_flat();

        void serialize_doc_pats(std::string fname);
        void deserialize_doc_pats(std::string fname);
        void serialize_doc_sets(std::string fname);
        void deserialize_doc_sets(std::string fname);
        void load_document_info();

        // The following two methods are not implemented
        void serialize_doc_rows();
        void deserialize_doc_rows();

        std::ifstream open_index_read();
        std::ofstream open_index_write();
        void read_index_header(std::ifstream& fin);
        void write_index_header(std::ofstream& fout);
        void write_basic_index_data(std::ofstream& fout);
        void read_basic_index_data(std::ifstream& fin);
        void write_overflow_tables(std::ofstream& fout);
        void read_overflow_tables(std::ifstream& fin);
        void write_counts_data(std::ofstream& fout);
        void read_counts_data(std::ifstream& fin);
        void write_main_table(std::ofstream& fout);
        void read_main_table(std::ifstream& fin, std::streamoff rlbwt_offset);
        void write_separators_thresholds(std::ofstream& fout);
        void read_separators_thresholds(std::ifstream& fin);
        void serialize();
        void deserialize();
        void output_ids();

#if BLOCKED_MODES
        void write_id_blocks(std::ofstream& fout);
        void read_id_blocks(std::ifstream& fin);
#endif

#if TALLY_MODES
        void write_tally_table(std::ofstream& fout);
        void read_tally_table(std::ifstream& fin);
#endif

        void serialize_sampled_SA();
        void deserialize_sampled_SA();

        void write_ftab();
        void read_ftab();
/*********  End of functions implemented in move_structure_io.cpp  *********/
/***************************************************************************/

	    KmerStatistics kmer_stats;
        friend class ReadProcessor;
    private:
        // Reference to output files for writing results
        OutputFiles* output_files;

        // The BWT file used for building the index
        std::ifstream bwt_file;

        // Sorted vector of the start offsets of each document.  
        std::vector<uint64_t> doc_offsets;
        // Species ID for each document
        std::vector<uint32_t> doc_ids;
        // Map from taxon id to compressed species index
        std::map<uint32_t, uint32_t> taxon_id_compress;
        // Compressed species index to taxon id
        std::vector<uint32_t> to_taxon_id;
        // log length of each species
        std::vector<double> log_lens;
        uint32_t num_docs;
        uint32_t num_species;
        uint64_t num_colors = 0;

        // Offset of run heads in the rlbwt.
        std::vector<uint64_t> run_offsets;

        // Document sets.
        std::vector<uint16_t> flat_colors;
        std::unordered_map<uint64_t, uint32_t> color_offset_to_id;
        std::vector<std::vector<uint16_t>> unique_doc_sets;
        std::vector<uint32_t> doc_set_inds;
        std::vector<MoveTally> doc_set_flat_inds;
        sdsl::bit_vector compressed;

        // Tree over documents
        std::vector<std::vector<uint16_t>> tree;
        std::vector<std::vector<uint16_t>> tree_doc_sets;
        std::vector<std::vector<uint16_t>> bin_lift;
        std::vector<uint16_t> t_in;
        std::vector<uint16_t> t_out;
        
        // Count of how much each doc set appears (by ID).
        std::vector<uint64_t> doc_set_cnts;

        // Document patterns (species that each row in BWT belongs to).
        std::vector<uint16_t> doc_pats;

        // Counts for multi classification
        std::vector<uint32_t> classify_cnts;
        std::vector<double> doc_scores;
        const double log4 = log(4);

        // Classifier object for binary classification
        Classifier *classifier;
    
        // Basic index configurations
        MoviOptions* movi_options;
        uint16_t nt_splitting;
        bool constant;

        // Vector of sampled SA entries. For experiment purposes.
        std::vector<uint64_t> sampled_SA_entries;

        // Basic index characteristics
        std::string bwt_string;
        uint64_t length;
        uint64_t r;
        uint64_t original_r;

        void compute_number_of_build_steps();
        uint64_t total_build_steps = 0;
        uint64_t current_build_step = 0;

        // The explicit values for the end bwt row
        uint64_t end_bwt_idx;
        uint64_t end_bwt_idx_thresholds[4];
        uint64_t end_bwt_idx_next_up[4];
        uint64_t end_bwt_idx_next_down[4];

        // The explicit thresholds for the separators
        std::vector<ThresholdsRow> separators_thresholds;
        std::unordered_map<uint64_t, uint64_t> separators_thresholds_map;

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
        std::vector<MoveRowColored> rlbwt_colored;
#if TALLY_MODES
        uint32_t tally_checkpoints;
        std::vector<std::vector<MoveTally>> tally_ids;
#endif

#if BLOCKED_MODES
        std::vector<std::vector<uint32_t>> id_blocks;
        uint64_t block_size = BLOCK_SIZE;
#endif

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

        // Used for gathering statistics
        uint64_t no_ftab;
        uint64_t all_initializations;
        std::unordered_map<uint32_t, uint32_t> repositions;
        std::unordered_map<uint32_t, uint32_t> ff_counts;
        std::unordered_map<uint64_t, uint64_t> run_lengths;

        // The following are only used for construction, not stored in the index
        std::vector<uint64_t> thresholds;
        std::vector<uint64_t> all_p;
        std::vector<char> heads;
        std::vector<char> original_run_heads;
        std::vector<std::unique_ptr<sdsl::bit_vector> > occs;
        std::vector<std::unique_ptr<sdsl::rank_support_v<> > > occs_rank;
        std::vector<uint64_t> heads_rank;
        std::vector<uint64_t> lens;
        std::vector<uint32_t> original_lens;

        sdsl::bit_vector bits;
        sdsl::rank_support_v<> rbits;
};

#endif