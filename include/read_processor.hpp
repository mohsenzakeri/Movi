#ifndef READ_PROCESSOR_HPP
#define READ_PROCESSOR_HPP

#include "move_structure.hpp"
#include "utils.hpp"
#include "batch_loader.hpp"

struct Strand {
    Strand() {}
    uint16_t st_length;
    std::string read_name;
    std::string read;
    MoveQuery mq;
    MoveInterval range;
    MoveInterval range_prev;

#if TALLY_MODE
    bool tally_state = false;
    uint64_t tally_b;
    uint64_t rows_until_tally;
    uint64_t next_check_point;
    uint64_t last_id;
    uint64_t run_id;
    bool id_found;
    uint8_t char_index;
    uint16_t tally_offset;
    int find_next_id_attempt;
#endif

    bool kmer_extension;
    bool finished;
    int32_t pos_on_r;
    uint32_t kmer_end;
    int32_t kmer_start;
    uint64_t length_processed;

    uint64_t idx;
    uint64_t offset;
    uint64_t match_len;
    uint64_t match_count;

    uint64_t ff_count;
    uint64_t scan_count;
    std::chrono::time_point<std::chrono::high_resolution_clock> t1;
    std::chrono::time_point<std::chrono::high_resolution_clock> t2;
    std::chrono::time_point<std::chrono::high_resolution_clock> t3;

    // Counts for classification
    std::vector<uint32_t> classify_cnts;
    uint32_t best_doc;
    uint32_t sum_matching_lengths;
    uint32_t second_best_doc;
    uint32_t colors_count;
};

class ReadProcessor {
    public:
        ReadProcessor(std::string reads_file_name, MoveStructure& mv_, int strands_, bool verbose_, bool reverse_);
        // void process_regular();
        uint64_t initialize_strands(std::vector<Strand>& processes, BatchLoader& reader);
        void process_latency_hiding(BatchLoader& reader);
        // void ziv_merhav_latency_hiding();
        // void backward_search_latency_hiding();
        void kmer_search_latency_hiding(uint32_t k, BatchLoader& reader);
        bool next_read(Strand& process, BatchLoader& reader);
        void write_mls(Strand& process);
        void compute_match_count(Strand& process);
        void write_count(Strand& process);
        void process_char(Strand& process);
        void end_process();
        bool backward_search(Strand& process, uint64_t end_pos);
        void reset_process(Strand& process, BatchLoader& reader);
        void reset_backward_search(Strand& process);
        void reset_kmer_search(Strand& process, BatchLoader& reader);
        void next_kmer_search(Strand& process);
        void next_kmer_search_negative_skip_all_heuristic(Strand& process, BatchLoader& reader);
        bool verify_kmer(Strand& process, uint64_t k);

#if TALLY_MODE
        void process_char_tally(Strand& process);
        void find_next_id(Strand& process);
        void count_rows_untill_tally(Strand& process);
        void find_tally_b(Strand& process);
        void process_latency_hiding_tally(BatchLoader& reader);
#endif
    private:
        MoveStructure& mv;
        int cache_line_size;
        int prefetch_step;
        gzFile fp;
        kseq_t *seq;
        int l;
        std::ofstream mls_file;
        std::ofstream matches_file;
        std::ofstream out_file;
        std::ofstream colors_file;
        int strands;
        uint32_t k;
        bool verbose = false;
        bool reverse = false;
        uint64_t read_processed;
        uint64_t total_kmer_count;
        uint64_t positive_kmer_count;
        uint64_t negative_kmer_count;
        uint64_t kmer_extension_count;
        uint64_t kmer_extension_stopped_count;
        uint64_t negative_kmer_extension_count;
        std::ofstream costs_file;
        std::ofstream scans_file;
        std::ofstream fastforwards_file;
        std::chrono::time_point<std::chrono::high_resolution_clock> t1;
};

#endif
