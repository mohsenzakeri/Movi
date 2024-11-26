#ifndef __READ_PROCESSOR__
#define __READ_PROCESSOR__

#define my_prefetch_r(address) __builtin_prefetch((void *)address, 0, 1)
#define my_prefetch_w(address) __builtin_prefetch((void *)address, 1, 2)

#include "move_structure.hpp"

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

struct Strand {
    Strand() {}
    uint16_t st_length;
    std::string read_name;
    std::string read;
    MoveQuery mq;
    MoveInterval range;
    MoveInterval range_prev;

#if MODE == 5 or MODE == 7
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
};

class ReadProcessor {
    public:
        ReadProcessor(std::string reads_file_name, MoveStructure& mv_, int strands_, bool verbose_, bool reverse_);
        // void process_regular();
        uint64_t initialize_strands(std::vector<Strand>& processes);
        void process_latency_hiding();
        // void ziv_merhav_latency_hiding();
        // void backward_search_latency_hiding();
        void kmer_search_latency_hiding(uint32_t k);
        bool next_read(Strand& process);
        void write_mls(Strand& process);
        void compute_match_count(Strand& process);
        void write_count(Strand& process);
        void process_char(Strand& process);
        bool backward_search(Strand& process, uint64_t end_pos);
        void reset_process(Strand& process);
        void reset_backward_search(Strand& process);
        void reset_kmer_search(Strand& process);
        void next_kmer_search(Strand& process);
        void next_kmer_search_negative_skip_all_heuristic(Strand& process);
        bool verify_kmer(Strand& process, uint64_t k);

#if MODE == 5 or MODE == 7
    void process_char_tally(Strand& process);
    void find_next_id(Strand& process);
    void count_rows_untill_tally(Strand& process);
    void find_tally_b(Strand& process);
    void process_latency_hiding_tally();
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