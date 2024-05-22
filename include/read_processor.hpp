#ifndef __READ_PROCESSOR__
#define __READ_PROCESSOR__

#define my_prefetch_r(address) __builtin_prefetch((void *)address, 0, 1)
#define my_prefetch_w(address) __builtin_prefetch((void *)address, 1, 2)

#include "move_structure.hpp"

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

struct Range {
    Range() {}
    Range& operator =(const Range& range) {
        run_start = range.run_start;
        offset_start = range.offset_start;
        run_end = range.run_end;
        offset_end = range.offset_end;
        return *this;
    }
    uint64_t run_start;
    uint64_t offset_start;
    uint64_t run_end;
    uint64_t offset_end;
};

struct Strand {
    Strand() {}
    uint16_t st_length;
    std::string read_name;
    std::string read;
    MoveQuery mq;
    Range range;
    Range range_prev;

    bool finished;
    int32_t pos_on_r;
    uint64_t length_processed;

    uint64_t idx;
    uint64_t offset;
    uint64_t match_len;

    uint64_t ff_count;
    uint64_t scan_count;
    std::chrono::time_point<std::chrono::high_resolution_clock> t1;
    std::chrono::time_point<std::chrono::high_resolution_clock> t2;
    std::chrono::time_point<std::chrono::high_resolution_clock> t3;
};

class ReadProcessor {
    public:
        ReadProcessor(std::string reads_file_name, MoveStructure& mv_, int strands_, bool query_pml, bool reverse_);
        // void process_regular();
        void process_latency_hiding(MoveStructure& mv);
        void backward_search_latency_hiding(MoveStructure& mv);
        bool next_read(Strand& process);
        void write_pmls(Strand& process, bool logs, bool write_stdout);
        void process_char(Strand& process, MoveStructure& mv);
        bool backward_search(Strand& process, MoveStructure& mv, uint64_t& match_count);
        void reset_process(Strand& process, MoveStructure& mv);
        void reset_backward_search(Strand& process, MoveStructure& mv);
    private:
        gzFile fp;
        kseq_t *seq;
        int l;
        std::ofstream pmls_file;
        std::ofstream matches_file;
        int strands;
        bool verbose = false;
        bool reverse = false;
        uint64_t read_processed;

        std::ofstream costs_file;
        std::ofstream scans_file;
        std::ofstream fastforwards_file;
        std::chrono::time_point<std::chrono::high_resolution_clock> t1;
};

#endif