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
    char* read_name;
    std::string read;
    MoveQuery mq;

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
        ReadProcessor(char* reads_file_name, MoveStructure& mv_, int strands_);
        // void process_regular();
        void process_latency_hiding(MoveStructure& mv);
        bool next_read(Strand& process);
        void write_pmls(Strand& process, bool logs);
        void process_char(Strand& process, MoveStructure& mv);
        void reset_process(Strand& process, MoveStructure& mv);
    private:
        // MoveStructure& mv;
        gzFile fp;
        kseq_t *seq;
        int l;
        std::ofstream pmls_file;
        int strands;
        bool verbose = false;
        uint64_t read_processed;

        std::ofstream costs_file;
        std::ofstream scans_file;
        std::ofstream fastforwards_file;
        std::chrono::time_point<std::chrono::high_resolution_clock> t1;
};

#endif