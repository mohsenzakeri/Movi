
#include <cstdint>
#include <zlib.h>
#include <stdio.h>
#include <cstdio>
#include <chrono>
#include <cstddef>
#include <unistd.h>
#include <sys/stat.h>

#include "kseq.h"
#include <sdsl/int_vector.hpp>
#include "cxxopts.hpp"

#include "move_structure.hpp"
#include "move_query.hpp"
#include "read_processor.hpp"
#include "movi_options.hpp"

int main(int argc, char** argv) {
    std::cerr << "\n\n**** The ftab with ftab-k = 12 has to be built.\n";
    std::cerr << "**** Only works with a compact index.\n\n\n";

    MoviOptions movi_options;
    movi_options.set_index_dir(std::string(argv[1]));
    movi_options.set_kmer();
    movi_options.set_no_header(true);

    MoveStructure mv_(&movi_options);
    mv_.deserialize();
    auto rlbwt = mv_.get_rlbwt();
    mv_.print_stats();

    movi_options.set_ftab_k(12);
    mv_.read_ftab();
    std::cerr<<"Ftab was read!\n";

    // Test 1, with initialization
    std::string R_ = "AATCGCACCGGACATCAAATGACAAAGAACTTATACGGTGGGACAAAGACTATGCTAAAA";
    MoveQuery mq(R_);
    uint64_t matched_len = 0;
    std::cerr << R_ << "\n";
    int match_right_end = 40;
    int32_t p = match_right_end;
    MoveBiInterval bi_init_test = mv_.initialize_bidirectional_search(mq, p, matched_len);
    int match_left_end = p;

    std::cerr << "\n\n**** Test 1, with initialization\n\n";
    std::cerr << R_.substr(match_left_end, matched_len - 1) << "\t"
              << bi_init_test.fw_interval << "\t" << bi_init_test.rc_interval << "\t"
              << bi_init_test.fw_interval.count(rlbwt) << "\n";
    while (match_right_end < R_.size() - 1 and  mv_.extend_right(R_[match_right_end + 1], bi_init_test)) {
        matched_len++;
        match_right_end++;
        std::cerr << match_left_end << " " << bi_init_test.match_len << " ";
        std::cerr << R_.substr(match_left_end, matched_len - 1) << "\t"
                  << bi_init_test.fw_interval << "\t" << bi_init_test.rc_interval << "\t"
                  << bi_init_test.fw_interval.count(rlbwt) << "\n";
    }

    while (match_left_end > 0 and mv_.extend_left(R_[match_left_end - 1], bi_init_test)) {        
        match_left_end--;
        matched_len++;
        std::cerr << match_left_end << " " << bi_init_test.match_len << " ";
        std::cerr << R_.substr(match_left_end, matched_len - 1) << "\t"
                  << bi_init_test.fw_interval << "\t" << bi_init_test.rc_interval << "\t"
                  << bi_init_test.fw_interval.count(rlbwt) << "\n";
    }



    // Test 2, without the initialization
    std::cerr << "\n\n**** Test 2, without the initialization\n\n";
    R_ = "GTACCGGGAC";
    MoveBiInterval bi_interval;
    bi_interval.fw_interval.run_start = 13;
    bi_interval.fw_interval.offset_start = 1;
    bi_interval.fw_interval.run_end = 15;
    bi_interval.fw_interval.offset_end = 1;

    bi_interval.rc_interval.run_start = 22;
    bi_interval.rc_interval.offset_start = 0;
    bi_interval.rc_interval.run_end = 23;
    bi_interval.rc_interval.offset_end = 3;
    std::cerr << "CC:\t" << bi_interval.fw_interval << "\t" << bi_interval.rc_interval << "\n";
    mv_.extend_right('G', bi_interval);
    std::cerr << "CCG:\t" << bi_interval.fw_interval << "\t" << bi_interval.rc_interval << "\n";
    mv_.extend_left('A', bi_interval);
    //std::cerr << "ACC:\t" << bi_interval.fw_interval << "\t" << bi_interval.rc_interval << "\n";
    std::cerr << "ACCG:\t" << bi_interval.fw_interval << "\t" << bi_interval.rc_interval << "\n";
    mv_.extend_right('G', bi_interval);
    std::cerr << "ACCGG:\t" << bi_interval.fw_interval << "\t" << bi_interval.rc_interval << "\n";
    mv_.extend_left('T', bi_interval);
    std::cerr << "TACCGG:\t" << bi_interval.fw_interval << "\t" << bi_interval.rc_interval << "\n";
    mv_.extend_left('G', bi_interval);
    std::cerr << "GTACCGG:\t" << bi_interval.fw_interval << "\t" << bi_interval.rc_interval << "\n";
    mv_.extend_right('G', bi_interval);
    std::cerr << "GTACCGGG:\t" << bi_interval.fw_interval << "\t" << bi_interval.rc_interval << "\n";
    mv_.extend_right('A', bi_interval);
    std::cerr << "GTACCGGGA:\t" << bi_interval.fw_interval << "\t" << bi_interval.rc_interval << "\n";
    mv_.extend_right('C', bi_interval);
    std::cerr << "GTACCGGGAC:\t" << bi_interval.fw_interval << "\t" << bi_interval.rc_interval << "\n";
    std::cerr << "---\n";
    if (mv_.extend_right('C', bi_interval))
        std::cerr << "GTACCGGGAC:\t" << bi_interval.fw_interval << "\t" << bi_interval.rc_interval << "\n";
    else
        std::cerr << "Extension not possible.\n";
}