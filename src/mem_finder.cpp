#include <sys/stat.h> 

#include "move_structure.hpp"
#include "utils.hpp"

// IMPORTANT: For the bidirectional search, the index has to be built with the reverse complement of the reference
uint64_t MoveStructure::query_mems(MoveQuery& mq) {
    int32_t min_mem_length = movi_options->get_min_mem_length();
    if (min_mem_length <= 1) return query_all_mems(mq);

    int32_t pos_on_r = 0; // MEM finding goes left to right
    std::string& query_seq = mq.query();
    #pragma omp atomic
    mem_stats.total_length += query_seq.length();

    size_t ftab_k = movi_options->get_ftab_k();
    // At the time, the initialization works if ftab with ftab_k exists only
    // And the multi-ftab strategy must be turned off
    movi_options->set_multi_ftab(false);
    
    uint64_t mems_found = 0;
    do {
        mems_found += query_mem_bml(mq, pos_on_r, min_mem_length, query_seq, ftab_k);
    } while (pos_on_r < query_seq.length());

    return mems_found;
}

// We search for a MEM in P[pos_on_r..m-1] starting at pos_on_r
// Returns true if a MEM is found, false otherwise, and updates pos_on_r
bool MoveStructure::query_mem_bml(MoveQuery& mq, int32_t& pos_on_r, int32_t& min_mem_length, std::string& query_seq, size_t& ftab_k) {
    // Reached end of query sequence
    if (pos_on_r + min_mem_length > query_seq.length()) {
        pos_on_r = query_seq.length();
        return false;
    }

    uint64_t match_len = 0;
    int32_t init_pos = pos_on_r + min_mem_length - 1;
    MoveBiInterval bi_interval;
    bi_interval = initialize_bidirectional_search(mq, init_pos, match_len);
    #pragma omp atomic
    mem_stats.ftab_bidirectional_extensions += match_len;
    if (match_len <= 1) {
        #pragma omp atomic
        ++mem_stats.ftab_bidirectional_fail;
    }
    // If min_mem_length is at least ftab_k, we skip using BML
    bool ftab_skip = match_len <= 1 && ftab_k <= min_mem_length;
    --init_pos; // Move to next character left of initial match range

    if (ftab_skip) {
        #pragma omp atomic
        ++mem_stats.ftab_skips;
        // If BML skip is known from ftab, find next left end with only backward extension (extend_left slow without ftab)
        MoveInterval fw_interval = bi_interval.fw_interval;
        for (size_t j = 0; j <= init_pos - pos_on_r; ++j) {
            if (!backward_search_step(query_seq[init_pos - j], fw_interval)) {
                pos_on_r = init_pos - j + 1;
                #pragma omp atomic
                ++mem_stats.bml_skips;
                #pragma omp atomic
                mem_stats.bml_chars_skipped += j - (init_pos - pos_on_r);
                return false;
            }
            #pragma omp atomic
            ++mem_stats.backward_extensions;
            ++match_len;
        }
        // Should never reach here
        std::cerr << "Extended past failed ftab\n";
        exit(0);
    }
    else {
        // Backward extension to find if left end permits a sufficiently long MEM, should equal pos_on_r if true
        for (size_t j = 0; j <= init_pos - pos_on_r; ++j) {
            if (!extend_left(query_seq[init_pos - j], bi_interval)) {
                pos_on_r = init_pos - j + 1;
                #pragma omp atomic
                ++mem_stats.bml_skips;
                #pragma omp atomic
                mem_stats.bml_chars_skipped += j - (init_pos - pos_on_r);
                return false;
            }
            #pragma omp atomic
            ++mem_stats.bidirectional_left_extensions;
            ++match_len;
        }
    }

    // Forward extension to find right end of MEM (exclusive)
    MoveInterval rc_interval = bi_interval.rc_interval; // We don't need to keep the fw_interval
    MoveInterval rc_interval_before_extension = rc_interval;
    size_t i;
    for (i = pos_on_r + min_mem_length; i < query_seq.length(); ++i) {
        rc_interval_before_extension = rc_interval;
        if (!forward_search_step(query_seq[i], rc_interval)) {
            rc_interval = rc_interval_before_extension;
            break;
        }
        #pragma omp atomic
        ++mem_stats.forward_extensions;
        ++match_len;
    }

    #pragma omp atomic
    ++mem_stats.mems;
    #pragma omp atomic
    mem_stats.total_counts += rc_interval.count(rlbwt);
    mq.add_mem(pos_on_r, i, rc_interval.count(rlbwt));

    // Backward extension to find next left end of candidate MEM (inclusive)
    // TODO: if this extension is at least L, we can skip the first BML step in next recursion
    size_t end_pos_on_r = i;
    size_t j;
    if (end_pos_on_r < query_seq.length()) {
        init_pos = end_pos_on_r;
        match_len = 0;
        MoveInterval fw_interval = initialize_backward_search(mq, init_pos, match_len);
        ++match_len; // TODO: this is still off by one for initialize_backward_search
        #pragma omp atomic
        mem_stats.ftab_backward_extensions += match_len;
        if (match_len <= 1) {
            #pragma omp atomic
            ++mem_stats.ftab_backward_fail;
        }
        --init_pos; // Move to next character left of ftab k-mer or character left of initial range if ftab fails
        for (i = 0; i <= init_pos - (pos_on_r + 1); ++i) {
            if (!backward_search_step(query_seq[init_pos - i], fw_interval)) {
                break;
            }
            #pragma omp atomic
            ++mem_stats.backward_extensions;
            ++match_len;
        }
    }

    pos_on_r = (init_pos - i) + 1;
    return true;
}

uint64_t MoveStructure::query_all_mems(MoveQuery& mq) {
    std::string& query_seq = mq.query();
    #pragma omp atomic
    mem_stats.total_length += query_seq.length();
    
    size_t ftab_k = movi_options->get_ftab_k();
    // At the time, the initialization works if ftab with ftab_k exists only
    // And the multi-ftab strategy must be turned off
    movi_options->set_multi_ftab(false);
    
    size_t s = 0; // start, inclusive
    size_t e = 0; // end, exclusive
    uint64_t match_len = 0; // bases matched so far
    int32_t init_pos = s; 
    MoveBiInterval bi_interval = initialize_bidirectional_search(mq, init_pos, match_len);
    #pragma omp atomic
    mem_stats.ftab_bidirectional_extensions += match_len;
    if (match_len <= 1) {
        #pragma omp atomic
        ++mem_stats.ftab_bidirectional_fail;
    }

    // Forward extension to find right end of MEM (exclusive)
    while (s < query_seq.length()) {
        // Extend MEM until mismatch
        MoveBiInterval bi_interval_before_extension = bi_interval;
        while (s + match_len < query_seq.length() && extend_right(query_seq[s + match_len], bi_interval)) {
            bi_interval_before_extension = bi_interval;
            ++match_len;
            #pragma omp atomic
            ++mem_stats.bidirectional_right_extensions;
        }
        e = s + match_len;
        #pragma omp atomic
        mem_stats.total_counts += bi_interval_before_extension.fw_interval.count(rlbwt);
        mq.add_mem(s, e, bi_interval_before_extension.fw_interval.count(rlbwt));

        // Backward extension to find next MEM start
        match_len = 0;
        if (e < query_seq.length()) {
            init_pos = e;
            bi_interval = initialize_bidirectional_search(mq, init_pos, match_len);
            #pragma omp atomic
            mem_stats.ftab_bidirectional_extensions += match_len;
            if (match_len <= 1) {
                #pragma omp atomic
                ++mem_stats.ftab_bidirectional_fail;
            }
            while (extend_left(query_seq[e - match_len], bi_interval)) {
                bi_interval_before_extension = bi_interval;
                ++match_len;
                #pragma omp atomic
                ++mem_stats.bidirectional_left_extensions;
            }
            bi_interval = bi_interval_before_extension;
        }
        s = e - match_len + 1;
    }

    return mq.get_mems().size();
}