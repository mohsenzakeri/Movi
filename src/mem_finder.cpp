#include <sys/stat.h> 

#include "move_structure.hpp"
#include "utils.hpp"

// IMPORTANT: For the bidirectional search, the index has to be built with the reverse complement of the reference
uint64_t MoveStructure::query_mems(MoveQuery& mq) {
    int32_t min_mem_length = movi_options->get_min_mem_length();
    // if (min_mem_length <= 1) return query_all_mems(mq);

    int32_t pos_on_r = 0;
    std::string& query_seq = mq.query();

    size_t ftab_k = movi_options->get_ftab_k();
    // At the time, the initialization works if ftab with ftab_k exists only
    // And the multi-ftab strategy must be turned off
    movi_options->set_multi_ftab(false);
    
    uint64_t mems_found = 0;

    do {
        // std::cerr << "Starting MEM finding at position " << pos_on_r << "\n";
        mems_found += query_mem_bml(mq, min_mem_length, query_seq, ftab_k, pos_on_r);
    } while (pos_on_r < query_seq.length());

    return mems_found;
}

// We search for a MEM in P[pos_on_r..m-1] starting at pos_on_r
// Returns true if a MEM is found, false otherwise
bool MoveStructure::query_mem_bml(MoveQuery& mq, int32_t& min_mem_length, std::string& query_seq, size_t& ftab_k, int32_t& pos_on_r) {
    // Reached end of query sequence
    if (pos_on_r + min_mem_length > query_seq.length()) {
        pos_on_r = query_seq.length();
        return false;
    }

    uint64_t match_len = 0;
    int32_t init_pos = pos_on_r + min_mem_length - 1;
    MoveBiInterval bi_interval;
    bi_interval = initialize_bidirectional_search(mq, init_pos, match_len);
    if (match_len > 0) {
        --init_pos;
    }
    else {
        std::cerr << "No match with ftab\n";
    }

    // Backward extension to find if left end permits a sufficiently long MEM, should equal pos_on_r if true
    for (size_t j = 0; j <= init_pos - pos_on_r; ++j) {
        if (!extend_left(query_seq[init_pos - j], bi_interval)) {
            pos_on_r = init_pos - j + 1;
            return false;
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
    }
    mq.add_mem(pos_on_r, i, rc_interval.count(rlbwt));

    // Backward extension to find next left end of candidate MEM (inclusive)
    // TODO: if this extension is at least L, we can skip the first BML step in next recursion
    size_t end_pos_on_r = i;
    size_t j;
    if (end_pos_on_r < query_seq.length()) {
        init_pos = end_pos_on_r;
        match_len = 0;
        MoveInterval fw_interval = initialize_backward_search(mq, init_pos, match_len);
         if (match_len > 0) {
            --init_pos;
        }
        else {
            std::cerr << "No match with ftab\n";
        }
        for (i = 0; i <= init_pos - (pos_on_r + 1); ++i) {
            if (!backward_search_step(query_seq[init_pos - i], fw_interval)) {
                break;
            }
        }
    }

    pos_on_r = (init_pos - i) + 1;
    return true;
}

uint64_t MoveStructure::query_all_mems(MoveQuery& mq) {
    std::string& query_seq = mq.query();
    
    size_t ftab_k = movi_options->get_ftab_k();
    // At the time, the initialization works if ftab with ftab_k exists only
    // And the multi-ftab strategy must be turned off
    movi_options->set_multi_ftab(false);
    
    size_t s = 0;
    uint64_t init_match_len = 0; // bases matched so far
    int32_t init_pos_on_r = s;
    MoveBiInterval bi_interval = initialize_bidirectional_search(mq, init_pos_on_r, init_match_len);
    // ftab gives no match, reinitialize to do entire forward extension to find right end of MEM
    if (init_match_len == 0 && ftab_k > 1) {
        bi_interval = MoveBiInterval();
        // TODO: MEM stats
    }

    while (s < query_seq.length()) {
        // Extend MEM until mismatch
        MoveBiInterval bi_interval_before_extension;
        do {
            bi_interval_before_extension = bi_interval;
        } while (extend_right(query_seq[s + bi_interval.match_len], bi_interval));
        size_t e = s + bi_interval_before_extension.match_len;
        // std::cerr << "Add MEM of length " << e - s << " from " << s << " to " << e << " with " << bi_interval_before_extension.fw_interval.count(rlbwt) << " matches\n";
        mq.add_mem(s, e, bi_interval_before_extension.fw_interval.count(rlbwt));

        // Backward extension to find next MEM start
        init_match_len = 0;
        init_pos_on_r = e;
        bi_interval = initialize_bidirectional_search(mq, init_pos_on_r, init_match_len);
        if (init_match_len == 0 && ftab_k > 1) {
            bi_interval = MoveBiInterval();
            // TODO: MEM stats
        }
        do {
            bi_interval_before_extension = bi_interval;
        } while (extend_left(query_seq[e - bi_interval.match_len], bi_interval));
        s = e - bi_interval_before_extension.match_len + 1;
        bi_interval = bi_interval_before_extension;
    }

    return mq.get_mems().size();
}