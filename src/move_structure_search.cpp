#include "move_structure.hpp"
#include "classifier.hpp"

void MoveStructure::update_interval(MoveInterval& interval, char next_char) {
    if (!check_alphabet(next_char)) {
        throw std::runtime_error(ERROR_MSG("[update interval] This should not happen! The character should have been checked before."));
    }
#if USE_NEXT_POINTERS
    // if (movi_options->is_debug())
    //     dbg << alphabet[rlbwt[interval.run_start].get_c()] << " " << alphabet[rlbwt[interval.run_end].get_c()] << " " << next_char << "\n";
    uint64_t read_alphabet_index = alphamap[static_cast<uint64_t>(next_char)];
    if ((interval.run_start <= interval.run_end) and (alphabet[rlbwt[interval.run_start].get_c()] != next_char)) {
        if (interval.run_start == 0) {
            // To check if this case ever happens. If not, we should get rid of this condition.
            // if (movi_options->is_debug())
            //     dbg << "run_start is 0 before updating the interval!\n";
            while ((interval.run_start <= interval.run_end) and (alphabet[rlbwt[interval.run_start].get_c()] != next_char)) {
                interval.run_start += 1;
                interval.offset_start = 0;
                if (interval.run_start >= r) {
                    break;
                }
            }
        } else {
            char rlbwt_char = alphabet[rlbwt[interval.run_start].get_c()];
            uint64_t alphabet_index = alphamap_3[alphamap[rlbwt_char]][read_alphabet_index];
            if (rlbwt[interval.run_start].get_next_down(alphabet_index) == std::numeric_limits<uint16_t>::max()) {
                interval.run_start = r;
            } else {
                uint64_t run_start_ = interval.run_start + rlbwt[interval.run_start].get_next_down(alphabet_index);
                interval.run_start = run_start_;
                interval.offset_start = 0;
            }
        }
    }
    if ((interval.run_end >= interval.run_start) and (alphabet[rlbwt[interval.run_end].get_c()] != next_char)) {
        char rlbwt_char = alphabet[rlbwt[interval.run_end].get_c()];
        uint64_t alphabet_index = alphamap_3[alphamap[rlbwt_char]][read_alphabet_index];
        if (rlbwt[interval.run_end].get_next_up(alphabet_index) == std::numeric_limits<uint16_t>::max()) {
            interval.run_end = r;
        } else {
            uint64_t run_end_ = interval.run_end - rlbwt[interval.run_end].get_next_up(alphabet_index);
            interval.run_end = run_end_;
            interval.offset_end = rlbwt[interval.run_end].get_n() - 1;
        }
    }
#else
    while (interval.run_start <= interval.run_end and get_char(interval.run_start) != next_char) { //  >= or >
        interval.run_start += 1;
        interval.offset_start = 0;
        if (interval.run_start >= r) {
            break;
        }
    }
    while (interval.run_end >= interval.run_start and get_char(interval.run_end) != next_char) { //  >= or >
        interval.run_end -= 1;
        interval.offset_end = rlbwt[interval.run_end].get_n() - 1;
        if (interval.run_end == 0) {
            break;
        }
    }
#endif

}

bool MoveStructure::extend_bidirectional(char c_, MoveInterval& fw_interval, MoveInterval& rc_interval) {
    MoveInterval fw_interval_before_extension = fw_interval;
    char c_comp = complement(c_);

    bool res = backward_search_step(c_, fw_interval);
    if (res) {
        // The alphabet is already checked to be legal (ACGT)
        uint64_t skip = 0;
        uint64_t current_run = fw_interval_before_extension.run_start;
        uint64_t current_offset = fw_interval_before_extension.offset_start;
        while (current_run <= fw_interval_before_extension.run_end ) {
            if (current_run != end_bwt_idx) {
                if (complement(get_char(current_run)) < c_comp) {
                    uint64_t char_count = current_run != fw_interval_before_extension.run_end ?
                                        get_n(current_run) - current_offset : fw_interval_before_extension.offset_end - current_offset + 1;
                    skip += char_count;
                }
            } else {
                skip += 1;
            }
            current_run += 1;
            current_offset = 0;
        }

        while (skip != 0) {
            int rows_after = get_n(rc_interval.run_start) - 1 - rc_interval.offset_start;
            if (rows_after >= skip) {
                rc_interval.offset_start += skip;
                skip = 0;
            } else {
                rc_interval.run_start += 1;
                rc_interval.offset_start = 0;
                skip -= rows_after + 1;
            }
        }
        // Compute the run end for the rc interval
        skip = fw_interval.count(rlbwt) - 1;
        rc_interval.run_end = rc_interval.run_start;
        rc_interval.offset_end = rc_interval.offset_start;
        while (skip != 0) {
            int rows_after = get_n(rc_interval.run_end) - 1 - rc_interval.offset_end;
            if (rows_after >= skip) {
                rc_interval.offset_end += skip;
                skip = 0;
            } else {
                rc_interval.run_end += 1;
                rc_interval.offset_end = 0;
                skip -= rows_after + 1;
            }
        }
        return true;
    } else {
        return false;
    }
}

bool MoveStructure::extend_left(char c, MoveBiInterval& bi_interval) {
    char c_ = c;
    if (extend_bidirectional(c_, bi_interval.fw_interval, bi_interval.rc_interval)) {
        bi_interval.match_len += 1;
        return true;
    } else {
        return false;
    }
}

bool MoveStructure::extend_right(char c, MoveBiInterval& bi_interval) {
    char c_ = complement(c);
    if (extend_bidirectional(c_, bi_interval.rc_interval, bi_interval.fw_interval)) {
        bi_interval.match_len += 1;
        return true;
    } else {
        return false;
    }
}

MoveBiInterval MoveStructure::backward_search_bidirectional(std::string& R, int32_t& pos_on_r, MoveBiInterval interval, int32_t max_length) {
    // If the pattern is found, the pos_on_r will be equal to 0 and the interval will be non-empty
    // Otherwise the interval corresponding to match from the end until and including the updated pos_on_r will be returned
    // The input interval is non-empty and corresponds to the interval that matches the read (R) at pos_on_r
    MoveBiInterval prev_interval = interval;
    int32_t pos_on_r_saved = pos_on_r;
    while (pos_on_r > 0 and !interval.fw_interval.is_empty()) {
        /*if (!check_alphabet(R[pos_on_r - 1])) {
            return interval;
        }*/
        prev_interval = interval;
        bool res = extend_left(R[pos_on_r - 1], interval);
        if (!interval.fw_interval.is_empty()) {
            pos_on_r -= 1;
        }

        // The following is only for the backward_search is called for "look ahead" for kmer skipping
        if (pos_on_r_saved - pos_on_r > max_length)
            break;
    }
    if (interval.fw_interval.is_empty()) {
        return prev_interval;
    } else {
        return interval;
    }
}

MoveInterval MoveStructure::backward_search(std::string& R, int32_t& pos_on_r, MoveInterval interval, int32_t max_length) {
    // If the pattern is found, the pos_on_r will be equal to 0 and the interval will be non-empty
    // Otherwise the interval corresponding to match from the end until and including the updated pos_on_r will be returned
    // The input interval is non-empty and corresponds to the interval that matches the read (R) at pos_on_r
    MoveInterval prev_interval = interval;
    // uint64_t match_len = 1; --- unused variable
    int32_t pos_on_r_saved = pos_on_r;
    while (pos_on_r > 0 and !interval.is_empty()) {
        /*if (!check_alphabet(R[pos_on_r - 1])) {
            return interval;
        }*/
        prev_interval = interval;
        backward_search_step(R, pos_on_r, interval);
        if (!interval.is_empty()) {
            pos_on_r -= 1;
        }
        /*update_interval(interval, R[pos_on_r - 1]);
        if (!interval.is_empty()) {
            LF_move(interval.offset_start, interval.run_start);
            LF_move(interval.offset_end, interval.run_end);
            pos_on_r -= 1;
        }*/

        // The following is only for the backward_search is called for "look ahead" for kmer skipping
        if (pos_on_r_saved - pos_on_r > max_length)
            break;
    }
    if (interval.is_empty()) {
        return prev_interval;
    } else {
        return interval;
    }
}

MoveInterval MoveStructure::try_ftab(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len, size_t ftab_k, bool rc) {
    auto& query_seq = mq.query();
    if (ftab_k > 1 and pos_on_r >= ftab_k - 1) {
        // uint64_t kmer_code = kmer_to_number(ftab_k, query_seq.substr(pos_on_r - ftab_k + 1, ftab_k), alphamap);
        uint64_t kmer_code = kmer_to_number(ftab_k, query_seq, pos_on_r - ftab_k + 1, alphamap, rc);
        if (kmer_code != std::numeric_limits<uint64_t>::max()) {
            auto& current_ftab = movi_options->is_multi_ftab() ? ftabs[ftab_k - 1] : ftab;
            if (!current_ftab[kmer_code].is_empty()) {
                // Add the skipped matching length, e.g., for zml computation
                if (movi_options->is_zml()) {
                    for (size_t i = 0; i < ftab_k - 1; i++) {
                        mq.add_ml(i, movi_options->is_stdout());
                    }
                }
                if (movi_options->is_zml() or movi_options->is_kmer()) {
                    match_len = ftab_k - 1;
                }
                pos_on_r = pos_on_r - ftab_k + 1;
                return current_ftab[kmer_code];
            }
        }
    }
    // If we reach here, we know the ftab could not be used, so we return an empty interval
    MoveInterval empty_interval;
    empty_interval.make_empty();
    return empty_interval;
}

MoveBiInterval MoveStructure::initialize_bidirectional_search(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len) {
    MoveBiInterval bi_interval;
    int32_t pos_on_r_before = pos_on_r;
    bi_interval.fw_interval = initialize_backward_search(mq, pos_on_r, match_len);
    // If the initialization was unsuccessfull, match_len of the bidirectional interval should be set to 0
    // Without the following condition, the match_len will be increamented in the following line
    if (match_len == 0) {
        bi_interval.match_len = match_len;
        return bi_interval;
    }
    // This is needed because of the way match_len is off by one
    match_len += 1;
    bi_interval.match_len = match_len;

    int32_t pos_on_r_rc = pos_on_r_before;
    uint64_t match_len_rc = 0;
    bi_interval.rc_interval = initialize_backward_search(mq, pos_on_r_rc, match_len_rc, true);

    // match_len_rc will be equal to match_len if both fw and rc are present int he reference
    if (match_len - 1 != match_len_rc) {
        throw std::runtime_error(ERROR_MSG("[initialize bidirectional search] The reverse complement might not be present in the reference."));
    }

    return bi_interval;
}

MoveInterval MoveStructure::initialize_backward_search(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len, bool rc) {
    all_initializations += 1;
    // Initialize assuming that the character at pos_on_r exists in the alphabet
    size_t ftab_k = movi_options->get_ftab_k();
    if (movi_options->is_multi_ftab()) {
        while (ftab_k > 1 and pos_on_r >= ftab_k - 1) {
            //std::cerr << ftab_k << "\n";
            MoveInterval ftab_res = try_ftab(mq, pos_on_r, match_len, ftab_k, rc);
            if (!ftab_res.is_empty())
                return ftab_res;
            ftab_k -= 2;
        }
    } else if (ftab_k > 1) {
        MoveInterval ftab_res = try_ftab(mq, pos_on_r, match_len, ftab_k, rc);
        if (!ftab_res.is_empty())
            return ftab_res;
    }
    // If we reach here, we know the ftab could not be used, so we initialize regularly
    if (movi_options->is_multi_ftab() and ftab_k < movi_options->get_ftab_k()) {
        ftab_k += 2;
    }
    if (pos_on_r >= ftab_k - 1) {
        no_ftab += 1;
    }
    auto& query_seq = mq.query();
    auto first_char_index = alphamap[rc ? complement(query_seq[pos_on_r]): query_seq[pos_on_r]] + 1;
    MoveInterval interval(
        first_runs[first_char_index],
        first_offsets[first_char_index],
        last_runs[first_char_index],
        last_offsets[first_char_index]
    );
    return interval;
}

bool MoveStructure::backward_search_step(char c, MoveInterval& interval) {
    if (!check_alphabet(c)) {
        interval.make_empty();
        return false;
    }

    update_interval(interval, c);
    if (!interval.is_empty()) {
        LF_move(interval.offset_start, interval.run_start);
        LF_move(interval.offset_end, interval.run_end);
        return true;
    } else {
        return false;
    }
}

uint64_t MoveStructure::backward_search_step(std::string& R, int32_t& pos_on_r, MoveInterval& interval) {
    // It is assumed that the interval represents a match until and including the position pos_on_r
    // Then we try to see if the match can be extended to position pos_on_r - 1
    // The interval becomes empty if the match cannot be extended
    // otherwise it is updated according to the character on pos_on_r - 1
    uint64_t ff_count = 0;
    if (pos_on_r <= 0) {
        throw std::runtime_error(ERROR_MSG("[backward search step] The backward search step never be called on position 0 of the read."));
    }

    if (!check_alphabet(R[pos_on_r - 1])) {
        interval.make_empty();
        return ff_count;
    }

    update_interval(interval, R[pos_on_r - 1]);
    if (!interval.is_empty()) {
        ff_count += LF_move(interval.offset_start, interval.run_start);
        ff_count += LF_move(interval.offset_end, interval.run_end);
    }

    return ff_count;
}

uint64_t MoveStructure::query_backward_search(MoveQuery& mq, int32_t& pos_on_r) {
    auto& query_seq = mq.query();
    // Check the special case of non-existing character at the end of the read
    // before initializing the interval based on that character
    if (!check_alphabet(query_seq[pos_on_r])) {
        pos_on_r += 1; // even the character at the end_pos was not found.
        return 0;
    }
    // Initial the interval by matching the character at the end of the read (pos_on_r)
    uint64_t not_used = 0;
    MoveInterval initial_interval = initialize_backward_search(mq, pos_on_r, not_used);
    return backward_search(query_seq, pos_on_r, initial_interval, std::numeric_limits<int32_t>::max()).count(rlbwt);
}

bool MoveStructure::look_ahead_ftab(MoveQuery& mq, uint32_t pos_on_r, int32_t& step) {
    size_t ftab_k = movi_options->get_ftab_k();
    size_t k = movi_options->get_k();
    auto& query_seq = mq.query();
    // int32_t pos_on_r_ahead = pos_on_r - static_cast<int32_t>(k/2);
    for (step = 0; step <= 19 ; step += 1) {
        int32_t pos_on_r_ahead = pos_on_r - k + ftab_k + step;
        uint64_t kmer_code = kmer_to_number(ftab_k, query_seq, pos_on_r_ahead - ftab_k, alphamap);
        if (kmer_code != std::numeric_limits<uint64_t>::max() and !ftab[kmer_code].is_empty()) {
            // return true;
        } else {
            return false;
        }
    }
    return true;
}

bool MoveStructure::look_ahead_backward_search(MoveQuery& mq, uint32_t pos_on_r, int32_t step) {
    size_t ftab_k = movi_options->get_ftab_k();
    size_t k = movi_options->get_k();
    auto& query_seq = mq.query();

    uint64_t match_len = 0;
    int32_t pos_on_r_ahead = pos_on_r - step;
    MoveInterval initial_interval = initialize_backward_search(mq, pos_on_r_ahead, match_len);
    backward_search(query_seq, pos_on_r_ahead, initial_interval, k - step - match_len);
    if (pos_on_r - pos_on_r_ahead >= k - 1) {
        return true;
    } else {
        return false;
    }
}