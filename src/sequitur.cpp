#include <sys/stat.h> 

#include "move_structure.hpp"
#include "utils.hpp"

std::ostringstream dbg;

// pos_on_r points to the furthest base of the kmer to be searched
// The search is from the right end of the kmer
// All kmers on the right side of the kmer_middle might be found in this call if they exist
uint64_t MoveStructure::query_kmers_from_bidirectional(MoveQuery& mq, int32_t& pos_on_r) {
    uint64_t kmers_found = 0;
    int32_t last_kmer_looked_at = pos_on_r;

    size_t ftab_k = movi_options->get_ftab_k();
    size_t k = movi_options->get_k();
    auto& query_seq = mq.query();

    // The middle point of the kmer
    int32_t kmer_middle = pos_on_r - k/2;
    // To remember the kmer which initiated the search
    int32_t pos_on_r_saved = pos_on_r;
    // The number of bases matched so far
    uint64_t match_len = 0;
    // The position of the left side of the kmer on the read
    int32_t kmer_left = pos_on_r - k + 1;
    // The position on the right side of the match to be found by ftab
    int32_t ftab_right = kmer_left + ftab_k - 1;

    // At the time, the initialization works if ftab with ftab_k exists only
    // And the multi-ftab strategy must be turned off
    movi_options->set_multi_ftab(false);

    MoveBiInterval bi_interval;
    int32_t ftab_right_initialize = ftab_right;

    bi_interval = initialize_bidirectional_search(mq, ftab_right_initialize, match_len);

    if (match_len == 0 and ftab_k > 1) {
        // If the ftab-k-mer at the end of the read is not found,
        // we can skip all the kmers until ftab_right - 1
        kmer_stats.initialize_skipped += 1;
        pos_on_r = ftab_right - 1;
        return 0;
    }

    // pos_on_r and pos_on_r_saved point to the beginning of the kmer for which the bidirectional initialization was performed above
    // Extend to right until the kmer is found and save the observed intervals beyond k/2
    std::vector<MoveBiInterval> partial_matches;
    partial_matches.resize(k);
    // kmer_right is the last position matched so far
    uint64_t kmer_right = ftab_right;
    int here = 0;
    while (kmer_right < pos_on_r_saved) {
        // The next k-mer we are looking at, ends at next_pos
        int next_pos = kmer_right + 1;

        if (movi_options->is_debug() and next_pos >= k)
            dbg << "\nkmer at " << next_pos << ": " << query_seq.substr(next_pos - k + 1, k) << " ";

        bool extend_right_res = extend_right(query_seq[next_pos], bi_interval);
        if (!extend_right_res) {
            here = 1;
            // The kmer was not found, we can skip kmers until ftab_right
            pos_on_r = kmer_right;
            last_kmer_looked_at = next_pos;

            // Printing all the kmers being skipped, because the extension to the right was not possible
            if (movi_options->is_debug()) {
                int j = next_pos + 1;
                while (j < pos_on_r_saved) {
                    if (j > k)
                        dbg << "\nkmer at " << j << ": " << query_seq.substr(j - k + 1, k) << "-";
                    j += 1;
                }
            }

            // It's important to break here, to avoid false extensions in the following iterations
            break;
        } else {
            match_len += 1;
            kmer_right = next_pos;
            bi_interval.match_len = match_len;
            // store the intervals for matches beyond half point of the kmer
            if (kmer_right > kmer_middle and kmer_right != pos_on_r) {
                bi_interval.match_len = match_len;
                // "kmer_right - kmer_left" is the index of the kmer from the left end
                partial_matches[kmer_right - kmer_left] = bi_interval;
            }
        }
    }

    if (kmer_right == pos_on_r_saved) {
        // The kmer at pos_on_r was found by k bidirectional backward search
        kmers_found += 1;
        kmer_stats.backward_search_empty += 1;
        // std::cerr << "The first kmer at " << pos_on_r << " was found.\n";
        pos_on_r = pos_on_r_saved - 1;
        // To avoid searching this kmer again in the next step
        kmer_right -= 1;

        if (movi_options->is_debug()) {
            if (pos_on_r_saved > k)
                dbg << "\nkmer at " << pos_on_r_saved << ": " << query_seq.substr(pos_on_r_saved - k + 1, k) << "\n";
            dbg << "1";
        }

    } else {
        if (movi_options->is_debug()) {
            if (pos_on_r_saved > k)
                dbg << "\nkmer at " << pos_on_r_saved << ": " << query_seq.substr(pos_on_r_saved - k + 1, k) << "-\n";
            dbg << "0";
        }
    }


    // Use partial matches for finding other overlapping kmers (until half point)
    if (kmer_right > kmer_middle) {
        if (movi_options->is_debug()) {
            for (int i = kmer_right + 1; i < pos_on_r_saved; i++) {
                dbg << "0";
            }
        }
        for (uint i = kmer_right; i > kmer_middle; i--) {
            // Set kmer_left_ext to be the last match position on the left end
            int32_t kmer_left_ext = kmer_left;
            auto& partial_match_interval = partial_matches[i - kmer_left];
            last_kmer_looked_at = i;
            while (partial_match_interval.match_len < k and kmer_left_ext > 0) {
                // We don't need to do bidirectional left extension here, simple backward search is enough
                bool res = backward_search_step(query_seq[kmer_left_ext - 1], partial_match_interval.fw_interval);
                // bool res = extend_left(query_seq[kmer_left_ext - 1], partial_match_interval);
                if (!res) {
                    // The current kmer is not present, move to the next partial match by breaking from the inner loop
                    break;
                } else {
                    kmer_left_ext -= 1;
                    partial_match_interval.match_len += 1;
                }
            }
            if (partial_match_interval.match_len >= k) {
                // The kmer was found by extending the partial match to left
                kmers_found += 1;
                kmer_stats.positive_skipped += 1;
                if (movi_options->is_debug()) {
                    std::cerr << "kmer at " << kmer_left_ext + k - 1 << " was found.\n";
                    dbg << "1";
                }
            } else {
                if (movi_options->is_debug()) {
                    dbg << "0";
                }
            }

            pos_on_r -= 1;
        }
        // At this point we have checked the presence of all the kmer beyond kmer_middle
    } else {
        // If we got here, we had to start the first while by breaking because of an unsuccessfull attempt to extend to the right
        if (here != 1)
            std::cerr << here << "\t" << last_kmer_looked_at << "\t" << pos_on_r << "\n";
        if (kmer_right != pos_on_r)
            std::cerr << "This should not happen: " << kmer_right << "\t" << pos_on_r << "\n";
        // pos_on_r should have already been assigned to be the last kmer_right
        // So we should never get here to do the assignment in practice
        pos_on_r = kmer_right;

        if (movi_options->is_debug()) {
            for (int i = kmer_middle + 1; i <= pos_on_r_saved - 1; i++) {
                dbg << "0";
            }
        }

    }
    return kmers_found;
}

uint64_t MoveStructure::query_kmers_from(MoveQuery& mq, int32_t& pos_on_r, bool single) {
    size_t ftab_k = movi_options->get_ftab_k();
    size_t k = movi_options->get_k();
    auto& query_seq = mq.query();
    int32_t pos_on_r_saved = pos_on_r;

    // An alternative strategy to look ahead for possible skipping
    // int32_t step = 0;
    // if (ftab_k > 1 and !look_ahead_ftab(mq, pos_on_r, step)) {
    //     kmer_stats.look_ahead_skipped += k - ftab_k - step;
    //     pos_on_r = pos_on_r - k + ftab_k + step - 1;
    //     pos_on_r_saved = pos_on_r;
    // }

    uint64_t match_len = 0;
    MoveInterval initial_interval;
    do {
        initial_interval = initialize_backward_search(mq, pos_on_r, match_len);
        if (match_len == 0 and ftab_k > 1) {
            kmer_stats.initialize_skipped += 1;
            pos_on_r -= 1;
            pos_on_r_saved = pos_on_r;
        }
    } while (match_len == 0 and pos_on_r >= k - 1 and ftab_k > 1);

    // I want to check how much slower it gets if turn off the positive skip:
    auto backward_search_result = backward_search(query_seq, pos_on_r, initial_interval, single ? k - match_len - 2 : std::numeric_limits<int32_t>::max());

    if (backward_search_result.is_empty()) {
        // We get here when there is an illegal character at pos_on_r, just skip the current position
        pos_on_r = pos_on_r_saved - 1;
        kmer_stats.backward_search_empty += 1;
        return 0;
    } else {
        if (pos_on_r_saved - pos_on_r >= k - 1) {
            // At leat one kmer was found, update the postion and return the count
            uint64_t kmers_found = pos_on_r_saved - pos_on_r - k + 2;
            kmer_stats.positive_skipped += kmers_found - 1;

            if (movi_options->is_debug()) {
                int32_t pos_on_r_ = pos_on_r_saved;
                auto backward_search_result_one_extra_base = backward_search(query_seq, pos_on_r_, initial_interval, k - match_len - 1);
                dbg << backward_search_result.count(rlbwt) << "----" <<  backward_search_result_one_extra_base.count(rlbwt) << "\n";
                dbg << backward_search_result << "\n";
                if (backward_search_result.count(rlbwt) == backward_search_result_one_extra_base.count(rlbwt)) {
                    // TODO: Use a counter to count the number of such incidents
                } else {
                }
            }

            pos_on_r = pos_on_r + k - 2;
	        return kmers_found;
        } else {
            // No kmer was found, update the postion
            kmer_stats.backward_search_failed += 1;
            pos_on_r = pos_on_r_saved - 1;
            return 0;
        }
    }
}

void MoveStructure::query_all_kmers(MoveQuery& mq, bool kmer_counts) {
    size_t ftab_k = movi_options->get_ftab_k();
    size_t k = movi_options->get_k();
    auto& query_seq = mq.query();
    int32_t pos_on_r = query_seq.length() - 1;

    // To handle a special case for k equal to 1
    if (k == 1) {
        uint64_t kmers_found = 0;
        while (pos_on_r >= 0) {
            kmers_found += check_alphabet(query_seq[pos_on_r]) ? 1 : 0;
            pos_on_r -= 1;
        }
        kmer_stats.positive_kmers += kmers_found;
        return;
    }

    while (!check_alphabet(query_seq[pos_on_r])) {
        pos_on_r -= 1; // Find the first position where the character is legal
    }


    int32_t step = k/3;
    // k - step has to be always greater than ftab-k
    if (k - step < ftab_k) {
        step = k - ftab_k - 1;
    }

    while (pos_on_r >= k - 1) {
        if (pos_on_r >= k -1 + step and !look_ahead_backward_search(mq, pos_on_r, step)) {
            kmer_stats.look_ahead_skipped += step + 1;
            pos_on_r = pos_on_r - step - 1;
        } else {
            if (kmer_counts) {
                if (pos_on_r == k - 1) {
                    kmer_stats.positive_kmers += query_kmers_from(mq, pos_on_r);
                } else {
                    kmer_stats.positive_kmers += query_kmers_from_bidirectional(mq, pos_on_r);

                    if (movi_options->is_debug()) {
                        dbg.str("");
                        dbg.clear();
                        int pos_on_r_before = pos_on_r;
                        int found_regular = 0;
                        dbg << "regular backward search:\n";
                        int k_m = pos_on_r - k/2;
                        dbg << pos_on_r << "\t" << k_m << "\n";
                        dbg << k/2 << "\n";
                        for (int j = pos_on_r; j > k_m; j--) {
                            dbg << " " << j << " ";
                            if (j >= k)
                                dbg << "\nkmer at " << j << ": " << query_seq.substr(j - k + 1, k) << " ";

                            int pos = j;
                            dbg << " pos:" << pos << " ";
                            int z = query_kmers_from(mq, pos, true);
                            if (z == 1 and pos == j - 1) {
                                dbg << "1";
                                found_regular += 1;
                            } else
                                dbg << "0";
                        }
                        dbg << "\n";
                        dbg << "bidirectional search:\n";
                        int pos = pos_on_r_before;
                        int found_bidirectional = query_kmers_from_bidirectional(mq, pos);
                        dbg << "\n";
                        kmer_stats.positive_kmers += found_bidirectional;
                        if (found_regular != found_bidirectional) {
                            std::cerr << "pos_on_r_before:" << pos_on_r_before << " pos_on_r:" << pos_on_r
                                    << " found_regular:" << found_regular << " found_bidirectional" << found_bidirectional << "\n";
                            std::cerr << dbg.str() << std::endl;
                        }
                    }
                }
            } else {
                kmer_stats.positive_kmers += query_kmers_from(mq, pos_on_r);
            }
        }

        while (!check_alphabet(query_seq[pos_on_r])) {
            pos_on_r -= 1; // Find the first position where the character is legal
        }
    }
}