#include "move_structure.hpp"

void MoveStructure::reconstruct_lf() {
    uint64_t offset = 0;
    uint64_t run_index = 0;
    uint64_t i = 0;
    uint64_t ff_count_tot = 0;

    uint64_t total_elapsed = 0;
    for (; run_index != end_bwt_idx; ) {
        auto begin = std::chrono::system_clock::now();
        ff_count_tot += LF_move(offset, run_index);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        total_elapsed += elapsed.count();

        if (i % 1000000 == 0)
            std::cerr << i << "\r";
        i += 1;
        // orig_string = rlbwt[run_index].get_c() + orig_string;
    }
    std::printf("Time measured for reconstructing the original text: %.3f seconds.\n", total_elapsed * 1e-9);

    std::cerr << "Finished reconstructing the original string.\n";
    std::cerr << "Total fast forward: " << ff_count_tot << "\n";
}

void MoveStructure::sequential_lf() {
    uint64_t line_index = 0;
    uint64_t row_index = 0;
    uint64_t ff_count_tot = 0;

    uint64_t total_elapsed = 0;
    for (uint64_t row_index = 0; row_index < r; row_index++) {
        auto& current = rlbwt[row_index];
        for (uint64_t j = 0; j < current.get_n(); j ++) {
            if (line_index % 1000000 == 0)
                std::cerr << line_index << "\r";
            uint64_t offset = j;
            uint64_t i = row_index;
            auto begin = std::chrono::system_clock::now();
            ff_count_tot += LF_move(offset, i);
            auto end = std::chrono::system_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            total_elapsed += elapsed.count();
            line_index += 1;
        }

    }
    std::printf("Time measured for LF-mapping of all the characters: %.3f seconds.\n", total_elapsed * 1e-9);

    std::cerr << line_index << "\n";
    std::cerr << "Finished performing LF-mapping for all the BWT characters.\n";
    std::cerr << "Total fast forward: " << ff_count_tot << "\n";
}

void MoveStructure::random_lf() {
    uint64_t ff_count_tot = 0;

    // generate the random order from 1 to length
    std::vector<uint64_t> random_order(length);
    for (uint64_t i = 0; i < length; i++) {
        random_order[i] = i;
    }
    std::srand(time(0));
    std::random_shuffle(random_order.begin(), random_order.end());

    // find the n and id for each random BWT row
    std::vector<uint64_t> n_to_id(length);
    std::vector<uint64_t> id_to_p(r);
    uint64_t current_n = 0;
    for (uint64_t i = 0; i < r; i++) {
        id_to_p[i] = current_n;
        for (uint64_t j = 0; j < get_n(i); j++) {
            n_to_id[current_n] = i;
            current_n += 1;
        }
    }


    uint64_t total_elapsed = 0;
    for (uint64_t i = 0; i < length; i++) {
        if (i % 1000000 == 0)
            std::cerr << i << "\r";

        // uint64_t n = std::rand() % r;
        uint64_t n = random_order[i];
        uint64_t id = n_to_id[n];
        uint64_t offset = n - id_to_p[id];
        auto begin = std::chrono::system_clock::now();
        ff_count_tot += LF_move(offset, id);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        total_elapsed += elapsed.count();
    }
    std::printf("Time measured for LF-mapping of all the characters in the random order: %.3f seconds.\n", total_elapsed * 1e-9);

    std::cerr << length << "\n";
    std::cerr << "Finished performing LF query for all the BWT characters.\n";
    std::cerr << "Total fast forward: " << ff_count_tot << "\n";
}

void MoveStructure::verify_lfs() {
    // This function only works if the LF function works correctly and
    // that depends on populating the occs bitvector correctly.
    // This function is only used for debugging purposes.
    uint64_t not_matched = 0;
    for (uint64_t i = 0; i < all_p.size(); i++) {
        std::uint64_t end_ = (i < all_p.size() - 1) ? all_p[i + 1] : length;
        for (uint64_t j = all_p[i]; j < end_; j++) {
            uint64_t offset_ = j - all_p[i];
            uint64_t idx_ = i;
            uint64_t lf = 0;
            if (i != end_bwt_idx) {
                uint64_t alphabet_index = rlbwt[idx_].get_c();
                lf = LF(j, alphabet_index);
            } else {
                std::cerr << "end_run = " << i << " len: " << rlbwt[i].get_n () << "\n";
            }
            LF_move(offset_, idx_);
            uint64_t lf_move = all_p[idx_] + offset_;
            if (lf != lf_move) {
                not_matched += 1;
                std::cerr << "j\t" << j << "\n";
                std::cerr << "idx\t" << i << "\n";
                std::cerr << "offset\t" << j - all_p[i] << "\n";
                std::cerr << "rlbwt[idx].get_id\t" << get_id(i) << "\n";
                std::cerr << "get_offset(i)\t" << get_offset(i) << "\n";
                for (uint64_t k = 0; k <= i; k++) {
                    std::cerr << rlbwt[k].get_n() << " ";
                }
                std::cerr << "\n\n";

                std::cerr << "lf\t" << lf << "\n";
                std::cerr << "lf_move\t" << lf_move << "\n";
                std::cerr << "idx_\t" << idx_ << "\n";
                std::cerr << "offset_\t" << offset_ << "\n";
                std::cerr << "all_p[idx_]\t" << all_p[idx_] << "\n";
                std::cerr << "\n\n\n";
            }
        }
    }
    if (not_matched == 0) {
        std::cerr << "All the LF_move operations are correct.\n";
    } else {
        original_r = 0;
        std::cerr << "There are " << not_matched << " LF_move operations that failed to match the true lf results.\n";
    }
}

void MoveStructure::verify_lf_loop() {
    // This function only works if the LF function works correctly and
    // It verifies that n LF_move operations loops back to the same BWT offeset.
    uint64_t idx = end_bwt_idx;
    uint64_t offset = 0;

    std::vector<std::vector<uint64_t>> bwt_offsets(r);
    for (uint64_t i = 0; i < r; i++) {
        bwt_offsets[i].resize(rlbwt[i].get_n(), 0);
    }

    std::cerr << "Verifying that all the LF_move operations loop back to the same BWT offset.\n";
    for (uint64_t i = 0; i < length; i++) {
        LF_move(offset, idx);
        bwt_offsets[idx][offset] = 1;
    }

    uint64_t visited_offsets = 0;
    for (uint64_t i = 0; i < r; i++) {
        for (uint64_t j = 0; j < rlbwt[i].get_n(); j++) {
            if (bwt_offsets[i][j] == 1) {
                visited_offsets += 1;
            }
        }
    }

    if (idx == end_bwt_idx and offset == 0 and visited_offsets == length) {
        std::cerr << "All the LF_move operations are correct.\n";
    } else {
        std::cerr << "\tlast idx: " << idx << " last offset: " << offset << "\n";
        std::cerr << "\tNumber of visited offsets: " << visited_offsets << "\n";
        std::cerr << "\tlength of the BWT: " << length << "\n";
        throw std::runtime_error(ERROR_MSG("[verify lf_loop] LF_move operations failed to loop back to the same BWT offset."));
    }
}

uint64_t MoveStructure::reposition_up(uint64_t idx, char c, uint64_t& scan_count) {
    if (idx == 0)
        return r;
    char row_c = alphabet[rlbwt[idx].get_c()];

    while (idx > 0 and row_c != c) {
        scan_count += 1;
        idx -= 1;
        row_c = alphabet[rlbwt[idx].get_c()];
    }
    /* if (logs) {
        if (repositions.find(scan_count) != repositions.end())
            repositions[scan_count] += 1;
        else
            repositions[scan_count] = 1;
    } */
    /* if (movi_options->is_verbose())
        std::cerr << "\t \t \t \t idx after the while in the reposition up" << idx << "\n";*/
    return (row_c == c) ? idx : r;
}

uint64_t MoveStructure::reposition_down(uint64_t idx, char c, uint64_t& scan_count) {
    if (idx == r - 1)
        return r;
    char row_c = alphabet[rlbwt[idx].get_c()];

    while (idx < r - 1 && row_c != c) {
        scan_count += 1;
        idx += 1;
        row_c = alphabet[rlbwt[idx].get_c()];
    }
    /* if (logs) {
        if (repositions.find(scan_count) != repositions.end())
            repositions[scan_count] += 1;
        else
            repositions[scan_count] = 1;
    } */
    /*if (movi_options->is_verbose())
        std::cerr << "\t \t \t \t idx after the while in the reposition down: " << idx << " " << c << " " << row_c << "\n";*/
    return (row_c == c) ? idx : r;
}

uint64_t MoveStructure::query_pml(MoveQuery& mq) {
    auto& R = mq.query();
    int32_t pos_on_r = R.length() - 1;
    uint64_t idx = r - 1; // or we can start from a random position in the rlbwt std::rand() % r
    uint64_t offset = get_n(idx) - 1;

    uint64_t match_len = 0;
    uint16_t ff_count = 0;
    uint64_t ff_count_tot = 0;
    uint64_t scan_count = 0;
    auto t1 = std::chrono::high_resolution_clock::now();

    if (movi_options->is_verbose()) {
        std::cerr << "beginning of the search \ton query: " << mq.query() << "\t";
        std::cerr << "and on BWT, idx(r-1): " << idx << " offset: " << offset << "\n";
    }

    // Multi-class classification
    if (movi_options->is_multi_classify()) {
        for (uint16_t i = 0; i < num_species; i++) {
            if (!movi_options->is_pvalue_scoring()) {
                classify_cnts[i] = 0;
            } else {
                doc_scores[i] = 0;
            }
        }
    }

    uint16_t best_doc = std::numeric_limits<uint16_t>::max(); // for multi-class classification
    uint16_t second_best_doc = std::numeric_limits<uint16_t>::max();
    uint64_t iteration_count = 0;
    uint32_t sum_matching_lengths = 0;
    while (pos_on_r > -1) {
        iteration_count += 1;
        if (movi_options->is_logs() and (iteration_count-1)%200 == 0) {
            t1 = std::chrono::high_resolution_clock::now();
        }

        if (movi_options->is_verbose())
            std::cerr << "Searching position " << pos_on_r << " of the read:\n";

        auto& row = rlbwt[idx];
        uint64_t row_idx = idx;
        char row_c = alphabet[row.get_c()];
        if (!check_alphabet(R[pos_on_r])) {
            // The character from the read does not exist in the reference
            match_len = 0;
            scan_count = 0;

            if (movi_options->is_verbose())
                std::cerr << "\t The character " << R[pos_on_r] << " does not exist.\n";
        } else if (row_c == R[pos_on_r]) {
            // Case 1
            match_len += 1;
            scan_count = 0;

            if (movi_options->is_verbose()) {
                std::cerr << "\t Cas1: It was a match. \n" << "\t Continue the search...\n";
                std::cerr << "\t match_len: " << match_len << "\n";
                std::cerr << "\t current_id: " << idx << "\t row.id: " << get_id(row_idx) << "\n"
                          << "\t row.get_n: " << get_n(row_idx) << " rlbwt[idx].get_n: " << get_n(get_id(row_idx)) << "\n"
                          << "\t offset: " << offset << "\t row.get_offset(): " << get_offset(row_idx) << "\n";
            }
        } else {
            // Case 2
            // Repositioning up or down (randomly or with thresholds)
            if (movi_options->is_verbose())
                std::cerr << "\t Case2: Not a match, looking for a match either up or down...\n";

            uint64_t idx_before_reposition = idx;
#if USE_THRESHOLDS
            bool up = movi_options->is_random_repositioning() ?
                               reposition_randomly(idx, offset, R[pos_on_r], scan_count) :
                               reposition_thresholds(idx, offset, R[pos_on_r], scan_count);
#else
            // When there is no threshold, reposition randomly
            bool up = reposition_randomly(idx, offset, R[pos_on_r], scan_count);
#endif

            match_len = 0;
            // scan_count = (!constant) ? std::abs((int)idx - (int)idx_before_reposition) : 0;

            char c = alphabet[rlbwt[idx].get_c()];

            if (movi_options->is_verbose())
                std::cerr << "\t up: " << up << " idx: " << idx << " c:" << c << "\n";

            // sanity check
            if (c == R[pos_on_r]) {
                // Observing a match after the reposition
                // The right match_len should be:
                // min(new_lcp, match_len + 1)
                // But we cannot compute lcp here
                offset = up ? get_n(idx) - 1 : 0;
                /* if (movi_options->is_verbose())
                    std::cerr << "\t idx: " << idx << " offset: " << offset << "\n"; */
            } else {
                std::cerr << "\t \t pos: " << pos_on_r << " r[pos]:" <<  R[pos_on_r] << " t[pointer]:" << c << "\n";
                std::cerr << "\t \t " << up << ", " << R[pos_on_r] << ", " << pos_on_r << "\n";
                std::cerr << "\t \t ";
                for (int k = 10; k > 0; --k)
                    std::cerr << alphabet[rlbwt[idx - k].get_c()] << "-";
                for (int k = 0; k < 10; k++)
                    std::cerr << alphabet[rlbwt[idx + k].get_c()] << "-";
                std::cerr << "\n";
                auto saved_idx = idx;

                movi_options->set_verbose(true);
#if USE_THRESHOLDS
                reposition_thresholds(saved_idx, offset, R[pos_on_r], scan_count);
#endif
                throw std::runtime_error(ERROR_MSG("[query pml] This should not happen!"));
            }
        }
    
        sum_matching_lengths += match_len;
        mq.add_ml(match_len, movi_options->is_stdout());
        if (movi_options->is_get_sa_entries()) {
            uint64_t sa_entry = get_SA_entries(idx, offset);
            mq.add_sa_entries(sa_entry);
        }
        pos_on_r -= 1;

        // LF step
        ff_count = LF_move(offset, idx);
        ff_count_tot += ff_count;
        if (movi_options->is_logs()) {
            if (iteration_count % 200 == 0) {
                auto t2 = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
                mq.add_cost(elapsed);
            }
            mq.add_fastforward(ff_count);
            mq.add_scan(scan_count);
        }

        if (movi_options->is_multi_classify()) {
            if (match_len >= movi_options->get_min_match_len()) {
                /*uint64_t full_ind = run_offsets[idx] + offset;
                uint16_t cur_doc = doc_pats[full_ind];
                classify_cnts[cur_doc]++;
                if (classify_cnts[cur_doc] >= classify_cnts[best_doc]) {
                    best_doc = cur_doc;
                }*/ 

                uint64_t color_id;
#if COLOR_MODE == 1
                color_id = static_cast<uint64_t>(rlbwt[idx].color_id);
                // Skip doc sets that weren't saved (thrown away by compression).
                if (color_id >= unique_doc_sets.size()) continue;
#else
                if (movi_options->is_doc_sets_vector_of_vectors()) {
                    color_id = static_cast<uint64_t>(doc_set_inds[idx]);
                    // Skip doc sets that weren't saved (thrown away by compression).
                    if (color_id >= unique_doc_sets.size()) continue;
                } else {
                    color_id = doc_set_flat_inds[idx].get();
                    // Skip doc sets that weren't saved (thrown away by compression).
                    if (color_id >= flat_colors.size()) continue;
                }
#endif
                std::span<uint16_t> cur_set;
                if (movi_options->is_doc_sets_vector_of_vectors()) {
                    std::vector<uint16_t> &cur_set_vec = unique_doc_sets[color_id];
                    cur_set = std::span<uint16_t>(cur_set_vec.data(), cur_set_vec.size());
                } else {
                    uint32_t cur_set_size = flat_colors[color_id];
                    cur_set = std::span<uint16_t>(flat_colors.data() + color_id + 1,
                                                  flat_colors.data() + color_id + 1 + cur_set_size);
                }

                for (int doc : cur_set) {
                    if (!movi_options->is_pvalue_scoring()) {
                        classify_cnts[doc]++;
                        if (doc != best_doc) {
                            if (best_doc == std::numeric_limits<uint16_t>::max() || classify_cnts[doc] > classify_cnts[best_doc]) {
                                second_best_doc = best_doc;
                                best_doc = doc;
                            } else if (second_best_doc == std::numeric_limits<uint16_t>::max() || classify_cnts[doc] > classify_cnts[second_best_doc]) {
                                second_best_doc = doc;
                            }
                        }
                    } else {
                        // p value strategy
                        double val = match_len - (log_lens[doc] / log4);
                        if (val >= 0) {
                            doc_scores[doc] += std::min(val, 1.);
                            if (doc != best_doc) {
                                if (best_doc == std::numeric_limits<uint16_t>::max() || doc_scores[doc] > doc_scores[best_doc]) {
                                    second_best_doc = best_doc;
                                    best_doc = doc;
                                } else if (second_best_doc == std::numeric_limits<uint16_t>::max() || doc_scores[doc] > doc_scores[second_best_doc]) {
                                    second_best_doc = doc;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (movi_options->is_multi_classify()) {

        output_files->out_file << mq.get_query_id() << ",";

        float PML_mean = static_cast<float>(sum_matching_lengths) / mq.query().length();
        if (PML_mean < UNCLASSIFIED_THRESHOLD || best_doc == std::numeric_limits<uint16_t>::max()) {
            // Not present
            output_files->out_file << "0,0\n";
        } else {
            if (second_best_doc == std::numeric_limits<uint16_t>::max()) {
                output_files->out_file << to_taxon_id[best_doc] << ",0";
            } else {
                float best_doc_cnt, second_best_doc_cnt, second_best_diff;
                if (movi_options->get_min_match_len() > 0) {
                    best_doc_cnt = classify_cnts[best_doc];
                    second_best_doc_cnt = classify_cnts[second_best_doc];
                    second_best_diff = (best_doc_cnt - second_best_doc_cnt);
                } else {
                    // p-value strategy
                    best_doc_cnt = doc_scores[best_doc];
                    second_best_doc_cnt = doc_scores[second_best_doc];
                    second_best_diff = (best_doc_cnt - second_best_doc_cnt);
                }

                if (second_best_diff < 0.05 * best_doc_cnt) {
                    output_files->out_file << to_taxon_id[best_doc] << "," << to_taxon_id[second_best_doc];
                } else {
                    output_files->out_file << to_taxon_id[best_doc] << ",0";
                }
            }
            output_files->out_file << "\n";
        }
    }

    return ff_count_tot;
}

#if USE_THRESHOLDS
bool MoveStructure::reposition_thresholds(uint64_t& idx, uint64_t offset, char r_char, uint64_t& scan_count) {
    // If offset is greather than or equal to the threshold, reposition down
    // otherwise, reposition up
    uint64_t saved_idx = idx;
    uint64_t alphabet_index = alphamap[static_cast<uint64_t>(r_char)];
    scan_count = 0;

    if (movi_options->is_verbose())
        std::cerr << "\t \t \t repositioning with thresholds ... \n";

    char rlbwt_char = alphabet[rlbwt[idx].get_c()];

    if (movi_options->is_verbose()) {
        std::cerr << "\t \t \t alphabet_index: " << alphabet_index << " r_char:" << r_char << " rlbwt_char:" << rlbwt_char << "\n";
        std::cerr << "\t \t \t idx:" << idx << "\n"
                  << "\t \t \t offset: " << offset << " threshold:" << get_thresholds(idx, alphamap_3[alphamap[rlbwt_char]][alphabet_index]) << "\n";
    }

    if (idx == end_bwt_idx) {
        if (movi_options->is_verbose()) std::cerr << "\t \t \t idx == end_bwt_idx"
                                << "\n\t \t \t idx: " << idx << " end_bwt_idx: " << end_bwt_idx << " alphabet_index: " << alphabet_index << "\n"
                                << "\t \t \t end_bwt_idx_thresholds[alphabet_index]: " << end_bwt_idx_thresholds[alphabet_index] << "\n";
        if (offset >= end_bwt_idx_thresholds[alphabet_index]) {
            if (movi_options->is_verbose())
                std::cerr << "\t \t \t Repositioning down with thresholds:\n";
#if USE_NEXT_POINTERS
            if (constant) {
                scan_count += 1;
                if (end_bwt_idx_next_down[alphabet_index] == std::numeric_limits<uint16_t>::max())
                    idx = r;
                else
                    idx = saved_idx + end_bwt_idx_next_down[alphabet_index];
            }
#endif

#if THRESHOLDS_WITHOUT_NEXTS
            idx = reposition_down(saved_idx, r_char, scan_count);
#endif
            if (r_char != alphabet[rlbwt[idx].get_c()])
                std::cerr << "1: " << r_char << " " << alphabet[rlbwt[idx].get_c()];
            return false;
        } else {
            if (movi_options->is_verbose())
                std::cerr << "\t \t \t Repositioning up with thresholds:\n";
#if USE_NEXT_POINTERS
            if (constant) {
                scan_count += 1;
                if (end_bwt_idx_next_up[alphabet_index] == std::numeric_limits<uint16_t>::max())
                    idx = r;
                else
                    idx = saved_idx - end_bwt_idx_next_up[alphabet_index];
            }
#endif

#if THRESHOLDS_WITHOUT_NEXTS
            idx = reposition_up(saved_idx, r_char, scan_count);
#endif
            if (r_char != alphabet[rlbwt[idx].get_c()])
                std::cerr << "2: " << r_char << " " << alphabet[rlbwt[idx].get_c()];
            return true;
        }
    }

    if (movi_options->is_verbose()) std::cerr << "\t \t \t rlbwt[idx].get_offset(): " << get_offset(idx)
                            << " get_thresholds(idx, alphabet_index): " << get_thresholds(idx, alphamap_3[alphamap[rlbwt_char]][alphabet_index]) 
                            << "\n\t \t \t idx:" << idx << "\n";

    alphabet_index = alphamap_3[alphamap[rlbwt_char]][alphabet_index];
    if (alphabet_index == 3)
        std::cerr << "error: alphamap_3 is not working in general - " 
                    << alphabet_index << "!\n";

    if (offset >= get_thresholds(idx, alphabet_index)) {
        if (movi_options->is_verbose())
            std::cerr << "\t \t \t Repositioning down with thresholds:\n";
#if USE_NEXT_POINTERS
        if (constant) {
            scan_count += 1;
            if (rlbwt[saved_idx].get_next_down(alphabet_index) == std::numeric_limits<uint16_t>::max())
                idx = r;
            else
                idx = saved_idx + rlbwt[saved_idx].get_next_down(alphabet_index);
        } else {
            throw std::runtime_error(ERROR_MSG("[reposition thresholds] MODE is set to " + std::to_string(MODE) + ", but the constant variable is false."));
        }
#endif
        auto tmp = idx;
#if THRESHOLDS_WITHOUT_NEXTS
        idx = reposition_down(saved_idx, r_char, scan_count);
#endif
        if (r_char != alphabet[rlbwt[idx].get_c()]) {
            std::cerr << "3: " << r_char << " " << alphabet[rlbwt[idx].get_c()] << "\n";
            std::cerr << "idx: " << idx << " saved_idx: " << saved_idx << " tmp: " << tmp << "\n";
            std::cerr << "offset: " << offset << "\n";
            std::cerr << "get_thresholds(saved_idx, alphabet_index): " << get_thresholds(saved_idx, alphabet_index) << "\n";
        }
        return false;
    } else {
        if (movi_options->is_verbose())
            std::cerr << "\t \t \t Repositioning up with thresholds:\n";
#if USE_NEXT_POINTERS
        scan_count += 1;
        if (constant) {
            if (rlbwt[saved_idx].get_next_up(alphabet_index) == std::numeric_limits<uint16_t>::max())
                idx = r;
            else
                idx = saved_idx - rlbwt[saved_idx].get_next_up(alphabet_index);
        }
#endif

#if THRESHOLDS_WITHOUT_NEXTS
        idx = reposition_up(saved_idx, r_char, scan_count);
#endif
        if (r_char != alphabet[rlbwt[idx].get_c()]) {
            std::cerr << "idx: " << idx << " saved_idx: " << saved_idx << "\n";
            std::cerr << "4: " << r_char << " " << alphabet[rlbwt[idx].get_c()] << "\n";
            std::cerr << "offset: " << offset << "\n";
            std::cerr << "get_thresholds(saved_idx, alphabet_index): " << get_thresholds(saved_idx, alphabet_index) << "\n";
        }
        return true;
    }

    // TODO: default return?

    if (r_char != alphabet[rlbwt[idx].get_c()])
        std::cerr << "5: " << r_char << " " << alphabet[rlbwt[idx].get_c()];
    return false;
}
#endif

bool MoveStructure::reposition_randomly(uint64_t& idx, uint64_t& offset, char r_char, uint64_t& scan_count) {
    uint64_t saved_idx = idx;
    thread_local ThreadRandom random_generator;
    // uint16_t reposition_direction = random_generator.get_random() % 2;
    uint16_t reposition_direction =  offset * 2 < get_n(idx) ? 1 : 0;
    bool up = false;
    scan_count = 0;
    if (movi_options->is_verbose())
        std::cerr << "idx before repositioning: " << idx << "\n";


    // Selecting the right direction for trivial cases at the beginning and end of the move table
    if (idx == r - 1) {
        reposition_direction = 1;
    }
    if (idx == 0) {
        reposition_direction = 0;
    }

    if ( (reposition_direction == 1 && idx > 0) or idx == r - 1) {

        // repositioning up
        up = true;

        if (movi_options->is_verbose())
            std::cerr << "Repositioning up randomly:\n";

        idx = reposition_up(saved_idx, r_char, scan_count);
        if (movi_options->is_verbose())
            std::cerr << "idx after repositioning up: " << idx << "\n";

        if (idx >= r) {
            if (movi_options->is_verbose())
                std::cerr << "Up didn't work, try repositioning down:\n";

            // repositioning down
            up = false;
            idx = reposition_down(saved_idx, r_char, scan_count);
            if (movi_options->is_verbose())
                std::cerr << "idx after repositioning down: " << idx << "\n";
            if (idx == r) {
                // TODO
                throw std::runtime_error(ERROR_MSG("[reposition randomly] Neither up or down repositioning works.\n" +
                                    "The character does not exist in the index.\n"));
                throw std::runtime_error(ERROR_MSG("[reposition randomly] Neither up or down repositioning works.\n"));
            }
        }
    } else {

        // repositioning down
        up = false;

        if (movi_options->is_verbose())
            std::cerr << "Repositioning down randomly:\n";

        idx = reposition_down(saved_idx, r_char, scan_count);
        if (movi_options->is_verbose())
            std::cerr << "idx after repositioning down: " << idx << "\n";

        if (idx >= r) {
            if (movi_options->is_verbose())
                std::cerr << "Down didn't work, try repositioning up:\n";

            // repositioning up
            up = true;
            idx = reposition_up(saved_idx, r_char, scan_count);
            if (movi_options->is_verbose())
                std::cerr << "idx after repositioning up: " << idx << "\n";
            if (idx == r) {
                // TODO
                throw std::runtime_error(ERROR_MSG("[reposition randomly] Neither up or down repositioning works.\n" +
                                                   "The character does not exist in the index.\n"));
            }
        }
    }

    // sanity check
    char c = alphabet[rlbwt[idx].get_c()];
    if (c != r_char or idx == r) {
        throw std::runtime_error(ERROR_MSG("[reposition randomly] This should never happen.""c: " + c + "\n" +
                                           "idx: " + std::to_string(idx) + "\n"));
    }

    return up;
}

uint64_t MoveStructure::query_zml(MoveQuery& mq) {
    auto& query_seq = mq.query();
    int32_t pos_on_r = query_seq.length() - 1;
    uint64_t match_len = 0;
    uint64_t ff_count_tot = 0;

    while (!check_alphabet(query_seq[pos_on_r]) and pos_on_r >= 0) {
        mq.add_ml(0, movi_options->is_stdout());
        pos_on_r -= 1;
    }

    if (pos_on_r < 0) {
        // Special case where no character in the read exists in the index.
        return 0;
    }

    // Multi-class classification
    if (movi_options->is_multi_classify()) {
        for (uint16_t i = 0; i < num_species; i++) {
            classify_cnts[i] = 0;
        }
    }

    MoveInterval interval = initialize_backward_search(mq, pos_on_r, match_len);
    while (pos_on_r > 0) {
        MoveInterval prev_interval = interval;
        ff_count_tot += backward_search_step(query_seq, pos_on_r, interval);
        if (!interval.is_empty()) {
            mq.add_ml(match_len, movi_options->is_stdout());
            pos_on_r -= 1;
            match_len += 1;
        } else {
            // Classification based on maximal matching
            if (movi_options->is_multi_classify() && match_len >= movi_options->get_min_match_len()) {
                if (prev_interval.run_start == prev_interval.run_end) {
                    for (uint64_t i = prev_interval.offset_start; i <= prev_interval.offset_end; i++) {
                        uint64_t full_ind = run_offsets[prev_interval.run_start] + i;
                        uint16_t cur_doc = doc_pats[full_ind];
                        classify_cnts[cur_doc] += match_len;
                    }   
                } else {
                    for (uint64_t i = prev_interval.offset_start; i < get_n(prev_interval.run_start); i++) {
                        uint64_t full_ind = run_offsets[prev_interval.run_start] + i;
                        uint16_t cur_doc = doc_pats[full_ind];
                        classify_cnts[cur_doc] += match_len;
                    }    
                    for (uint64_t r_ind = prev_interval.run_start + 1; r_ind < prev_interval.run_end; r_ind++) {
                        for (uint64_t i = 0; i < get_n(r_ind); i++) {
                            uint64_t full_ind = run_offsets[r_ind] + i;
                            uint16_t cur_doc = doc_pats[full_ind];
                            classify_cnts[cur_doc] += match_len;
                        }  
                    } 
                    for (uint64_t i = 0; i <= prev_interval.offset_end; i++) {
                        uint64_t full_ind = run_offsets[prev_interval.run_end] + i;
                        uint16_t cur_doc = doc_pats[full_ind];
                        classify_cnts[cur_doc] += match_len;
                    }           
                }
            }

            mq.add_ml(match_len, movi_options->is_stdout());
            pos_on_r -= 1;
            match_len = 0;
            while (!check_alphabet(query_seq[pos_on_r]) and pos_on_r > 0) {
                mq.add_ml(match_len, movi_options->is_stdout());
                pos_on_r -= 1;
            }
            // Special case where the character at position 0 of the read does not exist in the index.
            if (check_alphabet(query_seq[pos_on_r]))
                interval = initialize_backward_search(mq, pos_on_r, match_len);
        }
    }
    if (interval.is_empty()) {
        match_len = 0;
    }
    mq.add_ml(match_len, movi_options->is_stdout());

    // Document occuring the most is the genotype we think the query is from.
    if (movi_options->is_multi_classify()) {
        uint16_t best_doc = 0;
        for (uint16_t i = 1; i < num_species; i++) {
            if (classify_cnts[i] > classify_cnts[best_doc]) {
                best_doc = i;
            }
        }
        
        // Document occuring the most is the genotype we think the query is from.
        output_files->out_file << to_taxon_id[best_doc] << " ";
        for (uint16_t i = 0; i < num_species; i++) {
        //    output_files->out_file << classify_cnts[i] << " ";
        }
        output_files->out_file << "\n";
    }

    return ff_count_tot;
}