#include "read_processor.hpp"
#include <cpuid.h>

ReadProcessor::ReadProcessor(std::string reads_file_name, MoveStructure& mv_, int strands_ = 4, bool verbose_ = false, bool reverse_ = false) : mv(mv_) {

    int cpu_info[4];
    __cpuid(0x80000006, cpu_info[0], cpu_info[1], cpu_info[2], cpu_info[3]);
    cache_line_size = (cpu_info[2] & 0xFF);
    prefetch_step = cache_line_size/sizeof(MoveRow) - 1;

    verbose = verbose_;
    reverse = mv_.movi_options->is_reverse();
    // Solution for handling the stdin input: https://biowize.wordpress.com/2013/03/05/using-kseq-h-with-stdin/
    if (reads_file_name == "-") {
        FILE *instream = stdin;
        fp = gzdopen(fileno(instream), "r"); // STEP 2: open the file handler
    } else {
        fp = gzopen(reads_file_name.c_str(), "r"); // STEP 2: open the file handler
    }

    // seq = kseq_init(fp); // STEP 3: initialize seq
    // Open output files using the utility function
    open_output_files(*(mv_.movi_options), output_files);

    if (mv_.movi_options->is_classify()) {
        classifier.initialize_report_file(*mv_.movi_options);
    }

    total_kmer_count = 0;
    positive_kmer_count = 0;
    negative_kmer_count = 0;
    kmer_extension_count = 0;
    kmer_extension_stopped_count = 0;
    negative_kmer_extension_count = 0;
    read_processed = 0;
    strands = strands_;
}

bool ReadProcessor::next_read(Strand& process, BatchLoader& reader) {

    // l = kseq_read(seq); // STEP 4: read sequence
    bool valid_read = false;
    Read read_struct;
    valid_read = reader.grabNextRead(read_struct);

    if (valid_read) {
        #pragma omp atomic
        read_processed += 1;

        // process.st_length = seq->name.m;
        process.st_length = read_struct.id.length();
        // Read name is now stored in MoveQuery via set_query_id()
        // Create sequence string and optionally reverse it
        std::string sequence = std::string(read_struct.seq);
        if (reverse) {
            std::reverse(sequence.begin(), sequence.end());
        }
        process.match_len = 0;
        process.ff_count = 0;
        process.scan_count = 0;

        process.mq = MoveQuery(sequence);
        process.mq.set_query_id(read_struct.id);

        return false;
    } else {
        // std::cerr << "No more reads to process!\n";
        return true;
    }
}

void ReadProcessor::reset_process(Strand& process, BatchLoader& reader) {
    process.finished = next_read(process, reader);
    if (!process.finished) {
        process.length_processed = 0;
        process.pos_on_r = process.mq.query().length() - 1;
        process.idx = mv.r - 1;
        process.offset = mv.get_n(process.idx) - 1;

        // reset the multi-class classification variables
        process.best_doc = std::numeric_limits<uint16_t>::max();
        process.second_best_doc = std::numeric_limits<uint16_t>::max();
        process.sum_matching_lengths = 0;
        process.colors_count = 0;
        if (mv.movi_options->is_multi_classify()) {
            std::fill(process.classify_cnts.begin(), process.classify_cnts.end(), 0);
        }
#if TALLY_MODES
        // This is necessary
        process.tally_state = false;
        // These shouldn't matter
        process.id_found = true;
        process.tally_b = 0;
        process.rows_until_tally = 0;
        process.next_check_point = 0;
        process.run_id = 0;
        process.last_id = 0;
        process.char_index = 0;
        process.tally_offset = 0;
        process.find_next_id_attempt = 0;
#endif
    }
}

// This function is for computing PMLs
void ReadProcessor::process_char(Strand& process) {
    if (process.pos_on_r < process.mq.query().length() - 1) {
        // LF step
        // if (mv.logs)
        //     process.t3 = std::chrono::high_resolution_clock::now();
        process.ff_count = mv.LF_move(process.offset, process.idx);
        if (mv.movi_options->is_logs()) {
            // auto t4 = std::chrono::high_resolution_clock::now();
            auto t2 = std::chrono::high_resolution_clock::now();
            // auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - process.t1);
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
            // auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(process.t2 - process.t1) + 
            //                std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - process.t3);
            process.mq.add_cost(elapsed);
            process.mq.add_fastforward(process.ff_count);
            process.mq.add_scan(process.scan_count);
        }
    }
    if (mv.movi_options->is_logs()) {
        process.t1 = std::chrono::high_resolution_clock::now();
        t1 = process.t1;
    }

    if (mv.movi_options->is_multi_classify()) {
        if (process.match_len >= mv.movi_options->get_min_match_len()) {

            process.colors_count += 1;

            uint64_t color_id;
#if COLOR_MODE == 1
            color_id = static_cast<uint64_t>(mv.rlbwt[process.idx].color_id);
            // Skip doc sets that weren't saved (thrown away by compression).
            if (color_id >= mv.unique_doc_sets.size()) return;
#else
            if (mv.movi_options->is_doc_sets_vector_of_vectors()) {
                color_id = static_cast<uint64_t>(mv.doc_set_inds[process.idx]);
                // Skip doc sets that weren't saved (thrown away by compression).
                if (color_id >= mv.unique_doc_sets.size()) {
                    // std::cerr << "doc_set_inds[idx] >= unique_doc_sets.size()\n";
                    // std::cerr << "This should not happen when compression is not turned on.\n";
                    // std::cerr << "The compressed version of the prefetching mode is not supported yet.\n";
                    // exit(0);
                    return;
                }
            } else {
                color_id = mv.doc_set_flat_inds[process.idx].get();
                // Skip doc sets that weren't saved (thrown away by compression).
                if (color_id >= mv.flat_colors.size()) return;
            }
#endif
            std::span<uint16_t> cur_set;
            if (mv.movi_options->is_doc_sets_vector_of_vectors()) {
                std::vector<uint16_t> &cur_set_vec = mv.unique_doc_sets[color_id];
                cur_set = std::span<uint16_t>(cur_set_vec.data(), cur_set_vec.size());
            } else {
                uint32_t cur_set_size = mv.flat_colors[color_id];
                cur_set = std::span<uint16_t>(mv.flat_colors.data() + color_id + 1,
                                              mv.flat_colors.data() + color_id + 1 + cur_set_size);
            }

            for (int doc : cur_set) {
                process.classify_cnts[doc]++;
                if (doc != process.best_doc) {
                    if (process.best_doc == std::numeric_limits<uint16_t>::max() || process.classify_cnts[doc] > process.classify_cnts[process.best_doc]) {
                        process.second_best_doc = process.best_doc;
                        process.best_doc = doc;
                    } else if (process.second_best_doc == std::numeric_limits<uint16_t>::max() || process.classify_cnts[doc] > process.classify_cnts[process.second_best_doc]) {
                        process.second_best_doc = doc;
                    }
                }
            }

            if (mv.movi_options->is_report_color_ids()) {
                process.mq.add_color(mv.color_offset_to_id[color_id]);
            } else if (mv.movi_options->is_report_colors()) {
                process.mq.add_color(color_id);
            }

        } else {

            if (mv.movi_options->is_report_color_ids()) {
                process.mq.add_color(mv.num_colors);
            } else if (mv.movi_options->is_report_colors()) {
                process.mq.add_color(mv.flat_colors.size());
            }

        }
    }

    auto& row = mv.rlbwt[process.idx];
    uint64_t row_idx = process.idx;
    char row_c = mv.alphabet[row.get_c()];
    std::string& R = process.mq.query();
    process.scan_count = 0;
    if (mv.alphamap[static_cast<uint64_t>(R[process.pos_on_r])] == mv.alphamap.size()) {
        process.match_len = 0;
    } else if (row_c == R[process.pos_on_r]) {
        // Case 1
        process.match_len += 1;
    } else {
        // Case 2
        // Repositioning up or down (randomly or with thresholds)
        uint64_t idx_before_reposition = process.idx;

        bool up = false;
#if USE_THRESHOLDS
        up = mv.movi_options->is_random_repositioning() ?
                mv.reposition_randomly(process.idx, process.offset, R[process.pos_on_r], process.scan_count) :
                mv.reposition_thresholds(process.idx, process.offset, R[process.pos_on_r], process.scan_count);
#else
        // When there is no threshold, reposition randomly
        up = mv.reposition_randomly(process.idx, process.offset, R[process.pos_on_r], process.scan_count);
#endif
        process.match_len = 0;
        char c = mv.alphabet[mv.rlbwt[process.idx].get_c()];
        // sanity check
        if (c == R[process.pos_on_r]) {
            // Observing a match after the repositioning
            // The right match_len should be:
            // min(new_lcp, match_len + 1)
            // But we cannot compute lcp here
            process.offset = up ? mv.get_n(process.idx) - 1 : 0;
            if (verbose)
                std::cerr << "\t idx: " << process.idx << " offset: " << process.offset << "\n";
        } else {
            std::cerr << "\t \t This should not happen!\n";
        }
    }

    process.mq.add_ml(process.match_len, mv.movi_options->write_stdout_enabled());
    if (mv.movi_options->is_get_sa_entries()) {
        uint64_t sa_entry = mv.get_SA_entries(process.idx, process.offset);
        process.mq.add_sa_entries(sa_entry);
    }

    process.sum_matching_lengths += process.match_len;
    process.pos_on_r -= 1;

    // Check for early stopping if the read is unclassified
    if ( mv.movi_options->is_early_stop() and mv.movi_options->is_multi_classify() ) {
        if ( process.pos_on_r < process.mq.query().length() / 2 and process.pos_on_r % 100 == 0 ) {
            float PML_mean = static_cast<float>(process.sum_matching_lengths) / (process.mq.query().length() - process.pos_on_r);
            if (PML_mean < UNCLASSIFIED_THRESHOLD) {
                // setting the pos_on_r to -1 will stop the read processing
                process.pos_on_r = -1;
            }
        }
    }

    // if (mv.logs)
    //     process.t2 = std::chrono::high_resolution_clock::now();
    // LF step should happen here in the non-prefetch code
    // uint64_t ff_count = mv.LF_move(process.offset, process.idx);
    // process.ff_count_tot += ff_count;
}

#if TALLY_MODES
void ReadProcessor::process_char_tally(Strand& process) {
    if (process.tally_state == false) {

        if (process.pos_on_r < process.mq.query().length() - 1) {

            // First find the id if not found already
            if (process.id_found == false) {
                find_next_id(process);
                if (process.id_found == false) {
                    return;
                } else {
                    process.find_next_id_attempt = 0;
                }
            }

            // LF step
            process.ff_count = mv.LF_move(process.offset, process.idx, process.run_id);
        }

        // After LF is performed, calculate PML based on case1/case2
        auto& row = mv.rlbwt[process.idx];
        uint64_t row_idx = process.idx;
        char row_c = mv.alphabet[row.get_c()];
        std::string& R = process.mq.query();
        process.scan_count = 0;
        if (mv.alphamap[static_cast<uint64_t>(R[process.pos_on_r])] == mv.alphamap.size()) {
            process.match_len = 0;
        } else if (row_c == R[process.pos_on_r]) {
            // Case 1
            process.match_len += 1;
        } else {
            // Case 2
            // Repositioning up or down (randomly or with thresholds)
            uint64_t idx_before_reposition = process.idx;

            bool up = false;
#if USE_THRESHOLDS
            up = mv.movi_options->is_random_repositioning() ?
                    mv.reposition_randomly(process.idx, process.offset, R[process.pos_on_r], process.scan_count) :
                    mv.reposition_thresholds(process.idx, process.offset, R[process.pos_on_r], process.scan_count);
#else
            // When there is no threshold, reposition randomly
            up = mv.reposition_randomly(process.idx, process.offset, R[process.pos_on_r], process.scan_count);
#endif
            process.match_len = 0;
            char c = mv.alphabet[mv.rlbwt[process.idx].get_c()];
            // sanity check
            if (c == R[process.pos_on_r]) {
                // Observing a match after the repositioning
                // The right match_len should be:
                // min(new_lcp, match_len + 1)
                // But we cannot compute lcp here
                process.offset = up ? mv.get_n(process.idx) - 1 : 0;
                if (verbose)
                    std::cerr << "\t idx: " << process.idx << " offset: " << process.offset << "\n";
            } else {
                std::cerr << "\t \t This should not happen!\n";
            }
        }
        process.mq.add_ml(process.match_len, mv.movi_options->write_stdout_enabled());
        if (mv.movi_options->is_get_sa_entries()) {
            uint64_t sa_entry = mv.get_SA_entries(process.idx, process.offset);
            process.mq.add_sa_entries(sa_entry);
        }

        process.pos_on_r -= 1;

        // get_id is called at the beginning of the next LF
        // Now we should find out what is the next tally_b to prefetch it
        // Set tally_state to true, so that the other branch for id computation is called
        process.id_found = false;
        find_tally_b(process);
    } else {
        // The id for the run at the next_check_point
        process.run_id = mv.tally_ids[process.char_index][process.tally_b].get();

        count_rows_untill_tally(process);

        // The id stored at the checkpoint is actually the id for the current run
        // because there was no run with the same character between the current run and the next_check_point
        if (process.last_id == process.idx and mv.get_char(process.idx) != mv.get_char(process.next_check_point)) {
            // process.run_id is already set to the correct value above
            // process.run_id = mv.tally_ids[process.char_index][process.tally_b].get();

            process.id_found = true;
        }

        process.tally_state = false;
    }
}

void ReadProcessor::find_next_id(Strand& process) {

    // Sanity check, the offset in the destination run should be always smaller than the length of that run
    if (process.tally_offset >= mv.get_n(process.run_id)) {
        throw std::runtime_error(ERROR_MSG("[Read Processor - find next tally id] tally_offset: " +std::to_string(process.tally_offset) +
                                           " n: " + std::to_string(mv.get_n(process.run_id)) +
                                           "\nThe offset in the destination run should be always smaller than the length of that run."));
    }

    if (process.find_next_id_attempt == 0) {
        if (process.tally_offset >= process.rows_until_tally) {
            process.id_found = true;
            return;
        } else {
            process.rows_until_tally -= (process.tally_offset + 1);
            process.run_id -= 1;
            process.tally_offset = 0;
        }
    }

    int rows_visited = 0;
    while (process.rows_until_tally != 0 and rows_visited < mv.tally_checkpoints) {
        if (process.rows_until_tally >= mv.get_n(process.run_id)) {
            process.rows_until_tally -= mv.get_n(process.run_id);
            process.run_id -= 1;
            rows_visited += 1;
        } else {
            process.rows_until_tally = 0;
        }
    }
    // The id is not found yet, more prefetching is needed.
    if (process.rows_until_tally != 0) {
        process.find_next_id_attempt += 1;
        process.id_found = false;
        return;
    }

    process.id_found = true;
}

void ReadProcessor::count_rows_untill_tally(Strand& process) {

    process.rows_until_tally = 0;

    // // The id for the run at the next_check_point
    // process.run_id = mv.tally_ids[process.char_index][process.tally_b].get();

    uint64_t last_count = 0;
    process.last_id = mv.r;

    // Look for the rows between idx and the next_check_point
    // Count how many rows with the same character exists between
    // the current run head and the next run head for which we know the id
    // dbg << "from: " << idx << " to: " << next_check_point << "\n";
    for (uint64_t i = process.idx; i < process.next_check_point; i++) {
        // Only count the rows with the same character
        if (mv.get_char(i) == mv.get_char(process.idx)) {
            process.rows_until_tally += mv.get_n(i);
            process.last_id = i;
        }
    }

    // // The id stored at the checkpoint is actually the id for the current run
    // // because there was no run with the same character between the current run and the next_check_point
    // if (last_id == process.idx and mv.get_char(process.idx) != mv.get_char(process.next_check_point)) {
    //     // process.run_id is already set to the correct value above
    //     // process.run_id = mv.tally_ids[process.char_index][process.tally_b].get();

    //     process.tally_state = false;
    //     return;
    // }

    if (process.last_id == mv.r) {
        throw std::runtime_error(ERROR_MSG("[Read Processor - count rows until tally] last_id should never be equal to r = " + std::to_string(mv.r) + "\n"));
    }

    // We should know what will be the offset of the run head in the destination run (id)
    process.tally_offset = mv.get_offset(process.next_check_point);

    // At the checkpoint, one id is stored for each character
    // For the character of the checkpoint, the id of the checkpoint is stored
    // For other characters, the id of the last run before the checkpoint with that character is stored
    // So, we have to decrease the number of rows of the last run as they are not between the current run
    // and the run for which we know the id
    // Also, for other characters, offset should be stored based on that run (not the checkpoint run)
    if (mv.get_char(process.idx) != mv.get_char(process.next_check_point)) {
        process.rows_until_tally -= mv.get_n(process.last_id);
        process.tally_offset = mv.get_offset(process.last_id);
    }
}

void ReadProcessor::find_tally_b(Strand& process) {
    uint64_t id = 0;

    // The $ always goes to row/run 0
    if (process.idx == mv.end_bwt_idx) {
        process.run_id = 0;

        process.tally_state = false;
        process.id_found = true;
        return;
    }

    process.char_index = mv.rlbwt[process.idx].get_c();

    // The id of the last run is always stored at the last tally_id row
    if (process.idx == mv.r - 1) {
        uint64_t tally_ids_len = mv.tally_ids[process.char_index].size();
        process.run_id = mv.tally_ids[process.char_index][tally_ids_len - 1].get();

        process.tally_state = false;
        process.id_found = true;
        return;
    }

    uint64_t tally_a = process.idx / mv.tally_checkpoints;

    // If we are at a checkpoint, simply return the stored id
    if (process.idx % mv.tally_checkpoints == 0) {
        process.run_id = mv.tally_ids[process.char_index][tally_a].get();

        process.tally_state = false;
        process.id_found = true;
        return;
    }

    // find the next check point
    process.tally_b = tally_a + 1;

    process.next_check_point = process.tally_b * mv.tally_checkpoints;
    // If we are at the end, just look at the last run
    if (process.tally_b * mv.tally_checkpoints >= mv.r) {
        process.next_check_point = mv.r - 1;
    }

    process.tally_state = true;
}
#endif

void ReadProcessor::write_mls(Strand& process) {

    if (mv.movi_options->is_multi_classify()) {
        output_files.out_file << process.mq.get_query_id() << ",";
        // binary classification in the multi-class classification mode is handled at a different part of the code
        // if (mv.movi_options->is_classify() && !mv.classifier->is_present(process.mq.get_matching_lengths(), *mv.movi_options)) {
        float PML_mean = static_cast<float>(process.sum_matching_lengths) / process.mq.query().length();
        if (PML_mean < UNCLASSIFIED_THRESHOLD || process.best_doc == std::numeric_limits<uint16_t>::max()) {
            // Not present
            if (mv.movi_options->is_report_all()) {
                output_files.out_file << "0\n";
            } else {
                output_files.out_file << "0,0\n";
            }
        } else {
            if ((mv.movi_options->is_report_colors() or mv.movi_options->is_report_color_ids()) && mv.movi_options->write_output_allowed()) {
                // Writing the PMLs
                output_base_stats(DataType::match_length, mv.movi_options->write_stdout_enabled(), output_files.mls_file, process.mq);

                // Writing the colors
                output_base_stats(DataType::color, mv.movi_options->write_stdout_enabled(), output_files.colors_file, process.mq);
            }

            // If the second most occurring document is more than 95% of the most occurring one,
            // we report the other species as well and classify the read at a higher level.
            // out_file << mv.to_taxon_id[process.best_doc];
            if (mv.movi_options->is_report_all()) {

                if (mv.movi_options->get_min_score_frac() == 0) {
                    // For the min_diff_frac mode, we write the best document no matter what
                    // For the min_score_frac mode, the best document is outputed only if its score is high enough
                    output_files.out_file << mv.to_taxon_id[process.best_doc];
                }

                uint32_t output_document_count = 0;
                uint32_t best_doc_cnt = process.classify_cnts[process.best_doc];
                for (int i = 0; i < process.classify_cnts.size(); i++) {
                    if (mv.movi_options->get_min_score_frac() == 0) {

                        float diff_best = static_cast<float>(best_doc_cnt - process.classify_cnts[i]);

                        if (i!= process.best_doc and diff_best < mv.movi_options->get_min_diff_frac() * best_doc_cnt) {
                                output_files.out_file << "," << mv.to_taxon_id[i];
                        }
                    } else {

                        if (static_cast<float>(process.classify_cnts[i]) >= mv.movi_options->get_min_score_frac() * process.colors_count) {
                            output_files.out_file << "," << mv.to_taxon_id[i]; // << ":" << process.classify_cnts[i] << "/" << process.colors_count << "/" << process.mq.query().length();
                            output_document_count += 1;
                        }
                    }
                }

                if (mv.movi_options->get_min_score_frac() != 0 and output_document_count == 0) {
                    output_files.out_file << "0";
                }
            } else {
                if (process.second_best_doc == std::numeric_limits<uint16_t>::max()) {
                    output_files.out_file << mv.to_taxon_id[process.best_doc] << ",0";
                } else {

                    uint32_t best_doc_cnt = process.classify_cnts[process.best_doc];
                    uint32_t second_best_doc_cnt = process.classify_cnts[process.second_best_doc];
                    float second_best_diff = static_cast<float>(best_doc_cnt - second_best_doc_cnt);

                    if (second_best_diff < 0.05 * best_doc_cnt) {
                        output_files.out_file << mv.to_taxon_id[process.best_doc] << "," << mv.to_taxon_id[process.second_best_doc];
                    } else {
                        output_files.out_file << mv.to_taxon_id[process.best_doc] << ",0";
                    }
                }
            }
            output_files.out_file << "\n";
        }
    } else {

        if (mv.movi_options->is_classify()) {

            // Classification with 32 bits pmls works if the pml values are less than 2^16.
            std::vector<uint16_t> matching_lens_16(process.mq.get_matching_lengths().begin(), process.mq.get_matching_lengths().end());

            bool found = classifier.classify(process.mq.get_query_id(), matching_lens_16, *mv.movi_options);
            if (found and mv.movi_options->is_filter() && !mv.movi_options->is_no_output()) {
                output_read(process.mq);
            }
        }

        if (mv.movi_options->write_output_allowed()) {

            output_base_stats(DataType::match_length, mv.movi_options->write_stdout_enabled(), output_files.mls_file, process.mq);

            if (mv.movi_options->is_get_sa_entries()) {
                output_base_stats(DataType::sa_entry, mv.movi_options->write_stdout_enabled(), output_files.sa_entries_file, process.mq);
            }

            bool logs = mv.movi_options->is_logs();
            if (logs) {
                auto t2 = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
                process.mq.add_cost(elapsed);
                process.mq.add_fastforward(process.ff_count);
                process.mq.add_scan(process.scan_count);
                output_logs(output_files.costs_file, output_files.scans_file, output_files.fastforwards_file, process.mq);
            }

        }
    }
}

void ReadProcessor::write_count(Strand& process) {
    compute_match_count(process);
    if (mv.movi_options->write_output_allowed()) {
        auto& R = process.mq.query();
        output_counts(mv.movi_options->write_stdout_enabled(), output_files.matches_file, R.length(), process.pos_on_r, process.match_count, process.mq);
    }
}

void ReadProcessor::compute_match_count(Strand& process) {
    if (process.range.is_empty() and process.range_prev.is_empty())
        process.match_count = 0;
    else if (process.range.is_empty() and !process.range_prev.is_empty())
        process.match_count = process.range_prev.count(mv.rlbwt);
    else if (!process.range.is_empty())
        process.match_count = process.range.count(mv.rlbwt);
    else {
        throw std::runtime_error(ERROR_MSG("[Read Processor - compute match count] This should not happen, for match count computation.\n"));
    }
}

void ReadProcessor::end_process() {

    std::cerr << "no_ftab: " << mv.no_ftab << "\n";
    std::cerr << "all_initializations: " << mv.all_initializations << "\n";


    if (!mv.movi_options->is_filter()) {
        bool is_pml = mv.movi_options->is_pml();
        bool is_zml = mv.movi_options->is_zml();
        bool is_count = mv.movi_options->is_count();

        if (mv.movi_options->is_classify()) {
            classifier.close_report_file();
        }

        // Close output files using the utility function
        close_output_files(*(mv.movi_options), output_files);

        // kseq_destroy(seq); // STEP 5: destroy seq
        // std::cerr << "kseq destroyed!\n";
        // gzclose(fp); // STEP 6: close the file handler
        // std::cerr << "fp file closed!\n";
    }
}

void ReadProcessor::process_latency_hiding(BatchLoader& reader) {
    bool is_pml = mv.movi_options->is_pml();
    bool is_zml = mv.movi_options->is_zml();
    bool is_count = mv.movi_options->is_count();

    if ((is_pml and is_count) or (is_pml and is_zml) or (is_count and is_zml)) {
        throw std::runtime_error(ERROR_MSG("[Read Processor - process latency hiding] Multipe types of queries cannot be processed at the same time.\n"));
    }

    std::vector<Strand> processes;
    uint64_t finished_count = 0;
    #pragma omp critical
    {
        finished_count = initialize_strands(processes, reader);
        // std::cerr << strands << " processes are initiated.\n";
    }

    while (finished_count != strands) {
        for (uint64_t i = 0; i < strands; i++) {
            if (!processes[i].finished) {
                // 1: process next character -- doing fast forward
                bool backward_search_finished = false; // used for the count and zml queries
                if (is_pml) {
                    process_char(processes[i]);
                } else if (is_count) {
                    // backward_search_finished = backward_search(processes[i], processes[i].mq.query().length() - 1);
                    backward_search_finished = backward_search(processes[i], processes[i].kmer_end);
                } else if (is_zml) {
                    backward_search_finished = backward_search(processes[i], processes[i].kmer_end);
                }
                // 2: if the read is done -> Write the pmls and go to next read
                if ((is_pml and processes[i].pos_on_r <= -1) or
                    (is_count and backward_search_finished) or
                    (is_zml and (backward_search_finished and processes[i].pos_on_r > 0)) or
                    (is_zml and (processes[i].pos_on_r <= 0))) {
                    #pragma omp critical
                    {
                        if (read_processed % 1000 == 0)
                            std::cerr << read_processed << "\r";
                        if (is_pml) {
                            write_mls(processes[i]);
                            reset_process(processes[i], reader);
                        } else if (is_count) {
                            write_count(processes[i]);
                            reset_process(processes[i], reader);
                            reset_backward_search(processes[i]);
                        } else if (is_zml) {
                            if (backward_search_finished and processes[i].pos_on_r > 0) {
                                processes[i].mq.add_ml(processes[i].match_len, mv.movi_options->write_stdout_enabled());
                                if (mv.movi_options->is_get_sa_entries()) {
                                    uint64_t sa_entry = mv.get_SA_entries(processes[i].idx, processes[i].offset);
                                    processes[i].mq.add_sa_entries(sa_entry);
                                }

                                processes[i].pos_on_r -= 1;
                                reset_backward_search(processes[i]);
                                processes[i].kmer_end = processes[i].pos_on_r;
                                // continue;
                            } else if (processes[i].pos_on_r <= 0) {
                                processes[i].mq.add_ml(processes[i].match_len, mv.movi_options->write_stdout_enabled());
                                if (mv.movi_options->is_get_sa_entries()) {
                                    uint64_t sa_entry = mv.get_SA_entries(processes[i].idx, processes[i].offset);
                                    processes[i].mq.add_sa_entries(sa_entry);
                                }

                                write_mls(processes[i]);
                                reset_process(processes[i], reader);
                                reset_backward_search(processes[i]);
                                processes[i].kmer_end = processes[i].pos_on_r;
                            }
                        }
                        // 3: -- check if it was the last read in the file -> finished_count++
                        if (processes[i].finished) {
                            finished_count += 1;
                        }
                    }
                } else {
                    // 4: big jump with prefetch
                    if (is_pml) {
                        my_prefetch_r((void*)(&(mv.rlbwt[0]) + mv.get_id(processes[i].idx)));
                    } else if (is_count or is_zml) {
                        my_prefetch_r((void*)(&(mv.rlbwt[0]) + mv.get_id(processes[i].range.run_start)));
                        my_prefetch_r((void*)(&(mv.rlbwt[0]) + mv.get_id(processes[i].range.run_end)));
                    }
                }
            }
        }
    }
}

#if TALLY_MODES
void ReadProcessor::process_latency_hiding_tally(BatchLoader& reader) {
    bool is_pml = mv.movi_options->is_pml();

    std::vector<Strand> processes;
    uint64_t finished_count = 0;
    #pragma omp critical
    {
        finished_count = initialize_strands(processes, reader);
        // std::cerr << strands << " processes are initiated.\n";
    }

    while (finished_count != strands) {
        for (uint64_t i = 0; i < strands; i++) {
            if (!processes[i].finished) {
                if (is_pml) {
                    process_char_tally(processes[i]);
                }
                // 2: if the read is done -> Write the pmls and go to next read
                if (is_pml and processes[i].pos_on_r <= -1) {
                    #pragma omp critical
                    {
                        if (read_processed % 1000 == 0)
                            std::cerr << read_processed << "\r";

                        write_mls(processes[i]);
                        reset_process(processes[i], reader);
                    }
                }
                // 3: -- check if it was the last read in the file -> finished_count++
                if (processes[i].finished) {
                    finished_count += 1;
                } else {
                    // 4: prefetching
                    if (processes[i].tally_state) {
                        // prefetch tally
                        my_prefetch_r((void*)(&(mv.tally_ids[processes[i].char_index][0]) + processes[i].tally_b));
                        // prefetch following rows until the checkpoint
                        // Every prefetch loads 64 bytes which is about 20 move rows
                        for (uint64_t tally = processes[i].idx; tally <= processes[i].next_check_point; tally += prefetch_step)
                            my_prefetch_r((void*)(&(mv.rlbwt[0]) + tally));
                    } else {
                        // prefetch tally.id
                        for (uint64_t tally = 0; tally <= mv.tally_checkpoints; tally += prefetch_step)
                            my_prefetch_r((void*)(&(mv.rlbwt[0]) + processes[i].run_id - tally));
                    }
                }
            }
        }
    }
}
#endif

/* void ReadProcessor::ziv_merhav_latency_hiding() {
    std::vector<Strand> processes;
    for(int i = 0; i < strands; i++) processes.emplace_back(Strand());
    std::cerr << strands << " processes are created.\n";

    uint64_t finished_count = 0;
    for (uint64_t i = 0; i < strands; i++) {
        if (finished_count == 0) {
            // TODO:: update the reset function
            reset_process(processes[i]);
            reset_backward_search(processes[i]);
            processes[i].kmer_end = processes[i].pos_on_r;
            processes[i].match_len = 0;
        } else {
            processes[i].finished = true;
        }
        if (processes[i].finished) {
            std::cerr << "Warning: less than strands = " << strands << " reads.\n";
            finished_count += 1;
        }
    }
    std::cerr << strands << " processes are initiated.\n";

    uint64_t total_bs = 0;
    while (finished_count != strands) {
        for (uint64_t i = 0; i < strands; i++) {
            if (!processes[i].finished) {
                // 1: process next character
                bool backward_search_finished = backward_search(processes[i], processes[i].kmer_end);
                total_bs += 1;
                if (verbose)
                    std::cerr << backward_search_finished << " " << processes[i].kmer_end << " " << processes[i].pos_on_r << "\n";

                if (backward_search_finished and processes[i].pos_on_r > 0) {
                    reset_backward_search(processes[i]);
                    processes[i].kmer_end = processes[i].pos_on_r - 1;
                    processes[i].match_len = 0;
                } else if (processes[i].pos_on_r <= 0) {
                    processes[i].mq.add_ml(processes[i].match_len);
                    write_mls(processes[i]);
                    reset_process(processes[i]);
                    reset_backward_search(processes[i]);
                    processes[i].kmer_end = processes[i].pos_on_r;
                    processes[i].match_len = 0;
                    // 3: -- check if it was the last read in the file -> finished_count++
                    if (processes[i].finished) {
                        finished_count += 1;
                    }
                } else {
                    // 4: big jump with prefetch
                    if (!backward_search_finished) {
                        processes[i].mq.add_ml(processes[i].match_len);
                        processes[i].match_len += 1;
                    }
                    my_prefetch_r((void*)(&(mv.rlbwt[0]) + mv.rlbwt[processes[i].range.run_start].get_id()));
                    my_prefetch_r((void*)(&(mv.rlbwt[0]) + mv.rlbwt[processes[i].range.run_end].get_id()));
                }
            }
        }
    }

    std::cerr << "\n";
    std::cerr << "total_bs for ziv_merhav: " << total_bs << "\n";

    std::cerr << "The output file for the matching lengths closed.\n";

    kseq_destroy(seq); // STEP 5: destroy seq
    std::cerr << "kseq destroyed!\n";
    gzclose(fp); // STEP 6: close the file handler
    std::cerr << "fp file closed!\n";
} */

uint64_t ReadProcessor::initialize_strands(std::vector<Strand>& processes, BatchLoader& reader) {
    uint64_t finished_count = 0;
    uint64_t empty_strands = 0;
    for(int i = 0; i < strands; i++) {
        processes.emplace_back(Strand());
        if (mv.movi_options->is_multi_classify()) {
            processes[i].classify_cnts.resize(mv.get_num_species(), 0);
        }
    }

    // std::cerr << strands << " processes are created.\n";
    for (uint64_t i = 0; i < strands; i++) {
        if (finished_count == 0) {
            if (mv.movi_options->is_kmer()) {
                reset_kmer_search(processes[i], reader);
                if (!processes[i].finished) {
                    next_kmer_search(processes[i]);
                }
            } else {
                reset_process(processes[i], reader);
                if (!processes[i].finished)
                    reset_backward_search(processes[i]);
                /* if (mv.movi_options->is_zml()) {
                    processes[i].kmer_end = processes[i].pos_on_r;
                    //processes[i].match_len = 0;
                } */
            }
        } else {
            processes[i].finished = true;
        }
        if (processes[i].finished) {
            empty_strands += 1;
            finished_count += 1;
        }
    }

    if (empty_strands > 0) {
        // std::cerr << "Warning: there are fewer reads (" << strands - empty_strands << ") than the number of strands (" << strands << ").\n";
    }
    return finished_count;
}

void ReadProcessor::kmer_search_latency_hiding(uint32_t k_, BatchLoader& reader) {
    k = k_;

    std::vector<Strand> processes;
    uint64_t finished_count = initialize_strands(processes, reader);
    std::cerr << strands << " processes are initiated.\n";

    uint64_t total_bs = 0;
    // Assuming all the reads > k
    while (finished_count != strands) {
        for (uint64_t i = 0; i < strands; i++) {
            if (!processes[i].finished) {
                // 1: process next character
                bool backward_search_finished = backward_search(processes[i], processes[i].kmer_end);
                total_bs += 1;
                if (verbose)
                    std::cerr << backward_search_finished << " " << processes[i].kmer_start << " "
                              << processes[i].kmer_end << " " << processes[i].pos_on_r << "\n";
                // if ((processes[i].pos_on_r == processes[i].kmer_start and processes[i].kmer_start != 0) or (processes[i].pos_on_r == -1 and processes[i].kmer_start == 0)) {
                if (processes[i].pos_on_r == processes[i].kmer_start - 1) {
                    // 2: if the kmer is found
                    if (verbose)
                        std::cerr << processes[i].pos_on_r << " " << processes[i].kmer_start << "\n";
                    /*if (!verify_kmer(processes[i], k)) {
                        std::cerr << "kmer not verified!\n";
                        exit(0);
                    }*/
                    if (verbose)
                        std::cerr << "+ ";
                    if (processes[i].kmer_start >= 0) {
                        positive_kmer_count += 1;
                        if (processes[i].kmer_extension)
                            kmer_extension_count += 1;
                        else
                            processes[i].kmer_extension = true;
                        processes[i].kmer_end -= 1;
                        processes[i].kmer_start -= 1;
                        if (processes[i].kmer_start >= 0)
                            total_kmer_count += 1;
                        /////next_kmer_search(processes[i]);
                    } else {
                        reset_kmer_search(processes[i], reader);
                        next_kmer_search(processes[i]);
                        // 3: -- check if it was the last read in the file -> finished_count++
                        if (processes[i].finished) {
                            finished_count += 1;
                        }
                    }
                } else if (backward_search_finished and processes[i].kmer_extension) {
                    // 2: if the kmer was not found during an extension
                    if (processes[i].kmer_start >= 0)
                        total_kmer_count -= 1;
                    kmer_extension_stopped_count += 1;
                    processes[i].kmer_end += 1;
                    processes[i].kmer_start += 1;
                    /////////processes[i].pos_on_r = processes[i].kmer_end;
                    next_kmer_search(processes[i]);
                } else if (backward_search_finished) {
                    // 2: if the kmer was not found not during an extension
                    negative_kmer_count += 1;
                    if (processes[i].kmer_start >= 0) {
                        next_kmer_search(processes[i]);
                        // next_kmer_search_negative_skip_all_heuristic(processes[i], k);
                        if (processes[i].finished) {
                            finished_count += 1;
                        }
                    } else {
                        reset_kmer_search(processes[i], reader);
                        next_kmer_search(processes[i]);
                        // 3: -- check if it was the last read in the file -> finished_count++
                        if (processes[i].finished) {
                            finished_count += 1;
                        }
                    }
                } else if (processes[i].kmer_start < 0) {
                    if (verbose)
                        std::cerr << "- ";
                    reset_kmer_search(processes[i], reader);
                    next_kmer_search(processes[i]);
                    // 3: -- check if it was the last read in the file -> finished_count++
                    if (processes[i].finished) {
                        finished_count += 1;
                    }
                } else {
                    // 4: big jump with prefetch
                    my_prefetch_r((void*)(&(mv.rlbwt[0]) + mv.get_id(processes[i].range.run_start)));
                    my_prefetch_r((void*)(&(mv.rlbwt[0]) + mv.get_id(processes[i].range.run_end)));
                }
            }
        }
    }

    std::cerr << "\n";
    std::cerr << "total_bs: " << total_bs << "\n";
    std::cerr << "positive_kmer_count: " << positive_kmer_count << "\n";
    std::cerr << "negative_kmer_count: " << negative_kmer_count << "\n";
    std::cerr << "total_kmer_count: " << total_kmer_count << "\n";
    std::cerr << "kmer_extension_stopped_count: " << kmer_extension_stopped_count << "\n";
    std::cerr << "kmer_extension_count: " << kmer_extension_count << "\n";
    std::cerr << "negative_kmer_extension_count: " << negative_kmer_extension_count << "\n\n";

    kseq_destroy(seq); // STEP 5: destroy seq
    std::cerr << "kseq destroyed!\n";
    gzclose(fp); // STEP 6: close the file handler
    std::cerr << "fp file closed!\n";
}

void ReadProcessor::reset_backward_search(Strand& process) {
    // This function initialize the starting range for performing backward search
    // The initialization is based on the character at pos_on_r
    // If that character does not exist in the alphabet, it moves pos_on_r until finding a character that exists
    // If no character of the alphabet exists in the read, it makes the range empty
    process.match_count = 0;
    process.match_len = 0;
    std::string& query_seq = process.mq.query();
    if (mv.movi_options->is_count() and !mv.check_alphabet(query_seq[process.pos_on_r])) {
        process.range.make_empty();
        process.range_prev.make_empty();
        return;
    }
    while (!mv.check_alphabet(query_seq[process.pos_on_r]) and process.pos_on_r >= 0) {
        process.mq.add_ml(0, mv.movi_options->write_stdout_enabled());
        if (mv.movi_options->is_get_sa_entries()) {
            uint64_t sa_entry = mv.get_SA_entries(process.idx, process.offset);
            process.mq.add_sa_entries(sa_entry);
        }
        process.pos_on_r -= 1;
    }
    if (process.pos_on_r < 0) {
        // Special case where no character in the read exists in the index.
        process.range.make_empty();
        process.range_prev.make_empty();
        return;
    }
    process.range = mv.initialize_backward_search(process.mq, process.pos_on_r, process.match_len);
    process.kmer_end = process.pos_on_r;
    // Very expensive operation: process.match_count = process.range.count(mv.rlbwt);
}

void ReadProcessor::reset_kmer_search(Strand& process, BatchLoader& reader) {
    process.finished = next_read(process, reader);
    if (!process.finished) {
        process.length_processed = 0;
        process.kmer_start = process.mq.query().length() - k + 1;
        process.kmer_end = process.mq.query().length() - 1 + 1;
        process.pos_on_r = process.kmer_end;

        // Find the first position where the character is legal
        while (!mv.check_alphabet(process.mq.query()[process.pos_on_r])) {
            process.pos_on_r -= 1;
        }
    }
}

void ReadProcessor::next_kmer_search_negative_skip_all_heuristic(Strand& process, BatchLoader& reader) {
    process.kmer_extension = false;
    std::string& R = process.mq.query();
    process.kmer_start = process.pos_on_r - k;
    if (process.kmer_start >= 0) {
        total_kmer_count += (process.kmer_end - process.pos_on_r + 1);
        negative_kmer_count += (process.kmer_end - process.pos_on_r);
        negative_kmer_extension_count += (process.kmer_end - process.pos_on_r);
        process.kmer_end = process.pos_on_r - 1;
        process.pos_on_r = process.kmer_end;
        process.match_len = 0;
        process.range.run_start = mv.first_runs[mv.alphamap[R[process.pos_on_r]] + 1];
        process.range.offset_start = mv.first_offsets[mv.alphamap[R[process.pos_on_r]] + 1];
        process.range.run_end = mv.last_runs[mv.alphamap[R[process.pos_on_r]] + 1];
        process.range.offset_end = mv.last_offsets[mv.alphamap[R[process.pos_on_r]] + 1];
    } else {
        total_kmer_count += (process.kmer_end - k + 1);
        negative_kmer_count += (process.kmer_end - k + 1);
        negative_kmer_extension_count += (process.kmer_end - k + 1);
        reset_kmer_search(process, reader);
        next_kmer_search(process);
    }
}

void ReadProcessor::next_kmer_search(Strand& process) {
    process.kmer_extension = false;
    std::string& R = process.mq.query();
    process.kmer_start -=1;
    if (process.kmer_start >= 0)
        total_kmer_count += 1;
    process.kmer_end -=1;
    process.pos_on_r = process.kmer_end;
    process.match_len = 0;
    process.range.run_start = mv.first_runs[mv.alphamap[R[process.pos_on_r]] + 1];
    process.range.offset_start = mv.first_offsets[mv.alphamap[R[process.pos_on_r]] + 1];
    process.range.run_end = mv.last_runs[mv.alphamap[R[process.pos_on_r]] + 1];
    process.range.offset_end = mv.last_offsets[mv.alphamap[R[process.pos_on_r]] + 1];
}

bool ReadProcessor::verify_kmer(Strand& process, uint64_t k) {
    std::string& R = process.mq.query();
    std::string kmer = R.substr (process.kmer_start, k);
    int32_t pos_on_r = k - 1;
    MoveInterval interval(
        mv.first_runs[mv.alphamap[kmer[pos_on_r]] + 1],
        mv.first_offsets[mv.alphamap[kmer[pos_on_r]] + 1],
        mv.last_runs[mv.alphamap[kmer[pos_on_r]] + 1],
        mv.last_offsets[mv.alphamap[kmer[pos_on_r]] + 1]
    );
    auto match_count = mv.backward_search(kmer, pos_on_r, interval, std::numeric_limits<int32_t>::max()).count(mv.rlbwt);
    if (pos_on_r == 0 and match_count > 0) {
        return true;
    } else {
        return false;
    }
}

bool ReadProcessor::backward_search(Strand& process, uint64_t end_pos) {
    // This function performs one step backward search, suitable for queries with latency hiding
    // If it is the first step of the current search, it has to update the interval first,
    // Otherwise, it starts from the LF steps and continues until the next LF steps
    // The functions might stops before or after the LF step
    // The backward search is empty in two conditions:
    // 1) It is the count query and the last character on the read does not exist in the alphabet
    // 2) No character on the read exists in the alphabet
    // The first_iteration condition should be only true after the reset_backward_search function
    std::string& R = process.mq.query();
    if (verbose)
        std::cerr << "backward search begins:\n" << process.pos_on_r << " "
                    << R[process.pos_on_r] << " "
                    << mv.alphabet[mv.rlbwt[process.range.run_start].get_c()] << " "
                    << mv.alphabet[mv.rlbwt[process.range.run_end].get_c()] << "\n";
    bool first_iteration = process.pos_on_r == end_pos;
    if (first_iteration) {
        if (process.range.is_empty()) {
            process.pos_on_r += 1; // Why is this line necessary?
            return true;
        }

        if (!mv.check_alphabet(R[process.pos_on_r - 1])) {
            return true;
        }

        process.range_prev = process.range;
        mv.update_interval(process.range, R[process.pos_on_r - 1]);
    }

    if (!process.range.is_empty()) {
        mv.LF_move(process.range.offset_start, process.range.run_start);
        mv.LF_move(process.range.offset_end, process.range.run_end);
        // The last update_interval has been based on pos_on_r - 1
        // the range cannot be empty after the LF step if it was not empty before the LF steps
        // Therefore, we can safely say that the read at position pos_on_r - 1 is matched now.
        process.pos_on_r -= 1;
        if (mv.movi_options->is_zml()) {
            process.mq.add_ml(process.match_len, mv.movi_options->write_stdout_enabled());
            process.match_len += 1;
        }
        // To make the pos_on_r match the range currently represented after the two LF steps.
    } else {
        // Since the range has become empty after the update, we should count the previous position
        // Very expensive operation: process.match_count = process.range_prev.count(mv.rlbwt);
        return true;
    }

    if (process.pos_on_r <= 0) {
        // If we are at position 0 after the last LF step, we can finish the procedure
        // Very expensive operation: process.match_count = process.range.count(mv.rlbwt);
        if (process.pos_on_r < 0) {
            throw std::runtime_error(ERROR_MSG("[Read Processor - backward search] This should never happen pos_on_r: " +
                                               std::to_string(process.pos_on_r) + "\n"));
        }
        return true;
    }

    // The following steps is the preparations for the next LF steps
    if (!mv.check_alphabet(R[process.pos_on_r - 1])) {
        // If the next character does not belong, we can stop the backward search.
        // Very expensive operation: process.match_count = process.range.count(mv.rlbwt);
        return true;
    }
    // Store the current range as range_prev in case the range becomes empty after the update
    // If the range becomes empty, it means that the match cannot be extended.
    process.range_prev = process.range;
    if (verbose)
        std::cerr << "before: " << process.range.run_start << " " << process.range.run_end << " "
                    << static_cast<uint64_t>(mv.rlbwt[process.range.run_start].get_c()) << " "
                    << mv.alphabet[mv.rlbwt[process.range.run_end].get_c()] << "\n";
    mv.update_interval(process.range, R[process.pos_on_r - 1]);
    if (verbose)
        std::cerr << "after: " << process.range.run_start << " " << process.range.run_end << " "
                    << static_cast<uint64_t>(mv.rlbwt[process.range.run_start].get_c()) << " "
                    << mv.alphabet[mv.rlbwt[process.range.run_start].get_c()] << " "
                    << mv.alphabet[mv.rlbwt[process.range.run_end].get_c()] << "\n";

    return false;
}
