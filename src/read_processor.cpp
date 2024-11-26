#include "read_processor.hpp"

uint32_t alphamap_3_[4][4] = {{3, 0, 1, 2},
                             {0, 3, 1, 2},
                             {0, 1, 3, 2},
                             {0, 1, 2, 3}};

ReadProcessor::ReadProcessor(std::string reads_file_name, MoveStructure& mv_, int strands_ = 4, bool verbose_ = false, bool reverse_ = false) : mv(mv_) {
    verbose = verbose_;
    reverse = mv_.movi_options->is_reverse();
    // Solution for handling the stdin input: https://biowize.wordpress.com/2013/03/05/using-kseq-h-with-stdin/
    if (reads_file_name == "-") {
        FILE *instream = stdin;
        fp = gzdopen(fileno(instream), "r"); // STEP 2: open the file handler
    } else {
        fp = gzopen(reads_file_name.c_str(), "r"); // STEP 2: open the file handler
    }

    seq = kseq_init(fp); // STEP 3: initialize seq
    std::string index_type = mv_.index_type();
    if (!mv_.movi_options->is_stdout()) {
        if (mv_.movi_options->is_pml()) {
            std::string mls_file_name = reverse ? reads_file_name + "." + index_type + ".reverse.pml.bin" :
                                                reads_file_name + "." + index_type + ".pml.bin";
            mls_file = std::ofstream(mls_file_name, std::ios::out | std::ios::binary);
        } else if (mv_.movi_options->is_zml()) {
            std::string mls_file_name = reverse ? reads_file_name + "." + index_type + ".reverse.zml.bin" :
                                                reads_file_name + "." + index_type + ".zml.bin";
            mls_file = std::ofstream(mls_file_name, std::ios::out | std::ios::binary);
        } else if (mv_.movi_options->is_count()) {
            std::string matches_file_name = reads_file_name + "." + index_type + ".matches";
            matches_file = std::ofstream(matches_file_name);
        }
    }
    total_kmer_count = 0;
    positive_kmer_count = 0;
    negative_kmer_count = 0;
    kmer_extension_count = 0;
    kmer_extension_stopped_count = 0;
    negative_kmer_extension_count = 0;
    read_processed = 0;
    strands = strands_;
    if (mv_.movi_options->is_logs()) {
        costs_file = std::ofstream(reads_file_name + "." + index_type + ".costs");
        scans_file = std::ofstream(reads_file_name + "." + index_type + ".scans");
        fastforwards_file = std::ofstream(reads_file_name + "." + index_type + ".fastforwards");
    }
}

bool ReadProcessor::next_read(Strand& process) {
    if (read_processed % 1000 == 0)
        std::cerr << read_processed << "\r";
    read_processed += 1;
    l = kseq_read(seq); // STEP 4: read sequence
    if (l >= 0) {
        process.st_length = seq->name.m;
        process.read_name = seq->name.s;
        process.read = seq->seq.s;
        if (reverse) {
            std::reverse(process.read.begin(), process.read.end());
        }
        process.match_len = 0;
        process.ff_count = 0;
        process.scan_count = 0;
        process.mq = MoveQuery(process.read);
        return false;
    } else {
        return true;
    }
}

void ReadProcessor::reset_process(Strand& process) {
    process.finished = next_read(process);
    if (!process.finished) {
        process.length_processed = 0;
        process.pos_on_r = process.read.length() - 1;
        process.idx = mv.r - 1;
        process.offset = mv.get_n(process.idx) - 1;
    }
}

// This function is for computing PMLs
void ReadProcessor::process_char(Strand& process) {
    if (process.pos_on_r < process.read.length() - 1) {
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
        // Jumping up or down (randomly or with thresholds)
        uint64_t idx_before_jump = process.idx;
        bool up = false;
#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4 or MODE == 6
        up = mv.jump_thresholds(process.idx, process.offset, R[process.pos_on_r], process.scan_count);
#endif
#if MODE == 3
        up = mv.jump_randomly(process.idx, R[process.pos_on_r], process.scan_count);
#endif
        process.match_len = 0;
        char c = mv.alphabet[mv.rlbwt[process.idx].get_c()];
        // sanity check
        if (c == R[process.pos_on_r]) {
            // Observing a match after the jump
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
    process.mq.add_ml(process.match_len);
    process.pos_on_r -= 1;
    // if (mv.logs)
    //     process.t2 = std::chrono::high_resolution_clock::now();
    // LF step should happen here in the non-prefetch code
    // uint64_t ff_count = mv.LF_move(process.offset, process.idx);
    // process.ff_count_tot += ff_count;
}

void ReadProcessor::write_mls(Strand& process) {
    bool write_stdout = mv.movi_options->is_stdout();
    bool logs = mv.movi_options->is_logs();
    if (logs) {
        auto t2 = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
        process.mq.add_cost(elapsed);
        process.mq.add_fastforward(process.ff_count);
        process.mq.add_scan(process.scan_count);
    }
    if (write_stdout) {
        std::cout << ">" << process.read_name << " \n";
        auto& matching_lens = process.mq.get_matching_lengths();
        uint64_t mq_lens_size = matching_lens.size();
        for (int64_t i = mq_lens_size - 1; i >= 0; i--) {
            std::cout << matching_lens[i] << " ";
        }
        std::cout << "\n";
    } else {
        mls_file.write(reinterpret_cast<char*>(&process.st_length), sizeof(process.st_length));
        mls_file.write(reinterpret_cast<char*>(&process.read_name[0]), process.st_length);
        auto& ml_lens = process.mq.get_matching_lengths();
        uint64_t mq_lens_size = ml_lens.size();
        mls_file.write(reinterpret_cast<char*>(&mq_lens_size), sizeof(mq_lens_size));
        mls_file.write(reinterpret_cast<char*>(&ml_lens[0]), mq_lens_size * sizeof(ml_lens[0]));
    }

    if (logs) {
        costs_file << ">" << process.read_name << "\n";
        scans_file << ">" << process.read_name << "\n";
        fastforwards_file << ">" << process.read_name << "\n";
        for (auto& cost : process.mq.get_costs()) {
            costs_file << cost.count() << " ";
        }
        for (auto& scan: process.mq.get_scans()) {
            scans_file << scan << " ";
        }
        for (auto& fast_forward : process.mq.get_fastforwards()) {
            fastforwards_file << fast_forward << " ";
        }
        costs_file << "\n";
        scans_file << "\n";
        fastforwards_file << "\n";
    }
}

void ReadProcessor::write_count(Strand& process) {
    compute_match_count(process);
    auto& R = process.mq.query();
    bool write_stdout = mv.movi_options->is_stdout();
    if (write_stdout) {
        std::cout << process.read_name << "\t";
        std::cout << R.length() - process.pos_on_r << "/" << R.length() << "\t" << process.match_count << "\n";
    } else {
        matches_file << process.read_name << "\t";
        matches_file << R.length() - process.pos_on_r << "/" << R.length() << "\t" << process.match_count << "\n";
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
        std::cerr << "This should not happen, for match count computation.\n";
        exit(0);
    }
}

void ReadProcessor::process_latency_hiding() {
    bool is_pml = mv.movi_options->is_pml();
    bool is_zml = mv.movi_options->is_zml();
    bool is_count = mv.movi_options->is_count();

    if ((is_pml and is_count) or (is_pml and is_zml) or (is_count and is_zml)) {
        std::cerr << "Error parsing the input, multipe types of queries cannot be processed at the same time.\n";
        exit(0);
    }

    std::vector<Strand> processes;
    uint64_t finished_count = initialize_strands(processes);
    std::cerr << strands << " processes are initiated.\n";

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
                    if (is_pml) {
                        write_mls(processes[i]);
                        reset_process(processes[i]);
                    } else if (is_count) {
                        write_count(processes[i]);
                        reset_process(processes[i]);
                        reset_backward_search(processes[i]);
                    } else if (is_zml) {
                        if (backward_search_finished and processes[i].pos_on_r > 0) {
                            processes[i].mq.add_ml(processes[i].match_len);
                            processes[i].pos_on_r -= 1;
                            reset_backward_search(processes[i]);
                            processes[i].kmer_end = processes[i].pos_on_r;
                            continue;
                        } else if (processes[i].pos_on_r <= 0) {
                            processes[i].mq.add_ml(processes[i].match_len);
                            write_mls(processes[i]);
                            reset_process(processes[i]);
                            reset_backward_search(processes[i]);
                            processes[i].kmer_end = processes[i].pos_on_r;
                        }
                    }
                    // 3: -- check if it was the last read in the file -> finished_count++
                    if (processes[i].finished) {
                        finished_count += 1;
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

    if (!mv.movi_options->is_stdout()) {
        if (is_pml or is_zml) {
            mls_file.close();
            std::cerr << "Matching lengths file is closed!\n";
        } else if (is_count) {
            matches_file.close();
            std::cerr << "match_file file closed!\n";
        }
    }

    std::cerr << "no_ftab: " << mv.no_ftab << "\n";
    std::cerr << "all_initializations: " << mv.all_initializations << "\n";

    kseq_destroy(seq); // STEP 5: destroy seq
    std::cerr << "kseq destroyed!\n";
    gzclose(fp); // STEP 6: close the file handler
    std::cerr << "fp file closed!\n";

    if (mv.movi_options->is_logs()) {
        costs_file.close();
        scans_file.close();
        fastforwards_file.close();
    }
}

uint64_t ReadProcessor::initialize_strands(std::vector<Strand>& processes) {
    uint64_t finished_count = 0;
    for(int i = 0; i < strands; i++) processes.emplace_back(Strand());
    std::cerr << strands << " processes are created.\n";
    for (uint64_t i = 0; i < strands; i++) {
        if (finished_count == 0) {
            if (mv.movi_options->is_kmer()) {
                reset_kmer_search(processes[i]);
                if (!processes[i].finished) {
                    next_kmer_search(processes[i]);
                }
            } else {
                reset_process(processes[i]);
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
            std::cerr << "Warning: less than strands = " << strands << " reads.\n";
            finished_count += 1;
        }
    }
    return finished_count;
}

void ReadProcessor::kmer_search_latency_hiding(uint32_t k_) {
    k = k_;

    std::vector<Strand> processes;
    uint64_t finished_count = initialize_strands(processes);
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
                        reset_kmer_search(processes[i]);
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
                        reset_kmer_search(processes[i]);
                        next_kmer_search(processes[i]);
                        // 3: -- check if it was the last read in the file -> finished_count++
                        if (processes[i].finished) {
                            finished_count += 1;
                        }
                    }
                } else if (processes[i].kmer_start < 0) {
                    if (verbose)
                        std::cerr << "- ";
                    reset_kmer_search(processes[i]);
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
        process.mq.add_ml(0);
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

void ReadProcessor::reset_kmer_search(Strand& process) {
    process.finished = next_read(process);
    if (!process.finished) {
        process.length_processed = 0;
        process.kmer_start = process.read.length() - k + 1;
        process.kmer_end = process.read.length() - 1 + 1;
        process.pos_on_r = process.kmer_end;
    }
}

void ReadProcessor::next_kmer_search_negative_skip_all_heuristic(Strand& process) {
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
        reset_kmer_search(process);
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
            process.mq.add_ml(process.match_len);
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
            std::cerr << "This should never happen.\n";
            exit(0);
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