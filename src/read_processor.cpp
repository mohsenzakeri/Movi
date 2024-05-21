#include "read_processor.hpp"

uint32_t alphamap_3_[4][4] = {{3, 0, 1, 2},
                             {0, 3, 1, 2},
                             {0, 1, 3, 2},
                             {0, 1, 2, 3}};

ReadProcessor::ReadProcessor(std::string reads_file_name, MoveStructure& mv_, int strands_ = 4, bool query_pml = true, bool reverse_ = false) {
    // Solution for handling the stdin input: https://biowize.wordpress.com/2013/03/05/using-kseq-h-with-stdin/
    if (reads_file_name == "-") {
        FILE *instream = stdin;
        fp = gzdopen(fileno(instream), "r"); // STEP 2: open the file handler
    } else {
        fp = gzopen(reads_file_name.c_str(), "r"); // STEP 2: open the file handler
    }

    seq = kseq_init(fp); // STEP 3: initialize seq
    std::string index_type = mv_.index_type();
    if (query_pml) {
        std::string pmls_file_name = reads_file_name + "." + index_type + ".mpml.bin";
        pmls_file = std::ofstream(pmls_file_name, std::ios::out | std::ios::binary);
    } else {
        std::string matches_file_name = reads_file_name + "." + index_type + ".matches";
        matches_file = std::ofstream(matches_file_name);
    }

    read_processed = 0;
    strands = strands_;
    if (mv_.logs) {
        costs_file = std::ofstream(reads_file_name + "." + index_type + ".costs");
        scans_file = std::ofstream(reads_file_name + "." + index_type + ".scans");
        fastforwards_file = std::ofstream(reads_file_name + "." + index_type + ".fastforwards");
    }
    reverse = reverse_;
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
        if (reverse)
            std::reverse(process.read.begin(), process.read.end());
        process.match_len = 0;
        process.ff_count = 0;
        process.scan_count = 0;
        process.mq = MoveQuery(process.read);;
        return false;
    } else {
        return true;
    }
}

void ReadProcessor::reset_process(Strand& process, MoveStructure& mv) {
    process.finished = next_read(process);
    if (!process.finished) {
        process.length_processed = 0;
        process.pos_on_r = process.read.length() - 1;
        process.idx = mv.r - 1;
        process.offset = mv.get_n(process.idx) - 1;
    }
}

void ReadProcessor::process_char(Strand& process, MoveStructure& mv) {
    if (process.pos_on_r < process.read.length() - 1) {
        // LF step
        // if (mv.logs)
        //     process.t3 = std::chrono::high_resolution_clock::now();
        process.ff_count = mv.LF_move(process.offset, process.idx);
        if (mv.logs) {
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
    if (mv.logs) {
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
        bool up = mv.jump_thresholds(process.idx, process.offset, R[process.pos_on_r], process.scan_count);
        process.match_len = 0;
        char c = mv.alphabet[mv.rlbwt[process.idx].get_c_mm()];
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
    process.mq.add_pml(process.match_len);
    process.pos_on_r -= 1;
    // if (mv.logs)
    //     process.t2 = std::chrono::high_resolution_clock::now();
    // LF step should happen here in the non-prefetch code
    // uint64_t ff_count = mv.LF_move(process.offset, process.idx);
    // process.ff_count_tot += ff_count;
}

void ReadProcessor::write_pmls(Strand& process, bool logs, bool write_stdout) {
    if (logs) {
        auto t2 = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
        process.mq.add_cost(elapsed);
        process.mq.add_fastforward(process.ff_count);
        process.mq.add_scan(process.scan_count);
    }
    if (write_stdout) {
        std::cout << ">" << process.read_name << " \n";
        auto& pml_lens = process.mq.get_pml_lens();
        uint64_t mq_pml_lens_size = pml_lens.size();
        for (int64_t i = mq_pml_lens_size - 1; i >= 0; i--) {
            std::cout << pml_lens[i] << " ";
        }
        std::cout << "\n";
    } else {
        pmls_file.write(reinterpret_cast<char*>(&process.st_length), sizeof(process.st_length));
        pmls_file.write(reinterpret_cast<char*>(&process.read_name[0]), process.st_length);
        auto& pml_lens = process.mq.get_pml_lens();
        uint64_t mq_pml_lens_size = pml_lens.size();
        pmls_file.write(reinterpret_cast<char*>(&mq_pml_lens_size), sizeof(mq_pml_lens_size));
        pmls_file.write(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));
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

void ReadProcessor::process_latency_hiding(MoveStructure& mv) {
    std::vector<Strand> processes;
    for(int i = 0; i < strands; i++) processes.emplace_back(Strand());
    std::cerr << strands << " processes are created.\n";
    uint64_t fnished_count = 0;
    for (uint64_t i = 0; i < strands; i++) {
        if (fnished_count == 0) {
            reset_process(processes[i], mv);
            reset_backward_search(processes[i], mv);
        } else {
            processes[i].finished = true;
        }
        if (processes[i].finished) {
            std::cerr << "Warning: less than strands = " << strands << " reads.\n";
            fnished_count += 1;
        }
    }
    std::cerr << strands << " processes are initiated.\n";

    while (fnished_count != strands) {
        for (uint64_t i = 0; i < strands; i++) {
            if (!processes[i].finished) {
                // 1: process next character -- doing fast forward
                process_char(processes[i], mv);
                // 2: if the read is done -> Write the pmls and go to next read
                if (processes[i].pos_on_r <= -1) {
                    write_pmls(processes[i], mv.logs, mv.movi_options->is_stdout());
                    reset_process(processes[i], mv);
                    // 3: -- check if it was the last read in the file -> fnished_count++
                    if (processes[i].finished) {
                        fnished_count += 1;
                    }
                } else {
                    // 4: big jump with prefetch
                    my_prefetch_r((void*)(&(mv.rlbwt[0]) + mv.rlbwt[processes[i].idx].get_id()));
                }
            }
        }
    }

    pmls_file.close();
    std::cerr << "pmls file closed!\n";
    kseq_destroy(seq); // STEP 5: destroy seq
    std::cerr << "kseq destroyed!\n";
    gzclose(fp); // STEP 6: close the file handler
    std::cerr << "fp file closed!\n";

    if (mv.logs) {
        costs_file.close();
        scans_file.close();
        fastforwards_file.close();
    }
}

void ReadProcessor::backward_search_latency_hiding(MoveStructure& mv) {
    std::vector<Strand> processes;
    for(int i = 0; i < strands; i++) processes.emplace_back(Strand());
    std::cerr << strands << " processes are created.\n";

    uint64_t fnished_count = 0;
    for (uint64_t i = 0; i < strands; i++) {
        if (fnished_count == 0) {
            reset_process(processes[i], mv);
            reset_backward_search(processes[i], mv);
        } else {
            processes[i].finished = true;
        }
        if (processes[i].finished) {
            std::cerr << "Warning: less than strands = " << strands << " reads.\n";
            fnished_count += 1;
        }
    }
    std::cerr << strands << " processes are initiated.\n";

    while (fnished_count != strands) {
        for (uint64_t i = 0; i < strands; i++) {
            if (!processes[i].finished) {
                // 1: process next character -- doing fast forward
                uint64_t match_count = 0;
                bool backward_search_finished = backward_search(processes[i], mv, match_count);
                // 2: if the read is done -> Write the pmls and go to next read
                if (backward_search_finished) {
                    auto& R = processes[i].mq.query();
                    // matches_file << processes[i].read_name << (processes[i].pos_on_r == 0 ? "\tFound\t" : "\tNot-Found\t");
                    matches_file << processes[i].read_name << "\t";
                    if (processes[i].pos_on_r != 0) processes[i].pos_on_r += 1;
                    matches_file << R.length() - processes[i].pos_on_r << "/" << R.length() << "\t" << match_count << "\n";

                    reset_process(processes[i], mv);
                    reset_backward_search(processes[i], mv);
                    // 3: -- check if it was the last read in the file -> fnished_count++
                    if (processes[i].finished) {
                        fnished_count += 1;
                    }
                } else {
                    // 4: big jump with prefetch
                    my_prefetch_r((void*)(&(mv.rlbwt[0]) + mv.rlbwt[processes[i].range.run_start].get_id()));
                    my_prefetch_r((void*)(&(mv.rlbwt[0]) + mv.rlbwt[processes[i].range.run_end].get_id()));
                }
            }
        }
    }

    std::cerr << "pmls file closed!\n";
    kseq_destroy(seq); // STEP 5: destroy seq
    std::cerr << "kseq destroyed!\n";
    gzclose(fp); // STEP 6: close the file handler
    std::cerr << "fp file closed!\n";
}

void ReadProcessor::reset_backward_search(Strand& process, MoveStructure& mv) {
    std::string& R = process.mq.query();
    process.range.run_start = mv.first_runs[mv.alphamap[R[process.pos_on_r]] + 1];
    process.range.offset_start = mv.first_offsets[mv.alphamap[R[process.pos_on_r]] + 1];
    process.range.run_end = mv.last_runs[mv.alphamap[R[process.pos_on_r]] + 1];
    process.range.offset_end = mv.last_offsets[mv.alphamap[R[process.pos_on_r]] + 1];
}

bool ReadProcessor::backward_search(Strand& process, MoveStructure& mv, uint64_t& match_count) {
    std::string& R = process.mq.query();

    if (process.pos_on_r < R.length() - 1) {
        if (((process.range.run_start < process.range.run_end) or
            (process.range.run_start == process.range.run_end and process.range.offset_start <= process.range.offset_end)) and
            (mv.alphabet[mv.rlbwt[process.range.run_start].get_c()] == R[process.pos_on_r]) and
            (mv.alphabet[mv.rlbwt[process.range.run_end].get_c()] == R[process.pos_on_r])) {
            mv.LF_move(process.range.offset_start, process.range.run_start);
            mv.LF_move(process.range.offset_end, process.range.run_end);
            if (process.pos_on_r == 0) {
                if (process.range.run_start == process.range.run_end) {
                    match_count = process.range.offset_end - process.range.offset_start + 1;
                } else {
                    match_count = (mv.rlbwt[process.range.run_start].get_n() - process.range.offset_start) + (process.range.offset_end + 1);
                    for (uint64_t k = process.range.run_start + 1; k < process.range.run_end; k ++) {
                        match_count += mv.rlbwt[k].get_n();
                    }
                }
                return true;
            }
        } else {
            // The read was not found.
            if (process.range_prev.run_start == process.range_prev.run_end) {
                match_count = process.range_prev.offset_end - process.range_prev.offset_start + 1;
            } else {
                match_count = (mv.rlbwt[process.range_prev.run_start].get_n() - process.range_prev.offset_start) +
                                (process.range_prev.offset_end + 1);
                for (uint64_t k = process.range_prev.run_start + 1; k < process.range_prev.run_end; k ++) {
                    match_count += mv.rlbwt[k].get_n();
                }
            }
            return true;
        }
    }

    process.range_prev = process.range;
    process.pos_on_r -= 1;
    if (process.range.run_start == mv.end_bwt_idx or process.range.run_end == mv.end_bwt_idx or !mv.check_alphabet(R[process.pos_on_r])) {
        // The read was not found.
        if (process.range_prev.run_start == process.range_prev.run_end) {
            match_count = process.range_prev.offset_end - process.range_prev.offset_start + 1;
        } else {
            match_count = (mv.rlbwt[process.range_prev.run_start].get_n() - process.range_prev.offset_start) +
                            (process.range_prev.offset_end + 1);
            for (uint64_t k = process.range_prev.run_start + 1; k < process.range_prev.run_end; k ++) {
                match_count += mv.rlbwt[k].get_n();
            }
        }
        return true;
    }

#if MODE == 0
    while ((process.range.run_start < process.range.run_end) and (mv.alphabet[mv.rlbwt[process.range.run_start].get_c()] != R[process.pos_on_r])) {
        process.range.run_start += 1;
        process.range.offset_start = 0;
        if (process.range.run_start >= mv.r) {
            break;
        }
    }
    while ((process.range.run_end > process.range.run_start) and (mv.alphabet[mv.rlbwt[process.range.run_end].get_c()] != R[process.pos_on_r])) {
        process.range.run_end -= 1;
        process.range.offset_end = mv.rlbwt[process.range.run_end].get_n() - 1;
        if (process.range.run_end == 0) {
            break;
        }
    }
#endif
#if MODE == 1
    uint64_t read_alphabet_index = mv.alphamap[static_cast<uint64_t>(R[process.pos_on_r])];
    if ((process.range.run_start < process.range.run_end) and (mv.alphabet[mv.rlbwt[process.range.run_start].get_c()] != R[process.pos_on_r])) {
        if (process.range.run_start == 0) {
            while ((process.range.run_start < process.range.run_end) and (mv.alphabet[mv.rlbwt[process.range.run_start].get_c()] != R[process.pos_on_r])) {
                process.range.run_start += 1;
                process.range.offset_start = 0;
                if (process.range.run_start >= mv.r) {
                    break;
                }
            }
        } else {
            char rlbwt_char = mv.alphabet[mv.rlbwt[process.range.run_start].get_c()];
            uint64_t alphabet_index = alphamap_3_[mv.alphamap[rlbwt_char]][read_alphabet_index];
            if (mv.rlbwt[process.range.run_start].get_next_down(alphabet_index) == std::numeric_limits<uint16_t>::max()) {
                process.range.run_start = mv.r;
            } else {
                uint64_t run_start = process.range.run_start + mv.rlbwt[process.range.run_start].get_next_down(alphabet_index);
                if (run_start <= process.range.run_end) {
                    process.range.run_start = run_start;
                } else {
                    process.range.run_start = process.range.run_end;
                }
                process.range.offset_start = 0;
            }
        }
    }
    if ((process.range.run_end > process.range.run_start) and (mv.alphabet[mv.rlbwt[process.range.run_end].get_c()] != R[process.pos_on_r])) {
        char rlbwt_char = mv.alphabet[mv.rlbwt[process.range.run_end].get_c()];
        uint64_t alphabet_index = alphamap_3_[mv.alphamap[rlbwt_char]][read_alphabet_index];
        if (mv.rlbwt[process.range.run_end].get_next_up(alphabet_index) == std::numeric_limits<uint16_t>::max()) {
            process.range.run_end = mv.r;
        } else {
            uint64_t run_end = process.range.run_end - mv.rlbwt[process.range.run_end].get_next_up(alphabet_index);
            if (run_end >= process.range.run_start) {
                process.range.run_end = run_end;
            } else {
                process.range.run_end = process.range.run_start;
            }
            process.range.offset_end = mv.rlbwt[process.range.run_end].get_n() - 1;
        }
    }
#endif

    return false;
}