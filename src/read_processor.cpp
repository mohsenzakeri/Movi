#include "read_processor.hpp"

ReadProcessor::ReadProcessor(char* reads_file_name, MoveStructure& mv_, int strands_ = 4) {
    fp = gzopen(reads_file_name, "r"); // STEP 2: open the file handler
    seq = kseq_init(fp); // STEP 3: initialize seq
    // mv = mv_;
    std::string index_type = mv_.index_type();
    std::string pmls_file_name = static_cast<std::string>(reads_file_name) + "." + index_type + ".mpml.bin";
    pmls_file = std::ofstream(pmls_file_name, std::ios::out | std::ios::binary);
    read_processed = 0;
    strands = strands_;
}

bool ReadProcessor::next_read(Strand& process) {
    if (read_processed % 100 == 0)
        std::cerr << read_processed << "\r";
    read_processed += 1;
    l = kseq_read(seq); // STEP 4: read sequence
    if (l >= 0) {
        process.st_length = seq->name.m;
        process.read_name = seq->name.s;
        process.read = seq->seq.s;
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
        uint64_t ff_count = mv.LF_move(process.offset, process.idx);
        process.ff_count_tot += ff_count;
    }

    auto& row = mv.rlbwt[process.idx];
    uint64_t row_idx = process.idx;
    char row_c = mv.alphabet[row.get_c()];
    std::string& R = process.mq.query();
    if (mv.alphamap[static_cast<uint64_t>(R[process.pos_on_r])] == mv.alphamap.size()) {
        process.match_len = 0;
    } else if (row_c == R[process.pos_on_r]) {
        // Case 1
        process.match_len += 1;
    } else {
        // Case 2
        // Jumping up or down (randomly or with thresholds)
        uint64_t idx_before_jump = process.idx;
        uint64_t scan_count = 0; // not used
        bool up = mv.jump_thresholds(process.idx, process.offset, R[process.pos_on_r], scan_count);
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
                std::cerr<<"\t idx: " << process.idx << " offset: " << process.offset << "\n";
        } else {
            std::cerr << "\t \t This should not happen!\n";
        }
    }
    process.mq.add_pml(process.match_len);
    process.pos_on_r -= 1;

    // LF step
    // uint64_t ff_count = mv.LF_move(process.offset, process.idx);
    // process.ff_count_tot += ff_count;
    // uint64_t row_size = mv.rlbwt[process.idx].row_size();
    // my_prefetch_r((void*)(mv.rlbwt + mv.rlbwt[process.idx].get_id()));
}

void ReadProcessor::write_pmls(Strand& process) {
    pmls_file.write(reinterpret_cast<char*>(&process.st_length), sizeof(process.st_length));
    pmls_file.write(reinterpret_cast<char*>(&process.read_name[0]), process.st_length);
    auto& pml_lens = process.mq.get_pml_lens();
    uint64_t mq_pml_lens_size = pml_lens.size();
    pmls_file.write(reinterpret_cast<char*>(&mq_pml_lens_size), sizeof(mq_pml_lens_size));
    pmls_file.write(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));
}

void ReadProcessor::process_latency_hiding(MoveStructure& mv) {
    // std::string queries[strands];
    // uint64_t length_processed[strands];
    // bool finished[strands];

    std::vector<Strand> processes;
    for(int i = 0; i < strands; i++) processes.emplace_back(Strand());
    std::cerr << strands << " processes are created.\n";
    // Strand processes[strands];
    uint64_t fnished_count = 0;
    for (uint64_t i = 0; i < strands; i++) {
        reset_process(processes[i], mv);
        if (processes[i].finished) {
            std::cerr << "Warning: less than strands = " << strands << " reads.\n";
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
                    write_pmls(processes[i]);
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
    std::cerr<<"pmls file closed!\n";
    kseq_destroy(seq); // STEP 5: destroy seq
    std::cerr<<"kseq destroyed!\n";
    gzclose(fp); // STEP 6: close the file handler
    std::cerr<<"fp file closed!\n";
}
