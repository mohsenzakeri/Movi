#include "move_structure.hpp"

uint64_t pow2[ARR_SIZE];

MoveStructure::MoveStructure(MoviOptions* movi_options_) {
    movi_options = movi_options_;
    no_ftab = 0;
    all_initializations = 0;
}

MoveStructure::MoveStructure(MoviOptions* movi_options_, uint16_t nt_splitting_, bool constant_) {
    movi_options = movi_options_;
    nt_splitting = nt_splitting_;
    constant = constant_;
    no_ftab = 0;
    all_initializations = 0;

    MoviHeader header;
    print_index_version(header);
}

std::vector<MoveRow> MoveStructure::get_rlbwt() {
    return rlbwt;
}

// TODO: The following is useful for mmaping, there is a slowdown though
// MoveRow& MoveStructure::get_move_row(uint64_t idx) {
//     if (false) { // movi_options->is_mmap()) {
//         return rlbwt_view[idx];
//     } else {
//         return rlbwt[idx];
//     }
// }

uint64_t MoveStructure::get_SA_entries(uint64_t idx, uint64_t offset) {
    uint64_t abs_offset = all_p[idx] + offset;

    // Reaching the nearest SA sample entry
    uint64_t distance = 0;
    uint64_t SA_sample_rate = movi_options->get_SA_sample_rate();
    while (abs_offset % SA_sample_rate != 0) {
        LF_move(offset, idx);
        abs_offset = all_p[idx] + offset;
        distance += 1;
    }

    return sampled_SA_entries[abs_offset / SA_sample_rate] + distance;
}

uint32_t MoveStructure::compute_index(char row_char, char lookup_char) {
    uint32_t alpha_index = alphamap[lookup_char];
    return alpha_index;
    /*if (lookup_char < row_char)
        return alpha_index;
    else
        return alpha_index+1;*/
}

uint16_t MoveStructure::LF_move(uint64_t& offset, uint64_t& i, uint64_t id) {

    auto& row = rlbwt[i];
    auto idx = (id == std::numeric_limits<uint64_t>::max()) ? get_id(i) : id;
    if (idx >= r) {
        throw std::runtime_error(ERROR_MSG("[LF_move] This should not happen.\nidx:::" + std::to_string(idx) + " i:" + std::to_string(i)));
    }
    offset = get_offset(i) + offset;
    uint16_t ff_count = 0;

    if (idx < r - 1 && offset >= get_n(idx)) {
        uint64_t idx_ = fast_forward(offset, idx, 0);
        idx += idx_;
        if (idx_ >= std::numeric_limits<uint16_t>::max()) {
            throw std::runtime_error(ERROR_MSG("[LF_move] Number of fast forwards for a query was greater than 2^16: " +
                                     std::to_string(idx_) + "\noffset: " + std::to_string(offset) + "\nidx: " + std::to_string(idx)));
        }
        ff_count = static_cast<uint16_t>(idx_);
    }

    if (movi_options->is_logs()) {
        if (ff_counts.find(ff_count) != ff_counts.end())
            ff_counts[ff_count] += 1;
        else
            ff_counts[ff_count] = 1;
    }
    i = idx;
    return ff_count;
}

#define my_prefetch_rr(address) __builtin_prefetch((void *)address, 0, 1)

uint64_t MoveStructure::get_id(uint64_t idx) {
#if NO_EXTRA_TABLE
    return rlbwt[idx].get_id();
#endif
#if BLOCKED_MODES
    if (idx != end_bwt_idx) {
        uint64_t block_number = idx / block_size;
        return rlbwt[idx].get_id() + static_cast<uint64_t>(id_blocks[rlbwt[idx].get_c()][block_number]) + first_runs[rlbwt[idx].get_c() + 1];
    }
    else
        return rlbwt[idx].get_id();
#endif
#if TALLY_MODES
    uint64_t id = 0;
    // The $ always goes to row/run 0
    if (idx == end_bwt_idx) {
        return 0;
    }

    uint8_t char_index = rlbwt[idx].get_c();
    uint64_t tally_a = idx / tally_checkpoints;

    // The id of the last run is always stored at the last tally_id row
    if (idx == r - 1) {
        uint64_t tally_ids_len = tally_ids[char_index].size();
        return tally_ids[char_index][tally_ids_len - 1].get();
    }


    // If we are at a checkpoint, simply return the stored id
    if (idx % tally_checkpoints == 0) {
        id = tally_ids[char_index][tally_a].get();
        return tally_ids[char_index][tally_a].get();
    }

    // find the closest check point
    uint64_t tally_b = tally_a + 1;

    // Sanity Check to see if we fall off the tally vector
    // if (tally_b >= tally_ids[char_index].size()) {
    //     dbg << "ta: " << tally_a << " tb: " << tally_b <<
    //            " tsize: " << tally_ids[char_index].size() << "\n";
    // }

    uint64_t next_check_point = tally_b * tally_checkpoints;
    // If we are at the end, just look at the last run
    if (tally_b * tally_checkpoints >= r) {
        next_check_point = r - 1;
    }

    int64_t prev_check_point = tally_a * tally_checkpoints;

    // Find the closest checkpoint which is either before or after the current row
    // For the first run at 0, the previous tallies are not stored, so we have to look below
    // Since the above direction never helps, for now always set forward_direction to true
    bool forward_direciton = true;
    if (forward_direciton or (tally_a == 0 or tally_b - idx <= idx - tally_a)) {
        uint64_t rows_until_tally = 0;

        // The id for the run at the next_check_point
        id = tally_ids[char_index][tally_b].get();

        uint64_t last_count = 0;
        uint64_t last_id = r;

        // TODO: Could prefetching help with performance?
        // Prelimninary results were not hopeful
        // for (uint64_t i = 0; i < tally_checkpoints; i++) {
        //     my_prefetch_rr((void*)(&(rlbwt[0]) + idx + i));
        //     my_prefetch_rr((void*)(&(rlbwt[0]) + id - i));
        // }


        // Look for the rows between idx and the next_check_point
        // Count how many rows with the same character exists between
        // the current run head and the next run head for which we know the id
        // dbg << "from: " << idx << " to: " << next_check_point << "\n";
        for (uint64_t i = idx; i < next_check_point; i++) {
            // Only count the rows with the same character
            if (get_char(i) == get_char(idx)) {
                rows_until_tally += get_n(i);
                last_id = i;
            }
        }

        // The id stored at the checkpoint is actually the id for the current run
        // because there was no run with the same character between the current run and the next_check_point
        if (last_id == idx and get_char(idx) != get_char(next_check_point)) {
            return id;
        }
        if (last_id == r) {
            throw std::runtime_error(ERROR_MSG("[get id - tally mode] last_id should never be equal to r."));
        }

        // We should know what will be the offset of the run head in the destination run (id)
        uint16_t offset = get_offset(next_check_point);

        // At the checkpoint, one id is stored for each character
        // For the character of the checkpoint, the id of the checkpoint is stored
        // For other characters, the id of the last run before the checkpoint with that character is stored
        // So, we have to decrease the number of rows of the last run as they are not between the current run
        // and the run for which we know the id
        // Also, for other characters, offset should be stored based on that run (not the checkpoint run)
        if (get_char(idx) != get_char(next_check_point)) {
            rows_until_tally -= get_n(last_id);
            offset = get_offset(last_id);
        }

        // Sanity check, the offset in the destination run should be always smaller than the length of that run
        if (offset >= get_n(id)) {
            throw std::runtime_error(ERROR_MSG("[get id - tally mode] offset: " + std::to_string(offset) + " n: " + std::to_string(get_n(id)) +
                                                "\nThe offset in the destination run should be always smaller than the length of that run."));
        }
        if (offset >= rows_until_tally) {
            return id;
        } else {
            rows_until_tally -= (offset + 1);
            id -= 1;
        }

        while (rows_until_tally != 0) {
            if (rows_until_tally >= get_n(id)) {
                rows_until_tally -= get_n(id);
                id -= 1;
                my_prefetch_rr((void*)(&(rlbwt[0]) + id));
            } else {
                rows_until_tally = 0;
            }
        }
    } else {
        uint64_t rows_until_tally = 0;

        // The id for the run at the previous_check_point
        id = tally_ids[char_index][tally_a].get();

        uint64_t last_count = 0;
        uint64_t last_id = r;

        // Look for the rows between idx and the prev_check_point
        // Count how many rows with the same character exists between
        // the current run head and the prev run head for which we know the id
        // dbg << "from: " << idx << " to: " << prev_check_point << "\n";
        for (int64_t i = idx - 1; i >= prev_check_point; i--) {
            // Only count the rows with the same character
            if (get_char(i) == get_char(idx)) {
                rows_until_tally += get_n(i);
                last_id = i;
            }
        }

        // The last checkpoint is storing an id for a run above, so we have to count the
        // number of rows in that run too
        // We have to first find that run with the same character
        if (get_char(idx) != get_char(prev_check_point)) {
            last_id = prev_check_point;
            while (get_char(idx) != get_char(last_id)) {
                last_id -= 1;
            }
            rows_until_tally += get_n(last_id);
        }

        if (last_id == r) {
            throw std::runtime_error(ERROR_MSG("[get id - tally mode] last_id should never be equal to r.\nidx: " +
                                              std::to_string(idx) + " tally_a: " + std::to_string(tally_a) + "\n"));
        }

        uint16_t offset = get_offset(last_id);

        // Sanity check, the offset in the destination run should be always smaller than the length of that run
        if (offset >= get_n(id)) {
            throw std::runtime_error(ERROR_MSG("[get id - tally mode] offset: " + std::to_string(offset) + " n: " + std::to_string(get_n(id)) +
                                                "idx: " + std::to_string(idx) + " last_id: " + std::to_string(last_id) + " id: " + std::to_string(id) +
                                                " prev_check_point: " + std::to_string(prev_check_point) +
                                                "\nThe offset in the destination run should be always smaller than the length of that run."));
        }
        if ((get_n(id) - offset - 1) >= rows_until_tally) {
            return id;
        } else {
            rows_until_tally -= (get_n(id) - offset);
            id += 1;
        }

        while (rows_until_tally != 0) {
            if (rows_until_tally >= get_n(id)) {
                rows_until_tally -= get_n(id);
                id += 1;
                my_prefetch_rr((void*)(&(rlbwt[0]) + id));
            } else {
                rows_until_tally = 0;
            }
        }
    }

    return id;
#endif
}

char MoveStructure::get_char(uint64_t idx) {
    if (idx == end_bwt_idx)
        return '$';
    else
        return alphabet[rlbwt[idx].get_c()];
}

#if SPLIT_MAX_RUN
uint64_t MoveStructure::get_n(uint64_t idx) {
    return rlbwt[idx].get_n();
}

uint64_t MoveStructure::get_offset(uint64_t idx) {
    return rlbwt[idx].get_offset();
}
#endif

#if SPLIT_THRESHOLDS_TRUE
uint64_t MoveStructure::get_thresholds(uint64_t idx, uint32_t alphabet_index) {
    return rlbwt[idx].get_threshold(alphabet_index) == 0 ? 0 : get_n(idx);
}
#endif

#if SPLIT_THRESHOLDS_FALSE
uint64_t MoveStructure::get_n(uint64_t idx) {
    if (rlbwt[idx].is_overflow_n()) {
        return n_overflow[rlbwt[idx].get_n()];
    } else {
        return rlbwt[idx].get_n();
    }
}

uint64_t MoveStructure::get_offset(uint64_t idx) {
    if (rlbwt[idx].is_overflow_offset()) {
        return offset_overflow[rlbwt[idx].get_offset()];
    } else {
        return rlbwt[idx].get_offset();
    }
}

// for all the threshold related functions
uint64_t MoveStructure::get_thresholds(uint64_t idx, uint32_t alphabet_index) {
    if (rlbwt[idx].is_overflow_thresholds()) {
        return thresholds_overflow[get_rlbwt_thresholds(idx, alphabet_index)][alphabet_index];
    } else {
        return get_rlbwt_thresholds(idx, alphabet_index);
    }
}

uint16_t MoveStructure::get_rlbwt_thresholds(uint64_t idx, uint16_t i) {
    if (i >= alphabet.size() - 1) {
        throw std::runtime_error(ERROR_MSG("[get_rlbwt_thresholds] get_thresholds: " + std::to_string(i) +
                                 " is greater than or equal to " +std::to_string(alphabet.size() - 1)));
    }

    uint8_t status = rlbwt[idx].get_threshold_status(i);
    switch (status) {
        case 0: return 0; break;
        case 1: return rlbwt[idx].get_threshold(); break;
        case 3: return get_n(idx); break;
        default:
            throw std::runtime_error(ERROR_MSG("[get_rlbwt_thresholds] Undefined status for thresholds status: " + std::to_string(status)));
    }
}

void MoveStructure::set_rlbwt_thresholds(uint64_t idx, uint16_t i, uint16_t value) {
    if (i >= alphabet.size() - 1) {
        throw std::runtime_error(ERROR_MSG("[set_rlbwt_thresholds] get_thresholds: " + std::to_string(i) +
                                 " is greater than or equal to " + std::to_string(alphabet.size() - 1)));
    }

    uint8_t status = 0;
    if (value == 0) {
        status = 0;
    } else if (value == get_n(idx)) {
        status = 3;
    } else {
        status = 1;
        // [TODO] Not all the states where the multiple non-trivial thresholds exists are checked here
        if (rlbwt[idx].get_threshold() != value and rlbwt[idx].get_threshold() != 0 and
            rlbwt[idx].get_threshold() != get_n(idx) and !rlbwt[idx].is_overflow_thresholds()) {
            if (movi_options->is_debug()) {
                DEBUG_MSG("idx: " + std::to_string(idx) + " i: " + std::to_string(i) + " value: " + std::to_string(value));
                DEBUG_MSG(std::to_string(rlbwt[idx].get_threshold()) + " " + std::to_string(!rlbwt[i].is_overflow_thresholds()));
                DEBUG_MSG("There are more than 1 non-trivial threshold values.");
            }
            rlbwt[i].set_overflow_thresholds();
            return;
        }
        rlbwt[idx].set_threshold(value);
    }
    rlbwt[idx].set_threshold_status(i, status);
}
#endif // for all the threshold related functions

bool MoveStructure::check_alphabet(char& c) {
    if (use_separator()) {
        if (c == SEPARATOR) {
            return false;
        }
    }
    if (movi_options->ignore_illegal_chars_status() > 0) {
        if (alphamap[static_cast<uint64_t>(c)] == alphamap.size()) {
            thread_local ThreadRandom random_generator;
            c = movi_options->ignore_illegal_chars_status() == 1 ?  'A' : alphabet[ random_generator.get_random() % alphabet.size() ];
            return true;
        }
    }
    return alphamap[static_cast<uint64_t>(c)] != alphamap.size();
}

void MoveStructure::analyze_rows() {
    for (int i = 0; i < first_runs.size(); i++) {
        INFO_MSG(std::to_string(i) + "\t" + std::to_string(first_runs[i]) + "\t" + std::to_string(last_runs[i]));
    }
    std::vector<uint64_t> counts_length(16,0);
    std::vector<uint64_t> counts_offset(16,0);
    std::vector<uint64_t> counts_threshold0(16,0);
    std::vector<uint64_t> counts_threshold1(16,0);
    std::vector<uint64_t> counts_threshold2(16,0);
    uint64_t counter = 0;
    uint64_t split_thresholds = 0;
    uint64_t end_row = 0;
    for (uint64_t i = 0; i < r; i++) {
        end_row += get_n(i);
        std::cout << end_row << "\t" << get_char(i) << "\n";
        counter += 1;
        if (i % 1000000 == 0) PROGRESS_MSG(std::to_string(i) + "\t" + std::to_string(split_thresholds));
        for (int j = 0; j < 16; j ++) {
            if (get_n(i) >= std::pow(2,j + 1)) {
                counts_length[j] += 1;
            }
            if (get_offset(i) >= std::pow(2,j + 1)) {
                counts_offset[j] += 1;
            }
#if SPLIT_THRESHOLDS_FALSE
            if (get_thresholds(i, 0) >= std::pow(2,j + 1)) {
                counts_threshold0[j] += 1;
            }
            if (get_thresholds(i, 1) >= std::pow(2,j + 1)) {
                counts_threshold1[j] += 1;
            }
            if (get_thresholds(i, 2) >= std::pow(2,j + 1)) {
                counts_threshold2[j] += 1;
            }
#endif
        }
    }
    INFO_MSG("split_thresholds: " + std::to_string(split_thresholds));
    INFO_MSG("counter: " + std::to_string(counter));
    INFO_MSG("counts_length:");
    for (int j=0; j < 16; j++) {
        INFO_MSG(std::to_string(j + 1) + "\t" + std::to_string(counts_length[j]));
    }
    INFO_MSG("counts_offset:");
    for (int j=0; j < 16; j++) {
        INFO_MSG(std::to_string(j + 1) + "\t" + std::to_string(counts_offset[j]));
    }
    INFO_MSG("counts_threshold0:");
    for (int j=0; j < 16; j++) {
        INFO_MSG(std::to_string(j) + ": " + std::to_string(counts_threshold0[j]));
    }
    INFO_MSG("counts_threshold1:");
    for (int j=0; j < 16; j++) {
        INFO_MSG(std::to_string(j) + ": " + std::to_string(counts_threshold1[j]));
    }
    INFO_MSG("counts_threshold2:");
    for (int j=0; j < 16; j++) {
        INFO_MSG(std::to_string(j) + ": " + std::to_string(counts_threshold2[j]));
    }
}


void MoveStructure::print_basic_index_info() {
    INFO_MSG("Basic index characteristics:");
    INFO_MSG("\tn: " + std::to_string(length));
    INFO_MSG("\tr: " + std::to_string(r));
    if (original_r != 0) {
        INFO_MSG("\toriginal_r: " + std::to_string(original_r));
    }
    INFO_MSG("\tend_bwt_idx ($): " + std::to_string(end_bwt_idx));
}

void MoveStructure::print_stats() {
    print_basic_index_info();
    INFO_MSG("n/r:\t" + std::to_string(static_cast<double>(length)/r));

    for (int i = 0; i < first_runs.size(); i++) {
        INFO_MSG(std::string(1, i == 0 ? '$' : alphabet[i-1])
                  + "\t" + std::to_string(i)
                  + "\t" + std::to_string(first_runs[i]) + ":" + std::to_string(first_offsets[i])
                  + "\t" + std::to_string(last_runs[i])  + ":" + std::to_string(last_offsets[i]));
    }

    for (int i = 0; i < counts.size(); i++) {
        INFO_MSG("counts[" + std::to_string(i) + "]: " + std::to_string(counts[i]));
    }

    if (original_r != 0) {
        INFO_MSG("original_r: " + std::to_string(original_r));
        INFO_MSG("n/original_r: " + std::to_string(static_cast<double>(length)/original_r));
    }
    INFO_MSG("Size of the rlbwt table: " + std::to_string(sizeof(rlbwt[0]) * rlbwt.size() * (0.000000001)));
#if BLOCKED_MODES
    INFO_MSG("Size of the block table: " + std::to_string(sizeof(id_blocks[0][0]) * id_blocks[0].size() * id_blocks.size() * (0.000000001)));
#endif
#if TALLY_MODES
    INFO_MSG("Size of the tally table: " + std::to_string(sizeof(tally_ids[0][0]) * tally_ids[0].size() * tally_ids.size() * (0.000000001)));
#endif

    if (movi_options->is_output_ids()) {
        output_ids();
    }
}

uint64_t MoveStructure::LF(uint64_t row_number, uint64_t alphabet_index) {
    uint64_t lf = 0;
    lf += 1;
    for (uint64_t i = 0; i < alphabet_index; i++) {
        lf += counts[i];
    }
    auto& occ_rank = *occs_rank[alphabet_index];
    lf += static_cast<uint64_t>(occ_rank(row_number));
    return lf;
}

uint64_t MoveStructure::LF_heads(uint64_t run_number, uint64_t alphabet_index) {
    uint64_t lf = 0;
    lf += 1;
    for (uint64_t i = 0; i < alphabet_index; i++) {
        lf += counts[i];
    }
    lf += heads_rank[run_number];
    return lf;
}

uint64_t MoveStructure::fast_forward(uint64_t& offset, uint64_t idx, uint64_t x) {

    uint64_t idx_ = idx;
    if (movi_options->is_debug()) {
        DEBUG_MSG("\t\tFast forwarding:");
        DEBUG_MSG("\t\tidx: " + std::to_string(idx) + " offset: " + std::to_string(offset) + " n:" + std::to_string(get_n(idx)));
    }

    while (idx < r - 1 && offset >= get_n(idx)) {
        offset -= get_n(idx);
        idx += 1;
        if (movi_options->is_debug()) {
            DEBUG_MSG("\t\tff offset based: +" + std::to_string(idx - idx_));
        }
    }

    if (movi_options->is_debug()) {
        DEBUG_MSG("\t\tAfter fast forwarding: idx: " + std::to_string(idx) + " offset: " + std::to_string(offset) + " n:" + std::to_string(get_n(idx)));
    }

    return idx - idx_;
}

bool MoveStructure::use_separator() {
    if (movi_options->use_separators() or (alphabet[0] == SEPARATOR and alphabet.size() == 5)) {
        return true;
    }
    return false;
}