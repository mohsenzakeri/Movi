#include "move_structure.hpp"

void MoveStructure::build() {

    length = compute_length_from_bwt();
    std::cerr << GENERAL_MSG("The total length of the BWT (n): " + std::to_string(length));

    std::cerr << GENERAL_MSG("Building starts...");

    initialize_bits();

#if SPLIT_THRESHOLDS_TRUE
    // If thresholds splitting is enabled, bits will track the splitting points
    if (movi_options->is_thresholds()) {
        fill_bits_by_thresholds();
    }
#endif

    detect_move_row_boundaries();

    find_run_heads_information();

    if (alphabet.size() > 4) {
        std::cerr << "Warning: There are more than 4 characters in the reference:\n";
        for (uint64_t i = 0; i < alphabet.size(); i++) {
            std::cerr << "(" << i << ")" << "\t" << alphabet[i] << "\t" << counts[i] << "\n";
        }
        std::cerr << "\n";
    }

    build_move_rows();

#if USE_THRESHOLDS
    // compute the thresholds
    compute_thresholds();
#endif

#if USE_NEXT_POINTERS
    if (constant) {
        std::cerr << GENERAL_MSG("Computing the next ups and downs.");
        compute_nexts();
    }
#endif
    std::cerr << INFO_MSG("The move structure building is done.");
}

void MoveStructure::find_run_heads_information() {
    std::vector<uint64_t> bwt_character_counts(alphabet.size(), 0);

    all_p.resize(r + 1);
    all_p[r] = length;

    heads_rank.resize(r);
    heads_rank[0] = 0;

    uint64_t bwt_offset = 0;

    for (uint64_t i = 0; i < r; i++) {

        // Boundary check for bits vector
        if (bwt_offset >= bits.size()) {
            throw std::runtime_error(ERROR_MSG("[find run heads information] " + std::to_string(bwt_offset) +
                                               " offset is greater than the bits vector length (" + std::to_string(bits.size()) + ")"));
        }

        // The following line is important to tag the run boundaries in the main bitvector
        // Required for computing the id field of the rlbwt rows
        bits[bwt_offset] = 1;

        if (i>0 && i % 1000000 == 0) {
            std::cerr << "processed runs: " << i << "/" << r << "\t";
            std::cerr << "bwt_offset: " << bwt_offset << "\r";
        }

        uint64_t len = lens[i];
        if (heads[i] != END_CHARACTER) {
            bwt_character_counts[alphamap[static_cast<uint64_t>(heads[i])]] += len;
        } else{
            end_bwt_idx = i;
        }

        if (i < r - 1) {
            if (heads[i + 1] != END_CHARACTER) {
                uint64_t next_head_rank = bwt_character_counts[alphamap[static_cast<uint64_t>(heads[i + 1])]];
                heads_rank[i + 1] = next_head_rank;
            } else{
                heads_rank[i + 1] = 0;
            }
        }

        all_p[i] = bwt_offset;
        bwt_offset += len;
    }

    std::cerr << "\n";
}

void MoveStructure::add_detected_run(uint64_t scanned_bwt_length, uint64_t run_char,  uint16_t& run_length) {

    heads.push_back(run_char);
    lens.push_back(run_length);

    if (!movi_options->is_preprocessed()) {
        r += 1;

        // Boundary check for bits vector
        if (scanned_bwt_length + 1 >= bits.size()) {
            throw std::runtime_error(ERROR_MSG("[add_detected_run] " + std::to_string(scanned_bwt_length + 1) +
                                        " offset is greater than the bits vector length (" + std::to_string(bits.size()) + ")"));
        }

        bits[scanned_bwt_length + 1] = 1;
        run_length = 0;
    }
}

uint64_t MoveStructure::compute_length_from_bwt() {

    std::string bwt_filename = movi_options->get_ref_file() + std::string(".bwt");

    uint64_t total_length = 0;

    // The following finds the length of the BWT (n or end_pos)
    // We need to know n to initialize the main bitvector (bits)
    // The main bitvector is used to find the positions of the runs after any splitting
    if (movi_options->is_preprocessed()) {
        std::ifstream heads_file(bwt_filename + ".heads", std::ios::in | std::ios::binary);
        if (!heads_file.good()) {
            throw std::runtime_error(ERROR_MSG("[build] Failed to open the heads file: " + bwt_filename + ".heads"));
        }

        original_run_heads = std::vector<char>((std::istreambuf_iterator<char>(heads_file)), std::istreambuf_iterator<char>());

        original_r  = original_run_heads.size();
        std::cerr << "Number of BWT runs: " << original_r << "\n";

        original_lens.resize(original_r);

        std::ifstream len_file(bwt_filename + ".len", std::ios::in | std::ios::binary);
        if (!len_file.good()) {
            throw std::runtime_error(ERROR_MSG("[build] Failed to open the len file: " + bwt_filename + ".len"));
        }

        for (uint64_t i = 0; i < original_r; i++) {
            if (i>0 && i % 1000000 == 0)
                std::cerr << "Computing total length of the BWT: " << i << "/" << original_r << "\r";

            size_t len = 0;
            len_file.read(reinterpret_cast<char*>(&len), 5);
            original_lens[i] = len;
            total_length += len;
        }
        std::cerr << "\n";
    } else {
        bwt_file.open(bwt_filename);
        if (!bwt_file.good()) {
            throw std::runtime_error(ERROR_MSG("[build] Failed to open the BWT file: " + bwt_filename));
        }

        bwt_file.clear();

        bwt_file.seekg(0, std::ios_base::end);
        std::streampos end_pos = bwt_file.tellg();
        total_length = static_cast<uint64_t>(end_pos);
        if (movi_options->is_verbose())
            std::cerr << "total_length: " << total_length << "\n";

        bwt_file.seekg(0);
    }    

    return total_length;
}

void MoveStructure::initialize_bits() {
    if (nt_splitting) {
        std::string splitting_filename = movi_options->get_ref_file() + std::string(".d_col");
        std::ifstream splitting_file(splitting_filename);
        if (!splitting_file.good()) {
            throw std::runtime_error(ERROR_MSG("[build] Failed to open the d_col file: " + splitting_filename));
        }

        bits.load(splitting_file);
        std::cerr << "bits.size after loading the d_col file: " << bits.size() << "\n";
        rbits = sdsl::rank_support_v<>(&bits);
        std::cerr << "The main bit vector (bits) is loaded from the d_col file.\n";
    } else {
        bits = sdsl::bit_vector(length + 1, 0);
        bits[0] = 1;
        std::cerr << "The main bit vector (bits) is initialized.\n";
    }
}

void MoveStructure::detect_move_row_boundaries() {

    // Assume all the characters in the alphabet might exist in the BWT
    uint64_t all_chars_count = 256;

    // Map from the character ASCII code to the index in the alphabet
    alphamap.resize(all_chars_count);
    std::fill(alphamap.begin(), alphamap.end(), alphamap.size());

    // Counter for the number of occurrences of each character in the BWT
    std::vector<uint64_t> all_possible_chars(all_chars_count, 0);

    // Counters for the number of runs split by MAX_RUN_LENGTH and by thresholds
    uint64_t split_by_max_run = 0;
    uint64_t split_by_thresholds = 0;

    // Any mode with pre-built splitting bitvector (modes 1 and 4) does not work with the preprocessed build
    if (movi_options->is_preprocessed()) {

        uint64_t bwt_offset = 0;

#if SPLIT_THRESHOLDS_TRUE or SPLIT_ARRAY
        rbits = sdsl::rank_support_v<>(&bits);
#endif

        for (uint64_t i = 0; i < original_r; i++) {
            if (i>0 && i % 1000000 == 0)
                std::cerr << "Iterating over original_r: " << i << "/" << original_r << "\r";

            size_t len = static_cast<size_t>(original_lens[i]);

            if (original_run_heads[i] != END_CHARACTER) {
                all_possible_chars[static_cast<size_t>(original_run_heads[i])] += len;
            }

            size_t remaining_length = len;

            // First attempt to split the run by thresholds
            if ((nt_splitting or SPLIT_THRESHOLDS_TRUE) and
                    rbits(bwt_offset + len) != rbits(bwt_offset)) {

                uint64_t run_head = bwt_offset;

                for (uint64_t j = bwt_offset + 1; j < bwt_offset + len; j++) {

                    // Boundary check for bits vector
                    if (j >= bits.size()) {
                        throw std::runtime_error(ERROR_MSG("[build - finding the runs (preprocessed)] " + std::to_string(j) +
                                                " offset is greater than the bits vector length (" + std::to_string(bits.size()) + ")"));
                    }

                    if (bits[j]) {
                        size_t current_run_length = j - run_head;
                        remaining_length -= current_run_length;
#if SPLIT_THRESHOLDS_TRUE
                        split_by_thresholds += 1;
#endif

                        // Check if the run length is greater than MAX_RUN_LENGTH
                        while (current_run_length > MAX_RUN_LENGTH) {
                            split_by_max_run += 1;
                            uint16_t run_length = static_cast<uint16_t>(MAX_RUN_LENGTH);
                            add_detected_run(bwt_offset, original_run_heads[i], run_length);
                            current_run_length -= MAX_RUN_LENGTH;
                        }

                        // Anything remaining will be added as a separate run
                        if (current_run_length > 0) {
                            uint16_t run_length = static_cast<uint16_t>(current_run_length);
                            add_detected_run(bwt_offset, original_run_heads[i], run_length);
                        }
                        run_head = j;
                    }
                }

            }

            bwt_offset += len;

            // Split the remaining length by MAX_RUN_LENGTH
            while (remaining_length > MAX_RUN_LENGTH) {
                split_by_max_run += 1;
                uint16_t run_length = static_cast<uint16_t>(MAX_RUN_LENGTH);
                add_detected_run(bwt_offset, original_run_heads[i], run_length);
                remaining_length -= MAX_RUN_LENGTH;
            }

            // Anything remaining will be added as a separate run
            if (remaining_length > 0) {
                uint16_t run_length = static_cast<uint16_t>(remaining_length);
                add_detected_run(bwt_offset, original_run_heads[i], run_length);
            }

        }
        std::cerr << "\n";

        r = heads.size();

        build_alphabet(all_possible_chars);

    } else {

        // Reading the BWT from the file
        uint64_t current_char = bwt_file.get();
        uint64_t next_char = bwt_file.get();


        if (current_char != END_CHARACTER) {
            all_possible_chars[current_char] += 1;
        }

        all_p.push_back(0);
        heads.push_back(current_char);
        heads_rank.push_back(0);

        uint16_t run_length = 0;
        original_r = 1;
        r = 1;

        // Tracks how many characters in the BWT have been processed
        uint64_t scanned_bwt_length = 0;

        std::cerr << "Reading over the uncompressed BWT...\n";
        while (next_char != EOF) { // && current_char != 10
            run_length += 1;
            if (next_char != END_CHARACTER) {
                all_possible_chars[next_char] += 1;
            }

            if (original_r % 1000000 == 0) {
                std::cerr << "original_r: " << original_r << "\t";
                std::cerr << "r: " << r << "\t";
                std::cerr << "scanned_bwt_length: " << scanned_bwt_length << "/" << length << "\r";
            }

            // Boundary check for bits vector
            if (scanned_bwt_length + 1 >= bits.size()) {
                throw std::runtime_error(ERROR_MSG("[build - finding the runs] Finding the runs: " + std::to_string(scanned_bwt_length + 1) +
                                    " offset is greater than the bits vector length (" + std::to_string(bits.size()) + ")"));
            }

            // The first row is already set and accounted for, so we skip
            if (current_char != next_char) {
                // 1) A new run is detected if the next character is different
                original_r += 1;

                if (nt_splitting and !bits[scanned_bwt_length + 1]) {
                    throw std::runtime_error(ERROR_MSG("[build - finding the runs] There is something wrong with the splitting vector.\n" +
                                    "The run boundaries should have been set to 1 since a new character was detected."));
                }

                add_detected_run(scanned_bwt_length, next_char, run_length);

            } else if ((nt_splitting or SPLIT_THRESHOLDS_TRUE) and bits[scanned_bwt_length + 1] == 1) {
                // 2) A new run is detected based on a non-trivial threshold or the Nishimoto-Tabei splitting bitvector
                // The bit was already set by one of the threshold values or the splitting bitvector
                // So, we have found a new run, and reset the run length

                add_detected_run(scanned_bwt_length, next_char, run_length);

#if SPLIT_THRESHOLDS_TRUE
                split_by_thresholds += 1;
#endif
            } else if (run_length == MAX_RUN_LENGTH) {
                // 3) A new run is detected if the length of the run is greater than MAX_RUN_LENGTH

                add_detected_run(scanned_bwt_length, next_char, run_length);
                split_by_max_run += 1;
            }

            current_char = next_char;
            next_char = bwt_file.get();
            scanned_bwt_length += 1;
        }
        std::cerr << "\n";

        if (r != heads.size()) {
            throw std::runtime_error(ERROR_MSG("[build - finding the runs] r: " + std::to_string(r) + " heads.size(): " + std::to_string(heads.size()) +
                                    " The number of runs is not consistent."));
        }

        if (length != scanned_bwt_length + 1) {
            throw std::runtime_error(ERROR_MSG("[build - finding the runs] length: " + std::to_string(length) +
                                    " scanned_bwt_length: " + std::to_string(scanned_bwt_length) +
                                    " The length of the BWT is not consistent."));
        }

        build_alphabet(all_possible_chars);

    }

    // we don't need the original lengths anymore after all the runs are detected
    original_lens.clear();
    original_lens.shrink_to_fit();

    if (movi_options->is_thresholds()) {
        std::cerr << GENERAL_MSG("Number of runs added by thresholds splitting: " + std::to_string(split_by_thresholds));
    }
    std::cerr << "n: " << length << "\n";
    std::cerr << "r: " << r << "\n";
    std::cerr << "n/r:\t" << static_cast<double>(length)/r << "\n";
    std::cerr << "original_r: " << original_r << "\n";
}

void MoveStructure::build_alphabet(std::vector<uint64_t>& all_possible_chars) {
    uint64_t alphabet_index = 0;
    for (uint64_t i = 1; i < all_possible_chars.size(); i++) {
        if (all_possible_chars[i] != 0) {
            auto current_char = static_cast<unsigned char>(i);
            if (movi_options->is_verbose())
                std::cerr << "(" << i << ")"
                          << "\t" << current_char
                          << "\t index:" << alphabet_index
                          << "\t count:" << all_possible_chars[i] << "\n";

            alphabet.push_back(current_char);
            counts.push_back(all_possible_chars[i]);
            alphamap[i] = alphabet_index;
            alphabet_index += 1;
        }
    }
    std::cerr << GENERAL_MSG("All the characters are indexed.") << "\n";
}

void MoveStructure::build_move_rows() {
    if (movi_options->is_verbose() and bits.size() < 1000)
        std::cerr << "bits: " << bits << "\n";

    if (!nt_splitting) // The rank vector is already built if it was the splitting mode
        rbits = sdsl::rank_support_v<>(&bits);


    rlbwt.resize(r);

    if (movi_options->is_verbose() ) {
        std::cerr << "bits.size(): " << bits.size() << "\n";
        std::cerr << "rank_support_v<>(&bits)(bits.size()): "
                  << sdsl::rank_support_v<>(&bits)(bits.size()) << "\n";

        std::cerr << "r: " << r << "\t"
                  << "rlbwt.size(): " << rlbwt.size() << "\t"
                  << "all_p.size(): " << all_p.size() << "\n";
        std::cerr << "end_bwt_idx: " << end_bwt_idx << "\n\n";
    }

#if TALLY_MODES
    tally_ids.resize(alphabet.size());

    tally_checkpoints = movi_options->get_tally_checkpoints();
    uint64_t tally_ids_rows_count = r / tally_checkpoints + 2;

    std::vector<uint64_t> current_tally_ids;
    for (uint32_t alphabet_ind = 0; alphabet_ind < alphabet.size(); alphabet_ind++) {
        current_tally_ids.push_back(r);
        tally_ids[alphabet_ind].resize(tally_ids_rows_count);
    }
#endif

#if BLOCKED_MODES
    std::vector<uint64_t> raw_ids;
    raw_ids.resize(r);
#endif

    uint64_t offset = 0;
    uint64_t max_len = 0;

    // Building the move structure rows with O(r) loop
    for (uint64_t r_idx = 0; r_idx < r; r_idx++) {
        if (r_idx % 1000000 == 0)
            std::cerr << "Building Movi rows: " << r_idx << "/" << r << "\r";

        uint64_t lf  = 0;
        if (r_idx != end_bwt_idx) {
            uint64_t alphabet_index = alphamap[static_cast<uint64_t>(heads[r_idx])];
            // lf = LF(all_p[r_idx], alphabet_index);
            lf = LF_heads(r_idx, alphabet_index);
        } else {
            lf = 0;
        }

        uint64_t pp_id = rbits(lf) - 1;

        // Boundary check for bits vector
        if (lf >= bits.size()) {
            throw std::runtime_error(ERROR_MSG("[build rows] Building the move structure rows: " + std::to_string(lf) +
                                    " offset is greater than the bits vector length (" + std::to_string(bits.size()) + ")"));
        }

        if (bits[lf] == 1)
            pp_id += 1;

        // TODO: can we use lens instead? or use it to confirm the following computation
        uint64_t len = all_p[r_idx + 1] - all_p[r_idx];

        // check the boundaries before performing select
        if (pp_id >= r) {
            throw std::runtime_error(ERROR_MSG("[build rows] pp_id: " + std::to_string(pp_id) + " r: " + std::to_string(r) +
                                    " r_idx: " + std::to_string(r_idx) + " lf: " + std::to_string(lf)));
        }

        if (lf < all_p[pp_id]) {
            throw std::runtime_error(ERROR_MSG("[build rows] pp_id: " + std::to_string(pp_id) + " lf: " + std::to_string(lf) +
                                    " all_p[pp_id]: " + std::to_string(all_p[pp_id])));
        }

        offset = lf - all_p[pp_id];

        if (movi_options->is_verbose() and r_idx == 0) { // 0 or any run to be inspected
            std::cerr << "r_idx: " << r_idx
                        << " len: " << len
                        << " lf: " << lf
                        << " offset: " << offset
                        << " pp_id: " << pp_id
                        << "\n";
        }

#if TALLY_MODES
        rlbwt[r_idx].init(len, offset);

        if (r_idx != end_bwt_idx) {
            // Only update the tally corresponding to the current run's character
            uint16_t char_index = alphamap[heads[r_idx]];
            // The first tally id might have been set to the wrong value since no run with that character was observed
            // So, we fix the values once we the first run with that character
            if (current_tally_ids[char_index] == r) {
                uint64_t tally_ids_index = r_idx / tally_checkpoints;
                for (int tid = 0; tid <= tally_ids_index; tid++) {
                    tally_ids[char_index][tid].set_value(pp_id);
                }
            }
            current_tally_ids[char_index] = pp_id;
        }

        if (r_idx % tally_checkpoints == 0) {
            uint64_t tally_ids_index = r_idx / tally_checkpoints;
            // We reached a check_point, store the tally ids for all the characters
            for (uint32_t alphabet_ind = 0; alphabet_ind < alphabet.size(); alphabet_ind++) {
                tally_ids[alphabet_ind][tally_ids_index].set_value(current_tally_ids[alphabet_ind]);
            }
        }
#else


#if BLOCKED_MODES
        raw_ids[r_idx] = pp_id;
        pp_id = 0;
#endif

        // rlbwt[r_idx].init(bwt_row, len, lf, offset, pp_id);
        rlbwt[r_idx].init(len, offset, pp_id);
#endif

        // To take care of cases where length of the run
        // does not fit in uint16_t
#if SPLIT_MAX_RUN
        if (offset > MAX_RUN_LENGTH or len > MAX_RUN_LENGTH) {
            // Should not get here in the regular mode: MODE = 3
            throw std::runtime_error(ERROR_MSG("[build rows] The length or the offset are too large.\n" +
                                                "offset: " + std::to_string(offset) + "\tlength: " + std::to_string(length) + "\n"));
        }
#endif
#if SPLIT_THRESHOLDS_FALSE
        if (len > MAX_RUN_LENGTH) {
            throw std::runtime_error(ERROR_MSG("[build rows] This shouldn't happen anymore, because the run splitting should be applied in every mode now."));
            n_overflow.push_back(len);
            if (n_overflow.size() - 1 >= MAX_RUN_LENGTH) {
                throw std::runtime_error(ERROR_MSG("[build rows] The number of runs with overflow n is beyond " + std::to_string(MAX_RUN_LENGTH) + "! " + std::to_string(n_overflow.size() - 1)));
            }
            rlbwt[r_idx].set_n(n_overflow.size() - 1);
            rlbwt[r_idx].set_overflow_n();
        }
        if (offset > MAX_RUN_LENGTH) {
            throw std::runtime_error(ERROR_MSG("[build rows] This shouldn't happen anymore, because the run splitting should be applied in every mode now."));
            offset_overflow.push_back(offset);
            if (offset_overflow.size() - 1 >= MAX_RUN_LENGTH) {
                throw std::runtime_error(ERROR_MSG("[build rows] The number of runs with overflow offset is beyond " + std::to_string(MAX_RUN_LENGTH) + "! " + std::to_string(offset_overflow.size() - 1)));
            }
            rlbwt[r_idx].set_offset(offset_overflow.size() - 1);
            rlbwt[r_idx].set_overflow_offset();
        }
#endif
        if (len > max_len)
            max_len = len;
        if (movi_options->is_logs()) {
            if (run_lengths.find(len) != run_lengths.end())
                run_lengths[len] += 1;
            else
                run_lengths[len] = 1;
        }
        rlbwt[r_idx].set_c(heads[r_idx], alphamap);
    }
#if TALLY_MODES
    // Set the last tally_id using the current_tally_ids
    std::cout << "\n";
    for (uint32_t alphabet_ind = 0; alphabet_ind < alphabet.size(); alphabet_ind++) {
        tally_ids[alphabet_ind][tally_ids[alphabet_ind].size()-1].set_value(current_tally_ids[alphabet_ind]);
    }
#endif

    std::cerr << GENERAL_MSG("Max run length (after splitting if enabled): " + std::to_string(max_len));

    std::cerr << INFO_MSG("All the move rows are built.");

    find_base_interval_data();

#if BLOCKED_MODES
    compute_blocked_ids(raw_ids);
#endif
}

void MoveStructure::find_base_interval_data() {
    first_runs.push_back(0);
    first_offsets.push_back(0);
    last_runs.push_back(0);
    last_offsets.push_back(0);
    uint64_t char_count = 1;
    for (uint64_t i = 0; i < counts.size(); i++) {
        uint64_t last_run = last_runs.back();
        uint64_t last_offset = last_offsets.back();
        if (last_offset + 1 >= rlbwt[last_run].get_n()) {
            first_runs.push_back(last_run + 1);
            first_offsets.push_back(0);
        } else {
            first_runs.push_back(last_run);
            first_offsets.push_back(last_offset + 1);
        }
        char_count += counts[i];
        // char_count points to the position of the first row with character i in the BWT
        if (movi_options->is_verbose())
            std::cerr << i << " char_count: " << char_count << "\n";
        auto occ_rank = rbits(char_count);
        // occ_rank is the number of set bits until and including the last row with character i - 1
        if (movi_options->is_verbose())
            std::cerr << i << " occ_rank: " << occ_rank << "\n";
        // The index of the run is number of set bits - 1 (because we count from 0):
        last_runs.push_back(static_cast<uint64_t>(occ_rank - 1));
        if (movi_options->is_verbose())
            std::cerr << char_count << " " << all_p[last_runs.back()] << " " << get_n(last_runs.back()) <<  "\n";
        last_offsets.push_back(char_count - all_p[last_runs.back()] - 1);
    }
    if (movi_options->is_verbose()) {
        for (uint64_t i = 0; i < first_runs.size(); i++) {
            std::cerr << "<--- " << first_runs[i] << "\t" << first_offsets[i] << "\n";
            std::cerr << "    -\n    -\n    -\n";
            std::cerr << ">--- " << last_runs[i] << "\t" << last_offsets[i] << "\n";
        }
    }
}

void MoveStructure::fill_bits_by_thresholds() {
    for (int i = 0; i < thresholds.size(); i++) {

        // Boundary check for bits vector
        if (thresholds[i] >= bits.size()) {
            throw std::runtime_error(ERROR_MSG("[fill bits by thresholds] fill_bits_by_thresholds: " + std::to_string(thresholds[i]) +
                                               " offset is greater than the bits vector length (" + std::to_string(bits.size()) + ")"));
        }

        bits[thresholds[i]] = 1;
    }
    std::cerr << GENERAL_MSG("The bits vector is updated by thresholds.");
}


#if USE_THRESHOLDS
void MoveStructure::compute_thresholds() {
    // initialize the start threshold at the last row
    std::vector<uint64_t> alphabet_thresholds(alphabet.size(), length);

    uint64_t thr_i = original_r - 1;
    uint64_t run_p = 0;

    if (movi_options->is_verbose()) {
        std::cerr << "thresholds.size():" << thresholds.size() << " length: "
                  << length << " r: " << r <<  " original_r: " << original_r << "\n";
        std::cerr << "thresholds[r]: " << thresholds[original_r-1] << " r-1: "
                  << thresholds[original_r - 2] << " r-2: " << thresholds[original_r - 3] << "\n";
    }

    for (uint64_t i = rlbwt.size() - 1; i > 0; --i) {
        if (i % 1000000 == 0)
            std::cerr << "Updating the thresholds per run: " << (rlbwt.size() - i) << "/" << rlbwt.size() << "\r";

        char rlbwt_c = alphabet[rlbwt[i].get_c()];

        if (movi_options->is_verbose() and i >= rlbwt.size() - 10) {
            std::cerr << "i: " << i << "\n"
                << "rlbwt[i].get_offset(): " << get_offset(i) << "\n "
                << "get_n(i): " << get_n(i) << "\n"
                << "thresholds[thr_i]: " << thresholds[thr_i] << " "
                << "rlbwt_c: " << rlbwt_c << "\n";
        }

        std::vector<uint64_t> current_thresholds;
        current_thresholds.resize(alphabet.size() - 1);

        for (uint64_t j = 0; j < alphabet.size(); j++) {
            if (alphabet[j] == rlbwt_c) {
                if (thr_i >= thresholds.size()) {
                    throw std::runtime_error(ERROR_MSG("[compute thresholds] thr_i = " + std::to_string(thr_i) + " is out of bound:\n" +
                                                        " thresholds.size = " + std::to_string(thresholds.size()) +
                                                        "\nThe thresholds are not correct."));
                }
                alphabet_thresholds[j] = thresholds[thr_i];
            } else {
                if (alphamap_3[alphamap[rlbwt_c]][j] >= alphabet.size() - 1) {
                    throw std::runtime_error(ERROR_MSG("[compute thresholds] alphamap_3 is not working in general:\n"
                                "alphabet.size() - 1 = " + std::to_string(alphabet.size() - 1) + "\n"
                                "alphamap_3[alphamap[rlbwt_c]][j] = " + std::to_string(alphamap_3[alphamap[rlbwt_c]][j]) +
                                "\nThe alphabet map (alphamap_3) is not correct."));
                }

                if (alphabet_thresholds[j] >= all_p[i] + get_n(i)) {
                    // rlbwt[i].thresholds[j] = get_n(i);
                    if (i == end_bwt_idx) {
                        end_bwt_idx_thresholds[j] = get_n(i);
                        if (movi_options->is_verbose()) {
                            std::cerr << "condition 1: end_bwt_idx_thresholds[" << j << "]:" << end_bwt_idx_thresholds[j] << "\n";
                        }
                        continue;
                    }
#if SPLIT_THRESHOLDS_FALSE
                    set_rlbwt_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j], get_n(i));
#endif
#if SPLIT_THRESHOLDS_TRUE
                    rlbwt[i].set_threshold(alphamap_3[alphamap[rlbwt_c]][j], 1);
#endif
                    current_thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = get_n(i);
                } else if (alphabet_thresholds[j] <= all_p[i]) {
                    if (i == end_bwt_idx) {
                        end_bwt_idx_thresholds[j] = 0;
                        if (movi_options->is_verbose()) {
                            std::cerr << "condition 2: end_bwt_idx_thresholds[" << j << "]:" << end_bwt_idx_thresholds[j] << "\n";
                        }
                        continue;
                    }
#if SPLIT_THRESHOLDS_FALSE
                    set_rlbwt_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j], 0);
#endif
#if SPLIT_THRESHOLDS_TRUE
                    rlbwt[i].set_threshold(alphamap_3[alphamap[rlbwt_c]][j], 0);
#endif
                    current_thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = 0;
                } else {
                    if (i == end_bwt_idx) {
                        end_bwt_idx_thresholds[j] = alphabet_thresholds[j] - all_p[i];
                        if (movi_options->is_verbose()) {
                            std::cerr << "condition 3: end_bwt_idx_thresholds[" << j << "]:" << end_bwt_idx_thresholds[j] << "\n";
                        }
                        continue;
                    }
#if SPLIT_THRESHOLDS_FALSE
                    if (alphabet_thresholds[j] - all_p[i] >= std::numeric_limits<uint16_t>::max()) {
                        rlbwt[i].set_overflow_thresholds();
                    }
                    set_rlbwt_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j], alphabet_thresholds[j] - all_p[i]);
#endif
#if SPLIT_THRESHOLDS_TRUE
                    throw std::runtime_error(ERROR_MSG("[compute thresholds] j: " + std::to_string(j) + " i:" + std::to_string(i) +
                              " n:" + std::to_string(get_n(i)) + " at j:" + std::to_string(alphabet_thresholds[j]) +
                              " api:" + std::to_string(all_p[i]) +
                              "\nThis should never happen, since the runs are split at threshold boundaries."));
#endif
                    current_thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = alphabet_thresholds[j] - all_p[i];
                }
                // printing the values for last 10 runs to debug
                if (movi_options->is_verbose() and i >= rlbwt.size() - 10) {
                    std::cerr << "\t j: \t" << j << " "
                        << "alphabet[j]: " << alphabet[j] << "  "
                        << "alphamap_3[alphamap[rlbwt_c]][j]: " << alphamap_3[alphamap[rlbwt_c]][j] << " "
                        << "alphabet_thresholds[j]: " << alphabet_thresholds[j] << " "
#if SPLIT_THRESHOLDS_FALSE
                        << "rlbwt[i].thresholds[j]:" << get_rlbwt_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j]) << "\n";
#endif
#if SPLIT_THRESHOLDS_TRUE
                        << "rlbwt[i].thresholds[j]:" << get_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j]) << "\n";
#endif
                }
            }
        }
        if (i > 0 && (rlbwt[i].get_c() != rlbwt[i - 1].get_c()  || i == end_bwt_idx || i-1 == end_bwt_idx)) {
            thr_i--;
        }
#if SPLIT_THRESHOLDS_FALSE
        if (rlbwt[i].is_overflow_thresholds()) {
            if (thresholds_overflow.size() >= std::numeric_limits<uint16_t>::max()) {
                throw std::runtime_error(ERROR_MSG("[compute thresholds] Undefined behaviour: the number of runs with overflow thresholds is beyond uint16_t!" +
                                    std::to_string(thresholds_overflow.size()) + "\n"));
            }
            for (uint64_t k = 0; k < alphabet.size() - 1; k++) {
                set_rlbwt_thresholds(i, k, thresholds_overflow.size());
            }
            thresholds_overflow.push_back(current_thresholds);
        }
#endif
        run_p += get_n(i);
    }
    std::cerr << "\n";

    // since the thresholds for the first run was not calculated in the for
    for (uint64_t j = 0; j < alphabet.size() - 1; j++) {
#if SPLIT_THRESHOLDS_FALSE
        set_rlbwt_thresholds(0, j, 0);
#endif
#if SPLIT_THRESHOLDS_TRUE
        rlbwt[0].set_threshold(j, 0);
#endif
    }
}
#endif

#if BLOCKED_MODES
void MoveStructure::compute_blocked_ids(std::vector<uint64_t>& raw_ids) {

    block_size = BLOCK_SIZE;
    bool blocked_ids_computed = false;

    while (!blocked_ids_computed) {

        std::vector<uint32_t> block_start_id;
        block_start_id.resize(alphabet.size(), 0);

        uint64_t max_raw_id = 0;
        uint64_t max_blocked_id = 0;
        uint64_t block_count = 0;
        uint32_t max_diff = 0;

        id_blocks.resize(alphabet.size());

        uint64_t max_allowed_blocked_id = MAX_ALLOWED_BLOCKED_ID;
        bool small_block = false;

        for (uint64_t i = 0; i < rlbwt.size(); i++) {
            if (i % 1000000 == 0)
                std::cerr << "i: " << i << "\r";

            uint64_t block_boundary = (block_size) * block_count; // BLOCK_SIZE = 1 << 22 - 1

            if (i >= block_boundary) {

                if (movi_options->is_verbose())
                    std::cerr << "\n\n" << block_count << "\t";

                // A new block is being initiated
                block_count += 1;

                // Store the largest (last) observed id for each character in the last block as a check point for the new block
                for (uint64_t j = 0; j < alphabet.size(); j++) {

                    id_blocks[j].push_back(block_start_id[j]);

                    if (movi_options->is_verbose())
                        std::cerr << id_blocks[j][block_count - 1] << "\t";

                    if (block_count > 1)
                        max_diff = std::max(max_diff, id_blocks[j][block_count - 1] - id_blocks[j][block_count - 2] );
                }

                if (movi_options->is_verbose())
                    std::cerr << "\n";

            }

            // Check for potential errors in blocked id computation
            if (i != end_bwt_idx) {

                uint64_t id = raw_ids[i];
                max_raw_id = std::max(id, max_raw_id);

                uint64_t block_number = i / block_size;

                if (block_number != block_count - 1) {
                    throw std::runtime_error(ERROR_MSG("[compute blocked ids] The block calculation is incorrect.\n" +
                                    "block_count: " + std::to_string(block_count) + " block_number: " + std::to_string(block_number) + "\n"));
                }

                // First Calcuate the id with respect to the first run with the same character -> adjusted_id
                // The first entry of the first _runs stores the ids for the global end run, so +1 is required
                uint64_t adjusted_id = id - first_runs[rlbwt[i].get_c() + 1];
                // Calcuated the distance of the ajustedt_id from the check point of that character
                uint64_t blocked_id = adjusted_id - static_cast<uint64_t>(id_blocks[rlbwt[i].get_c()][block_number]);

                // Check for potential errors in blocked id computation
                if (blocked_id > max_allowed_blocked_id) {

                    std::cerr << "The number of bits in the runs are not enough for storing the blocked_id.\n";
                    std::cerr << "adjusted_id: " << adjusted_id << " id: " << id
                              << " first_runs[rlbwt[i].get_c() + 1]: " << first_runs[rlbwt[i].get_c() + 1] << "\n";
                    std::cerr << "id_blocks[rlbwt[i].get_c()][block_number]: " << id_blocks[rlbwt[i].get_c()][block_number] << "\n";
                    std::cerr << "rlbwt[i].get_c(): " << static_cast<uint32_t>(rlbwt[i].get_c()) << "\n";
                    std::cerr << "blocked_id: " << blocked_id << " max_blocked_id: " << max_blocked_id << "\n";
                    std::cerr << "max_allowed_blocked_id: " << max_allowed_blocked_id << "\n";
                    std::cerr << "block_count: " << block_count << " block_number: " << block_number << "\n";
                    std::cerr << "i: " << i << " block_size: " << block_size << "\n";

                    small_block = true;
                    break;
                }

                rlbwt[i].set_id(blocked_id);

                // Check if the blocked_id is stored correctly in the run
                if (rlbwt[i].get_id() != blocked_id) {
                    throw std::runtime_error(ERROR_MSG("[compute blocked ids] The id was not stored in the move row correctly.\n" +
                                    "rlbwt[i].get_id(): " + std::to_string(rlbwt[i].get_id()) + "\t" +
                                    "blocked_id: " + std::to_string(blocked_id) + "\n"));
                }

                max_blocked_id = std::max(blocked_id, max_blocked_id);

                // always holds the last id seen for each character
                block_start_id[rlbwt[i].get_c()] = static_cast<uint32_t>(id - first_runs[rlbwt[i].get_c() + 1]);

                if (id - first_runs[rlbwt[i].get_c() + 1] > std::numeric_limits<uint32_t>::max()) {
                    throw std::runtime_error(ERROR_MSG("[compute blocked ids] The block_start_id does not fit in uint32_t\n" +
                                    "id - first_runs[rlbwt[i].get_c() + 1]: " + std::to_string(id - first_runs[rlbwt[i].get_c() + 1]) + "\n"));

                }
            }

        }

        if (small_block) {
            blocked_ids_computed = false;
            block_size = block_size / 2;
            max_allowed_blocked_id = ((max_allowed_blocked_id + 1) / 2) - 1;

            for (auto& block: id_blocks) {
                block.clear();
                block.shrink_to_fit();
            }
            id_blocks.clear();
            id_blocks.shrink_to_fit();

            block_start_id.clear();
            block_start_id.shrink_to_fit();

        } else {
            blocked_ids_computed = true;

            std::cerr << "Blocked index:\n";
            std::cerr << "max raw id: " << max_raw_id << "\t max blocked id: " << max_blocked_id << "\n";
            std::cerr << "max allowed blocked id: " << max_allowed_blocked_id << "\n";
            std::cerr << "block_size: " << block_size << "\n";
            std::cerr << "Maximum distance between the check points: " << max_diff << "\n\n";
        }
    }
}
#endif


uint64_t scan_count;
#if USE_NEXT_POINTERS
void MoveStructure::compute_nexts() {
    for (uint64_t i = rlbwt.size() - 1; i > 0; --i) {
        if (i % 1000000 == 0)
            std::cerr << i << "\r";

        char rlbwt_c = alphabet[rlbwt[i].get_c()];
        for (uint64_t j = 0; j < alphabet.size(); j++) {
            if (i == end_bwt_idx) {
                auto idx = reposition_up(i, alphabet[j], scan_count);
                end_bwt_idx_next_up[j] = (idx == r) ? std::numeric_limits<uint16_t>::max() : i - idx;
                idx = reposition_down(i, alphabet[j], scan_count);
                end_bwt_idx_next_down[j] = (idx == r) ? std::numeric_limits<uint16_t>::max() : idx - i;
                continue;
            }
            if (alphabet[j] != rlbwt_c) {
                auto alphabet_idx = alphamap_3[alphamap[rlbwt_c]][j];
                auto idx = reposition_up(i, alphabet[j], scan_count);
                if (idx == r) {
                    rlbwt[i].set_next_up(alphabet_idx, std::numeric_limits<uint16_t>::max());
                } else {
                    if (i - idx > std::numeric_limits<uint16_t>::max())
                        std::cerr << "Warning - reposition up " << i - idx << " does not fit in 16 bits.\n";
                    rlbwt[i].set_next_up(alphabet_idx, i - idx);
                }

                idx = reposition_down(i, alphabet[j], scan_count);
                if (idx == r) {
                    rlbwt[i].set_next_down(alphabet_idx, std::numeric_limits<uint16_t>::max());
                } else {
                    if (idx - i > std::numeric_limits<uint16_t>::max())
                        std::cerr << "Warning - reposition down " << idx - i << " does not fit in 16 bits.\n";
                    rlbwt[i].set_next_down(alphabet_idx, idx - i);
                }
            }
        }
    }
}
#endif

void MoveStructure::build_ftab() {
    size_t ftab_k = movi_options->get_ftab_k();

    uint64_t ftab_size = std::pow(4, ftab_k);
    std::cerr << "ftab_size: " << ftab_size*sizeof(MoveInterval)*std::pow(10, -6) << " MB \n";
    ftab.clear();
    ftab.resize(ftab_size);
    std::cerr << "Number of ftab entries (ftab-k=" << ftab_k << "): " << ftab_size << "\n";
    for (uint64_t i = 0; i < ftab_size; i++) {
        std::string kmer = number_to_kmer(i, 2*(ftab_k), alphabet, alphamap);
        uint64_t kmer_code = kmer_to_number(ftab_k, kmer, 0, alphamap);
        if (kmer_code != i) {
            std::cerr << kmer << " " << i << " " << kmer_code << "\n";
        }
        int32_t pos_on_kmer = ftab_k - 1;
        MoveInterval interval(
            first_runs[alphamap[kmer[pos_on_kmer]] + 1],
            first_offsets[alphamap[kmer[pos_on_kmer]] + 1],
            last_runs[alphamap[kmer[pos_on_kmer]] + 1],
            last_offsets[alphamap[kmer[pos_on_kmer]] + 1]
        );
        MoveInterval kmer_interval = backward_search(kmer, pos_on_kmer, interval, std::numeric_limits<int32_t>::max());
        uint64_t match_count = kmer_interval.count(rlbwt);
        if (match_count >= 0 and pos_on_kmer == 0) {
            // std::cerr << kmer << " " << match_count << " " << kmer_interval << "\n";
            // if (!(ftab[i] == kmer_interval)) {
            //     std::cerr << "ftab[i] is different from kmer_interval: " << i << "\n";
            //     std::cerr << ftab[i] << "\n" << kmer_interval << "\n";
            //     exit(0);
            // }
            ftab[i] = kmer_interval;
        } else {
            // std::cerr << kmer << " not found!\n";
            // if (!ftab[i].is_empty()) {
            //     std::cerr << "ftab[i] is non-empty and different from kmer_interval: " << i << "\n";
            //     std::cerr << ftab[i].is_empty() << "\n" << match_count << " " << (pos_on_kmer == 0) << "\n";
            //     exit(0);
            // }
            ftab[i].make_empty();
        }
        // if (pos_on_kmer != 0) pos_on_kmer += 1;
        // std::cerr << kmer << "\n";
        // std::cerr << kmer.length() - pos_on_kmer << "/" << kmer.length() << "\t" << match_count << "\n";
    }
}

// Finds all SA entries in O(n).
void MoveStructure::find_sampled_SA_entries() {
    uint64_t tot_len = 0;
    all_p.resize(r);
    std::cerr << "r: " << r << "\n";
    std::cerr << "Finding the BWT offset of run starts (the p array)..\n";
    for (uint64_t i = 0; i < r; i++) {
        if (i % 1000000 == 0)
            std::cerr << i << "\r";
        all_p[i] = tot_len;
        tot_len += rlbwt[i].get_n();
    }
    std::cerr << "tot_len: " << tot_len << "\n";

    // Create a sampled SA
    uint64_t SA_sample_rate = movi_options->get_SA_sample_rate();
    uint64_t SA_sample_size = tot_len / SA_sample_rate + 1;
    sampled_SA_entries.resize(SA_sample_size);

    uint64_t offset = 0;
    uint64_t index = 0;
    uint64_t SA_val = tot_len;
    std::cerr << "Finding the sampled SA entries..\n";
    for (uint64_t i = 0; i < tot_len; i++) {
        if (i % 1000000 == 0)
            std::cerr << i << "\r";
        SA_val--;
        uint64_t row_ind = all_p[index] + offset;
        if (row_ind % SA_sample_rate == 0) {
            sampled_SA_entries[row_ind / SA_sample_rate] = SA_val;
        }
        LF_move(offset, index);
    }
    std::cerr << INFO_MSG("Finished building the sampled SA entries.");
}