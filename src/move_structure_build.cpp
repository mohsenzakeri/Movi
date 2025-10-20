#include "move_structure.hpp"

void MoveStructure::compute_number_of_build_steps() {
    current_build_step = 1;
    total_build_steps = 4;
#if USE_THRESHOLDS
    total_build_steps += 1;
#endif
#if BLOCKED_MODES
    total_build_steps += 1;
#endif
#if USE_NEXT_POINTERS
    total_build_steps += 1;
#endif
}

void MoveStructure::build() {
    if (movi_options->is_thresholds()) {
        std::string thr_filename = movi_options->get_ref_file() + std::string(".thr_pos");
        read_thresholds(thr_filename, thresholds);
    }

    INFO_MSG("Building starts...");
    compute_number_of_build_steps();

    length = compute_length_from_bwt();
    if (movi_options->is_verbose()) {
        INFO_MSG("The total length of the BWT (n): " + std::to_string(length));
    }

    initialize_bits();

#if SPLIT_THRESHOLDS_TRUE
    // If thresholds splitting is enabled, bits will track the splitting points
    if (movi_options->is_thresholds()) {
        fill_bits_by_thresholds();
    }
#endif

    detect_move_row_boundaries();

    find_run_heads_information();

    if (movi_options->use_separators() and alphabet.size() == 5) {
        // separators are present
        if (alphabet[0] != SEPARATOR) {
            throw std::runtime_error(ERROR_MSG("[build] The separator is not correct: " + std::to_string(alphabet[0]) + " != " + std::to_string(SEPARATOR)));
        }
    } else if (alphabet.size() > 4) {
        WARNING_MSG("There are more than 4 characters in the reference with no separators.");
        for (uint64_t i = 0; i < alphabet.size(); i++) {
            INFO_MSG("(" + std::to_string(i) + ")" + "\t" + std::to_string(alphabet[i]) + "\t" + std::to_string(counts[i]));
        }
    }

    build_move_rows();

#if USE_THRESHOLDS
    // compute the thresholds
    compute_thresholds();
#endif

#if USE_NEXT_POINTERS
    if (constant) {
        if (movi_options->is_verbose()) {
            INFO_MSG("Computing the next ups and downs..");
        }
        compute_nexts();
    }
#endif
    SUCCESS_MSG("The Movi index building is done.");
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

        if (i % 1000000 == 0 or i == r - 1) {
            print_progress_bar(i, r - 1, "Finding the run heads information", current_build_step, total_build_steps);
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
    PROGRESS_MSG("Successfully processed " + std::to_string(r) + " runs");
    current_build_step++;
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
        INFO_MSG("Number of BWT runs: " + std::to_string(original_r));
#if TALLY_MODES
        INFO_MSG("Checkpoint for the id field is stored every " + std::to_string(movi_options->get_tally_checkpoints()) + " move rows");
#endif

        original_lens.resize(original_r);

        std::ifstream len_file(bwt_filename + ".len", std::ios::in | std::ios::binary);
        if (!len_file.good()) {
            throw std::runtime_error(ERROR_MSG("[build] Failed to open the len file: " + bwt_filename + ".len"));
        }

        for (uint64_t i = 0; i < original_r; i++) {
            if (i % 1000000 == 0 or i == original_r - 1) {
                print_progress_bar(i, original_r - 1, "Computing total length of the BWT", current_build_step, total_build_steps);
            }
            size_t len = 0;
            len_file.read(reinterpret_cast<char*>(&len), 5);
            original_lens[i] = len;
            total_length += len;
        }
        PROGRESS_MSG("Successfully computed total BWT length: " + std::to_string(total_length));
        current_build_step++;
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
            INFO_MSG("total_length: " + std::to_string(total_length));

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
        INFO_MSG("bits.size after loading the d_col file: " + std::to_string(bits.size()));
        rbits = sdsl::rank_support_v<>(&bits);
        INFO_MSG("The main bit vector (bits) is loaded from the d_col file.");
    } else {
        bits = sdsl::bit_vector(length + 1, 0);
        bits[0] = 1;
        INFO_MSG("The main bit vector (bits) is initialized.");
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
            if (i % 1000000 == 0 or i == original_r - 1) {
                print_progress_bar(i, original_r - 1, "Detecting the move row boundaries", current_build_step, total_build_steps);
            }

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
        PROGRESS_MSG("Successfully detected " + std::to_string(r) + " move row boundaries");
        current_build_step++;

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

        INFO_MSG("Reading over the uncompressed BWT...");
        while (next_char != EOF) { // && current_char != 10
            run_length += 1;
            if (next_char != END_CHARACTER) {
                all_possible_chars[next_char] += 1;
            }

            if (original_r % 1000000 == 0) {
                print_progress_bar(scanned_bwt_length, length - 1, "Detecting the move row boundaries", current_build_step, total_build_steps);
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
        PROGRESS_MSG("Successfully detected " + std::to_string(r) + " move row boundaries");

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
        INFO_MSG("Number of runs added by thresholds splitting: " + std::to_string(split_by_thresholds));
    }
    INFO_MSG("Main statistics of the index:");
    INFO_MSG("\tn: " + std::to_string(length));
    INFO_MSG("\tr: " + std::to_string(r));
    INFO_MSG("\tn/r:\t" + std::to_string(static_cast<double>(length)/r));
    INFO_MSG("\toriginal_r: " + std::to_string(original_r));
}

void MoveStructure::build_alphabet(std::vector<uint64_t>& all_possible_chars) {
    uint64_t alphabet_index = 0;
    for (uint64_t i = 1; i < all_possible_chars.size(); i++) {
        if (all_possible_chars[i] != 0) {
            auto current_char = static_cast<unsigned char>(i);
            if (movi_options->is_verbose()) {
                INFO_MSG("(" + std::to_string(i) + ")"
                          + "\t" + std::to_string(current_char)
                          + "\t index:" + std::to_string(alphabet_index)
                          + "\t count:" + std::to_string(all_possible_chars[i]));
            }

            alphabet.push_back(current_char);
            counts.push_back(all_possible_chars[i]);
            alphamap[i] = alphabet_index;
            alphabet_index += 1;
        }
    }
    INFO_MSG("All the characters are indexed.");
}

void MoveStructure::build_move_rows() {
    if (movi_options->is_debug()) {
        DEBUG_MSG("Starting build_move_rows()");
    }

    if (movi_options->is_verbose() and bits.size() < 1000) {
        std::string bits_string = "";
        for (auto bit : bits)
            bits_string += (bit ? '1' : '0');
        INFO_MSG("bits: " + bits_string);
    }

    if (movi_options->is_debug()) {
        DEBUG_MSG("Building rank support vector for " + std::to_string(bits.size()) + " bits");
    }

    if (!nt_splitting) // The rank vector is already built if it was the splitting mode
        rbits = sdsl::rank_support_v<>(&bits);

    if (movi_options->is_debug()) {
        DEBUG_MSG("Rank support vector built successfully");
    }


    rlbwt.resize(r);

    if (movi_options->is_verbose() ) {
        INFO_MSG("bits.size(): " + std::to_string(bits.size()));
        INFO_MSG("rank_support_v<>(&bits)(bits.size()): "
                  + std::to_string(sdsl::rank_support_v<>(&bits)(bits.size())));

        INFO_MSG("r: " + std::to_string(r) + "\t"
                  + "rlbwt.size(): " + std::to_string(rlbwt.size()) + "\t"
                  + "all_p.size(): " + std::to_string(all_p.size()) + "\t"
                  + "end_bwt_idx: " + std::to_string(end_bwt_idx));
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
    if (movi_options->is_debug()) {
        DEBUG_MSG("Starting main loop for " + std::to_string(r) + " rows");
    }

    for (uint64_t r_idx = 0; r_idx < r; r_idx++) {
        if (r_idx % 1000000 == 0 or r_idx == r - 1) {
            print_progress_bar(r_idx, r - 1, "Building Movi rows", current_build_step, total_build_steps);
        }

        if (movi_options->is_debug()) {
            DEBUG_MSG("Processing row (r_idx=" + std::to_string(r_idx) + ")");
        }

        uint64_t lf  = 0;
        if (r_idx != end_bwt_idx) {
            uint64_t alphabet_index = alphamap[static_cast<uint64_t>(heads[r_idx])];
            // lf = LF(all_p[r_idx], alphabet_index);
            lf = LF_heads(r_idx, alphabet_index);
            if (movi_options->is_debug()) {
                DEBUG_MSG("LF_heads: " + std::to_string(lf));
            }
        } else {
            lf = 0;
        }

        uint64_t pp_id = rbits(lf) - 1;
        if (movi_options->is_debug()) {
            DEBUG_MSG("pp_id (id of lf of head of the run): " + std::to_string(pp_id));
        }

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

        if (movi_options->is_debug()) {
            DEBUG_MSG("r_idx: " + std::to_string(r_idx)
                        + " len: " + std::to_string(len)
                        + " lf: " + std::to_string(lf)
                        + " offset: " + std::to_string(offset)
                        + " pp_id: " + std::to_string(pp_id));
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
        if (movi_options->is_debug()) {
            DEBUG_MSG("rlbwt[r_idx].id: " + std::to_string(rlbwt[r_idx].get_id()));
        }
#endif

        // To take care of cases where length of the run
        // does not fit in uint16_t
#if SPLIT_MAX_RUN
        if (offset > MAX_RUN_LENGTH or len > MAX_RUN_LENGTH) {
            // Should not get here in the regular mode: MODE = 3
            throw std::runtime_error(ERROR_MSG("[build rows] The length or the offset are too large.\n" +
                                                "offset: " + std::to_string(offset) + "\tlength: " + std::to_string(length)));
        }
#endif
#if SPLIT_THRESHOLDS_FALSE
        if (movi_options->is_debug()) {
            DEBUG_MSG("rlbwt[r_idx].len: " + std::to_string(rlbwt[r_idx].get_n()));
        }

        if (len > MAX_RUN_LENGTH) {
            throw std::runtime_error(ERROR_MSG("[build rows] This shouldn't happen anymore, because the run splitting should be applied in every mode now."));
            n_overflow.push_back(len);
            if (n_overflow.size() - 1 >= MAX_RUN_LENGTH) {
                throw std::runtime_error(ERROR_MSG("[build rows] The number of runs with overflow n is beyond " + std::to_string(MAX_RUN_LENGTH) + "! " + std::to_string(n_overflow.size() - 1)));
            }
            rlbwt[r_idx].set_n(n_overflow.size() - 1);
            rlbwt[r_idx].set_overflow_n();
        }

        if (movi_options->is_debug()) {
            DEBUG_MSG("rlbwt[r_idx].is_overflow_n: " + std::to_string(rlbwt[r_idx].is_overflow_n()));
        }

        if (movi_options->is_debug()) {
            DEBUG_MSG("rlbwt[r_idx].offset: " + std::to_string(rlbwt[r_idx].get_offset()));
            DEBUG_MSG("offset: " + std::to_string(offset));
            DEBUG_MSG("MAX_RUN_LENGTH: " + std::to_string(MAX_RUN_LENGTH));
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

        if (movi_options->is_debug()) {
            DEBUG_MSG("offset: " + std::to_string(offset));
            DEBUG_MSG("rlbwt[r_idx].get_offset: " + std::to_string(rlbwt[r_idx].get_offset()));
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

        if (movi_options->is_debug()) {
            DEBUG_MSG("rlbwt[r_idx].set_c: " + std::to_string(rlbwt[r_idx].get_c()));
        }
    }
    PROGRESS_MSG("Successfully built " + std::to_string(r) + " move structure rows");
    current_build_step++;
#if TALLY_MODES
    // Set the last tally_id using the current_tally_ids
    for (uint32_t alphabet_ind = 0; alphabet_ind < alphabet.size(); alphabet_ind++) {
        tally_ids[alphabet_ind][tally_ids[alphabet_ind].size()-1].set_value(current_tally_ids[alphabet_ind]);
    }
#endif
    if (movi_options->is_verbose()) {
        INFO_MSG("Max run length (after splitting if enabled): " + std::to_string(max_len));
    }

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
            INFO_MSG("\t" + std::to_string(i) + " char_count: " + std::to_string(char_count));
        auto occ_rank = rbits(char_count);
        // occ_rank is the number of set bits until and including the last row with character i - 1
        if (movi_options->is_verbose())
            INFO_MSG("\t" + std::to_string(i) + " occ_rank: " + std::to_string(occ_rank));
        // The index of the run is number of set bits - 1 (because we count from 0):
        last_runs.push_back(static_cast<uint64_t>(occ_rank - 1));
        if (movi_options->is_verbose())
            INFO_MSG("\t" + std::to_string(char_count) + " " + std::to_string(all_p[last_runs.back()])
                        + " " + std::to_string(get_n(last_runs.back())));
        last_offsets.push_back(char_count - all_p[last_runs.back()] - 1);
    }
    if (movi_options->is_verbose()) {
        for (uint64_t i = 0; i < first_runs.size(); i++) {
            INFO_MSG("<--- " + std::to_string(first_runs[i]) + ":" + std::to_string(first_offsets[i]) + "\t---\t"
                        + std::to_string(last_runs[i]) + ":" + std::to_string(last_offsets[i]) + "--->" );
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
    INFO_MSG("The bits vector is updated by thresholds.");
}


#if USE_THRESHOLDS
void MoveStructure::debug_threshold_calculation(uint64_t i, uint64_t thr_i, uint64_t j, char rlbwt_c,
                                               const std::vector<uint64_t>& alphabet_thresholds) {
    if (!movi_options->is_debug() || i < rlbwt.size() - 10) return;

    DEBUG_MSG("i: " + std::to_string(i) + "\n"
              + "rlbwt[i].get_offset(): " + std::to_string(get_offset(i)) + "\n"
              + "get_n(i): " + std::to_string(get_n(i)) + "\n"
              + "thresholds[thr_i]: " + std::to_string(thresholds[thr_i]) + " "
              + "rlbwt_c: " + std::to_string(rlbwt_c));

    if (j < alphabet.size()) {
        DEBUG_MSG("\t j: \t" + std::to_string(j) + " "
                  + "alphabet[j]: " + std::to_string(alphabet[j]) + "  "
                  + "alphamap_3[alphamap[rlbwt_c]][j]: " + std::to_string(alphamap_3[alphamap[rlbwt_c]][j]) + " "
                  + "alphabet_thresholds[j]: " + std::to_string(alphabet_thresholds[j]) + " "
#if SPLIT_THRESHOLDS_FALSE
                  + "rlbwt[i].thresholds[j]:" + std::to_string(get_rlbwt_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j]))
#endif
#if SPLIT_THRESHOLDS_TRUE
                  + "rlbwt[i].thresholds[j]:" + std::to_string(get_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j]))
#endif
                  + "\n");
    }
}

void MoveStructure::set_threshold_for_one_character(uint64_t i, uint16_t threshold_index, uint64_t value,
                                                    uint16_t value_split, std::vector<uint64_t>& current_thresholds) {
    if (i == end_bwt_idx) {
        end_bwt_idx_thresholds[threshold_index] = value;
        return;
    }

    if (use_separator() and alphabet[rlbwt[i].get_c()] == SEPARATOR) {
        separators_thresholds.back().values[threshold_index] = value;
        return;
    }
#if SPLIT_THRESHOLDS_FALSE
    if (value > std::numeric_limits<uint16_t>::max()) {
        rlbwt[i].set_overflow_thresholds();
    } else {
        set_rlbwt_thresholds(i, threshold_index, static_cast<uint16_t>(value));
    }

    if (i != end_bwt_idx) {
        current_thresholds[threshold_index] = value;
    }
#endif
#if SPLIT_THRESHOLDS_TRUE
    if (value_split == std::numeric_limits<uint16_t>::max()) {
        throw std::runtime_error(ERROR_MSG("[set_thresholds_for_dna_rows]" +
                                           " This should never happen, since the runs are split at threshold boundaries.\n" +
                                           " threshold_index: " + std::to_string(threshold_index) + " i:" + std::to_string(i) +
                                           " value:" + std::to_string(value) + " all_p[i]:" + std::to_string(all_p[i])));
    }
    rlbwt[i].set_threshold(threshold_index, value_split);
#endif
}

void MoveStructure::compute_thresholds() {
    // initialize the start threshold at the last row
    std::vector<uint64_t> alphabet_thresholds(alphabet.size(), length);

    uint64_t thr_i = original_r - 1;

    // Iterate over all the move rows in reverse order
    for (uint64_t i = rlbwt.size() - 1; i > 0; --i) {
        if (i % 1000000 == 0 or i == rlbwt.size() - 1 or i == 1) {
            print_progress_bar(rlbwt.size() - 1 - i, rlbwt.size() - 2, "Updating the thresholds for move rows", current_build_step, total_build_steps);
        }

        if (thr_i >= thresholds.size()) {
            throw std::runtime_error(ERROR_MSG("[compute thresholds] thr_i = " + std::to_string(thr_i) + " is out of bound:\n" +
                                                " thresholds.size = " + std::to_string(thresholds.size())));
        }

        char rlbwt_c = alphabet[rlbwt[i].get_c()];
        // TODO: we are probably not looking at the right character here
        // char rlbwt_c = get_char(i);

        if (use_separator() and rlbwt_c == SEPARATOR) {
            separators_thresholds.push_back(ThresholdsRow{0, 0, 0, 0});
            separators_thresholds_map[i] = separators_thresholds.size() - 1;
        }

        // Used for non-split thresholds
        std::vector<uint64_t> current_thresholds;
        current_thresholds.resize(alphabet.size() - 1);
        if (use_separator()) {
            // We don't need to store a threshold for the separator character
            current_thresholds.resize(alphabet.size() - 2);
        }

        // For other move rows, store a threshold for each character other than the row's character
        for (uint64_t j = 0; j < alphabet.size(); j++) {
            if (alphabet[j] == rlbwt_c) {
                // We don't need to store a threshold for the row's character
                alphabet_thresholds[j] = thresholds[thr_i];
            } else {
                uint16_t threshold_index = i == end_bwt_idx ? j : static_cast<uint16_t>(alphamap_3[alphamap[rlbwt_c]][j]);
                if (use_separator()) {
                    if (j == 0) {
                        // No threshold for the separator character is stored
                        continue;
                    }
                    // We basically assume the main charaqcters of the alphabet start at index 1 so on
                    // Because the separator character is at the index 0
                    threshold_index = (i == end_bwt_idx or rlbwt_c == SEPARATOR) ?
                                            j - 1 : static_cast<uint16_t>(alphamap_3[alphamap[rlbwt_c] - 1][j - 1]);
                }

                uint64_t current_threshold = alphabet_thresholds[j];

                if (threshold_index >= alphabet.size() - 1 && i != end_bwt_idx && rlbwt_c != SEPARATOR) {
                    throw std::runtime_error(ERROR_MSG("[compute thresholds] alphamap_3 is not working in general:\n"
                                                       "alphabet.size() - 1 = " + std::to_string(alphabet.size() - 1) + "\n"
                                                       "alphamap_3[alphamap[rlbwt_c]][j] = " + std::to_string(threshold_index)));
                }

                uint64_t threshold_value;
                uint16_t threshold_value_split;
                if (current_threshold >= all_p[i] + get_n(i)) {
                    threshold_value = get_n(i);
                    threshold_value_split = 1;
                } else if (current_threshold <= all_p[i]) {
                    threshold_value = 0;
                    threshold_value_split = 0;
                } else {
                    threshold_value = current_threshold - all_p[i];
                    threshold_value_split = std::numeric_limits<uint16_t>::max();
                }

                set_threshold_for_one_character(i, threshold_index, threshold_value, threshold_value_split, current_thresholds);

                // printing the values for last 10 runs to debug
                debug_threshold_calculation(i, thr_i, j, rlbwt_c, alphabet_thresholds);
            }
        }

        if (i > 0 && (rlbwt[i].get_c() != rlbwt[i - 1].get_c()  || i == end_bwt_idx || i-1 == end_bwt_idx)) {
            thr_i--;
        }

#if SPLIT_THRESHOLDS_FALSE
        // Handle overflow thresholds
        if (rlbwt[i].is_overflow_thresholds()) {
            if (thresholds_overflow.size() >= std::numeric_limits<uint16_t>::max()) {
                throw std::runtime_error(ERROR_MSG("[compute thresholds] " +
                    "Undefined behaviour: the number of move rows with overflow thresholds is beyond uint16_t." +
                    std::to_string(thresholds_overflow.size())));
            }
            for (uint16_t j = 0; j < alphabet.size() - 1; j++) {
                set_rlbwt_thresholds(i, j, thresholds_overflow.size());
            }
            thresholds_overflow.push_back(current_thresholds);
        }
#endif
    }

    // since the thresholds for the first run was not calculated in the for
    if (!use_separator()) {
        for (uint16_t j = 0; j < alphabet.size() - 1; j++) {
#if SPLIT_THRESHOLDS_FALSE
            set_rlbwt_thresholds(0, j, 0);
#endif
#if SPLIT_THRESHOLDS_TRUE
            rlbwt[0].set_threshold(j, 0);
#endif
        }
    } else {
        if (alphabet[rlbwt[0].get_c()] == SEPARATOR) {
            separators_thresholds.push_back(ThresholdsRow{0, 0, 0, 0});
            separators_thresholds_map[0] = separators_thresholds.size() - 1;
        } else {
            for (uint16_t j = 0; j < alphabet.size() - 2; j++) {
#if SPLIT_THRESHOLDS_FALSE
                set_rlbwt_thresholds(0, j, 0);
#endif
#if SPLIT_THRESHOLDS_TRUE
                rlbwt[0].set_threshold(j, 0);
#endif
            }
        }
    }

    PROGRESS_MSG("Successfully updated the thresholds in all the rows.");
    current_build_step++;
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
            if (i % 1000000 == 0 or i == rlbwt.size() - 1)
                print_progress_bar(i, rlbwt.size() - 1, "Finding the blocked deltas", current_build_step, total_build_steps);

            uint64_t block_boundary = (block_size) * block_count; // BLOCK_SIZE = 1 << 22 - 1

            if (i >= block_boundary) {

                if (movi_options->is_verbose())
                    INFO_MSG("\n" + std::to_string(block_count) + "\t");

                // A new block is being initiated
                block_count += 1;

                // Store the largest (last) observed id for each character in the last block as a check point for the new block
                for (uint64_t j = 0; j < alphabet.size(); j++) {

                    id_blocks[j].push_back(block_start_id[j]);

                    if (movi_options->is_verbose())
                        INFO_MSG(std::to_string(id_blocks[j][block_count - 1]) + "\t");

                    if (block_count > 1)
                        max_diff = std::max(max_diff, id_blocks[j][block_count - 1] - id_blocks[j][block_count - 2] );
                }

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

                    INFO_MSG("The number of bits in the move rows are not enough for storing the blocked_id.");
                    INFO_MSG("adjusted_id: " + std::to_string(adjusted_id) + " id: " + std::to_string(id)
                              + " first_runs[rlbwt[i].get_c() + 1]: " + std::to_string(first_runs[rlbwt[i].get_c() + 1]));
                    INFO_MSG("id_blocks[rlbwt[i].get_c()][block_number]: " + std::to_string(id_blocks[rlbwt[i].get_c()][block_number]));
                    INFO_MSG("rlbwt[i].get_c(): " + std::to_string(static_cast<uint32_t>(rlbwt[i].get_c())));
                    INFO_MSG("blocked_id: " + std::to_string(blocked_id) + " max_blocked_id: " + std::to_string(max_blocked_id));
                    INFO_MSG("max_allowed_blocked_id: " + std::to_string(max_allowed_blocked_id));
                    INFO_MSG("block_count: " + std::to_string(block_count) + " block_number: " + std::to_string(block_number));
                    INFO_MSG("i: " + std::to_string(i) + " block_size: " + std::to_string(block_size));

                    small_block = true;
                    break;
                }

                rlbwt[i].set_id(blocked_id);

                // Check if the blocked_id is stored correctly in the run
                if (rlbwt[i].get_id() != blocked_id) {
                    throw std::runtime_error(ERROR_MSG("[compute blocked ids] The id was not stored in the move row correctly.\n" +
                                    "rlbwt[i].get_id(): " + std::to_string(rlbwt[i].get_id()) + "\t" +
                                    "blocked_id: " + std::to_string(blocked_id)));
                }

                max_blocked_id = std::max(blocked_id, max_blocked_id);

                // always holds the last id seen for each character
                block_start_id[rlbwt[i].get_c()] = static_cast<uint32_t>(id - first_runs[rlbwt[i].get_c() + 1]);

                if (id - first_runs[rlbwt[i].get_c() + 1] > std::numeric_limits<uint32_t>::max()) {
                    throw std::runtime_error(ERROR_MSG("[compute blocked ids] The block_start_id does not fit in uint32_t\n" +
                                    "id - first_runs[rlbwt[i].get_c() + 1]: " + std::to_string(id - first_runs[rlbwt[i].get_c() + 1])));

                }
            }

        }
        PROGRESS_MSG("Successfully computed blocked deltas for " + std::to_string(rlbwt.size()) + " move rows");
        current_build_step++;

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

            INFO_MSG("Blocked index statistics:");
            INFO_MSG("\tmax raw id: " + std::to_string(max_raw_id));
            INFO_MSG("\tmax blocked id: " + std::to_string(max_blocked_id));
            INFO_MSG("\tmax allowed blocked id: " + std::to_string(max_allowed_blocked_id));
            INFO_MSG("\tblock_size: " + std::to_string(block_size));
            INFO_MSG("\tMaximum distance between the check points: " + std::to_string(max_diff));
        }
    }
}
#endif


uint64_t scan_count;
#if USE_NEXT_POINTERS
void MoveStructure::compute_nexts() {
    for (uint64_t i = rlbwt.size() - 1; i > 0; --i) {
        if (i == rlbwt.size() - 1 or i % 1000000 == 0 or i == 1) {
            print_progress_bar(rlbwt.size() - 1 - i, rlbwt.size() - 2, "Computing the next ups and downs", current_build_step, total_build_steps);
        }
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
                        WARNING_MSG("reposition up " + std::to_string(i - idx) + " does not fit in 16 bits.");
                    rlbwt[i].set_next_up(alphabet_idx, i - idx);
                }

                idx = reposition_down(i, alphabet[j], scan_count);
                if (idx == r) {
                    rlbwt[i].set_next_down(alphabet_idx, std::numeric_limits<uint16_t>::max());
                } else {
                    if (idx - i > std::numeric_limits<uint16_t>::max())
                        WARNING_MSG("reposition down " + std::to_string(idx - i) + " does not fit in 16 bits.");
                    rlbwt[i].set_next_down(alphabet_idx, idx - i);
                }
            }
        }
    }
    PROGRESS_MSG("Successfully computed next ups and downs for " + std::to_string(rlbwt.size()) + " move rows");
    current_build_step++;
}
#endif

void MoveStructure::build_ftab() {
    size_t ftab_k = movi_options->get_ftab_k();

    uint64_t ftab_size = std::pow(4, ftab_k);
    INFO_MSG("ftab_size: " + std::to_string(ftab_size*sizeof(MoveInterval)*std::pow(10, -6)) + " MB");
    ftab.clear();
    ftab.resize(ftab_size);
    INFO_MSG("Number of ftab entries (ftab-k=" + std::to_string(ftab_k) + "): " + std::to_string(ftab_size));
    for (uint64_t i = 0; i < ftab_size; i++) {
        std::string kmer = number_to_kmer(i, 2*(ftab_k), alphabet, alphamap);
        uint64_t kmer_code = kmer_to_number(ftab_k, kmer, 0, alphamap);
        if (kmer_code != i) {
            INFO_MSG(kmer + " " + std::to_string(i) + " " + std::to_string(kmer_code));
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
            if (movi_options->is_debug()) {
                DEBUG_MSG(kmer + " " + std::to_string(match_count) + " " + kmer_interval.to_string());
            }
            // if (!(ftab[i] == kmer_interval)) {
            //     std::cerr << "ftab[i] is different from kmer_interval: " << i << "\n";
            //     std::cerr << ftab[i] << "\n" << kmer_interval << "\n";
            //     exit(0);
            // }
            ftab[i] = kmer_interval;
        } else {
            if (movi_options->is_debug()) {
                DEBUG_MSG(kmer + " not found!");
            }
            // if (!ftab[i].is_empty()) {
            //     std::cerr << "ftab[i] is non-empty and different from kmer_interval: " << i << "\n";
            //     std::cerr << ftab[i].is_empty() << "\n" << match_count << " " << (pos_on_kmer == 0) << "\n";
            //     exit(0);
            // }
            ftab[i].make_empty();
        }
        // if (pos_on_kmer != 0) pos_on_kmer += 1;
        if (movi_options->is_debug()) {
            DEBUG_MSG(kmer);
            DEBUG_MSG(kmer.length() - pos_on_kmer << "/" << kmer.length() << "\t" << match_count);
        }
    }
}

// Finds all SA entries in O(n).
void MoveStructure::find_sampled_SA_entries() {
    uint64_t tot_len = 0;
    all_p.resize(r);

    for (uint64_t i = 0; i < r; i++) {
        if (i % 1000000 == 0 or i == r - 1) {
            print_progress_bar(i, r - 1, "Finding the BWT offset of run heads");
        }
        all_p[i] = tot_len;
        tot_len += rlbwt[i].get_n();
    }
    PROGRESS_MSG("Successfully found BWT offsets for " + std::to_string(r) + " runs");

    if (movi_options->is_verbose()) {
        INFO_MSG("tot_len: " + std::to_string(tot_len));
    }

    // Create a sampled SA
    uint64_t SA_sample_rate = movi_options->get_SA_sample_rate();
    uint64_t SA_sample_size = tot_len / SA_sample_rate + 1;
    sampled_SA_entries.resize(SA_sample_size);

    uint64_t offset = 0;
    uint64_t index = 0;
    uint64_t SA_val = tot_len;
    INFO_MSG("Finding the sampled SA entries..");
    for (uint64_t i = 0; i < tot_len; i++) {
        if (i % 1000000 == 0 or i == tot_len - 1) {
            print_progress_bar(i, tot_len - 1, "Iterating over all of the BWT offsets");
        }
        SA_val--;
        uint64_t row_ind = all_p[index] + offset;
        if (row_ind % SA_sample_rate == 0) {
            sampled_SA_entries[row_ind / SA_sample_rate] = SA_val;
        }
        LF_move(offset, index);
    }
    PROGRESS_MSG("Successfully built " + std::to_string(tot_len) + " sampled SA entries");
}