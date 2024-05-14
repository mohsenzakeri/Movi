#include <sys/stat.h> 

#include "move_structure.hpp"

uint32_t alphamap_3[4][4] = {{3, 0, 1, 2},
                             {0, 3, 1, 2},
                             {0, 1, 3, 2},
                             {0, 1, 2, 3}};

void read_thresholds(std::string tmp_filename, sdsl::int_vector<>& thresholds) {
    int log_n = 100;

    struct stat filestat;
    FILE *fd;

    if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
        std::cerr <<("open() file " + tmp_filename + " failed");

    int fn = fileno(fd);
    if (fstat(fn, &filestat) < 0)
        std::cerr <<("stat() file " + tmp_filename + " failed");

    if (filestat.st_size % THRBYTES != 0)
        std::cerr <<("invilid file " + tmp_filename);

    size_t length_thr = filestat.st_size / THRBYTES;
    size_t threshold = 0;

    thresholds = sdsl::int_vector<>(length_thr, 0, log_n);

    size_t i = 0;
    for (i = 0; i < length_thr; ++i) {
        size_t threshold = 0;
        if ((fread(&threshold, THRBYTES, 1, fd)) != 1)
            std::cerr <<("fread() file " + tmp_filename + " failed");
        thresholds[i] = threshold;
    }
    std::cerr << "Finished reading " << i << " thresholds.\n";
}

MoveStructure::MoveStructure(MoviOptions* movi_options_) {
    movi_options = movi_options_;
    onebit = false;
    verbose = movi_options->is_verbose();
    logs = movi_options->is_logs();
}

MoveStructure::MoveStructure(std::string input_file_, bool onebit_, bool verbose_, bool logs_, uint16_t splitting_, bool constant_) {
    onebit = onebit_;
    verbose = verbose_;
    logs = logs_;
    splitting = splitting_;
    constant = constant_;

    if (!check_mode()) {
        std::cerr << "Your settings: \n"
                    << "onebit: " << onebit << "\n"
                    << "constant: " << constant << "\n"
                    << "splitting " << splitting << "\n";
        exit(0);
    }

    reconstructed = false;
    input_file = input_file_;

    std::string bwt_filename = input_file + std::string(".bwt");
    std::ifstream bwt_file(bwt_filename);

    std::string thr_filename = input_file + std::string(".thr_pos");
    read_thresholds(thr_filename, thresholds);

    build(bwt_file);
}

/*char MoveStructure::compute_char(uint64_t idx) {
    if (verbose) {
        std::cerr << idx << " bit1_begin: " << bit1_begin << "\n";
        std::cerr << idx << " bit1_after_eof: " << bit1_after_eof << "\n";
        std::cerr << "alphabet.size(): " << alphabet.size() << "\n";
    }
    char row_c;
    if (idx < eof_row) {
        row_c = idx%2 == 0 ? alphabet[bit1_begin] : alphabet[(1 + bit1_begin) % 2];
    } else if (idx == eof_row) {
        row_c = END_CHARACTER;
    } else{
        row_c = (idx - eof_row + 1)%2 == 0 ? alphabet[bit1_after_eof] : alphabet[(1 + bit1_after_eof) % 2];
    }
    return row_c;
}*/

std::string MoveStructure::index_type() {
    if (!onebit and !constant and splitting == 0) {
        return "default";
    } else if (constant and splitting != 0 and !onebit) {
        return "constant";
    }else if (splitting != 0 and !onebit) {
        return "splitting";
    } else if (onebit and splitting == 0) {
        return "onebit_nosplitting";
    } else if (onebit and splitting != 0) {
        return "onebit";
    } else {
        std::cerr << "Mode is not defined! \n";
        std::cerr << "onebit: " << onebit
                  << "\nconstant:" << constant
                  << "\nsplitting:" << splitting
                  << "\n";
    }
    exit(0);
}

bool MoveStructure::check_mode() {
#if MODE == 0
    if (onebit || constant) {
        std::cerr << "MODE is set to be 0: regular!\n";
        return false;
    }
#endif

#if MODE == 1
    if (onebit || !constant || !splitting) {
        std::cerr << "MODE is set to be 1: constant!\n";
        return false;
    }
#endif

#if MODE == 2
    if (!onebit) {
        std::cerr << "MODE is set to be 2: one-bit!\n";
        return false;
    }
#endif
    return true;
}

uint64_t MoveStructure::LF(uint64_t row_number) {
    uint64_t lf = 0;
    uint64_t alphabet_index = alphamap[static_cast<uint64_t>(bwt_string[row_number])];
    lf += 1;
    for (uint64_t i = 0; i < alphabet_index; i++) {
        lf += counts[i];
    }
    auto& occ_rank = *occs_rank[alphabet_index];
    lf += static_cast<uint64_t>(occ_rank(row_number));
    return lf;
}

/*std::string MoveStructure::reconstruct() {
    if (bwt_string == "") {
        for (uint32_t i = 0; i < r; i++) {
            for (uint32_t j = 0; j < get_n(i); j++)
                bwt_string += rlbwt[i].get_c();
        }
    }
    if (!reconstructed) {
        orig_string = "";
        orig_string += static_cast<char>(END_CHARACTER);
        for (uint64_t bwt_row = 0; bwt_row != end_bwt_row; bwt_row = LF(bwt_row)) {
            orig_string = bwt_string[bwt_row] + orig_string;
            
        }
        reconstructed = true;
    }
    return orig_string;
}*/

/*uint64_t MoveStructure::naive_lcp(uint64_t row1, uint64_t row2) {
    if (row1 >= length or row2 >= length)
        return 0;

    if (!reconstructed)
        orig_string = reconstruct();
    uint64_t lcp = 0;
    uint64_t sa1 = naive_sa(row1);
    uint64_t sa2 = naive_sa(row2);

    while (orig_string[sa1 + lcp] == orig_string[sa2 + lcp]) {
        lcp += 1;
    }

    return lcp;
}*/

/*uint64_t MoveStructure::naive_sa(uint64_t bwt_row) {
    uint64_t sa = 0;
    for (; bwt_row != end_bwt_row; bwt_row = LF(bwt_row)) {
        sa += 1;
    }
    return sa;
}*/

uint32_t MoveStructure::compute_index(char row_char, char lookup_char) {
    uint32_t alpha_index = alphamap[lookup_char];
    return alpha_index;
    /*if (lookup_char < row_char)
        return alpha_index;
    else
        return alpha_index+1;*/
}

uint16_t MoveStructure::LF_move(uint64_t& offset, uint64_t& i) {
    if (verbose) {
        std::cerr << "\t in LF:\n";
        std::cerr << "\t \t i: " << i << " offset: " << offset << "\n";
    }
    auto& row = rlbwt[i];
    auto idx = row.get_id();
    offset = get_offset(i) + offset;
    uint16_t ff_count = 0;
    if (verbose) {
        std::cerr << "\t \t i: " << i << " offset: " << offset << " idx: " << idx << "\n";
    }

    if (idx < r - 1 && offset >= get_n(idx)) {
        uint64_t idx_ = fast_forward(offset, idx, 0);
        idx += idx_;
        if (idx_ >= std::numeric_limits<uint16_t>::max()) {
            std::cerr << "Number of fast forwards for a query was greater than 2^16: " << idx_ << "\n";
            exit(0);
        }
        ff_count = static_cast<uint16_t>(idx_);
    }

    if (verbose) {
        std::cerr << "\t \t after fast forward:\n";
        std::cerr << "\t \t i: " << i << " offset: " << offset << " idx: " << idx << "\n";
    }

    if (logs) {
        if (ff_counts.find(ff_count) != ff_counts.end())
            ff_counts[ff_count] += 1;
        else
            ff_counts[ff_count] = 1;
    }
    i = idx;
    return ff_count;
}

std::string MoveStructure::reconstruct_lf() {
    orig_string = "";
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

        if (i % 10000 == 0)
            std::cerr << i << "\r";
        i += 1;
        // orig_string = rlbwt[run_index].get_c() + orig_string;
    }
    std::printf("Time measured for reconstructing the original text: %.3f seconds.\n", total_elapsed * 1e-9);

    std::cerr << "Finished reconstructing the original string.\n";
    std::cerr << "Total fast forward: " << ff_count_tot << "\n";
    return orig_string;
}

void MoveStructure::sequential_lf() {
    uint64_t line_index = 0;
    uint64_t row_index = 0;
    uint64_t ff_count_tot = 0;

    uint64_t total_elapsed = 0;
    for (uint64_t row_index = 0; row_index < r; row_index++) {
        auto& current = rlbwt[row_index];
        for (uint64_t j = 0; j < current.get_n(); j ++) {
            if (line_index % 10000 == 0)
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
        if (i % 10000 == 0)
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

uint64_t MoveStructure::get_n(uint64_t idx) {
    if (rlbwt[idx].is_overflow_n()) {
        return n_overflow[rlbwt[idx].get_n()];
    } else {
        return rlbwt[idx].get_n();
    }
}

uint64_t MoveStructure::get_n_ff(uint64_t idx) {
    if (rlbwt[idx].is_overflow_n_ff()) {
        return n_overflow[rlbwt[idx].get_n_ff()];
    } else {
        return rlbwt[idx].get_n_ff();
    }
}


uint64_t MoveStructure::get_offset(uint64_t idx) {
    if (rlbwt[idx].is_overflow_offset()) {
        return offset_overflow[rlbwt[idx].get_offset()];
    } else {
        return rlbwt[idx].get_offset();
    }
}

uint64_t MoveStructure::get_thresholds(uint64_t idx, uint32_t alphabet_index) {
    if (rlbwt[idx].is_overflow_thresholds()) {
        return thresholds_overflow[get_rlbwt_thresholds(idx, alphabet_index)][alphabet_index];
    } else {
        return get_rlbwt_thresholds(idx, alphabet_index);
    }
}

uint16_t MoveStructure::get_rlbwt_thresholds(uint64_t idx, uint16_t i) {
    if (i >= alphabet.size() - 1) {
        std::cerr << "get_thresholds: " << i << " is greater than or equal to " << alphabet.size() - 1 << "\n";
        exit(0);
    }

#if MODE == 0 || MODE == 1
    if (!onebit) {
        return rlbwt[idx].get_thresholds(i);
    }
#endif

#if MODE == 2
    return rlbwt[idx].get_threshold();
#endif
    std::cerr << "Undefined behavior!\n";
    exit(0);
}

void MoveStructure::set_rlbwt_thresholds(uint64_t idx, uint16_t i, uint16_t value) {
    if (i >= alphabet.size() - 1) {
        std::cerr << "get_thresholds: " << i << " is greater than or equal to " << alphabet.size() - 1 << "\n";
        exit(0);
    }

#if MODE == 0 || MODE == 1
    if (!onebit) {
        // rlbwt_thresholds[idx][i] = value;
        rlbwt[idx].set_thresholds(i, value);
    }
#endif

#if MODE == 2
    // rlbwt_1bit_thresholds[idx] = value;
    rlbwt[idx].set_threshold(value);
#endif
}

void MoveStructure::set_onebit() {
    onebit = true;
}

void MoveStructure::build_rlbwt(std::string bwt_filename) {
    // std::string bwt_filename = static_cast<std::string>(input_file) + std::string(".bwt");
    std::ifstream bwt_file(bwt_filename);    
    bwt_file.clear();
    bwt_file.seekg(0,std::ios_base::end);
    std::ios_base::streampos end_pos = bwt_file.tellg();
    if (verbose)
        std::cerr << "end_pos: " << end_pos << "\n";
    std::cerr << static_cast<uint64_t>(end_pos) << "\n";
    bwt_file.seekg(0);
    char current_char = bwt_file.get();
    char last_char = current_char;
    r = 0;
    size_t len = 0;
    
    std::ofstream len_file(bwt_filename + ".len", std::ios::out | std::ios::binary);
    std::ofstream heads_file(bwt_filename + ".heads");
    while (current_char != EOF) {
        if (r % 10000 == 0)
            std::cerr << r << "\r";
        if (current_char != last_char) {
            r += 1;
            // write output
            heads_file << last_char;
            len_file.write(reinterpret_cast<char*>(&len), 5);
            len = 0;
        } 
        len += 1;
        last_char = current_char;
        current_char = bwt_file.get();
        
    }

    // write output
    heads_file << last_char;
    len_file.write(reinterpret_cast<char*>(&len), 5);

    heads_file.close();
    len_file.close();
}

void MoveStructure::build(std::ifstream &bwt_file) {
    bwt_file.clear();
    bwt_file.seekg(0,std::ios_base::end);
    std::ios_base::streampos end_pos = bwt_file.tellg();
    if (verbose)
        std::cerr << "end_pos: " << end_pos << "\n";
    bwt_file.seekg(0);    
    std::cerr << "building.. \n";
    // Reading the BWT from the file
    bwt_string = "";
    bwt_string.resize(static_cast<uint64_t>(end_pos) + 1);
    uint64_t all_chars_count = 256;
    alphamap.resize(all_chars_count);
    std::fill(alphamap.begin(), alphamap.end(), alphamap.size());
    std::vector<uint64_t> all_chars(all_chars_count, 0);
    uint64_t current_char = bwt_file.get();
    r = 1;
    original_r = 1;
    // TODO Use a size based on the input size

    if (splitting) {

        std::string splitting_filename = input_file + std::string(".d_col");
        std::ifstream splitting_file(splitting_filename);

        bits.load(splitting_file);
        std::cerr << "bits.size after loading the d_col file: " << bits.size() << "\n";
    }
    else {
        if (verbose)
            std::cerr << "static_cast<uint64_t>(end_pos): " << static_cast<uint64_t>(end_pos) << "\n";
        bits = sdsl::bit_vector(static_cast<uint64_t>(end_pos) + 1, 0); // 5137858051
        bits[0] = 1;
    }



    std::cerr << "The main bit vector is built.\n";
    uint64_t bwt_curr_length = 0;
    while (current_char != EOF) { // && current_char != 10
        uint64_t current_char_ = static_cast<uint64_t>(current_char);
        // if (current_char != 'A' and current_char != 'C' and current_char != 'G' and current_char != 'T')
        //    std::cerr << "\ncurrent_char:" << current_char << "---" << static_cast<uint64_t>(current_char) << "---\n";
        if (original_r % 100000 == 0) {
            std::cerr << "original_r: " << original_r << "\t";
            std::cerr << "bwt_curr_length: " << bwt_curr_length << "\r";
        }
        if (bwt_curr_length > 0 && current_char != bwt_string[bwt_curr_length - 1]) {
            original_r += 1;
            if (!splitting) bits[bwt_curr_length] = 1;
        }
        if (splitting && bwt_curr_length > 0 && bits[bwt_curr_length]) {
            r += 1;
        }
        bwt_string[bwt_curr_length] = current_char;
        bwt_curr_length++;
        all_chars[current_char] += 1;

        current_char = bwt_file.get();
    }

    if (!splitting) r = original_r;

    std::cerr << "\n\nr: " << r << "\n";
    std::cerr << "original_r: " << original_r << "\n";
    length = bwt_curr_length; // bwt_string.length();
    std::cerr << "length: " << length << "\n";
    rlbwt.resize(r);
    /*if (!onebit)
        rlbwt_thresholds.resize(r);
    else
        rlbwt_1bit_thresholds.resize(r);*/
    //    rlbwt_chars.resize(r);
    if (verbose and bits.size() < 1000)
        std::cerr << "bits: " << bits << "\n";
    rbits = sdsl::rank_support_v<>(&bits);

    // Building the auxilary structures
    uint64_t alphabet_index = 0;
    // END_CHARACTER in the bwt created by pfp is 0
    for (uint64_t i = 1; i < all_chars_count; i++) {
        if (all_chars[i] != 0) {
            auto current_char = static_cast<unsigned char>(i);
            if (verbose)
                std::cerr << "i is " << i << "\t" << current_char 
                        << "\t" << all_chars[i] << " alphabet_index: " << alphabet_index << "\n";

            alphabet.push_back(current_char);
            counts.push_back(all_chars[i]);
            alphamap[i] = alphabet_index;
            alphabet_index += 1;

            sdsl::bit_vector* new_bit_vector = new sdsl::bit_vector(length, 0);
            occs.emplace_back(std::unique_ptr<sdsl::bit_vector>(new_bit_vector));
        }
    }

    if (alphabet.size() == 2) {
        std::cerr << "one bit alphabet version detected.\n";
        onebit = true;
        // bit1_begin = alphamap[bwt_string[0]];
    }
    if (onebit && alphabet.size() != 2) {
        std::cerr << "The fasta file has more than 2 characters in the alphabet!\n";
        exit(0);
    }
    if (alphabet.size() > 4) {
        std::cerr << "Warning: There are more than 4 characters, the index expexts only A, C, T and G in the reference.\n";
    }

    std::cerr << "All the characters are indexed.\n";

    for (uint64_t i = 0; i < length; i++) {
        if (i % 10000 == 0)
            std::cerr <<"length processed: " << i << "\r";
        if (static_cast<uint64_t>(bwt_string[i]) == END_CHARACTER)
            continue;
        auto& bit_vec = *occs[alphamap[static_cast<uint64_t>(bwt_string[i])]];
        bit_vec[i] = 1;
    }
    std::cerr << "\nAll the Occ bit vectors are built.\n";

    for (auto& occ: occs) {
        std::cerr << occs_rank.size() << "\r";
        if (verbose and (*occ).size() < 1000)
            std::cerr << *occ << "\n";
        occs_rank.emplace_back(std::unique_ptr<sdsl::rank_support_v<> >(new sdsl::rank_support_v<>(occ.get())));
    }
    if (verbose) {
        std::cerr << "size occs_rank:" << occs_rank.size() << "\n";
        std::cerr << "All Occ rank vectors are built.\n";
    }


    // Building the move structure rows
    uint64_t len = 0;
    uint64_t bwt_row = 0;
    uint64_t r_idx = 0;
    uint64_t offset = 0;
    uint64_t max_len = 0;
    if (verbose) {
        std::cerr << "bits.size(): " << bits.size() << "\n";
        std::cerr << "rank_support_v<>(&bits)(bits.size()): " << sdsl::rank_support_v<>(&bits)(bits.size()) << "\n";
    }
    sbits = sdsl::select_support_mcl<>(&bits);
    std::vector<uint64_t> all_p;
    for (uint64_t i = 0; i < length; i++) {
        if (i % 10000 == 0)
            std::cerr << i << "\r";
        if (bwt_string[i] == static_cast<unsigned char>(END_CHARACTER) ) {
            end_bwt_idx = r_idx;
        }
        if (i == length - 1 or bwt_string[i] != bwt_string[i+1] or bits[i+1]) {
            len += 1;
            uint64_t lf  = 0;
            if (bwt_string[i] != static_cast<unsigned char>(END_CHARACTER))
                lf = LF(bwt_row);
            else
                lf = 0;
            // bits[bwt_row] = 1;
            uint64_t pp_id = rbits(lf) - 1;
            if (bits[lf] == 1)
                pp_id += 1;

            // check the boundaries before performing select
            if (pp_id >= r) {
                std::cerr << "pp_id: " << pp_id << "r: " << r << "i: " << i << "bwt_row: " << bwt_row << "lf: " << lf << "\n";
                exit(0); // TODO: add error handling
            }
            if (lf < sbits(pp_id + 1)) {
                std::cerr << lf << " " << sbits(pp_id + 1);
                exit(0); // TODO: add error handling
            }
            offset = lf - sbits(pp_id + 1);

            if (verbose and r_idx == 0) // or any run to be inspected
                std::cerr << "r_idx: " << r_idx 
                          << " bwt_row: " << bwt_row
                          << " len: " << len
                          << " lf: " << lf 
                          << " offset: " << offset
                          << " pp_id: " << pp_id
                          << " sbits(pp_id): " << sbits(pp_id)
                          << " sbits(pp_id + 1): " << sbits(pp_id + 1)
                          << " sbits(pp_id - 1): " << sbits(pp_id - 1) << "\n";

            // rlbwt[r_idx].init(bwt_row, len, lf, offset, pp_id);
            rlbwt[r_idx].init(len, offset, pp_id);
            all_p.push_back(bwt_row);
            // To take care of cases where length of the run 
            // does not fit in uint16_t
            if (len >= std::numeric_limits<uint16_t>::max()) {
                n_overflow.push_back(len);
                if (n_overflow.size() - 1 >= std::numeric_limits<uint16_t>::max()) {
                    std::cerr << "Warning: the number of runs with overflow n is beyond uint16_t! " << n_overflow.size() - 1 << "\n";
                }
                rlbwt[r_idx].set_n(n_overflow.size() - 1);
                rlbwt[r_idx].set_overflow_n();
            }
            if (offset >= std::numeric_limits<uint16_t>::max()) {
                offset_overflow.push_back(offset);
                if (offset_overflow.size() - 1 >= std::numeric_limits<uint16_t>::max()) {
                    std::cerr << "Warning: the number of runs with overflow offset is beyond uint16_t! " << offset_overflow.size() - 1 << "\n";
                }
                rlbwt[r_idx].set_offset(offset_overflow.size() - 1);
                rlbwt[r_idx].set_overflow_offset();
            }

            if (len > max_len)
                max_len = len;
            if (logs) {
                if (run_lengths.find(len) != run_lengths.end())
                    run_lengths[len] += 1;
                else
                    run_lengths[len] = 1;
            }

            rlbwt[r_idx].set_c(bwt_string[i], alphamap);

            if (bwt_string[i] == static_cast<unsigned char>(END_CHARACTER)) {
                eof_row = r_idx;
                // bit1_after_eof = alphamap[bwt_string[i+1]];
            }

            bwt_row += len;
            len = 0;
            r_idx += 1;
        } else if (bwt_string[i] == bwt_string[i+1]) {
            len += 1;
        }
    }
    std::cerr << "All the move rows are built.\n";
    std::cerr << "Max run length: " << max_len << "\n";

    // compute the thresholds
    // initialize the start threshold at the last row
    std::vector<uint64_t> alphabet_thresholds(alphabet.size(), length);
    uint64_t thr_i = original_r - 1;
    uint64_t run_p = 0;
    if (verbose) {
        std::cerr << "thresholds.size():" << thresholds.size() << " length: " << length << " r: " << r <<  " original_r: " << original_r << "\n";
        std::cerr << "thresholds[r]: " << thresholds[original_r-1] << " r-1: " << thresholds[original_r - 2] << " r-2: " << thresholds[original_r - 3] << "\n";
    }
    for (uint64_t i = rlbwt.size() - 1; i > 0; --i) {
        if (i % 10000 == 0)
            std::cerr << "i: " << i << "\r";
        char rlbwt_c = alphabet[rlbwt[i].get_c()];
        /* if (thr_i != i) {
            std::cerr << "thr_i: " << thr_i << " i: " << i << "\n";
            exit(0);
        } */
        if (verbose and i >= rlbwt.size() - 10) 
            std::cerr << "i: " << i << "\n"
                << "rlbwt[i].get_offset(): " << get_offset(i) << "\n "
                << "get_n(i): " << get_n(i) << "\n"
                << "thresholds[i]: " << thresholds[i] << " "
                << "rlbwt_c: " << rlbwt_c << "\n";
        
        std::vector<uint64_t> current_thresholds;
        current_thresholds.resize(alphabet.size() - 1);
        for (uint64_t j = 0; j < alphabet.size(); j++) {
            if (alphabet[j] == rlbwt_c) {
                if (thr_i >= thresholds.size()) {
                    std::cerr << " thr_i = " << thr_i << " is out of bound:\n";
                    std::cerr << " thresholds.size = " << thresholds.size() << "\n";
                    exit(0); // TODO: add error handling
                }
                alphabet_thresholds[j] = thresholds[thr_i];
            } else {
                if (onebit and alphamap_3[alphamap[rlbwt_c]][j] != 0) {
                    std::cerr << "the alphamap_3 is not working for the one-bit alphabet:\n"
                                << "alphamap_3[alphamap[rlbwt_c]][j] = " << alphamap_3[alphamap[rlbwt_c]][j] << "\n";
                    exit(0); // TODO: add error handling
                }
                if (alphamap_3[alphamap[rlbwt_c]][j] >= alphabet.size() - 1) {
                    std::cerr << "alphamap_3 is not working in general:\n"
                                << "alphabet.size() - 1 = " << alphabet.size() - 1 << "\n"
                                << "alphamap_3[alphamap[rlbwt_c]][j] = " << alphamap_3[alphamap[rlbwt_c]][j] << "\n";
                    exit(0); // TODO: add error handling
                }

                if (alphabet_thresholds[j] >= all_p[i] + get_n(i)) {
                    // rlbwt[i].thresholds[j] = get_n(i);
                    if (rlbwt_c == END_CHARACTER) {
                        end_bwt_idx_thresholds[j] = get_n(i);
                        continue;
                    }
                    if (get_n(i) >= std::numeric_limits<uint16_t>::max()) {
                        rlbwt[i].set_overflow_thresholds();                            
                    }
                    set_rlbwt_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j], get_n(i));
                    current_thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = get_n(i);
                } else if (alphabet_thresholds[j] < all_p[i]) {
                    if (rlbwt_c == END_CHARACTER) {
                        end_bwt_idx_thresholds[j] = 0;
                        continue;
                    }
                    set_rlbwt_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j], 0);
                    current_thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = 0;
                } else {
                    if (rlbwt_c == END_CHARACTER) {
                        end_bwt_idx_thresholds[j] = alphabet_thresholds[j] - all_p[i];
                        continue;
                    }
                    if (alphabet_thresholds[j] - all_p[i] >= std::numeric_limits<uint16_t>::max()) {
                        rlbwt[i].set_overflow_thresholds();
                    }
                    set_rlbwt_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j], alphabet_thresholds[j] - all_p[i]);
                    current_thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = alphabet_thresholds[j] - all_p[i];
                }
                // printing the values for last 10 runs to debug
                if (verbose and i >= rlbwt.size() - 10) {
                    std::cerr << "\t j: \t" << j << " "
                        << "alphabet[j]: " << alphabet[j] << "  "
                        << "alphamap_3[alphamap[rlbwt_c]][j]: " << alphamap_3[alphamap[rlbwt_c]][j] << " "
                        << "alphabet_thresholds[j]: " << alphabet_thresholds[j] << " "
                        << "rlbwt[i].thresholds[j]:" << get_rlbwt_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j]) << "\n";
                }
            }
        }

        if (i > 0 && (rlbwt[i].get_c() != rlbwt[i - 1].get_c()  || i == end_bwt_idx || i-1 == end_bwt_idx)) {
            thr_i--;
        }

        if (rlbwt[i].is_overflow_thresholds()) {
            if (thresholds_overflow.size() >= std::numeric_limits<uint16_t>::max()) {
                std::cerr << "Undefined behaviour: the number of runs with overflow thresholds is beyond uint16_t! " 
                            << thresholds_overflow.size() << "\n";
                exit(0);
            }

            for (uint64_t k = 0; k < alphabet.size() - 1; k++) {
                set_rlbwt_thresholds(i, k, thresholds_overflow.size());
            }

            thresholds_overflow.push_back(current_thresholds);
        }
        run_p += get_n(i);
    } // end of the main for
    // since the thresholds for the first run was not calculated in the for
    for (uint64_t j = 0; j < alphabet.size() - 1; j++) {
        set_rlbwt_thresholds(0, j, 0);
    }
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
        auto occ_rank = rbits(char_count);
        last_runs.push_back(static_cast<uint64_t>(occ_rank - 1));
        last_offsets.push_back(char_count - all_p[last_runs.back()] - 1);
    }
    if (verbose) {
        for (uint64_t i = 0; i < first_runs.size(); i++) {
            std::cerr << "<--- " << first_runs[i] << "\t" << first_offsets[i] << "\n";
            std::cerr << "    -\n    -\n    -\n";
            std::cerr << ">--- " << last_runs[i] << "\t" << last_offsets[i] << "\n";
        }
    }

#if MODE == 1
    if (constant) {
        std::cerr << "Computing the next ups and downs.\n";
        compute_nexts();
    }
#endif
    std::cerr << "The move structure building is done.\n";
}

uint64_t scan_count;
#if MODE == 1
void MoveStructure::compute_nexts() {
    for (uint64_t i = rlbwt.size() - 1; i > 0; --i) {
        if (i % 100000 == 0)
            std::cerr << i << "\r";

        char rlbwt_c = alphabet[rlbwt[i].get_c()];
        for (uint64_t j = 0; j < alphabet.size(); j++) {
            if (i == end_bwt_idx) {
                auto idx = jump_up(i, alphabet[j], scan_count);
                end_bwt_idx_next_up[j] = (idx == r) ? std::numeric_limits<uint16_t>::max() : i - idx;
                idx = jump_down(i, alphabet[j], scan_count);
                end_bwt_idx_next_down[j] = (idx == r) ? std::numeric_limits<uint16_t>::max() : idx - i;
                continue;
            }
            if (alphabet[j] != rlbwt_c) {
                auto alphabet_idx = alphamap_3[alphamap[rlbwt_c]][j];
                auto idx = jump_up(i, alphabet[j], scan_count);
                if (idx == r) {
                    rlbwt[i].set_next_up(alphabet_idx, std::numeric_limits<uint16_t>::max());
                } else {
                    if (i - idx > std::numeric_limits<uint16_t>::max())
                        std::cerr << "Warning - jump up " << i - idx << " does not fit in 16 bits.\n";
                    rlbwt[i].set_next_up(alphabet_idx, i - idx);
                }

                idx = jump_down(i, alphabet[j], scan_count);
                if (idx == r) {
                    rlbwt[i].set_next_down(alphabet_idx, std::numeric_limits<uint16_t>::max());
                } else {
                    if (idx - i > std::numeric_limits<uint16_t>::max())
                        std::cerr << "Warning - jump down " << idx - i << " does not fit in 16 bits.\n";
                    rlbwt[i].set_next_down(alphabet_idx, idx - i);
                }
            }
        }
    }
}
#endif

uint64_t MoveStructure::fast_forward(uint64_t& offset, uint64_t idx, uint64_t x) {
    uint64_t idx_ = idx;
    if (verbose) {
        std::cerr << "\t \t fast forwarding:\n";
        std::cerr << " \t \t idx: " << idx << " offset: " << offset << " n:" << get_n_ff(idx) << "\n";
    }
    while (idx < r - 1 && offset >= get_n_ff(idx)) {
        offset -= get_n_ff(idx);
        idx += 1;
        if (verbose) std::cerr << "\t \t ff offset based: +" << idx - idx_ << "\n";
    }
    if (verbose) 
        std::cerr << " \t \t idx: " << idx << " offset: " << offset << " n:" << get_n_ff(idx) << "\n";
    return idx - idx_;
}

uint64_t MoveStructure::jump_up(uint64_t idx, char c, uint64_t& scan_count) {
    if (idx == 0)
        return r;
    char row_c = alphabet[rlbwt[idx].get_c_jj()];

    while (idx > 0 and row_c != c) {
        scan_count += 1;
        idx -= 1;
        row_c = alphabet[rlbwt[idx].get_c_jj()];
        // if (idx == 0) {
        //     break;
        // }
    }
    /* if (logs) {
        if (jumps.find(scan_count) != jumps.end())
            jumps[scan_count] += 1;
        else
            jumps[scan_count] = 1;
    } */
    if (verbose) 
        std::cerr << "\t \t \t \t idx after the while in the jump" << idx << "\n";
    return (row_c == c) ? idx : r;
}

uint64_t MoveStructure::jump_down(uint64_t idx, char c, uint64_t& scan_count) {
    if (idx == r - 1)
        return r;
    char row_c = alphabet[rlbwt[idx].get_c_jj()];

    while (idx < r - 1 && row_c != c) {
        scan_count += 1;
        idx += 1;
        row_c = alphabet[rlbwt[idx].get_c_jj()];
    }
    /* if (logs) {
        if (jumps.find(scan_count) != jumps.end())
            jumps[scan_count] += 1;
        else
            jumps[scan_count] = 1;
    } */
    if (verbose) 
        std::cerr << "\t \t \t \t idx after the while in the jump: " << idx << " " << c << " " << row_c << "\n";
    return (row_c == c) ? idx : r;
}

uint64_t MoveStructure::backward_search(std::string& R,  int32_t& pos_on_r) {
    if (!check_alphabet(R[pos_on_r])) {
        return 0;
    }
    uint64_t run_start = first_runs[alphamap[R[pos_on_r]] + 1];
    uint64_t offset_start = first_offsets[alphamap[R[pos_on_r]] + 1];
    uint64_t run_end = last_runs[alphamap[R[pos_on_r]] + 1];
    uint64_t offset_end = last_offsets[alphamap[R[pos_on_r]] + 1];

    // save the current range for reporting
    uint64_t run_start_prev = run_start;
    uint64_t offset_start_prev = offset_start;
    uint64_t run_end_prev = run_end;
    uint64_t offset_end_prev = offset_end;
    uint64_t scan_count = 0;
    while (pos_on_r > -1) {
        // save the current interval for reporting
        run_start_prev = run_start;
        offset_start_prev = offset_start;
        run_end_prev = run_end;
        offset_end_prev = offset_end;
        pos_on_r -= 1;
        if (run_start == end_bwt_idx or run_end == end_bwt_idx or !check_alphabet(R[pos_on_r])) {
            // std::cerr << "Not found\n";
            uint64_t match_count = 0;
            if (run_start_prev == run_end_prev) {
                match_count = offset_end_prev - offset_start_prev + 1;
            } else {
                match_count = (rlbwt[run_start_prev].get_n() - offset_start_prev) + (offset_end_prev + 1);
                for (uint64_t k = run_start_prev + 1; k < run_end_prev; k ++) {
                    match_count += rlbwt[k].get_n();
                }
            }
            // pos_on_r -= 1;
            return match_count;
        }
        if (verbose) {
            std::cerr << ">>> " << pos_on_r << ": " << run_start << "\t" << run_end << " " << offset_start << "\t" << offset_end << "\n";
            std::cerr << ">>> " << alphabet[rlbwt[run_start].get_c()] << " " << alphabet[rlbwt[run_end].get_c()] << " " << R[pos_on_r] << "\n";
        }
#if MODE == 0
        while ((run_start < run_end) and (alphabet[rlbwt[run_start].get_c()] != R[pos_on_r])) {
            run_start += 1;
            offset_start = 0;
            if (run_start >= r) {
                break;
            }
        }
        while ((run_end > run_start) and alphabet[rlbwt[run_end].get_c()] != R[pos_on_r]) {
            run_end -= 1;
            offset_end = rlbwt[run_end].get_n() - 1;
            if (run_end == 0) {
                break;
            }
        }
#endif
#if MODE == 1
        uint64_t read_alphabet_index = alphamap[static_cast<uint64_t>(R[pos_on_r])];
        if ((run_start < run_end) and (alphabet[rlbwt[run_start].get_c()] != R[pos_on_r])) {
            if (run_start == 0) {
                while ((run_start < run_end) and (alphabet[rlbwt[run_start].get_c()] != R[pos_on_r])) {
                    run_start += 1;
                    offset_start = 0;
                    if (run_start >= r) {
                        break;
                    }
                }
            } else {
                char rlbwt_char = alphabet[rlbwt[run_start].get_c_jj()];
                uint64_t alphabet_index = alphamap_3[alphamap[rlbwt_char]][read_alphabet_index];
                if (rlbwt[run_start].get_next_down(alphabet_index) == std::numeric_limits<uint16_t>::max()) {
                    run_start = r;
                } else {
                    uint64_t run_start_ = run_start + rlbwt[run_start].get_next_down(alphabet_index);
                    if (run_start_ <= run_end) {
                        run_start = run_start_;
                    } else {
                        run_start = run_end;
                    }
                    offset_start = 0;
                }
            }
        }
        if ((run_end > run_start) and (alphabet[rlbwt[run_end].get_c()] != R[pos_on_r])) {
            char rlbwt_char = alphabet[rlbwt[run_end].get_c_jj()];
            uint64_t alphabet_index = alphamap_3[alphamap[rlbwt_char]][read_alphabet_index];
            if (rlbwt[run_end].get_next_up(alphabet_index) == std::numeric_limits<uint16_t>::max()) {
                run_end = r;
            } else {
                uint64_t run_end_ = run_end - rlbwt[run_end].get_next_up(alphabet_index);
                if (run_end_ >= run_start) {
                    run_end = run_end_;
                } else {
                    run_end = run_start;
                }
                offset_end = rlbwt[run_end].get_n() - 1;
            }
        }
#endif
        if (verbose) {
          std::cerr << "<<< " << pos_on_r << ": " << run_start << "\t" << run_end << " " << offset_start << "\t" << offset_end << "\n";
          std::cerr << "<<< " << alphabet[rlbwt[run_start].get_c()] << " " << alphabet[rlbwt[run_end].get_c()] << " " << R[pos_on_r] << "\n";
        }
        if (((run_start < run_end) or (run_start == run_end and offset_start <= offset_end)) and
            (alphabet[rlbwt[run_start].get_c()] == R[pos_on_r] and alphabet[rlbwt[run_end].get_c()] == R[pos_on_r])) {
            LF_move(offset_start, run_start);
            LF_move(offset_end, run_end);
            if (pos_on_r == 0) {
                uint64_t match_count = 0;
                if (run_start == run_end) {
                    match_count = offset_end - offset_start + 1;
                } else {
                    match_count = (rlbwt[run_start].get_n() - offset_start) + (offset_end + 1);
                    for (uint64_t k = run_start + 1; k < run_end; k ++) {
                        match_count += rlbwt[k].get_n();
                    }
                }
                // pos_on_r -= 1;
                return match_count;
            }
        } else {
            // std::cerr << "Not found\n";
            uint64_t match_count = 0;
            if (run_start_prev == run_end_prev) {
                match_count = offset_end_prev - offset_start_prev + 1;
            } else {
                match_count = (rlbwt[run_start_prev].get_n() - offset_start_prev) + (offset_end_prev + 1);
                for (uint64_t k = run_start_prev + 1; k < run_end_prev; k ++) {
                    match_count += rlbwt[k].get_n();
                }
            }
            // pos_on_r -= 1;
            return match_count;
        }
        // pos_on_r -= 1;
    }
    std::cerr << "Should not get here!\n";
    return 0;
}

uint64_t MoveStructure::exact_matches(MoveQuery& mq) {
    std::string& R = mq.query();
    int32_t pos_on_r = R.length() - 1;
    uint64_t backward_search_count = 0;
    while (pos_on_r > -1) {
        backward_search_count += 1;
        int32_t pos_on_r_before = pos_on_r;
        uint64_t match_count = backward_search(R, pos_on_r);
        uint16_t match_len = pos_on_r_before - pos_on_r;
        mq.add_pml(match_len);
    }
    return backward_search_count;
}

uint64_t MoveStructure::query_pml(MoveQuery& mq, bool random) {
    if (random) {
        if (verbose)
            std::cerr << "Jumps are random - not with thresholds! \n";
        std::srand(time(0));
    }
    
    auto& R = mq.query();
    int32_t pos_on_r = R.length() - 1;
    uint64_t idx = r - 1; // std::rand() % r; // r - 1
    uint64_t offset = get_n(idx) - 1;

    uint16_t match_len = 0;
    uint16_t ff_count = 0;
    uint64_t ff_count_tot = 0;
    uint64_t scan_count = 0;
    auto t1 = std::chrono::high_resolution_clock::now();

    if (verbose) {
        std::cerr << "beginning of the search \ton query: " << mq.query() << "\t";
        std::cerr << "and on BWT, idx(r-1): " << idx << " offset: " << offset << "\n";
    }

    uint64_t iteration_count = 0;
    while (pos_on_r > -1) {
        iteration_count += 1;
        if (logs and (iteration_count-1)%200 == 0) {
            t1 = std::chrono::high_resolution_clock::now();
        }

        if (verbose)
            std::cerr << "Searching position " << pos_on_r << " of the read:\n";

        auto& row = rlbwt[idx];
        uint64_t row_idx = idx;
        char row_c = alphabet[row.get_c()];

        // if (alphamap[static_cast<uint64_t>(R[pos_on_r])] == alphamap.size()) {
        if (!check_alphabet(R[pos_on_r])) {
            // The character from the read does not exist in the reference
            match_len = 0;
            scan_count = 0;

            if (verbose)
                std::cerr << "\t The character " << R[pos_on_r] << " does not exist.\n";
        } else if (row_c == R[pos_on_r]) {
            // Case 1
            match_len += 1;
            scan_count = 0;

            if (verbose) {
                std::cerr << "\t Cas1: It was a match. \n" << "\t Continue the search...\n";
                std::cerr << "\t match_len: " << match_len << "\n";
                std::cerr << "\t current_id: " << idx << "\t row.id: " << row.get_id() << "\n" 
                          << "\t row.get_n: " << get_n(row_idx) << " rlbwt[idx].get_n: " << get_n(row.get_id()) << "\n"
                          << "\t offset: " << offset << "\t row.get_offset(): " << get_offset(row_idx) << "\n";
            }
        } else {
            // Case 2
            // Jumping up or down (randomly or with thresholds)
            if (verbose)
                std::cerr << "\t Case2: Not a match, looking for a match either up or down...\n";

            uint64_t idx_before_jump = idx;
            bool up = random ? jump_randomly(idx, R[pos_on_r], scan_count) : 
                               jump_thresholds(idx, offset, R[pos_on_r], scan_count);
            match_len = 0;
            // scan_count = (!constant) ? std::abs((int)idx - (int)idx_before_jump) : 0;
 
            char c = alphabet[rlbwt[idx].get_c_mm()];

            if (verbose)
                std::cerr << "\t up: " << up << " idx: " << idx << " c:" << c << "\n";

            // sanity check
            if (c == R[pos_on_r]) {
                // Observing a match after the jump
                // The right match_len should be:
                // min(new_lcp, match_len + 1)
                // But we cannot compute lcp here
                offset = up ? get_n(idx) - 1 : 0;
                if (verbose)
                    std::cerr << "\t idx: " << idx << " offset: " << offset << "\n";
            } else {
                std::cerr << "\t \t This should not happen!\n";
                std::cerr << "\t \t r[pos]:" <<  R[pos_on_r] << " t[pointer]:" << c << "\n";
                std::cerr << "\t \t " << up << ", " << onebit << ", " << R[pos_on_r] << ", " << pos_on_r << "\n";
                std::cerr << "\t \t ";
                for (int k = 10; k > 0; --k)
                    std::cerr << alphabet[rlbwt[idx - k].get_c()] << "-";
                for (int k = 0; k < 10; k++)
                    std::cerr << alphabet[rlbwt[idx + k].get_c()] << "-";
                std::cerr << "\n";
                auto saved_idx = idx;

                verbose = true;
                jump_thresholds(saved_idx, offset, R[pos_on_r], scan_count);
                exit(0);
            }
        }

        mq.add_pml(match_len);
        pos_on_r -= 1;

        // LF step
        ff_count = LF_move(offset, idx);
        ff_count_tot += ff_count;
        if (logs) {
            if (iteration_count % 200 == 0) {
                auto t2 = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
                mq.add_cost(elapsed);
            }
            mq.add_fastforward(ff_count);
            mq.add_scan(scan_count);
        }
    }
    return ff_count_tot;
}

bool MoveStructure::jump_thresholds(uint64_t& idx, uint64_t offset, char r_char, uint64_t& scan_count) {
    uint64_t saved_idx = idx;
    uint64_t alphabet_index = alphamap[static_cast<uint64_t>(r_char)];
    scan_count = 0;

    if (verbose)
        std::cerr << "\t \t \t jumping with thresholds ... \n";

    char rlbwt_char = alphabet[rlbwt[idx].get_c_jj()];

    if (verbose) {
        std::cerr << "\t \t \t alphabet_index: " << alphabet_index << " r_char:" << r_char << " rlbwt_char:" << rlbwt_char << "\n";
        std::cerr << "\t \t \t idx:" << idx << "\n"
                  << "\t \t \t offset: " << offset << " threshold:" << get_thresholds(idx, alphamap_3[alphamap[rlbwt_char]][alphabet_index]) << "\n";
    }

    if (idx == end_bwt_idx) {
        if (verbose) std::cerr << "\t \t \t idx == end_bwt_idx" 
                                << "\n\t \t \t idx: " << idx << " end_bwt_idx: " << end_bwt_idx << "\n";
        if (offset >= end_bwt_idx_thresholds[alphabet_index]) {
            if (verbose)
                std::cerr << "\t \t \t Jumping down with thresholds:\n";
#if MODE == 1
            if (constant) {
                scan_count += 1;
                if (end_bwt_idx_next_down[alphabet_index] == std::numeric_limits<uint16_t>::max())
                    idx = r;
                else
                    idx = saved_idx + end_bwt_idx_next_down[alphabet_index];
            }
#endif

#if MODE == 0 || MODE == 2
            idx = jump_down(saved_idx, r_char, scan_count);
#endif
            if (r_char != alphabet[rlbwt[idx].get_c_mm()])
                std::cerr << "1: " << r_char << " " << alphabet[rlbwt[idx].get_c_mm()];
            return false;
        } else {
            if (verbose)
                std::cerr << "\t \t \t Jumping up with thresholds:\n";
#if MODE == 1
            if (constant) {
                scan_count += 1;
                if (end_bwt_idx_next_up[alphabet_index] == std::numeric_limits<uint16_t>::max())
                    idx = r;
                else
                    idx = saved_idx - end_bwt_idx_next_up[alphabet_index];
            }
#endif

#if MODE == 0 || MODE == 2
            idx = jump_up(saved_idx, r_char, scan_count);
#endif
            if (r_char != alphabet[rlbwt[idx].get_c_mm()])
                std::cerr << "2: " << r_char << " " << alphabet[rlbwt[idx].get_c_mm()];
            return true;
        }
    }

    if (verbose) std::cerr << "\t \t \t rlbwt[idx].get_offset(): " << get_offset(idx) 
                            << " get_thresholds(idx, alphabet_index): " << get_thresholds(idx, alphamap_3[alphamap[rlbwt_char]][alphabet_index]) 
                            << "\n\t \t \t idx:" << idx << "\n";

    alphabet_index = alphamap_3[alphamap[rlbwt_char]][alphabet_index];
    if (onebit and alphabet_index != 0)
        std::cerr << "error: the alphamap_3 is not working for the one-bit alphabet - " 
                    << alphabet_index << "!\n";
    if (alphabet_index == 3)
        std::cerr << "error: alphamap_3 is not working in general - " 
                    << alphabet_index << "!\n";

    if (offset >= get_thresholds(idx, alphabet_index)) {
        if (verbose)
            std::cerr << "\t \t \t Jumping down with thresholds:\n";
#if MODE == 1
        if (constant) {
            scan_count += 1;
            if (rlbwt[saved_idx].get_next_down(alphabet_index) == std::numeric_limits<uint16_t>::max())
                idx = r;
            else
                idx = saved_idx + rlbwt[saved_idx].get_next_down(alphabet_index);
        }
#endif
        auto tmp = idx;
#if MODE == 0 || MODE == 2
        idx = jump_down(saved_idx, r_char, scan_count);
#endif
        if (r_char != alphabet[rlbwt[idx].get_c_mm()]) {
            std::cerr << "3: " << r_char << " " << alphabet[rlbwt[idx].get_c_mm()] << "\n";
            std::cerr << "idx: " << idx << " saved_idx: " << saved_idx << " tmp: " << tmp << "\n";
            std::cerr << "offset: " << offset << "\n";
            std::cerr << "get_thresholds(saved_idx, alphabet_index): " << get_thresholds(saved_idx, alphabet_index) << "\n";
        }
        return false;
    } else {
        if (verbose)
            std::cerr << "\t \t \t Jumping up with thresholds:\n";
#if MODE == 1
        scan_count += 1;
        if (constant) {
            if (rlbwt[saved_idx].get_next_up(alphabet_index) == std::numeric_limits<uint16_t>::max())
                idx = r;
            else
                idx = saved_idx - rlbwt[saved_idx].get_next_up(alphabet_index);
        }
#endif

#if MODE == 0 || MODE == 2
        idx = jump_up(saved_idx, r_char, scan_count);
#endif
        if (r_char != alphabet[rlbwt[idx].get_c_mm()]) {
            std::cerr << "idx: " << idx << " saved_idx: " << saved_idx << "\n";
            std::cerr << "4: " << r_char << " " << alphabet[rlbwt[idx].get_c_mm()] << "\n";
            std::cerr << "offset: " << offset << "\n";
            std::cerr << "get_thresholds(saved_idx, alphabet_index): " << get_thresholds(saved_idx, alphabet_index) << "\n";
        }
        return true;
    }

    // TODO: default return?

    if (r_char != alphabet[rlbwt[idx].get_c_mm()])
        std::cerr << "5: " << r_char << " " << alphabet[rlbwt[idx].get_c_mm()];
    return false;
}

bool MoveStructure::jump_randomly(uint64_t& idx, char r_char, uint64_t& scan_count) {
    uint64_t saved_idx = idx;
    uint64_t jump = std::rand() % 2; // To replace with ...
    bool up = false;
    scan_count = 0;
    if (verbose)
        std::cerr << "idx before jump: " << idx << "\n";

    if ( (jump == 1 && idx > 0) or idx == r - 1) {
        if (verbose)
            std::cerr << "Jumping up randomly:\n";

        // jumping up
        up = true;
        idx = jump_up(saved_idx, r_char, scan_count);
        if (verbose)
            std::cerr << "idx after jump: " << idx << "\n";
        char c = alphabet[rlbwt[idx].get_c()];
        if (c != r_char) {
            if (verbose)
                std::cerr << "Up didn't work, try jumping down:\n";

            // jump down
            up = false;
            idx = jump_down(saved_idx, r_char, scan_count);
            if (verbose)
                std::cerr << "idx after jump: " << idx << "\n";
        }
    } else {
        if (verbose)
            std::cerr << "Jumping down randomly:\n";

        // jumping down
        up = false;
        idx = jump_down(saved_idx, r_char, scan_count);
        if (verbose)
            std::cerr << "idx after jump: " << idx << "\n";
        char c = alphabet[rlbwt[idx].get_c()];
        if (c != r_char) {
            if (verbose)
                std::cerr << "Down didn't work, try jumping up:\n";

            // jump up
            up = true;
            idx = jump_up(saved_idx, r_char, scan_count);
            if (verbose)
                std::cerr << "idx after jump: " << idx << "\n";
        }
    }
    return up;
}

/*bool MoveStructure::jump_naive_lcp(uint64_t& idx, uint64_t pointer, char r_char, uint64_t& lcp) {
    uint64_t up_idx = jump_up(idx, r_char);
    uint64_t up_pointer = up_idx != r ? rlbwt[up_idx].get_p() + rlbwt[up_idx].get_n() - 1 : length;
    uint64_t down_idx = jump_down(idx, r_char);
    uint64_t down_pointer = down_idx != r ? rlbwt[down_idx].get_p() : length;
    uint64_t up_lcp = naive_lcp(pointer, up_pointer);
    uint64_t down_lcp = naive_lcp(pointer, down_pointer);
    if (verbose) {
        std::cerr << "down_idx: " << down_idx << " up_idx: " << up_idx << "\n";
        std::cerr << "down_pointer: " << down_pointer << " up_pointer: " << up_pointer << "\n";
        std::cerr << "down_lcp: " << down_lcp << " up_lcp: " << up_lcp << "\n";
    }

    if ( (up_idx != r) &&
         (up_lcp >= down_lcp || down_idx == r) ) {
        idx = up_idx;
        lcp = up_lcp;
        return true;
    } else if (down_idx != r)  {
        idx = down_idx;
        lcp = down_lcp;
        return false;
    } else {
        std::cerr << "should not happen during naive lcp jump!\n";
        exit(0);
    }
}*/

bool MoveStructure::check_alphabet(char c) {
    // for (char c_: alphabet) {
    //     if (c == c_)
    //         return true;
    // }
    // return false;
    return alphamap[static_cast<uint64_t>(c)] != alphamap.size();
}

void MoveStructure::serialize(std::string output_dir) {
    mkdir(output_dir.c_str(),0777);
    std::string fname = output_dir + "/movi_index.bin";
    if (onebit)
        fname = output_dir + "/movi_index_onebit.bin";
    std::ofstream fout(fname, std::ios::out | std::ios::binary);
    std::cerr << "length: " << length << " r: " << r << " end_bwt_idx: " << end_bwt_idx << "\n";
    fout.write(reinterpret_cast<char*>(&length), sizeof(length));
    fout.write(reinterpret_cast<char*>(&r), sizeof(r));
    fout.write(reinterpret_cast<char*>(&end_bwt_idx), sizeof(end_bwt_idx));
    fout.write(reinterpret_cast<char*>(&end_bwt_idx_thresholds[0]), 4*sizeof(end_bwt_idx_thresholds[0]));
    fout.write(reinterpret_cast<char*>(&end_bwt_idx_next_down[0]), 4*sizeof(end_bwt_idx_next_down[0]));
    fout.write(reinterpret_cast<char*>(&end_bwt_idx_next_up[0]), 4*sizeof(end_bwt_idx_next_up[0]));

    uint64_t alphamap_size = alphamap.size();
    fout.write(reinterpret_cast<char*>(&alphamap_size), sizeof(alphamap_size));
    fout.write(reinterpret_cast<char*>(&alphamap[0]), alphamap.size()*sizeof(alphamap[0]));

    uint64_t alphabet_size = alphabet.size();
    fout.write(reinterpret_cast<char*>(&alphabet_size), sizeof(alphabet_size));    
    fout.write(reinterpret_cast<char*>(&alphabet[0]), alphabet.size()*sizeof(alphabet[0]));

    fout.write(reinterpret_cast<char*>(&splitting), sizeof(splitting));
    fout.write(reinterpret_cast<char*>(&constant), sizeof(constant));
    fout.write(reinterpret_cast<char*>(&onebit), sizeof(onebit));
    fout.write(reinterpret_cast<char*>(&rlbwt[0]), rlbwt.size()*sizeof(rlbwt[0]));

    uint64_t n_overflow_size = n_overflow.size();
    fout.write(reinterpret_cast<char*>(&n_overflow_size), sizeof(n_overflow_size));
    fout.write(reinterpret_cast<char*>(&n_overflow[0]), n_overflow.size()*sizeof(uint64_t));
    uint64_t offset_overflow_size = offset_overflow.size();
    fout.write(reinterpret_cast<char*>(&offset_overflow_size), sizeof(offset_overflow_size));
    fout.write(reinterpret_cast<char*>(&offset_overflow[0]), offset_overflow.size()*sizeof(uint64_t));
    uint64_t thresholds_overflow_size = thresholds_overflow.size();
    fout.write(reinterpret_cast<char*>(&thresholds_overflow_size), sizeof(thresholds_overflow_size));
    for (uint32_t i = 0; i < thresholds_overflow_size; i++) {
        fout.write(reinterpret_cast<char*>(&thresholds_overflow[i][0]), (alphabet.size() - 1)*sizeof(thresholds_overflow[i][0]));
    }

    size_t orig_size = orig_string.size();
    fout.write(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));
    fout.write(reinterpret_cast<char*>(&reconstructed), sizeof(reconstructed));

    fout.write(reinterpret_cast<char*>(&eof_row), sizeof(eof_row));

    uint64_t counts_size = counts.size();
    fout.write(reinterpret_cast<char*>(&counts_size), sizeof(counts_size));
    fout.write(reinterpret_cast<char*>(&counts[0]), counts.size()*sizeof(counts[0]));

    uint64_t last_runs_size = last_runs.size();
    fout.write(reinterpret_cast<char*>(&last_runs_size), sizeof(last_runs_size));
    fout.write(reinterpret_cast<char*>(&last_runs[0]), last_runs.size()*sizeof(last_runs[0]));
    fout.write(reinterpret_cast<char*>(&last_offsets[0]), last_offsets.size()*sizeof(last_offsets[0]));
    fout.write(reinterpret_cast<char*>(&first_runs[0]), first_runs.size()*sizeof(first_runs[0]));
    fout.write(reinterpret_cast<char*>(&first_offsets[0]), first_offsets.size()*sizeof(first_offsets[0]));

    fout.close();
}

void MoveStructure::deserialize(std::string index_dir) {
    std::string fname = index_dir + "/movi_index.bin";
    if (onebit)
        fname = index_dir + "/movi_index_onebit.bin";
    std::ifstream fin(fname, std::ios::in | std::ios::binary);
    fin.seekg(0, std::ios::beg); 

    fin.read(reinterpret_cast<char*>(&length), sizeof(length));
    fin.read(reinterpret_cast<char*>(&r), sizeof(r));
    fin.read(reinterpret_cast<char*>(&end_bwt_idx), sizeof(end_bwt_idx));
    fin.read(reinterpret_cast<char*>(&end_bwt_idx_thresholds[0]), 4*sizeof(end_bwt_idx_thresholds[0]));
    fin.read(reinterpret_cast<char*>(&end_bwt_idx_next_down[0]), 4*sizeof(end_bwt_idx_next_down[0]));
    fin.read(reinterpret_cast<char*>(&end_bwt_idx_next_up[0]), 4*sizeof(end_bwt_idx_next_up[0]));

    std::cerr << "length: " << length << " r: " << r << " end_bwt_idx: " << end_bwt_idx << "\n";

    uint64_t alphamap_size;
    fin.read(reinterpret_cast<char*>(&alphamap_size), sizeof(alphamap_size));
    alphamap.resize(alphamap_size);
    fin.read(reinterpret_cast<char*>(&alphamap[0]), alphamap_size*sizeof(alphamap[0]));
    uint64_t alphabet_size;
    fin.read(reinterpret_cast<char*>(&alphabet_size), sizeof(alphabet_size));
    alphabet.resize(alphabet_size);
    fin.read(reinterpret_cast<char*>(&alphabet[0]), alphabet_size*sizeof(alphabet[0]));
    std::cerr << "alphabet_size: " << alphabet_size << "\n";
    if (alphabet.size() > 4) {
        std::cerr << "Warning: There are more than 4 characters, the index expects only A, C, T and G in the reference.\n";
    }

    fin.read(reinterpret_cast<char*>(&splitting), sizeof(splitting));
    fin.read(reinterpret_cast<char*>(&constant), sizeof(constant));
    fin.read(reinterpret_cast<char*>(&onebit), sizeof(onebit));

    if (!check_mode()) {
        std::cerr << "Your settings: \n"
                    << "onebit: " << onebit << "\n"
                    << "constant: " << constant << "\n"
                    << "splitting " << splitting << "\n";
        exit(0);
    }

    rlbwt.resize(r);
    fin.read(reinterpret_cast<char*>(&rlbwt[0]), r*sizeof(MoveRow));
    uint64_t n_overflow_size;
    fin.read(reinterpret_cast<char*>(&n_overflow_size), sizeof(n_overflow_size));
    n_overflow.resize(n_overflow_size);
    fin.read(reinterpret_cast<char*>(&n_overflow[0]), n_overflow_size*sizeof(uint64_t));
    uint64_t offset_overflow_size;
    fin.read(reinterpret_cast<char*>(&offset_overflow_size), sizeof(offset_overflow_size));
    offset_overflow.resize(offset_overflow_size);
    fin.read(reinterpret_cast<char*>(&offset_overflow[0]), offset_overflow_size*sizeof(uint64_t));
    uint64_t thresholds_overflow_size;
    fin.read(reinterpret_cast<char*>(&thresholds_overflow_size), sizeof(thresholds_overflow_size));
    thresholds_overflow.resize(thresholds_overflow_size);
    for (uint32_t i = 0; i < thresholds_overflow_size; i++) {
        thresholds_overflow[i].resize(alphabet.size() - 1);
        fin.read(reinterpret_cast<char*>(&thresholds_overflow[i][0]), (alphabet.size() - 1)*sizeof(uint64_t));
    }

    std::cerr << "All the move rows are read.\n";

    size_t orig_size;
    fin.read(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));

    fin.read(reinterpret_cast<char*>(&reconstructed), sizeof(reconstructed));
    reconstructed = false;

    fin.read(reinterpret_cast<char*>(&eof_row), sizeof(eof_row));

    uint64_t counts_size = 0;
    fin.read(reinterpret_cast<char*>(&counts_size), sizeof(counts_size));
    counts.resize(counts_size);
    fin.read(reinterpret_cast<char*>(&counts[0]), counts_size*sizeof(uint64_t));

    uint64_t last_runs_size = 0;
    fin.read(reinterpret_cast<char*>(&last_runs_size), sizeof(last_runs_size));
    last_runs.resize(last_runs_size);
    fin.read(reinterpret_cast<char*>(&last_runs[0]), last_runs_size*sizeof(uint64_t));
    last_offsets.resize(last_runs_size);
    fin.read(reinterpret_cast<char*>(&last_offsets[0]), last_runs_size*sizeof(uint64_t));
    first_runs.resize(last_runs_size);
    fin.read(reinterpret_cast<char*>(&first_runs[0]), last_runs_size*sizeof(uint64_t));
    first_offsets.resize(last_runs_size);
    fin.read(reinterpret_cast<char*>(&first_offsets[0]), last_runs_size*sizeof(uint64_t));
    fin.close();
}

void MoveStructure::analyze_rows() {
    std::vector<uint64_t> counts(16,0);
    for (uint64_t i = 0; i < r; i++) {
        // std::cerr << i << " " << rlbwt[i].get_id() << " " << alphabet[rlbwt[i].get_c()] << " " << get_n(i) << " " << get_offset(i) << ":\t";
        // for (int  j = 0; j < alphabet.size() - 1; j ++)
        //     std::cerr << get_thresholds(i, j) << " ";
        // std::cerr << "\n";
        if (i%100000 == 0) std::cerr << i << "\r";
        for (int j = 0; j < 16; j ++) {
            if (get_n(i) >= std::pow(2,j)) {
                counts[j] += 1;
            }
        }
        if (i%10000 == 0) std::cerr << i << "\r";
        if (((get_thresholds(i, 0) != 0 and get_thresholds(i, 0) != get_n(i)) and
             (get_thresholds(i, 1) != 0 and get_thresholds(i, 1) != get_n(i)) and
             (get_thresholds(i, 2) != 0 and get_thresholds(i, 2) != get_n(i)) and
             (get_thresholds(i, 0) != get_thresholds(i, 1) or get_thresholds(i, 1) != get_thresholds(i, 2))) or
            ((get_thresholds(i,0) != 0 and get_thresholds(i,0) != get_n(i) and get_thresholds(i,1) != 0 and get_thresholds(i,1) != get_n(i) and get_thresholds(i, 0) != get_thresholds(i, 1)) or
             (get_thresholds(i,0) != 0 and get_thresholds(i,0) != get_n(i) and get_thresholds(i,2) != 0 and get_thresholds(i,2) != get_n(i) and get_thresholds(i, 0) != get_thresholds(i, 2)) or
             (get_thresholds(i,2) != 0 and get_thresholds(i,2) != get_n(i) and get_thresholds(i,1) != 0 and get_thresholds(i,1) != get_n(i) and get_thresholds(i, 2) != get_thresholds(i, 1)))
            ) {
            if (get_n(i) >= 256) {
                std::cerr << i << " " << rlbwt[i].get_id() << " " << alphabet[rlbwt[i].get_c()] << " " << get_n(i) << " " << get_offset(i) << ":\t";
                for (int  j = 0; j < alphabet.size() - 1; j ++)
                    std::cerr << get_thresholds(i, j) << " ";
                std::cerr << "\n";
            }
        }
    }
    for (int j=0; j < 16; j++) {
        std::cerr << j << ": " << counts[j] << "\n";
    }
}

void MoveStructure::print_stats() {
    std::cerr << "n: " << length << "\n";
    std::cerr << "r: " << r << "\n";
    std::cerr << "n/r: " << length/r << "\n";
    if (original_r != 0) {
        std::cerr << "original_r: " << original_r << "\n";
        std::cerr << "n/original_r: " << length/original_r << "\n";
    }
}
