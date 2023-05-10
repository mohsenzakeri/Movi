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
        std::cerr<<("open() file " + tmp_filename + " failed");

    int fn = fileno(fd);
    if (fstat(fn, &filestat) < 0)
        std::cerr<<("stat() file " + tmp_filename + " failed");

    if (filestat.st_size % THRBYTES != 0)
        std::cerr<<("invilid file " + tmp_filename);

    size_t length_thr = filestat.st_size / THRBYTES;
    size_t threshold = 0;

    thresholds = sdsl::int_vector<>(length_thr, 0, log_n);

    size_t i = 0;
    for (i = 0; i < length_thr; ++i) {
        size_t threshold = 0;
        if ((fread(&threshold, THRBYTES, 1, fd)) != 1)
            std::cerr<<("fread() file " + tmp_filename + " failed");
        thresholds[i] = threshold;
    }
    std::cerr << "Finished reading " << i << " thresholds.\n";
}

MoveStructure::MoveStructure(bool verbose_, bool logs_) {
    verbose = verbose_;
    logs = logs_;
}

MoveStructure::MoveStructure(char* input_file, bool bit1_, bool verbose_, bool logs_) {
    verbose = verbose_;
    logs = logs_;
    bit1 = bit1_;
    reconstructed = false;

    std::string bwt_filename = input_file + std::string(".bwt");
    std::ifstream bwt_file(bwt_filename);

    std::string thr_filename = input_file + std::string(".thr_pos");
    read_thresholds(thr_filename, thresholds);

    build(bwt_file);
}

char MoveStructure::compute_char(uint64_t idx) {
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

/*std::string MoveStructure::reconstruct_move() {
    orig_string = "";
    uint64_t bwt_index = 0;
    uint64_t run_index = 0;
    uint64_t i = 0;
    uint64_t ff_count_tot = 0;
    for (; bwt_index != end_bwt_row; ff_count_tot += LF_move(bwt_index, run_index)) {
        if (i % 10000 == 0)
            std::cerr<< i << "\r";
        i += 1;
        // orig_string = rlbwt[run_index].get_c() + orig_string;
    }
    std::cerr << i << " " << bwt_index << " " << run_index << "\n";
    std::cerr << length << "\n";
    std::cerr << "Finished reconstructing the original string.\n";
    std::cerr << "Total fast forward: " << ff_count_tot << "\n";
    return orig_string;
}*/

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

/*uint64_t MoveStructure::LF_move(uint64_t& pointer, uint64_t& i) {
    auto& row = rlbwt[i];
    auto idx = row.get_id();
    pointer = row.get_pp() + (pointer - row.get_p());
    uint64_t ff_count = 0;

    if (idx < r - 1 && pointer >= rlbwt[idx].get_p() + get_n(idx)) {
        uint64_t idx_ = fast_forward(pointer, idx);
        idx += idx_;
        ff_count += idx_;
    }

    if (logs) {
        if (ff_counts.find(ff_count) != ff_counts.end())
            ff_counts[ff_count] += 1;
        else
            ff_counts[ff_count] = 1;
    }
    i = idx;
    return ff_count;
}*/

/*void MoveStructure::all_lf_test() { // std::ifstream &bwt_file
    // bwt_file.clear();
    // bwt_file.seekg(0);
    // char  current_char = bwt_file.get();
    // char last_char = current_char;
    uint64_t line_index = 0;
    uint64_t row_index = 0;
    uint64_t ff_count_tot = 0;
    // while (current_char != EOF) { // && current_char != 10
    //     if (line_index % 10000 == 0)
    //         std::cerr<< line_index << "\r";

    //     uint64_t pointer = line_index;
    //     uint64_t i = row_index;
    //     ff_count_tot += LF_move(pointer, i);

    //     last_char = current_char;
    //     current_char = bwt_file.get();
    //     line_index += 1;
    //     if (current_char != last_char) {
    //         row_index += 1;
    //     }
    // }

    for (uint64_t row_index = 0; row_index < r; row_index++) {
        auto& current = rlbwt[row_index];
        for (uint64_t j = 0; j < current.get_n(); j ++) {
            if (line_index % 10000 == 0)
                std::cerr<< line_index << "\r";
            uint64_t pointer = line_index;
            uint64_t i = row_index;
            ff_count_tot += LF_move(pointer, i);
            line_index += 1;
        }

    }
    std::cerr<< line_index << "\r";
    std::cerr<< "Finished performing LF query for all the BWT characters.\n";
    std::cerr << "Total fast forward: " << ff_count_tot << "\n";
}*/

/*uint64_t MoveStructure::random_lf_test() {
    std::srand(time(0));
    uint64_t ff_count_tot = 0;
    for (uint64_t i = 0; i < length; i++) {
        if (i % 10000 == 0)
            std::cerr<< i << "\r";

        uint64_t idx = std::rand() % r;
        auto& row = rlbwt[idx];
        uint64_t pointer = row.get_p() + (std::rand() % row.get_n());
        
        ff_count_tot += LF_move(pointer, idx);
    }
    std::cerr << length << "\r";
    std::cerr << "Finished performing LF query for all the BWT characters.\n";
    std::cerr << "Total fast forward: " << ff_count_tot << "\n";
    return ff_count_tot;
}*/

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

uint64_t MoveStructure::get_thresholds(uint64_t idx, uint32_t alphabet_index) {
    if (rlbwt[idx].is_overflow_thresholds()) {
        return thresholds_overflow[rlbwt[idx].thresholds[alphabet_index]][alphabet_index];
    } else {
        return rlbwt[idx].thresholds[alphabet_index];
    }
}
void MoveStructure::build(std::ifstream &bwt_file) {
    bwt_file.clear();
    bwt_file.seekg(0,std::ios_base::end);
    std::ios_base::streampos end_pos = bwt_file.tellg();
    std::cerr << "end_pos: " << end_pos << "\n";
    std::cerr << static_cast<uint64_t>(end_pos) << "\n";
    bwt_file.seekg(0);    
    std::cerr<<"building.. \n";
    // Reading the BWT from the file
    bwt_string = "";
    uint64_t all_chars_count = 256;
    alphamap.resize(all_chars_count);
    std::fill(alphamap.begin(), alphamap.end(), alphamap.size());
    std::vector<uint64_t> all_chars(all_chars_count, 0);
    uint64_t current_char = bwt_file.get();
    r = 1;
    // TODO Use a size based on the input size

    bits = sdsl::bit_vector(static_cast<uint64_t>(end_pos) + 1, 0); // 5137858051
    bits[0] = 1;
    std::cerr<<"bit vector is built!\n";
    while (current_char != EOF) { // && current_char != 10
        if (r % 10000 == 0)
            std::cerr<< r << "\r";
        if (bwt_string.length() > 0 && current_char != bwt_string.back()) {
            r += 1;
            bits[bwt_string.length()] = 1;
        }
        bwt_string += current_char;
        all_chars[current_char] += 1;

        current_char = bwt_file.get();
    }
    std::cerr<< "r: " << r << "\n";
    length = bwt_string.length();
    std::cerr<<"length: " << length << "\n";
    rlbwt.resize(r);
    // if (!bit1)
    //    rlbwt_chars.resize(r);
    if (verbose and bits.size() < 1000)
        std::cerr<<"bits: " << bits << "\n";
    rbits = sdsl::rank_support_v<>(&bits);

    // Building the auxilary structures
    uint64_t alphabet_index = 0;
    // END_CHARACTER in the bwt created by pfp is 0
    for (uint64_t i = 1; i < all_chars_count; i++) {
        if (all_chars[i] != 0) {
            auto current_char = static_cast<unsigned char>(i);
            if (verbose)
                std::cerr<< "i is " << i << "\t" << current_char 
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
        bit1 = true;
        bit1_begin = alphamap[bwt_string[0]];
    }

    std::cerr << "All the characters are indexed.\n";

    for (uint64_t i = 0; i < length - 1; i++) {
        if (i % 10000 == 0)
            std::cerr<< i << "\r";
        if (static_cast<uint64_t>(bwt_string[i]) == END_CHARACTER)
            continue;
	    auto& bit_vec = *occs[alphamap[static_cast<uint64_t>(bwt_string[i])]];
        bit_vec[i] = 1;
    }
    std::cerr<< length << "\n";
    std::cerr<<"All Occ bit vectors are built.\n";

    for (auto& occ: occs) {
        std::cerr<< occs_rank.size() << "\r";
        if (verbose and (*occ).size() < 1000)
            std::cerr<< *occ << "\n";
        occs_rank.emplace_back(std::unique_ptr<sdsl::rank_support_v<> >(new sdsl::rank_support_v<>(occ.get())));
    }
    std::cerr<< occs_rank.size() << "\n";
    std::cerr<<"All Occ rank vectors are built.\n";


    // Building the move structure rows
    uint64_t len = 0;
    uint64_t bwt_row = 0;
    uint64_t r_idx = 0;
    uint64_t offset = 0;
    uint64_t max_len = 0;
    std::cerr<< "bits.size(): " << bits.size() << "\n";
    std::cerr<< "rank_support_v<>(&bits)(bits.size()): " << sdsl::rank_support_v<>(&bits)(bits.size()) << "\n";
    sbits = sdsl::select_support_mcl<>(&bits);
    for (uint64_t i = 0; i < length; i++) {
        if (i % 10000 == 0)
            std::cerr<< i << "\r";
        if (bwt_string[i] == static_cast<unsigned char>(END_CHARACTER) ) {
            end_bwt_idx = r_idx;
        }

        if (i == length - 1 or bwt_string[i] != bwt_string[i+1]) {
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
            

            if (pp_id == 0) {
                offset = 0;
            } else {
                // check the boundaries before performing select
                if (pp_id >= r) {
                    std::cerr << "pp_id: " << pp_id << "r: " << r << "i: " << i << "bwt_row: " << bwt_row << "lf: " << lf << "\n";
                    exit(0);
                }
                if (lf < sbits(pp_id + 1)) {
                    std::cerr << lf << " " << sbits(pp_id + 1);
                    exit(0);
                }

                offset = lf - sbits(pp_id + 1);
            }
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
            // To take care of cases where length of the run 
            // does not fit in uint16_t
            if (len >= std::numeric_limits<uint16_t>::max()) {
                n_overflow.push_back(len);
                if (n_overflow.size() - 1 >= std::numeric_limits<uint16_t>::max()) {
                    std::cerr << "Warning: the number of runs with overflow n is beyond uint16_t! " << n_overflow.size() - 1 << "\n";
                }
                rlbwt[r_idx].set_n(n_overflow.size() - 1);
                // std::cerr << r_idx << " " << len << " is_overflow_n: " << rlbwt[r_idx].is_overflow_n() << "\n";
                rlbwt[r_idx].set_overflow_n();
                // std::cerr << r_idx << " " << len << " is_overflow_n: " << rlbwt[r_idx].is_overflow_n() << "\n";
            }
            if (offset >= std::numeric_limits<uint16_t>::max()) {
                offset_overflow.push_back(offset);
                if (offset_overflow.size() - 1 >= std::numeric_limits<uint16_t>::max()) {
                    std::cerr << "Warning: the number of runs with overflow offset is beyond uint16_t! " << offset_overflow.size() - 1 << "\n";
                }
                rlbwt[r_idx].set_offset(offset_overflow.size() - 1);
                // std::cerr << r_idx << " " << len << " is_overflow_offset: " << rlbwt[r_idx].is_overflow_offset() << "\n";
                rlbwt[r_idx].set_overflow_offset();
                // std::cerr << r_idx << " " << len << " is_overflow_offset: " << rlbwt[r_idx].is_overflow_offset() << "\n";
            }

            if (len > max_len)
                max_len = len;
            if (logs) {
                if (run_lengths.find(len) != run_lengths.end())
                    run_lengths[len] += 1;
                else
                    run_lengths[len] = 1;
            }

            if (!bit1){
                rlbwt[r_idx].set_c(bwt_string[i], alphamap);
            }

            if (bwt_string[i] == static_cast<unsigned char>(END_CHARACTER)) {
                eof_row = r_idx;
                bit1_after_eof = alphamap[bwt_string[i+1]];
            }

            bwt_row += len;
            len = 0;
            r_idx += 1;
        } else if (bwt_string[i] == bwt_string[i+1]) {
            len += 1;
        }
    }
    std::cerr<<"All the move rows are built!\n";
    std::cerr<<"Max len: " << max_len << "\n";

    // compute the thresholds
    uint64_t alphabet_thresholds[4];
    // initialize the start threshold at the last row
    for (uint64_t j = 0; j < 4; j++)
        alphabet_thresholds[j] = length;
    if (bit1) {
        /*for (uint32_t i = 0; i < rlbwt.size() - 1; ++i) {
            if (thresholds[i + 1] >= rlbwt[i].get_p() + get_n(i)) {
                rlbwt[i].threshold_1bit = get_n(i);
            } else if (thresholds[i + 1] < rlbwt[i].get_p()) {
                rlbwt[i].threshold_1bit = 0;
            } else {
                rlbwt[i].threshold_1bit = thresholds[i + 1] - rlbwt[i].get_p();
            }
            // rlbwt[i].threshold_1bit = thresholds[i + 1];
        }
        rlbwt[r - 1].threshold_1bit = rlbwt[r - 1].get_n();*/
    } else {
        uint64_t run_p = 0;
        for (uint64_t i = rlbwt.size() - 1; i > 0; --i) {
            if (i % 10000 == 0)
                std::cerr<<"i: " << i << "\r";
            char rlbwt_c = bit1 ? compute_char(i) : alphabet[rlbwt[i].get_c()];
            if (verbose and i >= rlbwt.size() - 10) 
                std::cerr << "i: " << i << "\n"
                    << "rlbwt[i].get_offset(): " << get_offset(i) << "\n "
                    << "get_n(i): " << get_n(i) << "\n"
                    << "thresholds[i]: " << thresholds[i] << " "
                    << "rlbwt_c: " << rlbwt_c << "\n";
            
            std::vector<uint64_t> current_thresholds;
            current_thresholds.resize(3);
            for (uint64_t j = 0; j < alphabet.size(); j++) {
                if (alphabet[j] == rlbwt_c) {
                    alphabet_thresholds[j] = thresholds[i];
                } else {
	                if (alphabet_thresholds[j] >= run_p + get_n(i)) {
                        // rlbwt[i].thresholds[j] = get_n(i);
                        if (rlbwt_c == END_CHARACTER) {
                            end_bwt_idx_thresholds[j] = get_n(i);
                            continue;
                        }
                        if (get_n(i) >= std::numeric_limits<uint16_t>::max()) {
                            rlbwt[i].set_overflow_thresholds();                            
                        }
                        rlbwt[i].thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = get_n(i);
                        current_thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = get_n(i);
                        if (alphamap_3[alphamap[rlbwt_c]][j] == 3) std::cerr << "error: " << alphamap_3[alphamap[rlbwt_c]][j] << "\n";
                    } else if (alphabet_thresholds[j] < run_p) {
                        // rlbwt[i].thresholds[j] = 0;
                        if (rlbwt_c == END_CHARACTER) {
                            end_bwt_idx_thresholds[j] = 0;
                            continue;
			            }
                        rlbwt[i].thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = 0;
                        current_thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = 0;
                        if (alphamap_3[alphamap[rlbwt_c]][j] == 3) std::cerr << "error: " << alphamap_3[alphamap[rlbwt_c]][j] << "\n";
                    } else {
                        // rlbwt[i].thresholds[j] = alphabet_thresholds[j] - run_p;
                        if (rlbwt_c == END_CHARACTER) {
                            end_bwt_idx_thresholds[j] = alphabet_thresholds[j] - run_p;
                            continue;
			            }
                        if (alphabet_thresholds[j] - run_p >= std::numeric_limits<uint16_t>::max()) {
                            rlbwt[i].set_overflow_thresholds();
                        }
                        rlbwt[i].thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = alphabet_thresholds[j] - run_p;
                        current_thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = alphabet_thresholds[j] - run_p;
                        if (alphamap_3[alphamap[rlbwt_c]][j] == 3) std::cerr << "error: " << alphamap_3[alphamap[rlbwt_c]][j] << "\n";
                    }

                    if (verbose and i >= rlbwt.size() - 10)
                        std::cerr << "\t j: \t" << j << " "
                            << "alphabet[j]: " << alphabet[j] << "  "
                            << "alphamap_3[alphamap[rlbwt_c]][j]: " << alphamap_3[alphamap[rlbwt_c]][j] << " "
                            << "alphabet_thresholds[j]: " << alphabet_thresholds[j] << " "
                            << "rlbwt[i].thresholds[j]:" << rlbwt[i].thresholds[alphamap_3[alphamap[rlbwt_c]][j]] << "\n";

                    // rlbwt[i].thresholds[j] = alphabet_thresholds[j];
                }
            }
            if (rlbwt[i].is_overflow_thresholds()) {
                rlbwt[i].thresholds[0] = thresholds_overflow.size();
                rlbwt[i].thresholds[1] = thresholds_overflow.size();
                rlbwt[i].thresholds[2] = thresholds_overflow.size();
                if (thresholds_overflow.size() >= std::numeric_limits<uint16_t>::max())
                    std::cerr << "Warning: the number of runs with overflow thresholds is beyond uint16_t! " << thresholds_overflow.size() << "\n";

                thresholds_overflow.push_back(current_thresholds);
            }
            run_p += get_n(i);
        }
        for (uint64_t j = 0; j < alphabet.size() - 1; j++) {
            rlbwt[0].thresholds[j] = 0;
        }
    }
    std::cerr<< length << "\n";
    // std::cerr<<"Computing the next ups and downs.\n";
    std::cerr<< "The move structure building is done.\n";
}

/*uint64_t MoveStructure::fast_forward(uint64_t pointer, uint64_t idx) {
    uint64_t idx_ = idx;
    if (verbose) 
        std::cerr << " \t \t pointer: " << pointer << " p + n:" <<  rlbwt[idx].get_p() + get_n(idx) << "\n";
    while (idx < r - 1 && pointer >= rlbwt[idx].get_p() + get_n(idx)) {
        idx += 1;
        if (verbose) std::cerr << "\t \t ff: +" << idx - idx_ << "\n";
    }
    return idx - idx_;
}*/

uint64_t MoveStructure::fast_forward(uint64_t& offset, uint64_t idx, uint64_t x) {
    uint64_t idx_ = idx;
    if (verbose) 
        std::cerr << " \t \t offset: " << offset << " n:" << get_n(idx) << "\n";
    while (idx < r - 1 && offset >= get_n(idx)) {
        offset -= get_n(idx);
        idx += 1;
        if (verbose) std::cerr << "\t \t ff offset based: +" << idx - idx_ << "\n";
    }
    return idx - idx_;
}

uint64_t MoveStructure::jump_up(uint64_t idx, char c) {
    if (idx == 0)
        return r;
    char row_c = bit1 ? compute_char(idx) : alphabet[rlbwt[idx].get_c()];
    uint32_t jump_count = 0;
    while (idx > 0 and row_c != c) {
        jump_count += 1;
        idx -= 1;
        row_c = bit1 ? compute_char(idx) : alphabet[rlbwt[idx].get_c()];
    }
    if (logs) {
        if (jumps.find(jump_count) != jumps.end())
            jumps[jump_count] += 1;
        else
            jumps[jump_count] = 1;
    }
    if (verbose) 
        std::cerr << "\t \t \t \t idx after the while in the jump" << idx << "\n";
    return (row_c == c) ? idx : r;
}

uint64_t MoveStructure::jump_down(uint64_t idx, char c) {
    if (idx == r - 1)
        return r;
    char row_c = bit1 ? compute_char(idx) : alphabet[rlbwt[idx].get_c()];
    uint32_t jump_count = 0;
    while (idx < r - 1 && row_c != c) {
        jump_count += 1;
        idx += 1;
        row_c = bit1 ? compute_char(idx) : alphabet[rlbwt[idx].get_c()];
    }
    if (logs) {
        if (jumps.find(jump_count) != jumps.end())
            jumps[jump_count] += 1;
        else
            jumps[jump_count] = 1;
    }
    if (verbose) 
        std::cerr << "\t \t \t \t idx after the while in the jump: " << idx << " " << c << " " << row_c << "\n";
    return (row_c == c) ? idx : r;
}

uint64_t MoveStructure::query_ms(MoveQuery& mq, bool random) {
    std::srand(time(0));
    std::string R = mq.query();
    int32_t pos_on_r = R.length() - 1;
    uint64_t idx = r - 1; // std::rand() % r; // r - 1
    if (verbose) std::cerr<< "Begin search from idx = " << idx << "\n";
    // uint64_t pointer = rlbwt[idx].get_p();
    uint64_t offset = get_n(idx) - 1;
    uint64_t match_len = 0;

    if (verbose)
        std::cerr << "beginning of the search:\n query: " << mq.query() << "\n";

    if (verbose)
        std::cerr << "idx(r-1): " << idx << " offset: " << offset << "\n";

    uint64_t ff_count = 0;
    while (pos_on_r > -1) {
        if (idx == r) std::cerr << idx << "\n";
        if (verbose)
            std::cerr<< "Searching position " << pos_on_r << " of the read:\n";

        auto& row = rlbwt[idx];
        uint64_t row_idx = idx;
        char row_c = bit1 ? compute_char(idx) : alphabet[row.get_c()];

        if (alphamap[static_cast<uint64_t>(R[pos_on_r])] == alphamap.size()) { // not to use map
            // The character from the read does not exist in the reference
            match_len = 0;
            mq.add_ms(match_len);
            pos_on_r -= 1;

            if (verbose)
                std::cerr<< "\t The character " << R[pos_on_r] << " does not exist.\n";
        } else if (row_c == R[pos_on_r]) {
            if (verbose)
                std::cerr<< "\t Cas1: It was a match. \n" << "\t Continue the search...\n";

            // Case 1
            match_len += 1;
            if (verbose)
                std::cerr<<"\t match: " << match_len << "\n";
            mq.add_ms(match_len - 1);
            pos_on_r -= 1;

            // idx = row.id;
            idx = row.get_id();
            // offset based: pointer = row.get_pp() + (pointer - row.get_p());
            offset = get_offset(row_idx) + offset;
            if (verbose)
                std::cerr << "\t row.id: " << row.get_id() << " row.get_n: " << get_n(row_idx) << "\n"
                          << "\t rlbwt[idx].get_n: " << get_n(idx) << "\n"
                          << "\t idx: " << idx << "\n" // " pointer: " << pointer << " pointer-p: " << pointer - rlbwt[idx].get_p() << "\n"
                          << "\t offset: " << offset << "\t row.get_offset(): " << get_offset(row_idx) << "\n";

            // if (idx < r - 1 && pointer >= rlbwt[idx].get_p() + get_n(idx)) {
            if (idx < r - 1 && offset >= get_n(idx)) {
                if (verbose)
                    std::cerr<<"\t fast forwarding: " << idx << "\n";
                // uint64_t idx_ = fast_forward(pointer, idx);
                uint64_t idx__ = fast_forward(offset, idx, 0);
                // if (idx_ != idx__) { std::cerr<< "\t ff doesn't match" << idx_ << " " << idx__ << "\n";}
                ff_count += idx__;
                idx += idx__;
            }

        } else {
            if (verbose)
                std::cerr<< "\t Case2: Not a match, looking for a match either up or down...\n";

            // Case 2
            // Jumping randomly up or down or with naive lcp computation
            uint64_t lcp = 0;
            bool up = random ? jump_randomly(idx, R[pos_on_r]) : 
                               jump_thresholds(idx, offset, R[pos_on_r]);
            //                   jump_naive_lcp(idx, pointer, R[pos_on_r], lcp);
            char c = bit1 ? compute_char(idx) : alphabet[rlbwt[idx].get_c()];
            if (verbose)
                std::cerr<< "\t up: " << up << " lcp: " << lcp << " idx: " << idx << " c:" << c << "\n";

            auto saved_idx = idx;
            // sanity check
            if (c == R[pos_on_r]) {
                // Observing a match after the jump
                // The right match_len should be:
                // min(new_lcp, match_len + 1)
                // But we cannot compute lcp here
                if (up) {
                    // offset based: pointer = rlbwt[idx].get_p() + get_n(idx) - 1;
                    offset = get_n(idx) - 1;
                } else {
                    // offset based: pointer = rlbwt[idx].get_p();
                    offset = 0;
                }
                match_len = random ? 0 : std::min(match_len, lcp);
                if (verbose)
                    std::cerr<<"\t idx: " << " offset: " << offset << "\n";
            } else {
                std::cerr << "\t \t This should not happen!\n";
                std::cerr << "\t \t r[pos]:" <<  R[pos_on_r] << " t[pointer]:" << c << "\n";
                std::cerr << "\t \t " << up << ", " << bit1 << ", " << R[pos_on_r] << ", " << pos_on_r << "\n";
                std::cerr << "\t \t ";
                for (int k = 10; k > 0; --k)
                    std::cerr << alphabet[rlbwt[idx - k].get_c()] << "-";
                for (int k = 0; k < 10; k++)
                    std::cerr << alphabet[rlbwt[idx + k].get_c()] << "-";
                std::cerr<<"\n";
                /*char c_1 = bit1 ? compute_char(saved_idx) : rlbwt_chars[saved_idx];
                char c_2 = bit1 ? compute_char(idx) : rlbwt[idx].get_c();
                std::cerr<<rlbwt[saved_idx].get_p() << "\n";
                std::cerr<<saved_idx << " - " << c_1 << " - " << idx << " - " << c_2 << "\n";
                std::cerr << "saved: " << rlbwt[saved_idx].get_p() << "\n";
                for (uint32_t k = 0; k < 4; k++) {
                    std::cerr<< alphabet[k] << " " << rlbwt[saved_idx].thresholds[k] << " ";
                }
                std::cerr<<"\n";
                std::cerr << "new: " << rlbwt[idx].get_p() << "\n";
                for (uint32_t k = 0; k < 4; k++) {
                    std::cerr<< alphabet[k] << " " << rlbwt[idx].thresholds[k] << " ";
                }
                std::cerr<<"\n";*/
                
                /*for (uint32_t k = idx - 15; k < idx+1; k++) {
                    std::cerr<< k << " " << rlbwt_chars[k] << "\n";
                    for (uint32_t l = 0; l < 4; l++)
                        std::cerr<< alphabet[l] << "\t\t" << rlbwt[k].thresholds[l] << "\t\t";
                    std::cerr<<"\n";
                }*/
                verbose = true;
                jump_thresholds(saved_idx, offset, R[pos_on_r]);
                exit(0);
            }
        }
    }
    return ff_count;
}

bool MoveStructure::jump_thresholds(uint64_t& idx, uint64_t offset, char r_char) {
    uint64_t saved_idx = idx;
    uint64_t alphabet_index = alphamap[static_cast<uint64_t>(r_char)];
    if (verbose)
        std::cerr<<"\t \t \t jumping with thresholds ... \n";
    char rlbwt_char = bit1 ? compute_char(idx) : alphabet[rlbwt[idx].get_c()];
    if (verbose)
        std::cerr<<"\t \t \t alphabet_index: " << alphabet_index << " r_char:" << r_char << " rlbwt_char:" << rlbwt_char << "\n";
    // if (r_char > rlbwt_char and rlbwt_char != static_cast<unsigned char>(END_CHARACTER))
    //    alphabet_index -= 1;
    if (verbose)
        std::cerr << "\t \t \t idx:" << idx << "\n" // " pointer:" << pointer << " pointer-p: " << pointer - rlbwt[idx].get_p() << "\n"
                  << "\t \t \t offset: " << offset << " threshold:" << get_thresholds(idx, alphabet_index) << "\n";

    if (!bit1) {
	    if (idx == end_bwt_idx) {
            if (verbose) std::cerr << "\t \t \t idx == end_bwt_idx" 
                                   << "\n\t \t \t idx: " << idx << " end_bwt_idx: " << end_bwt_idx << "\n";
            // if (pointer >= rlbwt[idx].get_p() + end_bwt_idx_thresholds[alphabet_index]) { // and idx != r-1) {
            if (offset >= end_bwt_idx_thresholds[alphabet_index]) { // and idx != r-1) {    
                idx = jump_down(saved_idx, r_char);
	            return false;
            } else {
                idx = jump_up(saved_idx, r_char);
	            return true;
            }
        }
        if (verbose) std::cerr << "\t \t \t rlbwt[idx].get_offset(): " << get_offset(idx) 
                               << " get_thresholds(idx, alphabet_index): " << get_thresholds(idx, alphabet_index) 
                               << "\n\t \t \t idx:" << idx << "\n";
        alphabet_index = alphamap_3[alphamap[rlbwt_char]][alphabet_index];

	    // if (pointer >= rlbwt[idx].get_p() + get_thresholds(idx, alphabet_index)) { // and idx != r-1) {
        // if (offset >= get_thresholds(idx, alphabet_index)) { // and idx != r-1) {
        if (offset >= get_thresholds(idx, alphabet_index)) {
            if (verbose)
                std::cerr<< "\t \t \t Jumping down with thresholds:\n";
            idx = jump_down(saved_idx, r_char);
            return false;
        } else {
            if (verbose)
                std::cerr<< "\t \t \t Jumping up with thresholds:\n";
            idx = jump_up(saved_idx, r_char);
            return true;
        }
    } else {
        /*if (pointer >= rlbwt[idx].get_p() + rlbwt[idx].threshold_1bit and idx != r-1) {
            if (verbose)
                std::cerr<< "Jumping down with thresholds:\n";
            idx = jump_down(saved_idx, r_char);
            return false;
        } else {
            if (verbose)
                std::cerr<< "Jumping up with thresholds:\n";
            idx = jump_up(saved_idx, r_char);
            return true;
        }*/
    }
    // TODO: default return?
    return false;
}

bool MoveStructure::jump_randomly(uint64_t& idx, char r_char) {
    uint64_t saved_idx = idx;
    uint64_t jump = std::rand() % 2; // To replace with ...
    bool up = false;

    if (verbose)
        std::cerr<<"idx before jump: " << idx << "\n";

    if ( (jump == 1 && idx > 0) or idx == r - 1) {
        if (verbose)
            std::cerr<< "Jumping up randomly:\n";

        // jumping up
        up = true;
        idx = jump_up(saved_idx, r_char);
        if (verbose)
            std::cerr<<"idx after jump: " << idx << "\n";
        char c = bit1 ? compute_char(idx) : alphabet[rlbwt[idx].get_c()];
        if (c != r_char) {
            if (verbose)
                std::cerr<< "Up didn't work, try jumping down:\n";

            // jump down
            up = false;
            idx = jump_down(saved_idx, r_char);
            if (verbose)
                std::cerr<<"idx after jump: " << idx << "\n";
        }
    } else {
        if (verbose)
            std::cerr<< "Jumping down randomly:\n";

        // jumping down
        up = false;
        idx = jump_down(saved_idx, r_char);
        if (verbose)
            std::cerr<<"idx after jump: " << idx << "\n";
        char c = bit1 ? compute_char(idx) : alphabet[rlbwt[idx].get_c()];
        if (c != r_char) {
            if (verbose)
                std::cerr<< "Down didn't work, try jumping up:\n";

            // jump up
            up = true;
            idx = jump_up(saved_idx, r_char);
            if (verbose)
                std::cerr<<"idx after jump: " << idx << "\n";
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
        std::cerr<< "down_idx: " << down_idx << " up_idx: " << up_idx << "\n";
        std::cerr<< "down_pointer: " << down_pointer << " up_pointer: " << up_pointer << "\n";
        std::cerr<< "down_lcp: " << down_lcp << " up_lcp: " << up_lcp << "\n";
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
        std::cerr<< "should not happen during naive lcp jump!\n";
        exit(0);
    }
}*/

void MoveStructure::serialize(char* output_dir) {
    mkdir(output_dir,0777);
    std::string fname = static_cast<std::string>(output_dir) + "/rlbwt.bin";
    std::ofstream fout(fname, std::ios::out | std::ios::binary);
    std::cerr<< "length: " << length << " r: " << r << " end_bwt_idx: " << end_bwt_idx << "\n";
    fout.write(reinterpret_cast<char*>(&length), sizeof(length));
    fout.write(reinterpret_cast<char*>(&r), sizeof(r));
    fout.write(reinterpret_cast<char*>(&end_bwt_idx), sizeof(end_bwt_idx));
    fout.write(reinterpret_cast<char*>(&end_bwt_idx_thresholds[0]), 4*sizeof(end_bwt_idx_thresholds[0]));

    uint64_t alphamap_size = alphamap.size();
    fout.write(reinterpret_cast<char*>(&alphamap_size), sizeof(alphamap_size));
    fout.write(reinterpret_cast<char*>(&alphamap[0]), alphamap.size()*sizeof(alphamap[0]));

    uint64_t alphabet_size = alphabet.size();
    fout.write(reinterpret_cast<char*>(&alphabet_size), sizeof(alphabet_size));    
    fout.write(reinterpret_cast<char*>(&alphabet[0]), alphabet.size()*sizeof(alphabet[0]));

    fout.write(reinterpret_cast<char*>(&bit1), sizeof(bit1));
    std::cerr<< "sizeof(rlbwt[0]): " << sizeof(rlbwt[0]) << "\n";
    fout.write(reinterpret_cast<char*>(&rlbwt[0]), rlbwt.size()*sizeof(rlbwt[0]));
    /*for (uint32_t i = 0; i < r; i++){
        fout.write(reinterpret_cast<char*>(&rlbwt[i].p), sizeof(rlbwt[i].p));
        fout.write(reinterpret_cast<char*>(&rlbwt[i].n), sizeof(rlbwt[i].n));
        fout.write(reinterpret_cast<char*>(&rlbwt[i].pp), sizeof(rlbwt[i].pp));
        fout.write(reinterpret_cast<char*>(&rlbwt[i].id), sizeof(rlbwt[i].id));
        fout.write(reinterpret_cast<char*>(&rlbwt[i].overflow_bits), sizeof(rlbwt[i].overflow_bits));
        if (bit1) {
            fout.write(reinterpret_cast<char*>(&rlbwt[i].threshold_1bit), sizeof(rlbwt[i].threshold_1bit));
        } else {
            for (uint32_t j = 0; j < alphabet.size() - 1; j ++) {
                fout.write(reinterpret_cast<char*>(&(rlbwt[i].thresholds[j])), sizeof(rlbwt[i].thresholds[j]));
            }
        }
    }*/
    // if (!bit1)
    //    fout.write(reinterpret_cast<char*>(&rlbwt_chars[0]), rlbwt_chars.size()*sizeof(rlbwt_chars[0]));
    uint64_t n_overflow_size = n_overflow.size();
    fout.write(reinterpret_cast<char*>(&n_overflow_size), sizeof(n_overflow_size));
    fout.write(reinterpret_cast<char*>(&n_overflow[0]), n_overflow.size()*sizeof(uint64_t));
    uint64_t offset_overflow_size = offset_overflow.size();
    fout.write(reinterpret_cast<char*>(&offset_overflow_size), sizeof(offset_overflow_size));
    fout.write(reinterpret_cast<char*>(&offset_overflow[0]), offset_overflow.size()*sizeof(uint64_t));
    uint64_t thresholds_overflow_size = thresholds_overflow.size();
    fout.write(reinterpret_cast<char*>(&thresholds_overflow_size), sizeof(thresholds_overflow_size));
    for (uint32_t i = 0; i < thresholds_overflow_size; i++)
        fout.write(reinterpret_cast<char*>(&thresholds_overflow[i][0]), 3*sizeof(thresholds_overflow[i][0]));

    size_t orig_size = orig_string.size();
    fout.write(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));
    fout.write(reinterpret_cast<char*>(&reconstructed), sizeof(reconstructed));

    fout.write(reinterpret_cast<char*>(&eof_row), sizeof(eof_row));
    fout.write(reinterpret_cast<char*>(&bit1_begin), sizeof(bit1_begin));
    fout.write(reinterpret_cast<char*>(&bit1_after_eof), sizeof(bit1_after_eof));

    fout.close();
}

void MoveStructure::deserialize(char* index_dir) {
    std::cerr << "verbose: " << verbose << "\n";
    std::string fname = static_cast<std::string>(index_dir) + "/rlbwt.bin";
    std::ifstream fin(fname, std::ios::in | std::ios::binary);
    fin.seekg(0, std::ios::beg); 
    std::cerr<< "length: " << length << " r: " << r << " end_bwt_idx: " << end_bwt_idx << "\n";
    fin.read(reinterpret_cast<char*>(&length), sizeof(length));
    fin.read(reinterpret_cast<char*>(&r), sizeof(r));
    fin.read(reinterpret_cast<char*>(&end_bwt_idx), sizeof(end_bwt_idx));
    fin.read(reinterpret_cast<char*>(&end_bwt_idx_thresholds[0]), 4*sizeof(end_bwt_idx_thresholds[0]));

    std::cerr<< "length: " << length << " r: " << r << " end_bwt_idx: " << end_bwt_idx << "\n";

    uint64_t alphamap_size;
    fin.read(reinterpret_cast<char*>(&alphamap_size), sizeof(alphamap_size));
    alphamap.resize(alphamap_size);
    fin.read(reinterpret_cast<char*>(&alphamap[0]), alphamap_size*sizeof(alphamap[0]));
    std::cerr<<"alphamap_size: " << alphamap_size << "\n";
    uint64_t alphabet_size;
    fin.read(reinterpret_cast<char*>(&alphabet_size), sizeof(alphabet_size));
    alphabet.resize(alphabet_size);
    fin.read(reinterpret_cast<char*>(&alphabet[0]), alphabet_size*sizeof(alphabet[0]));
    std::cerr<<"alphabet_size: " << alphabet_size << "\n";
    rlbwt.resize(r);
    fin.read(reinterpret_cast<char*>(&bit1), sizeof(bit1));
    fin.read(reinterpret_cast<char*>(&rlbwt[0]), r*sizeof(MoveRow));
    /*for (uint32_t i = 0; i < r; i++){
        fin.read(reinterpret_cast<char*>(&rlbwt[i].p), sizeof(rlbwt[i].p));
        fin.read(reinterpret_cast<char*>(&rlbwt[i].n), sizeof(rlbwt[i].n));
        fin.read(reinterpret_cast<char*>(&rlbwt[i].pp), sizeof(rlbwt[i].pp));
        fin.read(reinterpret_cast<char*>(&rlbwt[i].id), sizeof(rlbwt[i].id));
        fin.read(reinterpret_cast<char*>(&rlbwt[i].overflow_bits), sizeof(rlbwt[i].overflow_bits));
        if (bit1) {
            fin.read(reinterpret_cast<char*>(&rlbwt[i].threshold_1bit), sizeof(rlbwt[i].threshold_1bit));
        } else {
            for (uint32_t j = 0; j < alphabet.size() - 1; j ++){
                fin.read(reinterpret_cast<char*>(&(rlbwt[i].thresholds[j])), sizeof(uint16_t));
            }
        }
    }*/
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
        thresholds_overflow[i].resize(3);
        fin.read(reinterpret_cast<char*>(&thresholds_overflow[i][0]), 3*sizeof(uint64_t));
    }

    /*if (!bit1) {
        rlbwt_chars.resize(r);
        fin.read(reinterpret_cast<char*>(&rlbwt_chars[0]), r*sizeof(char));
    }*/
    std::cerr << "All the move rows are read.\n";

    /*bwt_string.resize(length);
    fin.read(reinterpret_cast<char*>(&bwt_string[0]), length);*/

    size_t orig_size;
    fin.read(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));

    /*orig_string.resize(orig_size);
    fin.read(reinterpret_cast<char*>(&orig_string[0]), orig_size);*/

    fin.read(reinterpret_cast<char*>(&reconstructed), sizeof(reconstructed));
    reconstructed = false;

    fin.read(reinterpret_cast<char*>(&eof_row), sizeof(eof_row));
    fin.read(reinterpret_cast<char*>(&bit1_begin), sizeof(bit1_begin));
    fin.read(reinterpret_cast<char*>(&bit1_after_eof), sizeof(bit1_after_eof));

    fin.close();
}
