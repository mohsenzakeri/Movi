#include <sys/stat.h> 

#include "move_structure.hpp"

MoveStructure::MoveStructure(char* input_file, bool bit1_, bool verbose_) {
    verbose = verbose_;
    bit1 = bit1_;
    reconstructed = false;
    std::ifstream bwt_file(input_file);
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
        row_c = idx%2 == 0 ? alphabet[bit1_begin] : alphabet[1 + (bit1_begin % 2)];
    } else {
        row_c = (idx - eof_row)%2 == 1 ? alphabet[bit1_after_eof] : alphabet[1 + (bit1_after_eof % 2)];
    }
    return row_c;
}

uint64_t MoveStructure::LF(uint64_t row_number) {
    uint64_t lf = 0;
    uint64_t alphabet_index = alphamap[static_cast<uint64_t>(bwt_string[row_number])];
    for (uint64_t i = 0; i < alphabet_index; i++) {
        lf += counts[i];
    }
    auto& occ_rank = *occs_rank[alphabet_index];
    lf += static_cast<uint64_t>(occ_rank(row_number));
    return lf;
}

std::string MoveStructure::reconstruct() {
    if (!reconstructed) {
        orig_string = "";
        orig_string += static_cast<char>(END_CHARACTER);
        for (uint64_t bwt_row = 0; bwt_row != end_bwt_row; bwt_row = LF(bwt_row)) {
            orig_string = bwt_string[bwt_row] + orig_string;
        }
        reconstructed = true;
    }
    return orig_string;
}

uint64_t MoveStructure::naive_lcp(uint64_t row1, uint64_t row2) {
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
}

uint64_t MoveStructure::naive_sa(uint64_t bwt_row) {
    uint64_t sa = 0;
    for (; bwt_row != end_bwt_row; bwt_row = LF(bwt_row)) {
        sa += 1;
    }
    return sa;
}

void MoveStructure::build(std::ifstream &bwt_file) {
    bwt_file.clear();
    bwt_file.seekg(0);

    // Reading the BWT from the file
    bwt_string = "";
    uint64_t all_chars_count = 256;
    alphamap.resize(all_chars_count);
    std::fill(alphamap.begin(), alphamap.end(), alphamap.size());
    std::vector<uint64_t> all_chars(all_chars_count, 0);
    uint64_t current_char = bwt_file.get();
    r = 1;
    // TODO Use a size based on the input size
    bits = sdsl::bit_vector(100000000000, 0);
    bits[0] = 1;
    while (current_char != EOF && current_char != 10) {
        if (r % 10000 == 0)
            std::cerr<< r << "\r";
        if (bwt_string.length() > 0 && current_char != bwt_string.back()) {
            r += 1;
            bits[bwt_string.length() + 1] = 1;
        }
        bwt_string += current_char;
        all_chars[current_char] += 1;

        current_char = bwt_file.get();
    }
    std::cerr<< "r: " << r << "\n";
    length = bwt_string.length();
    std::cerr<<"length: " << length << "\n";
    rlbwt.resize(r);
    if (!bit1)
        rlbwt_chars.resize(r);
    if (verbose and bits.size() < 1000)
        std::cerr<<"bits: " << bits << "\n";
    rbits = sdsl::rank_support_v<>(&bits);

    // Building the auxilary structures
    uint64_t alphabet_index = 0;
    for (uint64_t i = 0; i < all_chars_count; i++) {
        std::cerr<< i << "\r";
        if (all_chars[i] != 0) {
            auto current_char = static_cast<unsigned char>(i);
            std::cerr<< "i is " << i << "\t" << current_char 
                        << "\t" << all_chars[i] << "\n";

            alphabet.push_back(current_char);
            counts.push_back(all_chars[i]);
            alphamap[i] = alphabet_index;
            alphabet_index += 1;

            sdsl::bit_vector* new_bit_vector = new sdsl::bit_vector(length, 0);
	        occs.push_back(new_bit_vector);
        }
    }

    if (alphabet.size() == 3) {
        bit1 = true;
        bit1_begin = alphamap[bwt_string[0]];

    }
    std::cerr << "All the characters are indexed.\n";

    for (uint64_t i = 0; i < length; i++) {
        if (i % 10000 == 0)
            std::cerr<< i << "\r";
	    auto& bit_vec = *occs[alphamap[static_cast<uint64_t>(bwt_string[i])]];
        bit_vec[i] = 1;
    }
    std::cerr<< length << "\n";
    std::cerr<<"All Occ bit vectors are built.\n";

    for (auto& occ: occs) {
        std::cerr<< occs_rank.size() << "\r";
        if (verbose and (*occ).size() < 1000)
            std::cerr<< *occ << "\n";
        occs_rank.push_back(new sdsl::rank_support_v<>(occ));
    }
    std::cerr<< occs_rank.size() << "\n";
    std::cerr<<"All Occ rank vectors are built.\n";


    // Building the move structure rows
    uint16_t len = 0;
    uint64_t offset = 0;
    uint64_t r_idx = 0;
    bits = sdsl::bit_vector(length, 0);
    for (uint64_t i = 0; i < length; i++) {
        if (i % 10000 == 0)
            std::cerr<< i << "\r";
        if (bwt_string[i] == static_cast<char>(END_CHARACTER) ) {
            end_bwt_row = i;
        }

        if (i == length - 1 or bwt_string[i] != bwt_string[i+1]) {
            len += 1;
            uint64_t lf  = 0;
            if (alphamap[static_cast<uint64_t>(bwt_string[i])] != 0)
                lf = LF(offset);
            // bits[offset] = 1;
            uint64_t pp_id = rbits(lf);
            if (bits[lf] == 1)
                pp_id += 1;
            rlbwt[r_idx].init(offset, len, lf, pp_id);

            if (!bit1)
                rlbwt_chars[r_idx] = bwt_string[i];
            if (bwt_string[i] == alphabet[0]) {
                eof_row = rlbwt.size();
                bit1_after_eof = alphamap[bwt_string[i+1]];
            }

            offset += len;
            len = 0;
            r_idx += 1;
        } else if (bwt_string[i] == bwt_string[i+1]) {
            len += 1;
        }
    }

    std::cerr<< length << "\n";
    std::cerr<< "The first pass of move structure building is done.\n";
}

uint64_t MoveStructure::fast_forward(uint64_t pointer, uint64_t idx) {
    if (verbose)
        std::cerr << idx << " + " << rlbwt[idx].get_p() << " + " << rlbwt[idx].get_n() << "\n";
    while (idx < r - 1 && pointer >= rlbwt[idx].get_p() + rlbwt[idx].get_n()) 
        idx += 1;
    return idx;
}

uint64_t MoveStructure::jump_up(uint64_t idx, char c) {
    if (idx == 0)
        return r;
    char row_c = bit1 ? compute_char(idx) : rlbwt_chars[idx];
    while (idx > 0 and row_c != c) {
        idx -= 1;
        row_c = bit1 ? compute_char(idx) : rlbwt_chars[idx];
    }
    if (verbose) 
        std::cerr << "idx after the while in the jump" << idx << "\n";
    return (row_c == c) ? idx : r;
}

uint64_t MoveStructure::jump_down(uint64_t idx, char c) {
    if (idx == r - 1)
        return r;
    char row_c = bit1 ? compute_char(idx) : rlbwt_chars[idx];
    while (idx < r - 1 && row_c != c) {
        idx += 1;
        row_c = bit1 ? compute_char(idx) : rlbwt_chars[idx];
    }
    if (verbose) 
        std::cerr << "idx after the while in the jump: " << idx << " " << c << " " << row_c << "\n";
    return (row_c == c) ? idx : r;
}

void MoveStructure::query_ms(MoveQuery& mq, bool random) {
    std::srand(time(0));
    std::string R = mq.query();
    int32_t pos_on_r = R.length() - 1;

    uint64_t idx = std::rand() % r;
    if (verbose) std::cerr<< "Begin search from idx = " << idx << "\n";
    // uint64_t idx = std::rand() % r;
    uint64_t pointer = rlbwt[idx].get_p();
    uint64_t match_len = 0;
    
    if (verbose) {
        std::cerr<< "beginning of the search:\n";
        std::cerr<< "query: " << mq.query() << "\n";
    }
    if (verbose)
        std::cerr<< "idx(r-1): " << idx << " pointer: " << pointer << "\n";
    while (pos_on_r > -1) {
        if (idx == r) std::cerr << idx << "\n";
        if (verbose)
            std::cerr<< "Searching position " << pos_on_r << " of the read:\n";

        auto& row = rlbwt[idx];
        char row_c = bit1 ? compute_char(idx) : rlbwt_chars[idx];
            
        if (alphamap[static_cast<uint64_t>(R[pos_on_r])] == alphamap.size()) { // not to use map
            // The character from the read does not exist in the reference
            match_len = 0;
            mq.add_ms(match_len);
            pos_on_r -= 1;

            if (verbose)
                std::cerr<< "The character " << R[pos_on_r] << " does not exist.\n";
        } else if (row_c == R[pos_on_r]) {
            if (verbose)
                std::cerr<< "It was a match. \n" << "Continue the search...\n";

            // Case 1
            match_len += 1;
            if (verbose)
                std::cerr<<"match: " << match_len << "\n";
            mq.add_ms(match_len);
            pos_on_r -= 1;

            // idx = row.id;
            idx = row.get_id();
            // pointer = row.pp + (pointer - row.p);
            pointer = row.get_pp() + (pointer - row.get_p());
            if (verbose)
                std::cerr<<"Case 1 idx: " << idx << " pointer: " << pointer << "\n";
            idx = fast_forward(pointer, idx);
            if (verbose)
                std::cerr<<"fast forwarding: " << idx << "\n";
        } else {
            if (verbose)
                std::cerr<< "Not a match, looking for a match either up or down...\n";

            // Case 2
            // Jumping randomly up or down or with naive lcp computation
            uint64_t lcp = 0;
            bool up = random ? jump_randomly(idx, R[pos_on_r]) : 
                               jump_naive_lcp(idx, pointer, R[pos_on_r], lcp);
            if (verbose)
                std::cerr<< "up: " << up << " lcp: " << lcp << " idx: " << idx << "\n";

            // sanity check
            char c = bit1 ? compute_char(idx) : rlbwt_chars[idx];
            if (c == R[pos_on_r]) {
                // Observing a match after the jump
                // The right match_len should be:
                // min(new_lcp, match_len + 1)
                // But we cannot compute lcp here
                if (up)
                    pointer = rlbwt[idx].get_p() + rlbwt[idx].get_n() - 1;
                else
                    pointer = rlbwt[idx].get_p();
                match_len = random ? 0 : std::min(match_len, lcp);
                if (verbose)
                    std::cerr<<"Case 2 idx: " << idx << " pointer: " << pointer << "\n";
            } else {
                std::cerr << "This should not happen!\n";
                exit(0);
            }
        }
    }
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
        char c = bit1 ? compute_char(idx) : rlbwt_chars[idx];
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
        char c = bit1 ? compute_char(idx) : rlbwt_chars[idx];
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

bool MoveStructure::jump_naive_lcp(uint64_t& idx, uint64_t pointer, char r_char, uint64_t& lcp) {
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
}

void MoveStructure::seralize(char* output_dir) {
    mkdir(output_dir,0777);
    std::string fname = static_cast<std::string>(output_dir) + "/rlbwt.bin";
    std::ofstream fout(fname, std::ios::out | std::ios::binary);
    std::cerr<< "length: " << length << " r: " << r << " end_bwt_row: " << end_bwt_row << "\n";
    fout.write(reinterpret_cast<char*>(&length), sizeof(length));
    fout.write(reinterpret_cast<char*>(&r), sizeof(r));
    fout.write(reinterpret_cast<char*>(&end_bwt_row), sizeof(end_bwt_row));

    uint64_t alphamap_size = alphamap.size();
    fout.write(reinterpret_cast<char*>(&alphamap_size), sizeof(alphamap_size));
    fout.write(reinterpret_cast<char*>(&alphamap[0]), alphamap.size()*sizeof(alphamap[0]));

    uint64_t alphabet_size = alphabet.size();
    fout.write(reinterpret_cast<char*>(&alphabet_size), sizeof(alphabet_size));    
    fout.write(reinterpret_cast<char*>(&alphabet[0]), alphabet.size()*sizeof(alphabet[0]));

    fout.write(reinterpret_cast<char*>(&bit1), sizeof(bit1));
    fout.write(reinterpret_cast<char*>(&rlbwt[0]), rlbwt.size()*sizeof(rlbwt[0]));
    if (!bit1) {
        fout.write(reinterpret_cast<char*>(&rlbwt_chars[0]), rlbwt_chars.size()*sizeof(rlbwt_chars[0]));
    }
    fout.write(reinterpret_cast<char*>(&bwt_string[0]), length);
    size_t orig_size = orig_string.size();
    fout.write(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));
    fout.write(reinterpret_cast<char*>(&orig_string[0]), orig_size);
    fout.write(reinterpret_cast<char*>(&reconstructed), sizeof(reconstructed));

    fout.write(reinterpret_cast<char*>(&eof_row), sizeof(eof_row));
    fout.write(reinterpret_cast<char*>(&bit1_begin), sizeof(bit1_begin));
    fout.write(reinterpret_cast<char*>(&bit1_after_eof), sizeof(bit1_after_eof));

    fout.close();
}

void MoveStructure::deseralize(char* index_dir) {
    std::string fname = static_cast<std::string>(index_dir) + "/rlbwt.bin";
    std::ifstream fin(fname, std::ios::in | std::ios::binary);
    uint64_t a, b, c;
    fin.seekg(0, std::ios::beg); 
    fin.read(reinterpret_cast<char*>(&length), sizeof(length));
    fin.read(reinterpret_cast<char*>(&r), sizeof(r));
    fin.read(reinterpret_cast<char*>(&end_bwt_row), sizeof(end_bwt_row));
    std::cerr<< "length: " << length << " r: " << r << " end_bwt_row: " << end_bwt_row << "\n";

    uint64_t alphamap_size;
    fin.read(reinterpret_cast<char*>(&alphamap_size), sizeof(alphamap_size));
    alphamap.resize(alphamap_size);
    fin.read(reinterpret_cast<char*>(&alphamap[0]), alphamap_size*sizeof(alphamap[0]));

    uint64_t alphabet_size;
    fin.read(reinterpret_cast<char*>(&alphabet_size), sizeof(alphabet_size));
    alphabet.resize(alphabet_size);
    fin.read(reinterpret_cast<char*>(&alphabet[0]), alphabet_size*sizeof(alphabet[0]));

    rlbwt.resize(r);
    fin.read(reinterpret_cast<char*>(&bit1), sizeof(bit1));
    fin.read(reinterpret_cast<char*>(&rlbwt[0]), r*sizeof(MoveRow));
    if (!bit1) {
        rlbwt_chars.resize(r);
        fin.read(reinterpret_cast<char*>(&rlbwt_chars[0]), r*sizeof(char));
    }
    std::cerr << "All the move rows are read.\n";

    bwt_string.resize(length);
    fin.read(reinterpret_cast<char*>(&bwt_string[0]), length);
    size_t orig_size;
    fin.read(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));
    orig_string.resize(orig_size);
    fin.read(reinterpret_cast<char*>(&orig_string[0]), orig_size);
    fin.read(reinterpret_cast<char*>(&reconstructed), sizeof(reconstructed));
    reconstructed = false;

    fin.read(reinterpret_cast<char*>(&eof_row), sizeof(eof_row));
    fin.read(reinterpret_cast<char*>(&bit1_begin), sizeof(bit1_begin));
    fin.read(reinterpret_cast<char*>(&bit1_after_eof), sizeof(bit1_after_eof));

    fin.close();
}