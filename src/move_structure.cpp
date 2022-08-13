#include <sys/stat.h> 

#include "move_structure.hpp"

MoveStructure::MoveStructure(char* input_file, bool verbose_) {
    verbose = verbose_;
    reconstructed = false;
    std::ifstream bwt_file(input_file);
    build(bwt_file);
}

uint32_t MoveStructure::LF(uint32_t row_number) {
    uint32_t lf = 0;
    uint32_t alphabet_index = alphamap[static_cast<uint32_t>(bwt_string[row_number])];
    for (uint32_t i = 0; i < alphabet_index; i++) {
        lf += counts[i];
    }
    auto& occ_rank = *occs_rank[alphabet_index];
    lf += occ_rank(row_number);
    return lf;
}

std::string MoveStructure::reconstruct() {
    if (!reconstructed) {
        orig_string = "";
        orig_string += static_cast<char>(END_CHARACTER);
        for (uint32_t bwt_row = 0; bwt_row != end_bwt_row; bwt_row = LF(bwt_row)) {
            orig_string = bwt_string[bwt_row] + orig_string;
        }
        reconstructed = true;
    }
    return orig_string;
}

uint32_t MoveStructure::naive_lcp(uint32_t row1, uint32_t row2) {
    if (row1 >= length or row2 >= length)
        return 0;

    if (!reconstructed)
        orig_string = reconstruct();
    uint32_t lcp = 0;
    uint32_t sa1 = naive_sa(row1);
    uint32_t sa2 = naive_sa(row2);

    while (orig_string[sa1 + lcp] == orig_string[sa2 + lcp]) {
        lcp += 1;
    }

    return lcp;
}

uint32_t MoveStructure::naive_sa(uint32_t bwt_row) {
    uint32_t sa = 0;
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
    uint32_t all_chars_count = 256;
    alphamap.resize(all_chars_count);
    std::fill(alphamap.begin(), alphamap.end(), alphamap.size());
    std::vector<uint32_t> all_chars(all_chars_count, 0);
    uint32_t current_char = bwt_file.get();
    r = 1;
    while (current_char != EOF && current_char != 10) {
        if (bwt_string.length() > 0 && current_char != bwt_string.back())
            r += 1;
        bwt_string += current_char;
        all_chars[static_cast<uint32_t>(current_char)] += 1;

        current_char = bwt_file.get();
    }
    std::cerr<< "r: " << r << "\n";
    length = bwt_string.length();
    std::cerr<<"length: " << length << "\n";
    rlbwt.resize(r);

    // Building the auxilary structures
    uint32_t alphabet_index = 0;
    for (uint32_t i = 0; i < all_chars_count; i++) {
        if (all_chars[i] != 0) {
            auto current_char = static_cast<unsigned char>(i);
            if (verbose)
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

    for (uint32_t i = 0; i < length; i++) {
        auto x = alphamap[static_cast<uint32_t>(bwt_string[i])];
        if (x >= occs.size())
            std::cerr<<occs.size() << " " << x << "\n";
	    auto& bit_vec = *occs[alphamap[static_cast<uint32_t>(bwt_string[i])]];
        bit_vec[i] = 1;
    }

    for (auto& occ: occs) {
        if (verbose)
            std::cerr<< *occ << "\n";
        occs_rank.push_back(new sdsl::rank_support_v<>(occ));
    }

    // Building the move structure rows
    uint32_t len = 0;
    uint32_t offset = 0;
    uint32_t r_idx = 0;
    bits = sdsl::bit_vector(length, 0);
    for (uint32_t i = 0; i < length; i++) {
        if (bwt_string[i] == static_cast<char>(END_CHARACTER) ) {
            end_bwt_row = i;
        }

        if (i == length - 1 or bwt_string[i] != bwt_string[i+1]) {
            len += 1;
            uint32_t lf  = 0;
            if (alphamap[static_cast<uint32_t>(bwt_string[i])] != 0)
                lf = LF(offset);
            bits[offset] = 1;
            rlbwt[r_idx].init(offset, len, lf, i, bwt_string[i]);

            offset += len;
            len = 0;
            r_idx += 1;
        } else if (bwt_string[i] == bwt_string[i+1]) {
            len += 1;
        }
    }    

    if (verbose)
        std::cerr<<"bits: " << bits << "\n";
    rbits = sdsl::rank_support_v<>(&bits);
    uint32_t idx = 0;
    for (auto& row: rlbwt) {
        uint32_t set_bit = row.pp;
        // Walk back untill the first preceeding 1
        while (set_bit > 0 && bits[set_bit] == 0)
            set_bit -= 1;
        row.id = rbits(set_bit);
        if (verbose)
            std::cerr<< idx << " " << row.c << " offset: " << row.p << " len: " << 
                    row.n << " lf_pp: " << row.pp << " pp_id: " << row.id << "\n";
        idx += 1;
    }

}

uint32_t MoveStructure::fast_forward(uint32_t pointer, uint32_t idx) {
    if (verbose)
        std::cerr << idx << " + " << rlbwt[idx].p << " + " << rlbwt[idx].n << "\n";
    while (idx < r - 1 && pointer >= rlbwt[idx].p + rlbwt[idx].n)
        idx += 1;
    return idx;
}

uint32_t MoveStructure::jump_up(uint32_t idx, char c) {
    if (idx == 0)
        return r;
    while (idx > 0 and rlbwt[idx].c != c) {
        idx -= 1;
    }
    return (rlbwt[idx].c == c) ? idx : r;
}

uint32_t MoveStructure::jump_down(uint32_t idx, char c) {
    if (idx == r - 1)
        return r;
    while (idx < r - 1 && rlbwt[idx].c != c) {
        idx += 1;
    }
    return (rlbwt[idx].c == c) ? idx : r;
}

void MoveStructure::query_ms(MoveQuery& mq, bool random) {
    std::srand(time(0));
    std::string R = mq.query();
    int32_t pos_on_r = R.length() - 1;

    uint32_t idx = std::rand() % r;
    if (verbose) std::cerr<< "Begin search from idx = " << idx << "\n";
    // uint32_t idx = std::rand() % r;
    uint32_t pointer = rlbwt[idx].p;
    uint32_t match_len = 0;
    
    std::cerr<< "beginning of the search:\n";
    std::cerr<< "query: " << mq.query() << "\n";
    if (verbose)
        std::cerr<< "idx(r-1): " << idx << " pointer: " << pointer << "\n";
    while (pos_on_r > -1) {
        if (verbose)
            std::cerr<< "Searching position " << pos_on_r << " of the read:\n";

        auto& row = rlbwt[idx];
        if (alphamap[static_cast<uint32_t>(R[pos_on_r])] == alphamap.size()) { // not to use map
            // The character from the read does not exist in the reference
            match_len = 0;
            mq.add_ms(match_len);
            pos_on_r -= 1;

            if (verbose)
                std::cerr<< "The character " << R[pos_on_r] << " does not exist.\n";
        } else if (row.c == R[pos_on_r]) {
            if (verbose)
                std::cerr<< "It was a match. \n" << "Continue the search...\n";

            // Case 1
            match_len += 1;
            if (verbose)
                std::cerr<<"match: " << match_len << "\n";
            mq.add_ms(match_len);
            pos_on_r -= 1;

            idx = row.id;
            pointer = row.pp + (pointer - row.p);
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
            uint32_t lcp = 0;
            bool up = random ? jump_randomly(idx, R[pos_on_r]) : 
                               jump_naive_lcp(idx, pointer, R[pos_on_r], lcp);
            if (verbose)
                std::cerr<< "up: " << up << " lcp: " << lcp << " idx: " << idx << "\n";

            // sanity check
            if (rlbwt[idx].c == R[pos_on_r]) {
                // Observing a match after the jump
                // The right match_len should be:
                // min(new_lcp, match_len + 1)
                // But we cannot compute lcp here
                if (up)
                    pointer = rlbwt[idx].p + rlbwt[idx].n - 1;
                else
                    pointer = rlbwt[idx].p;
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

bool MoveStructure::jump_randomly(uint32_t& idx, char r_char) {
    uint32_t saved_idx = idx;
    uint32_t jump = std::rand() % 2;
    bool up = false;

    if ( (jump == 1 && idx > 0) or idx == r - 1) {
        if (verbose)
            std::cerr<< "Jumping up randomly:\n";

        // jumping up
        up = true;
        idx = jump_up(saved_idx, r_char);
        if (verbose)
            std::cerr<<"idx after jump: " << idx << "\n";
        if (rlbwt[idx].c != r_char) {
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
        if (rlbwt[idx].c != r_char) {
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

bool MoveStructure::jump_naive_lcp(uint32_t& idx, uint32_t pointer, char r_char, uint32_t& lcp) {
    uint32_t up_idx = jump_up(idx, r_char);
    uint32_t up_pointer = up_idx != r ? rlbwt[up_idx].p + rlbwt[up_idx].n - 1 : length;;
    uint32_t down_idx = jump_down(idx, r_char);
    uint32_t down_pointer = down_idx != r ? rlbwt[down_idx].p : length;
    uint32_t up_lcp = naive_lcp(pointer, up_pointer);
    uint32_t down_lcp = naive_lcp(pointer, down_pointer);
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

    uint32_t alphamap_size = alphamap.size();
    fout.write(reinterpret_cast<char*>(&alphamap_size), sizeof(alphamap_size));
    for (uint32_t alphamap_item : alphamap) {
        fout.write(reinterpret_cast<char*>(&alphamap_item), sizeof(alphamap_item));
    }

    for (uint32_t i = 0; i < r; i++) {
	    unsigned char curr_char = rlbwt[i].c;
        fout.write(reinterpret_cast<char*>(&curr_char), sizeof(curr_char));
        uint32_t curr_p = rlbwt[i].p;
        fout.write(reinterpret_cast<char*>(&curr_p), sizeof(curr_p));
        uint32_t curr_n = rlbwt[i].n;
        fout.write(reinterpret_cast<char*>(&curr_n), sizeof(curr_n));
        uint32_t curr_pp = rlbwt[i].pp;
        fout.write(reinterpret_cast<char*>(&curr_pp), sizeof(curr_pp));
        uint32_t curr_id = rlbwt[i].id;
        fout.write(reinterpret_cast<char*>(&curr_id), sizeof(curr_id));
    }

    fout.write(reinterpret_cast<char*>(&bwt_string[0]), length);
    size_t orig_size = orig_string.size();
    fout.write(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));
    fout.write(reinterpret_cast<char*>(&orig_string[0]), orig_size);
    fout.write(reinterpret_cast<char*>(&reconstructed), sizeof(reconstructed));

    fout.close();
}

void MoveStructure::deseralize(char* index_dir) {
    std::string fname = static_cast<std::string>(index_dir) + "/rlbwt.bin";
    std::ifstream fin(fname, std::ios::in | std::ios::binary);
    uint32_t a, b, c;
    fin.seekg(0, std::ios::beg); 
    fin.read(reinterpret_cast<char*>(&length), sizeof(length));
    fin.read(reinterpret_cast<char*>(&r), sizeof(r));
    fin.read(reinterpret_cast<char*>(&end_bwt_row), sizeof(end_bwt_row));
    std::cerr<< "length: " << length << " r: " << r << " end_bwt_row: " << end_bwt_row << "\n";
    uint32_t alphamap_size;
    fin.read(reinterpret_cast<char*>(&alphamap_size), sizeof(alphamap_size));
    alphamap.resize(alphamap_size);
    for (uint32_t i = 0; i < alphamap_size; i++) {
        uint32_t alphamap_item;
        fin.read(reinterpret_cast<char*>(&alphamap_item), sizeof(alphamap_item));
        alphamap[i] = alphamap_item;
    }

    unsigned char curr_char;
    uint32_t curr_p, curr_n, curr_pp, curr_id;
    for(uint32_t i = 0; i < r; i++) {
        if (i%10000 == 0) std::cerr<< i << "\r";
        fin.read(reinterpret_cast<char*>(&curr_char), sizeof(curr_char));
        fin.read(reinterpret_cast<char*>(&curr_p), sizeof(curr_p));
        fin.read(reinterpret_cast<char*>(&curr_n), sizeof(curr_n));
        fin.read(reinterpret_cast<char*>(&curr_pp), sizeof(curr_pp));
        fin.read(reinterpret_cast<char*>(&curr_id), sizeof(curr_id));
        move_row mr(curr_p, curr_n, curr_pp, curr_id, curr_char);
        rlbwt.push_back(mr);
        for (uint32_t j = 0; j < curr_n; j++)
            bwt_string += curr_char;
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

    fin.close();
}