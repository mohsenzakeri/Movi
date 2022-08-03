#include <fstream>

#include "move_structure.hpp"

MoveStructure::MoveStructure(char* input_file) {
    std::ifstream bwt_file(input_file);
    build(bwt_file);
}

uint32_t MoveStructure::LF(uint32_t row_number) {
    uint32_t lf = 0;
    uint32_t alphabet_index = alphamap[bwt_string[row_number]];
    for (uint32_t i = 0; i < alphabet_index; i++) {
        lf += counts[i];
    }
    auto& occ_rank = *occs_rank[alphabet_index];
    lf += occ_rank(row_number);    
    return lf;
}

void MoveStructure::build(std::ifstream &bwt_file) {
    bwt_file.clear();
    bwt_file.seekg(0);

    // Reading the BWT from the file
    bwt_string = "";
    uint32_t all_chars_count = 256;
    std::vector<uint32_t> all_chars(all_chars_count, 0);
    uint32_t current_char = bwt_file.get();
    while (current_char != EOF)
    {
        bwt_string += current_char;

        all_chars[static_cast<uint32_t>(current_char)] += 1;
        current_char = bwt_file.get();
    }
    length = bwt_string.length();
    std::cerr<<"length: " << length << "\n";

    // Building auxilary structures
    uint32_t alphabet_index = 0;
    for (uint32_t i = 0; i < all_chars_count; i++) {
        if (all_chars[i] != 0) {
            auto current_char = static_cast<unsigned char>(i);
            std::cerr<< i << "\t" << current_char << "\t" << all_chars[i] << "\n";

            alphabet.push_back(current_char);
            counts.push_back(all_chars[i]);
            alphamap[current_char] = alphabet_index;
            alphabet_index += 1;

            sdsl::bit_vector* new_bit_vector = new sdsl::bit_vector(length, 0);
	        occs.push_back(new_bit_vector);
        }
    }
    for (uint32_t i = 0; i < length; i++) {
	    auto& bit_vec = *occs[alphamap[bwt_string[i]]];
        bit_vec[i] = 1;
    }
    for (auto& occ: occs) {
        occs_rank.push_back(new sdsl::rank_support_v<>(occ));
    }

    // Building the move structure rows
    uint32_t len = 0;
    uint32_t offset = 0;
    sdsl::bit_vector bits(length, 0);
    for (uint32_t i = 0; i < length; i++) {
        if (i == length - 1 or bwt_string[i] != bwt_string[i+1]) {
            len += 1;
            uint32_t lf  = 0;
            if (alphamap[bwt_string[i]] != 0)
                lf = LF(offset);
            move_row* row = new move_row(offset, len, lf, i, bwt_string[i]);
            bits[offset] = 1;
            rlbwt.push_back(row);
            offset += len;
            len = 0;
        }
        else if (bwt_string[i] == bwt_string[i+1])
            len += 1;
    }

    sdsl::rank_support_v<> rbits(&bits);
    for (auto& row: rlbwt) {
        uint32_t set_bit = row->pp;
        while (set_bit > 0 and bits[set_bit] == 0)
            set_bit -= 1;
        if (set_bit == 0 and bits[set_bit] == 0)
            row->id = 0;
        else
            row->id = rbits(set_bit) - 1;
        // std::cerr<< row->c << " offset: " << row->p << " len: " << row->n << 
        //                     " lf_pp: " << row->pp << " pp_id: " << row->id << "\n";
    }

}