#include "utils.hpp"

uint32_t alphamap_3[4][4] = {{3, 0, 1, 2},
                             {0, 3, 1, 2},
                             {0, 1, 3, 2},
                             {0, 1, 2, 3}};

std::string program() {
#if MODE == 0
    return "large"; // old name: "default"
#endif
#if MODE == 1
    return "constant";
#endif
#if MODE == 4 // Like the default constant mode, but without the pointers to the neighbors with the other characters
    return "split";
#endif
#if MODE == 3
    return "regular"; // old name: compact
#endif
#if MODE == 6
    return "regular-thresholds"; // old name: compact-thresholds
#endif
#if MODE == 2
    return "blocked";
#endif
#if MODE == 8
    return "blocked-thresholds";
#endif
#if MODE == 5
    return "tally";
#endif
#if MODE == 7
    return "tally-thresholds";
#endif
    std::cerr << "The mode is not defined.";
    exit(0);
}

char complement(char c) {
    // # is the separator, complement(#) = #
    char c_comp = c == '#' ? '#' : (c == 'A' ? 'T' : ( c == 'C' ? 'G' : (c == 'G' ? 'C' : 'A')));
    return c_comp;
}

std::string reverse_complement(std::string& fw) {
    std::string rc = "";
    rc.resize(fw.size());
    for (int i = fw.size() - 1;  i >= 0; i--) {
        rc[fw.size() - i - 1] = complement(fw[i]);
    }
    return rc;
}

std::string reverse_complement_from_pos(MoveQuery& mq_fw, int32_t pos_on_r, uint64_t match_len) {
    std::string rc = "";
    rc.resize(match_len);
    for (int i = pos_on_r + match_len - 1;  i >= pos_on_r; i--) {
        rc[pos_on_r + match_len - 1 - i] = complement(mq_fw.query()[i]);
    }
    return rc;
}

std::string number_to_kmer(size_t j, size_t m, std::vector<unsigned char>& alphabet, std::vector<uint64_t>& alphamap) {
    std::string kmer = "";
    for (int i = m - 2; i >= 0; i -= 2) {
        size_t pair = (j >> i) & 0b11;
        // std::cerr << i << " " << pair << " ";
        kmer += alphabet[pair + alphamap['A']];
    }
    return kmer;
}

uint64_t kmer_to_number(size_t k, std::string& r, int32_t pos, std::vector<uint64_t>& alphamap, bool rc) {
    if (r.length() < k) {
        std::cerr << "The k does not match the kmer length!\n";
        exit(0);
    }

    uint64_t res = 0;
    for (size_t i = 0; i < k; i++) {
        if (alphamap[static_cast<uint64_t>(r[pos + i])] == alphamap.size()) {
            return std::numeric_limits<uint64_t>::max();
        }
        if (rc) {
            uint64_t char_code = alphamap[complement(r[pos + i])] - alphamap['A'];
            res = (char_code << ((i)*2)) | res;
        } else {
            uint64_t char_code = alphamap[r[pos + i]] - alphamap['A'];
            res = (char_code << ((k-i-1)*2)) | res;
        }
    }
    return res;
}

uint8_t F_char(std::vector<uint64_t>& first_runs, uint64_t run) {
    for (uint8_t i = first_runs.size() - 1; i > 0; i--) {
        if (run >= first_runs[i]) {
            return i - 1;
        }
    }
    std::cerr << "Undefined behaviour.\n";
    exit(0);
}