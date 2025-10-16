#include <iostream>
#include <fstream>
#include <zlib.h>
#include <cstdint>
#include <vector>

#include "kseq.h"

#include "commons.hpp"

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

// Reads fasta file. Returns the total length in number of base pairs.
uint64_t read_fasta(const char* file_name,
                    std::ofstream& clean_fasta,
                    bool input_type,
                    bool rc,
                    bool separators,
                    std::vector<uint64_t>& doc_lengths) {
    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(file_name, "r"); // STEP 2: open the file handler
    seq = kseq_init(fp); // STEP 3: initialize seq
    uint64_t line = 0;
    uint64_t total_length = 0;

    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        line += 1;
        std::string seq_rc = "";
        seq_rc.resize(seq->seq.l);
        total_length += rc ? 2*seq->seq.l : seq->seq.l;

        for (size_t i = 0; i < seq->seq.l; i++) {
            auto c = seq->seq.s[i];
            switch (c) {
                case 'a': seq->seq.s[i] = 'A'; break;
                case 'c': seq->seq.s[i] = 'C'; break;
                case 'g': seq->seq.s[i] = 'G'; break;
                case 't': seq->seq.s[i] = 'T'; break;
                default: break;
            }
            // if (c != 65 && c != 67 && c != 71 && c != 84) seq->seq.s[i] = 'A';
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') seq->seq.s[i] = 'A';
            if (rc) {
                switch (seq->seq.s[i]) {
                    case 'A': seq_rc[seq->seq.l - 1 - i] = 'T'; break;
                    case 'C': seq_rc[seq->seq.l - 1 - i] = 'G'; break;
                    case 'G': seq_rc[seq->seq.l - 1 - i] = 'C'; break;
                    case 'T': seq_rc[seq->seq.l - 1 - i] = 'A'; break;
                    default: WARNING_MSG("The alphabet still includes non-ACTG after cleaning!"); exit(0);
                }
            }
        }

        if (separators) {
            // Adding the separators (%)
            clean_fasta << '>' << seq->name.s << '\n' << seq->seq.s << SEPARATOR << '\n';
            if (rc)
                clean_fasta << '>' << seq->name.s << "_rev_comp" << '\n' << seq_rc << SEPARATOR << '\n';
            total_length += rc ? 2 : 1;
        } else {
            clean_fasta << '>' << seq->name.s << '\n' << seq->seq.s << '\n';
            if (rc)
                clean_fasta << '>' << seq->name.s << "_rev_comp" << '\n' << seq_rc << '\n';
        }

        if (!input_type) {
            doc_lengths.push_back(total_length);
            total_length = 0;
        }
    }

    if (input_type) {
        doc_lengths.push_back(total_length);
    }
    kseq_destroy(seq);
    gzclose(fp);
    return total_length;
}

int main(int argc, char* argv[]) {
    if (argc < 2 or std::string(argv[1]) == "--help") {
        INFO_MSG("Usage: movi-prepare-ref <fasta_file> <clean_fasta_file> [list] [separators]");
        return 0;
    }
    // Fasta/q reader from http://lh3lh3.users.sourceforge.net/parsefastq.shtml
    bool input_type = (argc > 3 && std::string(argv[3]) == "list") || (argc > 4 && std::string(argv[4]) == "list");
    bool separators = (argc > 3 && std::string(argv[3]) == "separators") || (argc > 4 && std::string(argv[4]) == "separators");
    // TODO: rc flag settting should be fixed
    bool rc = (argc > 4 and std::string(argv[4]) == "fw") ? false : true;

    INFO_MSG("Reverse complement: " + std::to_string(rc));
    std::ofstream clean_fasta(static_cast<std::string>(argv[2]));
    // Length of each document
    std::vector<uint64_t> doc_lengths;
    if (input_type) {
        std::ifstream list_file(static_cast<std::string>(argv[1]));
        std::string fasta_file = "";
        while (std::getline(list_file, fasta_file)) {
            INFO_MSG(fasta_file);
            uint64_t length = read_fasta(fasta_file.data(), clean_fasta, input_type, rc, separators, doc_lengths);
        }
    } else {
        INFO_MSG("Input fasta file: " + static_cast<std::string>(argv[1]));
        INFO_MSG("Output fasta file: " + static_cast<std::string>(argv[2]));
        uint64_t length = read_fasta(argv[1], clean_fasta, input_type, rc, separators, doc_lengths);
    }
    clean_fasta.close();

    if (separators) {
        INFO_MSG("Separators between genomes are added");
    }
    INFO_MSG("The clean fasta with the reverse complement is stored at " + static_cast<std::string>(argv[2]));

    // doc_offsets file stores the start index of the next run (last index of current run + 1)
    std::ofstream doc_offsets(static_cast<std::string>(argv[2]) + ".doc_offsets");
    uint64_t cur_ind = 0;
    for (int i = 0; i < doc_lengths.size(); i++) {
        cur_ind += doc_lengths[i];
        doc_offsets << cur_ind << "\n";
    }
    doc_offsets.close();
    
    return 0;
}
