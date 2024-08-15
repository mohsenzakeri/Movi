#include <iostream>
#include <fstream>
#include <zlib.h>
#include <cstdint>

#include "kseq.h"

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

void read_fasta(const char* file_name, std::ofstream& clean_fasta, bool rc, bool kmer_mode) {
    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(file_name, "r"); // STEP 2: open the file handler
    seq = kseq_init(fp); // STEP 3: initialize seq
    uint64_t line = 0;

    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        line += 1;
        std::string seq_rc = "";
        seq_rc.resize(seq->seq.l);
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
                    default: std::cerr << "The alphabet still includes non-ACTG after cleaning!\n"; exit(0);
                }
            }
        }
        if (kmer_mode) {
            // Adding the separators (#) for the kmer mode
            clean_fasta << '>' << seq->name.s << '\n' << seq->seq.s << '#' << '\n';
            if (rc)
                clean_fasta << '>' << seq->name.s << "_rev_comp" << '\n' << seq_rc << '#' << '\n';

        } else {
            clean_fasta << '>' << seq->name.s << '\n' << seq->seq.s << '\n';
            if (rc)
                clean_fasta << '>' << seq->name.s << "_rev_comp" << '\n' << seq_rc << '\n';
        }
    }

    kseq_destroy(seq);
    gzclose(fp);
}

int main(int argc, char* argv[]) {
    // Fasta/q reader from http://lh3lh3.users.sourceforge.net/parsefastq.shtml
    bool input_type = (argc > 3 and std::string(argv[3]) == "list");
    bool kmer_mode = (argc > 3 and std::string(argv[3]) == "kmer");
    bool rc = (argc > 4 and std::string(argv[4]) == "fw") ? false : true;
    std::cerr << rc << "\n";
    std::ofstream clean_fasta(static_cast<std::string>(argv[2]));
    if (input_type) {
        std::ifstream list_file(static_cast<std::string>(argv[1]));
        std::string fasta_file = "";
        while (std::getline(list_file, fasta_file)) {
            std::cerr << fasta_file << "\n";
            read_fasta(fasta_file.data(), clean_fasta, rc, kmer_mode);
        }
    } else {
        std::cerr << argv[1] << "\n";
        std::cerr << argv[2] << "\n";
        read_fasta(argv[1], clean_fasta, rc, kmer_mode);
    }
    clean_fasta.close();

    std::cerr << "The clean fasta with the reverse complement is stored at " << static_cast<std::string>(argv[2]) << "\n";
    return 0;
}
