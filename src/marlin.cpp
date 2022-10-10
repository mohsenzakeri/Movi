
#include <zlib.h>
#include <stdio.h>
#include <chrono>

#include "kseq.h"

#include "move_structure.hpp"
#include "move_query.hpp"

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

int main(int argc, char* argv[]) {
    std::string command = argv[1];
    if (command == "build") {
        std::cerr<<"The move structure is being built.\n";
        bool mode = argv[2] == "1bit" ? true : false;
        bool verbose = (argc > 5 and std::string(argv[5]) == "verbose");
        MoveStructure mv_(argv[3], mode, verbose);
        std::cerr<<"The move structure is successfully built!\n";
        // mv_.reconstruct();
        // std::cerr<<"The original string is reconstructed.\n";
        // std::cerr<<"The original string is:\n" << mv_.reconstruct() << "\n";
        mv_.seralize(argv[4]);
        std::cerr<<"The move structure is successfully stored at " << argv[4] << "\n";
    } else if (command == "query") {
        bool verbose = (argc > 4 and std::string(argv[4]) == "verbose");
        MoveStructure mv_(verbose);
        auto begin = std::chrono::system_clock::now();
        mv_.deseralize(argv[2]);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::printf("Time measured for loading the index: %.3f seconds.\n", elapsed.count() * 1e-9);
        std::cerr << "The move structure is read from the file successfully.\n";
        // std::cerr<<"The original string is: " << mv_.reconstruct() << "\n";
        // std::string query = argv[3];

        // Fasta/q reader from http://lh3lh3.users.sourceforge.net/parsefastq.shtml
        gzFile fp;
        kseq_t *seq;
        int l;
        fp = gzopen(argv[3], "r"); // STEP 2: open the file handler
        seq = kseq_init(fp); // STEP 3: initialize seq
        std::ofstream pmls_file(static_cast<std::string>(argv[3]) + ".mpml");
        while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
            /*printf("name: %s\n", seq->name.s);
            if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            printf("seq: %s\n", seq->seq.s);
            if (seq->qual.l) printf("qual: %s\n", seq->qual.s); */

            std::string query_seq = seq->seq.s;
            MoveQuery mq(query_seq);
            bool random_jump = true;
            mv_.query_ms(mq, random_jump);
            pmls_file << seq->name.s << "\n";
            pmls_file << mq <<"\n";
        }
        pmls_file.close();
        // printf("return value: %d\n", l);
        kseq_destroy(seq); // STEP 5: destroy seq
        gzclose(fp); // STEP 6: close the file handler
    }
}