#include <cstdint>
#include <zlib.h>
#include <stdio.h>
#include <cstdio>
#include <chrono>
#include <cstddef>
#include <unistd.h>
#include <sys/stat.h>

#include "kseq.h"
#include <sdsl/int_vector.hpp>

#include "move_structure.hpp"
#include "move_query.hpp"

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::string command = argv[1];
    if (command == "build") {
        std::cerr<<"The move structure is being built.\n";
        bool mode = std::string(argv[2]) == "1bit" ? true : false;
        uint16_t splitting = std::string(argv[2]) == "split" ? 5 : 0;
        std::cerr << "splitting: " << splitting << "\n";
        std::cerr << "mode: " << mode << "\n";
        bool verbose = (argc > 5 and std::string(argv[5]) == "verbose");
        bool logs = (argc > 5 and std::string(argv[5]) == "logs");
        MoveStructure mv_(argv[3], mode, verbose, logs, splitting);
        std::cerr<<"The move structure is successfully built!\n";
        // mv_.reconstruct();
        // std::cerr<<"The original string is reconstructed.\n";
        // std::cerr<<"The original string is:\n" << mv_.R() << "\n";
        mv_.serialize(argv[4]);
        std::cerr<<"The move structure is successfully stored at " << argv[4] << "\n";
        if (logs) {
            std::ofstream rl_file(static_cast<std::string>(argv[4]) + "/run_lengths");
            for (auto& run_length : mv_.run_lengths) {
                rl_file <<run_length.first << "\t" << run_length.second << "\n";
            }
            rl_file.close();
        }
    } else if (command == "query") {
        bool verbose = (argc > 4 and std::string(argv[4]) == "verbose");
        bool logs = (argc > 4 and std::string(argv[4]) == "logs");
        bool reverse_query = (argc > 4 and std::string(argv[4]) == "reverse");
        std::cerr << verbose << " " << logs << "\n";
        MoveStructure mv_(verbose, logs); 
        auto begin = std::chrono::system_clock::now();
        mv_.deserialize(argv[2]);
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
        // std::ofstream pmls_file(static_cast<std::string>(argv[3]) + ".mpml");
        std::ofstream pmls_file(static_cast<std::string>(argv[3]) + ".mpml.bin", std::ios::out | std::ios::binary);
        uint64_t all_ff_count = 0;
        while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
            /*printf("name: %s\n", seq->name.s);
            if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            printf("seq: %s\n", seq->seq.s);
            if (seq->qual.l) printf("qual: %s\n", seq->qual.s); */

            std::string query_seq = seq->seq.s;
            // reverse for the null reads
            if (reverse_query)
                std::reverse(query_seq.begin(), query_seq.end());
            MoveQuery mq(query_seq);
            bool random_jump = false;
            // std::cerr << seq->name.s << "\n";
            all_ff_count += mv_.query_ms(mq, random_jump);

            /* pmls_file << ">" << seq->name.s << "\n";
            for (int64_t i = mq.ms_lens.size() - 1; i >= 0; i--) {
                pmls_file << mq.ms_lens[i] << " ";
            }
            pmls_file << "\n"; */
            uint16_t st_length = seq->name.m;
            pmls_file.write(reinterpret_cast<char*>(&st_length), sizeof(st_length));
            pmls_file.write(reinterpret_cast<char*>(&seq->name.s), st_length);
            uint64_t mq_ms_lens_size = mq.ms_lens.size();
            pmls_file.write(reinterpret_cast<char*>(&mq_ms_lens_size), sizeof(mq_ms_lens_size));
            pmls_file.write(reinterpret_cast<char*>(&mq.ms_lens[0]), mq_ms_lens_size * sizeof(mq.ms_lens[0]));
        }
        pmls_file.close();
        std::cerr<<"pmls file closed!\n";
        // printf("return value: %d\n", l);
        kseq_destroy(seq); // STEP 5: destroy seq
        std::cerr<<"kseq destroyed!\n";
        gzclose(fp); // STEP 6: close the file handler
        std::cerr<<"fp file closed!\n";
        std::cerr<<"all fast forward counts: " << all_ff_count << "\n";
        if (logs) {
            std::ofstream jumps_file(static_cast<std::string>(argv[3]) + ".jumps");
            for (auto& jump : mv_.jumps) {
                jumps_file <<jump.first << "\t" << jump.second << "\n";
            }
            jumps_file.close();
        }
    } else if (command == "rlbwt") {
        std::cerr<<"The run and len files are being built.\n";
        bool verbose = (argc > 3 and std::string(argv[3]) == "verbose");
        bool logs = (argc > 3 and std::string(argv[3]) == "logs");
        MoveStructure mv_(verbose, logs);
        mv_.build_rlbwt(argv[2]);
    } 
    /*else if (command == "LF") {
        bool verbose = (argc > 4 and std::string(argv[4]) == "verbose");
        MoveStructure mv_(verbose);
        auto begin = std::chrono::system_clock::now();
        mv_.deserialize(argv[2]);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::printf("Time measured for loading the index: %.3f seconds.\n", elapsed.count() * 1e-9);
        std::cerr << "The move structure is read from the file successfully.\n";

        // std::string bwt_filename = argv[3] + std::string(".bwt");
        // std::cerr << bwt_filename << "\n";
        // std::ifstream bwt_file(bwt_filename);
        mv_.all_lf_test();
        
	    // std::ofstream ff_counts_file(static_cast<std::string>(argv[3]) + ".ff_counts");
        // for (auto& ff_count : mv_.ff_counts) {
        //     ff_counts_file <<ff_count.first << "\t" << ff_count.second << "\n";
        // }
        // ff_counts_file.close();
    } else if (command == "randomLF") {
        bool verbose = (argc > 3 and std::string(argv[3]) == "verbose");
        MoveStructure mv_(verbose);
        auto begin = std::chrono::system_clock::now();
        mv_.deserialize(argv[2]);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::printf("Time measured for loading the ind`ex: %.3f seconds.\n", elapsed.count() * 1e-9);
        std::cerr << "The move structure is read from the file successfully.\n";

        // mv_.random_lf_test();
        std::cerr << mv_.random_lf_test() << "\n";
    } else if (command == "reconstruct") {
        bool verbose = (argc > 3 and std::string(argv[3]) == "verbose");
        MoveStructure mv_(verbose);
        auto begin = std::chrono::system_clock::now();
        mv_.deserialize(argv[2]);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::printf("Time measured for loading the index: %.3f seconds.\n", elapsed.count() * 1e-9);
        std::cerr << "The move structure is read from the file successfully.\n";

        mv_.reconstruct_move();
    }*/
}
