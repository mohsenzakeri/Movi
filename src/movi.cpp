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

        bool onebit = std::string(argv[2]) == "onebit" ? true : false;
        uint16_t splitting = std::string(argv[2]) == "split" ? 5 : 0;
        bool constant = std::string(argv[2]) == "constant" ? true : false;
        if (constant) {
            splitting = 5;
        }
        std::cerr << "splitting: " << splitting << "\n";
        std::cerr << "onebit: " << onebit << "\n";
        std::cerr << "constant: " << constant << "\n";

        bool verbose = (argc > 5 and std::string(argv[5]) == "verbose");
        bool logs = (argc > 5 and std::string(argv[5]) == "logs");

        MoveStructure mv_(argv[3], onebit, verbose, logs, splitting, constant);
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
    } else if (command == "query" or command == "query-onebit") {
        bool verbose = (argc > 4 and std::string(argv[4]) == "verbose");
        bool logs = (argc > 4 and std::string(argv[4]) == "logs");
        bool reverse_query = (argc > 4 and std::string(argv[4]) == "reverse");
        std::cerr << verbose << " " << logs << "\n";
        MoveStructure mv_(verbose, logs); 
        if (command == "query-onebit")
            mv_.set_onebit();
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
        std::ofstream costs_file;
        std::ofstream scans_file;
        std::ofstream fastforwards_file;
        if (logs) {
            costs_file = std::ofstream(static_cast<std::string>(argv[3]) + ".costs");
            scans_file = std::ofstream(static_cast<std::string>(argv[3]) + ".scans");
            fastforwards_file = std::ofstream(static_cast<std::string>(argv[3]) + ".fastforwards");
        }
        std::ofstream pmls_file(static_cast<std::string>(argv[3]) + ".mpml.bin", std::ios::out | std::ios::binary);
        uint64_t all_ff_count = 0;
        // uint64_t query_pml_tot_time = 0;
        // uint64_t iteration_tot_time = 0;
        while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
            /*printf("name: %s\n", seq->name.s);
            if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            printf("seq: %s\n", seq->seq.s);
            if (seq->qual.l) printf("qual: %s\n", seq->qual.s); */

            std::string query_seq = seq->seq.s;
            // reverse for the null reads
            if (reverse_query)
                std::reverse(query_seq.begin(), query_seq.end());
            // auto t1 = std::chrono::high_resolution_clock::now();    
            MoveQuery mq(query_seq);
            bool random_jump = false;
            // std::cerr << seq->name.s << "\n";
            all_ff_count += mv_.query_pml(mq, random_jump);
            // auto t2 = std::chrono::high_resolution_clock::now();
            // query_pml_tot_time += static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

            /* pmls_file << ">" << seq->name.s << "\n";
            for (int64_t i = mq.pml_lens.size() - 1; i >= 0; i--) {
                pmls_file << mq.pml_lens[i] << " ";
            }
            pmls_file << "\n"; */
            uint16_t st_length = seq->name.m;
            pmls_file.write(reinterpret_cast<char*>(&st_length), sizeof(st_length));
            pmls_file.write(reinterpret_cast<char*>(&seq->name.s[0]), st_length);
            auto& pml_lens = mq.get_pml_lens();
            uint64_t mq_pml_lens_size = pml_lens.size();
            pmls_file.write(reinterpret_cast<char*>(&mq_pml_lens_size), sizeof(mq_pml_lens_size));
            pmls_file.write(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));
            if (logs) {
                costs_file << ">" << seq->name.s << "\n";
                scans_file << ">" << seq->name.s << "\n";
                fastforwards_file << ">" << seq->name.s << "\n";
                for (auto& cost : mq.get_costs()) {
                    costs_file << cost.count() << " ";
                    // iteration_tot_time += static_cast<uint64_t>(mq.costs[i].count());
                }
                for (auto& scan: mq.get_scans()) {
                    scans_file << scan << " ";
                }
                for (auto& fast_forward : mq.get_fastforwards()) {
                    fastforwards_file << fast_forward << " ";
                }
                costs_file << "\n";
                scans_file << "\n";
                fastforwards_file << "\n";
            }
        }
        std::cerr<<"all fast forward counts: " << all_ff_count << "\n";
        // std::cerr << "query_pml_tot_time: " << query_pml_tot_time << "\n";
        // std::cerr << "iteration_tot_time: " << iteration_tot_time << "\n";
        if (logs) {
            costs_file.close();
            scans_file.close();
            fastforwards_file.close();
        }
        pmls_file.close();
        std::cerr<<"pmls file closed!\n";
        // printf("return value: %d\n", l);
        kseq_destroy(seq); // STEP 5: destroy seq
        std::cerr<<"kseq destroyed!\n";
        gzclose(fp); // STEP 6: close the file handler
        std::cerr<<"fp file closed!\n";

        /* if (logs) {
            std::ofstream jumps_file(static_cast<std::string>(argv[3]) + ".jumps");
            for (auto& jump : mv_.jumps) {
                jumps_file <<jump.first << "\t" << jump.second << "\n";
            }
            jumps_file.close();
        } */
    } else if (command == "rlbwt") {
        std::cerr<<"The run and len files are being built.\n";
        bool verbose = (argc > 3 and std::string(argv[3]) == "verbose");
        bool logs = (argc > 3 and std::string(argv[3]) == "logs");
        MoveStructure mv_(verbose, logs);
        mv_.build_rlbwt(argv[2]);
    } else if (command == "view") {
        std::cerr << command << "\n";
        std::cerr << argv[2] << "\n";

        std::string fname = static_cast<std::string>(argv[2]);
        std::ifstream pmls_file(fname, std::ios::in | std::ios::binary);
        pmls_file.seekg(0, std::ios::beg);
        while (true) {
            uint16_t st_length = 0;
            pmls_file.read(reinterpret_cast<char*>(&st_length), sizeof(st_length));
            if (pmls_file.eof()) break;

            std::string read_name;
            read_name.resize(st_length);
            pmls_file.read(reinterpret_cast<char*>(&read_name[0]), st_length);
            read_name.erase(std::find(read_name.begin(), read_name.end(), '\0'), read_name.end());
            std::cout << ">" << read_name << " \n";
            uint64_t mq_pml_lens_size = 0;
            pmls_file.read(reinterpret_cast<char*>(&mq_pml_lens_size), sizeof(mq_pml_lens_size));
            std::vector<uint16_t> pml_lens;
            pml_lens.resize(mq_pml_lens_size);
            pmls_file.read(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));
            for (int64_t i = mq_pml_lens_size - 1; i >= 0; i--) {
                std::cout << pml_lens[i] << " ";
            }
            std::cout << "\n";
        }
    } else if (command == "stats") {
        std::cerr << command << "\n";
        std::cerr << argv[2] << "\n";
        bool verbose = (argc > 3 and std::string(argv[3]) == "verbose");
        bool logs = (argc > 3 and std::string(argv[3]) == "logs");
        MoveStructure mv_(verbose, logs);
        mv_.deserialize(argv[2]);
        mv_.print_stats();
    }
    else if (command == "LF" or command == "randomLF" or command == "reconstruct") {
        bool verbose = (argc > 3 and std::string(argv[3]) == "verbose");
        bool logs = (argc > 3 and std::string(argv[3]) == "logs");
        MoveStructure mv_(verbose, logs);

        auto begin = std::chrono::system_clock::now();
        mv_.deserialize(argv[2]);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::printf("Time measured for loading the index: %.3f seconds.\n", elapsed.count() * 1e-9);
        std::cerr << "The move structure is read from the file successfully.\n";

        // std::string bwt_filename = argv[3] + std::string(".bwt");
        // std::cerr << bwt_filename << "\n";
        // std::ifstream bwt_file(bwt_filename);
        if (command == "LF")
            mv_.all_lf_test();
        else if (command == "randomLF")
            mv_.random_lf_test();
        else
            mv_.reconstruct_move();

        if (logs) {
            std::ofstream ff_counts_file(static_cast<std::string>(argv[3]) + ".ff_counts");
            for (auto& ff_count : mv_.ff_counts) {
                ff_counts_file <<ff_count.first << "\t" << ff_count.second << "\n";
            }
            ff_counts_file.close();
        }
    }
}
