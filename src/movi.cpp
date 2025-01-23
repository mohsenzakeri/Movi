#include <cstdint>
#include <stdio.h>
#include <cstdio>
#include <chrono>
#include <cstddef>
#include <unistd.h>
#include <sys/stat.h>

#include <omp.h>
#include <sdsl/int_vector.hpp>
#include "cxxopts.hpp"

#include "utils.hpp"
#include "move_structure.hpp"
#include "move_query.hpp"
#include "read_processor.hpp"
#include "movi_options.hpp"
#include "movi_parser.hpp"
#include "classifier.hpp"
#include "batch_loader.hpp"

void build_ftab(MoveStructure& mv_, MoviOptions& movi_options) {
    if (movi_options.is_multi_ftab() and movi_options.get_ftab_k() > 1) {
        int max_ftab = movi_options.get_ftab_k();
        for (int i = 2; i <= max_ftab; i++) {
            movi_options.set_ftab_k(i);
            mv_.compute_ftab();
            mv_.write_ftab();
            std::cerr << "The ftab table for k = " << i << " is built and stored in the index directory.\n";
        }
    } else if (movi_options.get_ftab_k() > 1) {
        mv_.compute_ftab();
        mv_.write_ftab();
    }
}

void query(MoveStructure& mv_, MoviOptions& movi_options) {
    if (movi_options.get_ftab_k() != 0) {
        mv_.read_ftab();
        std::cerr<<"Ftab was read!\n";
    }

    omp_set_num_threads(movi_options.get_threads());
    omp_set_nested(0);

    if (!movi_options.no_prefetch()) {
        ReadProcessor rp(movi_options.get_read_file(), mv_, movi_options.get_strands(), movi_options.is_verbose(), movi_options.is_reverse());

        std::ifstream input_file (movi_options.get_read_file().c_str());

#pragma omp parallel
        {
        BatchLoader reader;

            // Iterates over batches of data until none left
            while (true) {
                bool valid_batch = true;
                #pragma omp critical // one reader at a time
                {
                    valid_batch = reader.loadBatch(input_file, 1000, 4*movi_options.get_strands());
                }
                if (!valid_batch) {
                    break;
                }

                if (movi_options.is_pml() or movi_options.is_zml() or movi_options.is_count()) {


#if TALLY_MODE
                    rp.process_latency_hiding_tally(reader);
#else
                    rp.process_latency_hiding(reader);
#endif
                } else if (movi_options.is_kmer()) {
                    rp.kmer_search_latency_hiding(movi_options.get_k(), reader);
                }
            }
        }
        rp.end_process();

    } else {
        // gzFile fp;
        // int l;
        // kseq_t* seq = open_kseq(fp, movi_options.get_read_file());
        std::ofstream costs_file;
        std::ofstream scans_file;
        std::ofstream fastforwards_file;
        std::ofstream report_file;
        std::string index_type = program();
        if (movi_options.is_logs()) {
            costs_file = std::ofstream(movi_options.get_read_file() + "." + index_type + ".costs");
            scans_file = std::ofstream(movi_options.get_read_file() + "." + index_type + ".scans");
            fastforwards_file = std::ofstream(movi_options.get_read_file() + "." + index_type + ".fastforwards");
        }

        Classifier classifier;
        if (movi_options.is_classify()) {
            classifier.initialize_report_file(movi_options);
        }

        uint64_t total_ff_count = 0;

        std::ofstream mls_file;
        std::ofstream count_file;
        if (!movi_options.is_stdout()) {
            if (movi_options.is_pml() or movi_options.is_zml())
                mls_file = std::ofstream(movi_options.get_read_file() + "." + index_type + "." + query_type(movi_options) + ".bin", std::ios::out | std::ios::binary);
            else if (movi_options.is_count())
                count_file = std::ofstream(movi_options.get_read_file() + "." + index_type + ".matches");
        }

        std::ifstream input_file (movi_options.get_read_file().c_str());
        uint64_t read_processed = 0;

        #pragma omp parallel
        {
            BatchLoader reader;

            // Iterates over batches of data until none left
            while (true) {
                bool valid_batch = true;
                #pragma omp critical // one reader at a time
                {
                    valid_batch = reader.loadBatch(input_file, 1000, 1);
                }
                if (!valid_batch) break;


                Read read_struct;
                bool valid_read = false;

                // while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
                // Iterates over reads in a single batch
                while (true) {
                    valid_read = reader.grabNextRead(read_struct);
                    if (!valid_read) break;

                    #pragma omp atomic
                    read_processed += 1 ;
                    #pragma omp critical
                    {
                        if (read_processed % 1000 == 0)
                            std::cerr << read_processed << "\r";
                    }

                    // std::string query_seq = seq->seq.s;
                    std::string query_seq = std::string(read_struct.seq);
                    MoveQuery mq;
                    if (movi_options.is_pml() or movi_options.is_zml()) {
                        if (movi_options.is_reverse())
                            std::reverse(query_seq.begin(), query_seq.end());
                        mq = MoveQuery(query_seq);
                        bool random_jump = false;
                        if (movi_options.is_pml()) {
                            total_ff_count += mv_.query_pml(mq, random_jump);
                        } else if (movi_options.is_zml()) {
                            total_ff_count += mv_.query_zml(mq);
                        }

                        std::vector<uint16_t> matching_lens = mq.get_matching_lengths();
                        if (movi_options.is_classify()) {
                            // classifier.classify(seq->name.s, matching_lens, movi_options);
                            classifier.classify(read_struct.id, matching_lens, movi_options);
                        }
                        #pragma omp critical
                        {
                            output_matching_lengths(movi_options.is_stdout(), mls_file, read_struct.id, matching_lens);
                        }
                    } else if (movi_options.is_count()) {
                        int32_t pos_on_r = query_seq.length() - 1;
                        // uint64_t match_count = mv_.backward_search(query_seq, pos_on_r);
                        // if (pos_on_r != 0) pos_on_r += 1;
                        mq = MoveQuery(query_seq);
                        uint64_t match_count = mv_.query_backward_search(mq, pos_on_r);

                        #pragma omp critical
                        {
                            output_counts(movi_options.is_stdout(), count_file, read_struct.id, query_seq.length(), pos_on_r, match_count);
                        }
                    } else if (movi_options.is_kmer()) {
                        mq = MoveQuery(query_seq);
                        mv_.query_all_kmers(mq, movi_options.is_kmer_count());
                    }

                    if (movi_options.is_logs()) {
                        #pragma omp critical
                        {
                            output_logs(costs_file, scans_file, fastforwards_file, read_struct.id, mq);
                        }
                    }
                }
            }
        }
        
        if (movi_options.is_pml() or movi_options.is_zml()) {
            std::cerr << "all fast forward counts: " << total_ff_count << "\n";
            if (!movi_options.is_stdout()) {
                mls_file.close();
            }
            std::cerr << "The output file for the matching lengths closed.\n";
        } else if (movi_options.is_count()) {
            if (!movi_options.is_stdout()) {
                count_file.close();
            }
            std::cerr << "The count file is closed.\n";
        } else if (movi_options.is_kmer()) {
            mv_.kmer_stats.print(movi_options.is_kmer_count());
        }
        if (movi_options.is_classify()) {
            classifier.close_report_file();
        }
        if (movi_options.is_logs()) {
            costs_file.close();
            scans_file.close();
            fastforwards_file.close();
        }
        // close_kseq(seq, fp);
    }
}

void view(MoviOptions& movi_options) {
    std::ifstream mls_file(movi_options.get_mls_file(), std::ios::in | std::ios::binary);
    mls_file.seekg(0, std::ios::beg);


    Classifier classifier;
    if (movi_options.is_classify()) {
        classifier.initialize_report_file(movi_options);
    }

    while (true) {
        uint16_t st_length = 0;
        mls_file.read(reinterpret_cast<char*>(&st_length), sizeof(st_length));
        if (mls_file.eof()) break;

        std::string read_name;
        read_name.resize(st_length);
        mls_file.read(reinterpret_cast<char*>(&read_name[0]), st_length);
        read_name.erase(std::find(read_name.begin(), read_name.end(), '\0'), read_name.end());
        std::cout << ">" << read_name << " \n";
        uint64_t mq_pml_lens_size = 0;
        mls_file.read(reinterpret_cast<char*>(&mq_pml_lens_size), sizeof(mq_pml_lens_size));
        std::vector<uint16_t> pml_lens;
        pml_lens.resize(mq_pml_lens_size);
        mls_file.read(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));
        for (int64_t i = mq_pml_lens_size - 1; i >= 0; i--) {
            std::cout << pml_lens[i] << " ";
        }
        std::cout << "\n";

        if (movi_options.is_classify()) {
            classifier.classify(read_name, pml_lens, movi_options);
        }

    }

    if (movi_options.is_classify()) {
        classifier.close_report_file();
    }
}

int main(int argc, char** argv) {
    try {

        MoviOptions movi_options;
        if (!parse_command(argc, argv, movi_options)) {
            return 0;
        }

        if (movi_options.is_stdout()) {
            // Disable sync for faster I/O
            std::ios_base::sync_with_stdio(false);

            // Untie std::cin from std::cout for better performance
            std::cin.tie(nullptr);

            // Define a buffer of 1 MB
            constexpr size_t BUFFER_SIZE = 1024 * 1024; // 1 MB
            char buffer[BUFFER_SIZE];
            // Set the custom buffer for std::cout
            std::cout.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);
        }

        std::string command = movi_options.get_command();
        if (command == "build") {
            MoveStructure mv_(&movi_options, SPLIT_ARRAY, CONSTANT_MODE);
            if (movi_options.if_verify()) {
                std::cerr << "Verifying the LF_move results...\n";
                mv_.verify_lfs();
            }
            mv_.serialize();
            build_ftab(mv_, movi_options);
            std::cerr << "The move structure is successfully stored at " << movi_options.get_index_dir() << "\n";
            if (movi_options.is_output_ids()) {
                mv_.print_ids();
            }
            Classifier classifier;
            // generate pml null database
            movi_options.set_pml();
            movi_options.set_generate_null_reads(true);
            classifier.generate_null_statistics(mv_, movi_options);
            // generate zml null database
            movi_options.set_zml();
            movi_options.set_generate_null_reads(false); // do not regenerate the null reads
            classifier.generate_null_statistics(mv_, movi_options);
        } else if (command == "query") {
            MoveStructure mv_(&movi_options);
            auto begin = std::chrono::system_clock::now();
            mv_.deserialize();
            auto end = std::chrono::system_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            std::fprintf(stderr, "Time measured for loading the index: %.3f seconds.\n", elapsed.count() * 1e-9);
            begin = std::chrono::system_clock::now();
            query(mv_, movi_options);
            end = std::chrono::system_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            std::fprintf(stderr, "Time measured for processing the reads: %.3f seconds.\n", elapsed.count() * 1e-9);
        } else if (command == "view") {
            view(movi_options);
        } else if (command == "rlbwt") {
            std::cerr << "The run and len files are being built.\n";
            MoveStructure mv_(&movi_options);
            mv_.build_rlbwt();
        } else if (command == "LF") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();
            std::cerr << "The move structure is read from the file successfully.\n";
            if (movi_options.get_LF_type() == "sequential")
                mv_.sequential_lf();
            else if (movi_options.get_LF_type() == "random")
                mv_.random_lf();
            else if (movi_options.get_LF_type() == "reconstruct")
                mv_.reconstruct_lf();
        } else if (command == "inspect") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();
            mv_.print_stats();
            // mv_.compute_run_lcs();
            // mv_.analyze_rows();
        } else if (command == "ftab") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();
            build_ftab(mv_, movi_options);
        } else if (command == "null") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();
            Classifier classifier;
            classifier.generate_null_statistics(mv_, movi_options);
        } else {
            const std::string message = "Invalid action: \"" + command + "\"";
            throw std::runtime_error(message);
        }

        return 0;

    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}