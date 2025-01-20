#include <cstdint>
#include <stdio.h>
#include <cstdio>
#include <chrono>
#include <cstddef>
#include <unistd.h>
#include <sys/stat.h>

#include <sdsl/int_vector.hpp>
#include "cxxopts.hpp"

#include "classifier.hpp"
#include "utils.hpp"
#include "move_structure.hpp"
#include "move_query.hpp"
#include "read_processor.hpp"
#include "movi_options.hpp"
#include "movi_parser.hpp"

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

    if (!movi_options.no_prefetch()) {
        ReadProcessor rp(movi_options.get_read_file(), mv_, movi_options.get_strands(), movi_options.is_verbose(), movi_options.is_reverse());
        if (movi_options.is_pml() or movi_options.is_zml() or movi_options.is_count()) {
#if TALLY_MODE
            rp.process_latency_hiding_tally();
#else
            rp.process_latency_hiding();
#endif
        } else if (movi_options.is_kmer()) {
            rp.kmer_search_latency_hiding(movi_options.get_k());
        }
    } else {
        gzFile fp;
        int l;
        kseq_t* seq = open_kseq(fp, movi_options.get_read_file());
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

        uint64_t read_processed = 0;
        while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
            if (read_processed % 1000 == 0)
                std::cerr << read_processed << "\r";
            read_processed += 1 ;
            std::string query_seq = seq->seq.s;
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
                if (movi_options.is_stdout()) {
                    std::cout << ">" << seq->name.s << " \n";
                    auto& pml_lens = mq.get_matching_lengths();
                    uint64_t mq_pml_lens_size = pml_lens.size();
                    for (int64_t i = mq_pml_lens_size - 1; i >= 0; i--) {
                        std::cout << pml_lens[i] << " ";
                    }
                    std::cout << "\n";
                } else {
                    uint16_t st_length = seq->name.m;
                    mls_file.write(reinterpret_cast<char*>(&st_length), sizeof(st_length));
                    mls_file.write(reinterpret_cast<char*>(&seq->name.s[0]), st_length);
                    auto& pml_lens = mq.get_matching_lengths();
                    uint64_t mq_pml_lens_size = pml_lens.size();
                    mls_file.write(reinterpret_cast<char*>(&mq_pml_lens_size), sizeof(mq_pml_lens_size));
                    mls_file.write(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));
                }
                if (movi_options.is_classify()) {
                    std::vector<uint16_t> matching_lens = mq.get_matching_lengths();
                    classifier.classify(seq->name.s, matching_lens, movi_options);
                }
            } else if (movi_options.is_count()) {
                int32_t pos_on_r = query_seq.length() - 1;
                // uint64_t match_count = mv_.backward_search(query_seq, pos_on_r);
                // if (pos_on_r != 0) pos_on_r += 1;
                mq = MoveQuery(query_seq);
                uint64_t match_count = mv_.query_backward_search(mq, pos_on_r);
                if (movi_options.is_stdout()) {
                    std::cout << seq->name.s << "\t";
                    std::cout << query_seq.length() - pos_on_r << "/" << query_seq.length() << "\t" << match_count << "\n";
                } else {
                    count_file << seq->name.s << "\t";
                    count_file << query_seq.length() - pos_on_r << "/" << query_seq.length() << "\t" << match_count << "\n";
                }
            } else if (movi_options.is_kmer()) {
                mq = MoveQuery(query_seq);
                mv_.query_all_kmers(mq, movi_options.is_kmer_count());
            }

            if (movi_options.is_logs()) {
                costs_file << ">" << seq->name.s << "\n";
                scans_file << ">" << seq->name.s << "\n";
                fastforwards_file << ">" << seq->name.s << "\n";
                for (auto& cost : mq.get_costs()) {
                    costs_file << cost.count() << " ";
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
        close_kseq(seq, fp);
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