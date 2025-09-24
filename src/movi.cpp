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

// Function to handle PML/ZML processing for a single read
uint64_t handle_pml_zml(MoveQuery& mq, MoviOptions& movi_options,
                        MoveStructure& mv_, OutputFiles& output_files, Classifier& classifier) {
    uint64_t ff_count = 0;

    if (movi_options.is_pml()) {
        ff_count += mv_.query_pml(mq);
    } else if (movi_options.is_zml()) {
        ff_count += mv_.query_zml(mq);
    }

    #pragma omp critical
    {
        if (movi_options.is_classify()) {
            std::vector<uint16_t> matching_lens;

            // TODO: Classify for the large pml lens is not supported yet.
            if (movi_options.is_small_pml_lens()) {
                for (uint32_t i = 0;  i < mq.get_matching_lengths().size(); i++) {
                    matching_lens.push_back(static_cast<uint16_t>(mq.get_matching_lengths()[i]));
                }
            }

            bool found = classifier.classify(mq.get_query_id(), matching_lens, movi_options);

            if (found and movi_options.is_filter() && !movi_options.is_no_output()) {
                output_read(mq);
            }
        }

        if (movi_options.write_output_allowed()) {
            output_base_stats(DataType::match_length, movi_options.write_stdout_enabled(), output_files.mls_file, mq);

            if (movi_options.is_get_sa_entries()) {
                output_base_stats(DataType::sa_entry, movi_options.write_stdout_enabled(), output_files.sa_entries_file, mq);
            }
        }
    }

    return ff_count;
}

// Function to handle count processing for a single read
void handle_count(MoveQuery& mq, MoviOptions& movi_options,
                  MoveStructure& mv_, OutputFiles& output_files) {
    int32_t pos_on_r = mq.query().length() - 1;
    // uint64_t match_count = mv_.backward_search(query_seq, pos_on_r);
    // if (pos_on_r != 0) pos_on_r += 1;
    uint64_t match_count = mv_.query_backward_search(mq, pos_on_r);

    if (movi_options.write_output_allowed()) {
        #pragma omp critical
        {
            output_counts(movi_options.write_stdout_enabled(), output_files.matches_file, mq.query().length(), pos_on_r, match_count, mq);
        }
    }
}

// Function to handle k-mer processing for a single read
void handle_kmer(MoveQuery& mq, MoviOptions& movi_options,
                 MoveStructure& mv_, OutputFiles& output_files) {
    mv_.query_all_kmers(mq, movi_options.is_kmer_count());

    if (!movi_options.is_kmer_count()) {
        if (movi_options.write_output_allowed()) {
            #pragma omp critical
            {
                output_kmers(movi_options.write_stdout_enabled(), output_files.kmer_file, mq.query().length() - movi_options.get_k() + 1, mq);
            }
        }
    }
}

// Helper function to setup input file stream
void setup_input_file(std::ifstream& input_file, const std::string& read_file) {
    if (read_file == "-") {
        input_file.copyfmt(std::cin);
        input_file.clear(std::cin.rdstate());
        input_file.basic_ios<char>::rdbuf(std::cin.rdbuf());
    } else {
        input_file.open(read_file.c_str());
    }
}

// Function to load color table/document sets
void load_color_table(MoveStructure& mv_, MoviOptions& movi_options) {
    auto begin = std::chrono::system_clock::now();

    if (movi_options.is_full_color()) {
        mv_.fill_run_offsets();
        std::string fname = movi_options.get_index_dir() + "/doc_pats.bin";
        mv_.deserialize_doc_pats(fname);
    } else {
        if (movi_options.is_doc_sets_vector_of_vectors()) {
            if (!movi_options.is_freq_compressed() and !movi_options.is_tree_compressed()) {
                std::string fname = movi_options.get_index_dir() + "/doc_sets.bin";
                mv_.deserialize_doc_sets(fname);
            } else if (movi_options.is_freq_compressed()) {
                std::string fname = movi_options.get_index_dir() + "/compress_doc_sets.bin";
                mv_.deserialize_doc_sets(fname);
            } else if (movi_options.is_tree_compressed()) {
                std::string fname = movi_options.get_index_dir() + "/tree_doc_sets.bin";
                mv_.deserialize_doc_sets(fname);
            }
        } else {
            mv_.deserialize_doc_sets_flat();
        }
        mv_.load_document_info();
    }
    mv_.initialize_classify_cnts();

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::fprintf(stderr, "Time measured for loading the document sets: %.3f seconds.\n", elapsed.count() * 1e-9);
}

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

void color(MoveStructure& mv_, MoviOptions& movi_options) {
    mv_.load_document_info();

    auto begin = std::chrono::system_clock::now();
    if (movi_options.is_full_color()) {
        // Build document patterns (full information)
        mv_.fill_run_offsets();
        mv_.build_doc_pats();
        std::cerr << "Done building document info for each BWT row" << std::endl;
        mv_.serialize_doc_pats(movi_options.get_index_dir() + "/doc_pats.bin");

        mv_.build_doc_sets();
        std::cerr << "Done building document sets" << std::endl;
        mv_.serialize_doc_sets(movi_options.get_index_dir() + "/doc_sets.bin");
    } else {
        if (!movi_options.is_compressed()) {
            mv_.fill_run_offsets();

            std::string doc_pats_name = movi_options.get_index_dir() + "/doc_pats.bin";
            std::ifstream doc_pats_file(doc_pats_name);
            if (doc_pats_file.good()) {
                mv_.deserialize_doc_pats(doc_pats_name);
            } else {
                std::cerr << "Doc patterns are not available, building... \n";

                auto begin = std::chrono::system_clock::now();
                mv_.build_doc_pats();
                auto end = std::chrono::system_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
                std::printf("Time measured for building the document patterns: %.3f seconds.\n", elapsed.count() * 1e-9);
                mv_.serialize_doc_pats(movi_options.get_index_dir() + "/doc_pats.bin");
            }
            std::cerr << "Done reading document pattern information" << std::endl;
            mv_.build_doc_sets();
            std::cerr << "Done building document sets" << std::endl;
            if (movi_options.is_doc_sets_vector_of_vectors()) {
                mv_.serialize_doc_sets(movi_options.get_index_dir() + "/doc_sets.bin");
            } else {
                mv_.flat_and_serialize_colors_vectors();
            }
        } else {
            mv_.deserialize_doc_sets(movi_options.get_index_dir() + "/doc_sets.bin");
            mv_.compress_doc_sets();
            mv_.serialize_doc_sets(movi_options.get_index_dir() + "/compress_doc_sets.bin");

            //mv_.build_doc_set_similarities();
            //mv_.build_tree_doc_sets();
            //mv_.serialize_doc_sets(movi_options.get_index_dir() + "/tree_doc_sets.bin");
        }
    }

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::printf("Time measured for building colors: %.3f seconds.\n", elapsed.count() * 1e-9);
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

        std::ifstream input_file;
        if (movi_options.get_read_file() == "-") {
            input_file.copyfmt(std::cin);
            input_file.clear(std::cin.rdstate());
            input_file.basic_ios<char>::rdbuf(std::cin.rdbuf());
        } else {
            input_file.open(movi_options.get_read_file().c_str());
        }

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

        OutputFiles output_files;

        // Open output files using the utility function
        open_output_files(movi_options, output_files);

        Classifier classifier;
        mv_.set_classifier(&classifier);
        mv_.set_output_files(&output_files);
        if (movi_options.is_classify() or movi_options.is_filter()) {
            classifier.initialize_report_file(movi_options);
        }

        uint64_t total_ff_count = 0;

        std::ifstream input_file;
        if (movi_options.get_read_file() == "-") {
            input_file.copyfmt(std::cin);
            input_file.clear(std::cin.rdstate());
            input_file.basic_ios<char>::rdbuf(std::cin.rdbuf());
        } else {
            input_file.open(movi_options.get_read_file().c_str());
        }

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

                    if (movi_options.is_reverse())
                        std::reverse(query_seq.begin(), query_seq.end());

                    MoveQuery mq = MoveQuery(query_seq);
                    mq.set_query_id(read_struct.id);

                    if (movi_options.is_pml() or movi_options.is_zml()) {
                        total_ff_count += handle_pml_zml(mq, movi_options, mv_, output_files, classifier);

                    } else if (movi_options.is_count()) {

                        handle_count(mq, movi_options, mv_, output_files);

                    } else if (movi_options.is_kmer()) {

                        handle_kmer(mq, movi_options, mv_, output_files);

                    }

                    if (movi_options.is_logs()) {
                        #pragma omp critical
                        {
                            output_logs(output_files.costs_file, output_files.scans_file, output_files.fastforwards_file, mq);
                        }
                    }
                }
            }
        }

        if (!movi_options.is_no_output() and !movi_options.is_filter()) {
            if (movi_options.is_pml() or movi_options.is_zml()) {
                std::cerr << "all fast forward counts: " << total_ff_count << "\n";
            } else if (movi_options.is_kmer()) {
                mv_.kmer_stats.print(movi_options.is_kmer_count());
            }

            if (movi_options.is_classify()) {
                classifier.close_report_file();
            }

        if (movi_options.is_classify()) {
            classifier.close_report_file();
        }

        // Close output files using the utility function
        close_output_files(movi_options, output_files);
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

        // TODO: There is a lot of duplicate code here to handle 16 and 32 bits pml variations.
        if (movi_options.is_small_pml_lens()) {
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

        } else if (movi_options.is_large_pml_lens()) {

            std::vector<uint64_t> pml_lens;
            pml_lens.resize(mq_pml_lens_size);
            mls_file.read(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));

            for (int64_t i = mq_pml_lens_size - 1; i >= 0; i--) {
                std::cout << pml_lens[i] << " ";
            }

            std::cout << "\n";

            // TODO: used for color offsets (instad of color ids)

        } else {
            std::vector<uint32_t> pml_lens;
            pml_lens.resize(mq_pml_lens_size);
            mls_file.read(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));

            for (int64_t i = mq_pml_lens_size - 1; i >= 0; i--) {
                std::cout << pml_lens[i] << " ";
            }

            std::cout << "\n";

            if (movi_options.is_classify()) {
                // Classification with 32 bits pmls works if the pml values are less than 2^16.
                std::vector<uint16_t> matching_lens_16(pml_lens.begin(), pml_lens.end());
                classifier.classify(read_name, matching_lens_16, movi_options);
            }
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
            return 1;
        }

        if (movi_options.is_stdout() and !movi_options.is_no_output()) {
            // Disable sync for faster I/O
            std::ios_base::sync_with_stdio(false);

            // Untie std::cin from std::cout for better performance
            std::cin.tie(nullptr);

            // Define a buffer of 1 MB
            constexpr size_t BUFFER_SIZE = 1024 * 1024; // 1 MB
            char buffer[BUFFER_SIZE];
            // Set the custom buffer for std::cout
            std::cout.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

            // Set custom buffer for stdout (printf uses stdout)
            if (setvbuf(stdout, buffer, _IOFBF, BUFFER_SIZE) != 0) {
                perror("Failed to set buffer for stdout");
                return 1;
            }
        }

        std::string command = movi_options.get_command();
        if (command == "build") {
            MoveStructure mv_(&movi_options, SPLIT_ARRAY, CONSTANT_MODE);
            if (movi_options.is_verify()) {
                std::cerr << "Verifying the LF_move results...\n";
                mv_.verify_lf_loop();
            }
            mv_.serialize();
            build_ftab(mv_, movi_options);
            std::cerr << "The move structure is successfully stored at " << movi_options.get_index_dir() << "\n\n";
            if (movi_options.is_output_ids()) {
                mv_.print_ids();
            }

            std::cerr << "Generating the null statistics...\n\n";
            Classifier classifier;

            // generate pml null database
            std::cerr << "With PML:\n";
            movi_options.set_pml();
            movi_options.set_generate_null_reads(true);
            classifier.generate_null_statistics(mv_, movi_options);

            // generate zml null database
            std::cerr << "\nWith ZML:\n";
            movi_options.set_zml();
            movi_options.set_generate_null_reads(false); // do not regenerate the null reads
            classifier.generate_null_statistics(mv_, movi_options);

            if (movi_options.is_color()) {
                color(mv_, movi_options);
            }

        } else if (command == "build-SA") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();
            mv_.find_sampled_SA_entries();
            mv_.serialize_sampled_SA();
        } else if (command == "color") {
            MoveStructure mv_(&movi_options);
            auto begin = std::chrono::system_clock::now();
            mv_.deserialize();
            auto end = std::chrono::system_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            std::printf("Time measured for loading the index: %.3f seconds.\n", elapsed.count() * 1e-9);

            color(mv_, movi_options);
        } else if (command == "query") {
            MoveStructure mv_(&movi_options);

            auto begin = std::chrono::system_clock::now();
            mv_.deserialize();
            auto end = std::chrono::system_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            std::fprintf(stderr, "Time measured for loading the index: %.3f seconds.\n", elapsed.count() * 1e-9);

            if (movi_options.is_get_sa_entries()) {
                begin = std::chrono::system_clock::now();
                mv_.deserialize_sampled_SA();
                end = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
                std::fprintf(stderr, "Time measured for deserializing the sampled SA entries: %.3f seconds.\n", elapsed.count() * 1e-9);
            }

            if (movi_options.is_multi_classify()) {
                load_color_table(mv_, movi_options);

                if (movi_options.is_report_colors() or movi_options.is_report_color_ids()) {
                    // Copmuter color ids for the output
                    mv_.compute_color_ids_from_flat();
                }
            }

            begin = std::chrono::system_clock::now();
            query(mv_, movi_options);
            end = std::chrono::system_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

            std::fprintf(stderr, "Time measured for processing the reads: %.3f seconds.\n", elapsed.count() * 1e-9);

            // Avoid taking too long for dealloction of large data structures at the end of the program
            std::quick_exit(0);

        } else if (command == "view") {
            view(movi_options);
        } else if (command == "rlbwt") {
            std::cerr << "The run and len files are being built.\n";
            MoveStructure mv_(&movi_options);
            mv_.build_rlbwt();
        } else if (command == "color-move-rows") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();

            load_color_table(mv_, movi_options);

            mv_.add_colors_to_rlbwt();

            mv_.serialize();

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
            if (movi_options.is_flat_color_vectors()) {
                std::string fname = movi_options.get_index_dir() + "/doc_sets.bin";
                mv_.deserialize_doc_sets(fname);
                mv_.load_document_info();
                std::cerr << "The color table is read successfully.\n";
                mv_.flat_and_serialize_colors_vectors();
                std::cerr << "The flat color table is serialized successfully in the index directory (doc_sets_flat.bin).\n";
            }
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