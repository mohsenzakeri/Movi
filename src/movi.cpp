#include <cstdint>
#include <stdio.h>
#include <cstdio>
#include <chrono>
#include <cstddef>
#include <unistd.h>
#include <sys/stat.h>

#include <omp.h>
#include "sdsl_wrapper.hpp"
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
            std::vector<uint16_t> matching_lens_16(mq.get_matching_lengths().begin(), mq.get_matching_lengths().end());
            // Classification with 32 bits pmls works if the pml values are less than 2^16.

            bool found = classifier.classify(mq.get_query_id(), matching_lens_16, movi_options);

            if (movi_options.is_filter() && !movi_options.is_no_output()) {
                if (found && !movi_options.is_invert()) {
                    output_read(mq);
                } else if (!found && movi_options.is_invert()) {
                    output_read(mq);
                }
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

// Function to handle mem processing for a single read
void handle_mem(MoveQuery& mq, MoviOptions& movi_options,
                MoveStructure& mv_, OutputFiles& output_files) {
    mv_.query_mems(mq);
    if (movi_options.write_output_allowed()) {
        #pragma omp critical
        {
            output_mems(movi_options.write_stdout_enabled(), output_files.mems_file, mq);
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
        // This check is already done in the parser too.
        if (!std::filesystem::exists(read_file)) {
            throw std::runtime_error(ERROR_MSG("The input file " + read_file + " does not exist."));
        }
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
    TIMING_MSG(elapsed, "loading the document sets");
}

void build_ftab(MoveStructure& mv_, MoviOptions& movi_options) {
    if (movi_options.is_multi_ftab() and movi_options.get_ftab_k() > 1) {
        int max_ftab = movi_options.get_ftab_k();
        for (int i = 2; i <= max_ftab; i++) {
            movi_options.set_ftab_k(i);
            mv_.build_ftab();
            mv_.write_ftab();
            SUCCESS_MSG("The ftab table for k = " + std::to_string(i) + " is built and stored in the index directory.");
        }
    } else if (movi_options.get_ftab_k() > 1) {
        mv_.build_ftab();
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
        mv_.serialize_doc_pats(movi_options.get_index_dir() + "/doc_pats.bin");

        mv_.build_doc_sets();
        SUCCESS_MSG("Done building document sets.");
        mv_.serialize_doc_sets(movi_options.get_index_dir() + "/doc_sets.bin");
    } else {
        if (!movi_options.is_compressed()) {
            mv_.fill_run_offsets();

            std::string doc_pats_name = movi_options.get_index_dir() + "/doc_pats.bin";
            std::ifstream doc_pats_file(doc_pats_name);
            if (doc_pats_file.good()) {
                mv_.deserialize_doc_pats(doc_pats_name);
                INFO_MSG("Done reading document pattern information");
            } else {
                INFO_MSG("Doc patterns are not available, building...");

                auto begin = std::chrono::system_clock::now();
                mv_.build_doc_pats();
                auto end = std::chrono::system_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
                TIMING_MSG(elapsed, "building the document patterns");
                mv_.serialize_doc_pats(movi_options.get_index_dir() + "/doc_pats.bin");
            }
            mv_.build_doc_sets();
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
    TIMING_MSG(elapsed, "building colors");
}

void query(MoveStructure& mv_, MoviOptions& movi_options) {

    if (movi_options.get_ftab_k() != 0) {
        mv_.read_ftab();
        INFO_MSG("Ftab was read!");
    }

#if TALLY_MODES
    if (movi_options.is_zml() or movi_options.is_count()) {
        //TODO: Implement tally modes for zml and count queries.
        throw std::runtime_error(ERROR_MSG("This query type is not suppported by the sampled modes yet."));
    }
#endif

    if (movi_options.is_mem()) {
        if (!movi_options.no_prefetch()) {
            movi_options.set_prefetch(false);
            WARNING_MSG("MEM finding does not support prefetching. Continuing with prefetching disabled.");
        }
        if (movi_options.get_ftab_k() == 0) {
            throw std::runtime_error(ERROR_MSG("MEM finding requires ftab. Please build the ftab using the ./movi ftab --ftab-k <k>, then pass --ftab-k <k> to the query step."));
        } else {
            if (movi_options.get_min_mem_length() > movi_options.get_ftab_k()) {
                WARNING_MSG("Setting minimum MEM (length " + std::to_string(movi_options.get_min_mem_length()) + ") greater than ftab k (" + std::to_string(movi_options.get_ftab_k()) + ") causes a slower MEM search.");
            }
        }
    }

    omp_set_num_threads(movi_options.get_threads());
    omp_set_nested(0);

    // This is for the no-prefetch mode (the prefetch mode has its own classifier)
    Classifier classifier;
    mv_.set_classifier(&classifier);
    if (movi_options.is_classify()) {
        classifier.initialize_report_file(movi_options);
    }

    std::ifstream input_file;
    setup_input_file(input_file, movi_options.get_read_file());

    OutputFiles output_files;
    open_output_files(movi_options, output_files);
    mv_.set_output_files(&output_files);

    uint64_t total_ff_count = 0;

    auto begin = std::chrono::system_clock::now();

    if (!movi_options.no_prefetch()) {

        ReadProcessor rp(mv_, movi_options.get_strands(), movi_options.is_verbose(), movi_options.is_reverse(), output_files, classifier);

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


#if TALLY_MODES
                    rp.process_latency_hiding_tally(reader);
#else
                    rp.process_latency_hiding(reader);
#endif
                } else if (movi_options.is_kmer()) {
                    rp.kmer_search_latency_hiding(movi_options.get_k(), reader);
                }
            }
        }

        // TODO: total ff_count is not correct in the prefetch mode.
        total_ff_count += rp.get_total_ff_count();
        rp.end_process();

        SUCCESS_MSG(format_number_with_commas(rp.get_read_processed()) + " reads are processed.");

    } else {
        if (!movi_options.is_kmer()) {
            // For kmer queries, latency hiding is disabled by default.
            INFO_MSG("Latency hiding is disabled...");
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

                // Iterates over reads in a single batch
                while (true) {
                    valid_read = reader.grabNextRead(read_struct);
                    if (!valid_read) break;

                    #pragma omp atomic
                    read_processed += 1 ;

                    #pragma omp critical
                    {
                        if (read_processed % 1000 == 0) {
                            QUERY_PROGRESS_MSG("Number of reads processed: " + format_number_with_commas(read_processed));
                        }
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


                    } else if (movi_options.is_mem()) {
                        handle_mem(mq, movi_options, mv_, output_files);
                    }

                    if (movi_options.is_logs()) {
                        if (movi_options.write_output_allowed()) {
                            #pragma omp critical
                            {
                                output_logs(output_files.costs_file, output_files.scans_file, output_files.fastforwards_file, mq);
                            }
                        }
                    }
                }
            }
        }
        SUCCESS_MSG(format_number_with_commas(read_processed) + " reads are processed.");

    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    TIMING_MSG(elapsed, "processing the reads");

    if (movi_options.write_output_allowed()) {
        print_query_stats(movi_options, total_ff_count, mv_);
    }

    if (movi_options.is_classify()) {
        classifier.close_report_file();
    }

    close_output_files(movi_options, output_files);
}

void view(MoviOptions& movi_options) {
    std::ifstream mls_file(movi_options.get_mls_file(), std::ios::in | std::ios::binary);
    if (!mls_file.good()) {
        throw std::runtime_error(ERROR_MSG("Failed to open the MLS file: " + movi_options.get_mls_file()));
    }

    mls_file.seekg(0, std::ios::beg);

    Classifier classifier;
    if (movi_options.is_classify()) {
        classifier.initialize_report_file(movi_options);
    }

    uint8_t entry_size = 32;
    if (!movi_options.is_no_header()) {
        // Read BPF header
        BPFHeader header;
        mls_file.read(reinterpret_cast<char*>(&header), sizeof(header));
        if (header.magic != BPF_MAGIC) {
            throw std::runtime_error("Invalid BPF header.");
        }
        if (header.version != BPF_VERSION_MAJOR) {
            throw std::runtime_error("Invalid BPF version.");
        }
        entry_size = header.entry_size;
    } else {
        if (movi_options.is_small_pml_lens()) {
            entry_size = 16;
        } else if (movi_options.is_large_pml_lens()) {
            entry_size = 64;
        }
    }

    while (true) {
        uint16_t st_length = 0;
        mls_file.read(reinterpret_cast<char*>(&st_length), sizeof(st_length));
        if (mls_file.eof()) break;

        std::string read_name;
        read_name.resize(st_length);
        mls_file.read(reinterpret_cast<char*>(&read_name[0]), st_length);
        read_name.erase(std::find(read_name.begin(), read_name.end(), '\0'), read_name.end());
        std::cout << ">" << read_name << "\n";
        uint64_t mq_pml_lens_size = 0;
        mls_file.read(reinterpret_cast<char*>(&mq_pml_lens_size), sizeof(mq_pml_lens_size));

        // TODO: There is a lot of duplicate code here to handle 16 and 32 bits pml variations.
        if (entry_size == 16) {
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

        } else if (entry_size == 64) {

            std::vector<uint64_t> pml_lens;
            pml_lens.resize(mq_pml_lens_size);
            mls_file.read(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));

            for (int64_t i = mq_pml_lens_size - 1; i >= 0; i--) {
                std::cout << pml_lens[i] << " ";
            }

            std::cout << "\n";

            // TODO: used for color offsets (instad of color ids)

        } else if (entry_size == 32) {
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
        } else {
            throw std::runtime_error("Invalid BPF entry size.");
        }

    }

    if (movi_options.is_classify()) {
        classifier.close_report_file();
    }
}

void build_rlbwt(MoviOptions& movi_options) {
    INFO_MSG("The run and len files are being built.");

    std::ifstream bwt_file(movi_options.get_bwt_file());
    if (!bwt_file.good()) {
        throw std::runtime_error(ERROR_MSG("[build_rlbwt] Failed to open the BWT file: " + movi_options.get_bwt_file()));
    }

    bwt_file.clear();
    bwt_file.seekg(0,std::ios_base::end);
    std::streampos end_pos = bwt_file.tellg();

    if (movi_options.is_verbose()) {
        INFO_MSG("end_pos: " + std::to_string(end_pos));
    }

    bwt_file.seekg(0);
    char current_char = bwt_file.get();
    char last_char = current_char;
    uint64_t r = 0;
    size_t len = 0;

    std::ofstream len_file(movi_options.get_bwt_file() + ".len", std::ios::out | std::ios::binary);
    if (!len_file.good()) {
        throw std::runtime_error(ERROR_MSG("[build_rlbwt] Failed to open the length file: " + movi_options.get_bwt_file() + ".len"));
    }

    std::ofstream heads_file(movi_options.get_bwt_file() + ".heads");
    if (!heads_file.good()) {

        throw std::runtime_error(ERROR_MSG("[build_rlbwt] Failed to open the heads file: " + movi_options.get_bwt_file() + ".heads"));
    }

    while (current_char != EOF) {
        if (r % 1000000 == 0)
            PROGRESS_MSG(std::to_string(r));
        if (current_char != last_char) {
            r += 1;
            // write output
            heads_file << last_char;
            len_file.write(reinterpret_cast<char*>(&len), 5);
            len = 0;
        }
        len += 1;
        last_char = current_char;
        current_char = bwt_file.get();
    }

    // write output
    heads_file << last_char;
    len_file.write(reinterpret_cast<char*>(&len), 5);

    heads_file.close();
    len_file.close();
}

int main(int argc, char** argv) {

    try {

        MoviOptions movi_options;

        if (!parse_command(argc, argv, movi_options)) {
            return 1;
        }

        // If validate flags is set, just return 0, no execution is needed.
        if (movi_options.is_validate_flags()) {
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

            // Set custom buffer for stdout (printf uses stdout)
            if (setvbuf(stdout, buffer, _IOFBF, BUFFER_SIZE) != 0) {
                perror("Failed to set buffer for stdout");
                return 1;
            }
        }

        std::string command = movi_options.get_command();

        if (command == "build") {
            if (movi_options.use_separators()) {
                if (!SUPPORTS_SEPARATORS) {
                    // TODO: Fully support separators for large, split, and constant indexes
                    throw std::runtime_error(ERROR_MSG("[build] Separators are not supported for the " + program() + " index."));
                }
            }

            MoveStructure mv_(&movi_options, SPLIT_ARRAY, CONSTANT_INDEX);

            mv_.build();

            if (movi_options.is_verify()) {
                INFO_MSG("Verifying the LF_move results...");
                mv_.verify_lf_loop();
            }
            mv_.serialize();
            build_ftab(mv_, movi_options);
            SUCCESS_MSG("The Movi index is successfully stored at " + movi_options.get_index_dir());
            if (movi_options.is_output_ids()) {
                mv_.output_ids();
            }

            INFO_MSG("Generating the null statistics...");
            Classifier classifier;

            // generate pml null database
            movi_options.set_pml();
            movi_options.set_generate_null_reads(true);
            classifier.generate_null_statistics(mv_, movi_options);
            INFO_MSG("Successfully generated null statistics with PML");

            // generate zml null database
            movi_options.set_zml();
            movi_options.set_generate_null_reads(false); // do not regenerate the null reads
            classifier.generate_null_statistics(mv_, movi_options);
            INFO_MSG("Successfully generated null statistics with ZML");

            if (movi_options.is_color()) {
                color(mv_, movi_options);
            }

        } else if (command == "build-SA") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();
            mv_.find_sampled_SA_entries();
            mv_.serialize_sampled_SA();
            SUCCESS_MSG("Successfully stored sampled SA entries at " + movi_options.get_index_dir());
        } else if (command == "color") {
            MoveStructure mv_(&movi_options);
            auto begin = std::chrono::system_clock::now();
            mv_.deserialize();
            auto end = std::chrono::system_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            TIMING_MSG(elapsed, "loading the index");

            color(mv_, movi_options);
        } else if (command == "query") {
            // Check if the input file exists
            if (movi_options.get_read_file() != "-" and !std::filesystem::exists(movi_options.get_read_file())) {
                throw std::runtime_error(ERROR_MSG("The input file " + movi_options.get_read_file() + " does not exist."));
            }

            MoveStructure mv_(&movi_options);

            auto begin = std::chrono::system_clock::now();
            mv_.deserialize();
            auto end = std::chrono::system_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            TIMING_MSG(elapsed, "loading the index");

            if (movi_options.is_get_sa_entries()) {
                begin = std::chrono::system_clock::now();
                mv_.deserialize_sampled_SA();
                end = std::chrono::system_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
                TIMING_MSG(elapsed, "loading the sampled SA entries");
            }

            if (movi_options.is_multi_classify()) {
                load_color_table(mv_, movi_options);

                if (movi_options.is_report_colors() or movi_options.is_report_color_ids()) {
                    // Copmuter color ids for the output
                    mv_.compute_color_ids_from_flat();
                }
            }

            query(mv_, movi_options);

            // Avoid taking too long for dealloction of large data structures at the end of the program
            std::quick_exit(0);

        } else if (command == "view") {
            view(movi_options);
        } else if (command == "rlbwt") {
            build_rlbwt(movi_options);
        } else if (command == "color-move-rows") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();

            load_color_table(mv_, movi_options);

            mv_.add_colors_to_rlbwt();

            mv_.serialize();

        } else if (command == "LF") {
            MoveStructure mv_(&movi_options);
            mv_.deserialize();
            INFO_MSG("The Movi index is read from the file successfully.");
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
                INFO_MSG("The color table is read successfully.");
                mv_.flat_and_serialize_colors_vectors();
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
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
