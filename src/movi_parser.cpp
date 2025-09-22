#include "movi_parser.hpp"

bool parse_command(int argc, char** argv, MoviOptions& movi_options, bool supress_messages) {
    // movi_options.print_options();

    // cxxopts::Options options("movi-" + program(), "Please use the following format:");
    cxxopts::Options options("movi", "");

    options.add_options()
        ("command", "Command to execute", cxxopts::value<std::string>())
        ("type", "Which index type should be built or used.", cxxopts::value<std::string>())
        ("h,help", "Print help")
        ("no-header", "Header information in not stored")
        ("default-blocks", "The block size should not be read from the index (for blocked indexes)")
        ("d,dbg", "Enable debug mode")
        ("v,verbose", "Enable verbose mode")
        ("l,logs", "Enable logs");

    auto buildOptions = options.add_options("build")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("f,fasta", "Reference file", cxxopts::value<std::string>())
        ("l,list", "List of fasta files, only works with 'movi' binary", cxxopts::value<std::string>())
        ("color", "Add colors to the index for multi-class classification")
        ("color-vectors", "Build a vector of vectors for colors (builds \"ref.fa.doc_sets.bin\")")
        ("thresholds", "Store the threshold values in the compact mode by splitting the runs at threshold boundaries")
        ("preprocessed", "The BWT is preprocessed into heads and lens files")
        ("verify", "Verify if all the LF_move operations are correct")
        ("output-ids", "Output the adjusted ids of all the runs to ids.* files, one file per character")
        ("ftab-k", "The length of the ftab kmer", cxxopts::value<uint32_t>())
        ("tally", "Sample id at every tally runs", cxxopts::value<uint32_t>())
        ("multi-ftab", "Use ftabs with smaller k values if the largest one fails")
        ("keep", "Keep the extra files after the build step")
        ("skip-prepare", "Skip the prepare_ref step")
        ("skip-pfp", "Skip the pfp_thresholds step")
        ("skip-rlbwt", "Skip the rlbwt step -- constant and split index")
        ("skip-r-permute", "Skip the r-permute step -- constant and split index");

    auto buildSAOptions = options.add_options("build-SA")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("sample-rate", "The sample rate for storing the sampled SA (default: 100)", cxxopts::value<uint64_t>());

    auto queryOptions = options.add_options("query")
        ("pml", "Compute the pseudo-matching lengths (PMLs)")
        ("rpml", "Compute the pseudo-matching lengths using random repositioning (RPMLs)")
        ("zml", "Compute the Ziv-Merhav cross parsing length (ZMLs)")
        ("count", "Compute the count queries")
        ("kmer", "Search all the kmers")
        ("kmer-count", "Find the count of every kmer")
        ("sa-entries", "Find the SA entries for each read")
        ("classify", "Classify the reads")
        ("filter", "Filter the reads based on the matching lengths, output the filtered reads to stdout")
        ("multi-classify", "Multi-class classification with PMLs")
        ("early-stop", "Early stop the read processing for unclassified reads")
        ("report-all", "Report all the taxon ids for each read (default: min-diff-frac = 0.05), not available in no-prefetch mode")
        ("report-colors", "Report the colors (as offsets in the flat color table) for each PML")
        ("report-color-ids", "Report the color ids for each PML")
        ("min-diff-frac", "Report all the taxon ids which have score close to the best document by this fraction", cxxopts::value<float>())
        ("min-score-frac", "Report all the taxon ids which have score >= min-score-frac * read-length. If set, min-diff-frac mode is turned off.", cxxopts::value<float>())
        ("min-len", "Minimum matching length for classification (only consider PMLs >= min-len), default is 1", cxxopts::value<uint8_t>())
        ("pvalue-scoring", "Use p-value scoring for classification")
        ("full", "Use full coloring information to compute pseudo-matching lengths (PMLs)")
        ("compress", "Use compressed document sets for classification")
        ("freq-compress", "Use frequency compressed document sets for classification")
        ("tree-compress", "Use tree compressed document sets for classification")
        ("color-move-rows", "Color the move rows, query is not performed")
        ("color-vectors", "Use vector of vectors for colors (requires \"ref.fa.doc_sets.bin\")")
        ("bin-width", "The width of the bin used for classification", cxxopts::value<uint32_t>())
        ("reverse", "Use the reverse (not reverse complement) of the reads to perform queries")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("r,read", "fasta/fastq Read file for query", cxxopts::value<std::string>())
        ("mmap", "Use memory mapping to read the index")
        ("n,no-prefetch", "Disable prefetching for query")
        ("k,k-length", "The length of the kmer", cxxopts::value<uint32_t>())
        ("ftab-k", "The length of the ftba kmer", cxxopts::value<uint32_t>())
        ("multi-ftab", "Use ftabs with smaller k values if the largest one fails")
        ("s,strands", "Number of strands for query", cxxopts::value<int>())
        ("t,threads", "Number of threads for query", cxxopts::value<int>())
        ("o,out-file", "Output file if computing PMLs for classification", cxxopts::value<std::string>())
        ("stdout", "Write the output to stdout, writes the matching lengths by default, or the report of classification if --classify is passed")
        ("no-output", "Do not write any output, ignores other options about the output")
        ("ignore-illegal-chars", "In the case of illegal characters (i.e., non-ACGT for genomic data), substitute the character with \'A\'(1) or a random character from the alphabet (2).", cxxopts::value<int>());

    auto colorOptions = options.add_options("color")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("color-vectors", "Build a vector of vectors for colors (builds \"ref.fa.doc_sets.bin\")")
        ("compress", "Whether or not we compress doc sets (only keep most frequent few)")
        ("full", "Whether or not to store all document information (or just the sets for each run)")
        ("t,threads", "Number of threads for query", cxxopts::value<int>());

    auto viewOptions = options.add_options("view")
        ("mls-file", "The matching lengths (PML or ZML) file in the binary format", cxxopts::value<std::string>())
        ("small-pml", "Read the binary file with PMLs stored as uint16_t.")
        ("large-pml", "Read the binary file with PMLs stored as uint64_t.");

    auto rlbwtOptions = options.add_options("rlbwt")
        ("bwt-file", "BWT file", cxxopts::value<std::string>());

    auto LFOptions = options.add_options("LF")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("lf-type", "type of the LF query: \"reconstruct\", \"sequential\", or \"random\"", cxxopts::value<std::string>());

    auto inspectOptions = options.add_options("inspect")
        ("output-ids", "Output the adjusted ids of all the runs to ids.* files, one file per character")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("flat-color-vectors", "Flat and serialize the colors vectors");

    auto ftabOptions = options.add_options("ftab")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("ftab-k", "The length of the ftab kmer", cxxopts::value<uint32_t>())
        ("multi-ftab", "Use ftabs with smaller k values if the largest one fails");

    auto nullOptions = options.add_options("null")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("gen-reads", "Generate null reads");

    options.parse_positional({ "command" });
    std::vector<std::string> help_groups;

    try {
        auto result = options.parse(argc, argv);

        // Set global logs flags first

        if (result.count("no-header")) {
            // For hadnling old indexes
            movi_options.set_no_header(true);
        }

        if (result.count("default-block")) {
            // For hadnling old indexes
            movi_options.set_adjusted_block(false);
        }

        if (result.count("verbose")) {
            movi_options.set_verbose(true);
        }

        if (result.count("logs")) {
            movi_options.set_logs(true);
        }

        if (result.count("dbg")) {
            movi_options.set_debug(true);
        }

        // Set command-specific flags
        if (result.count("command")) {
            std::string command = result["command"].as<std::string>();
            movi_options.set_command(command);
            help_groups.push_back("");
            help_groups.push_back(command);

            if (command == "build") {
                if (result.count("index") == 1 and result.count("fasta") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    movi_options.set_ref_file(result["fasta"].as<std::string>());
                    if (result.count("color")) { movi_options.set_color(true); }
                    if (result.count("color-vectors") == 1) { movi_options.set_doc_sets_vector_of_vectors(true); }
                    if (result.count("ftab-k") >= 1) { movi_options.set_ftab_k(static_cast<uint32_t>(result["ftab-k"].as<uint32_t>())); }
                    if (result.count("multi-ftab") >= 1) { movi_options.set_multi_ftab(true); }
                    if (result.count("output-ids") >= 1) { movi_options.set_output_ids(true); }
                    if (result.count("verify")) {
                        movi_options.set_verify(true);
                    }
                    if (result.count("preprocessed")) {
                        movi_options.set_preprocessed(true);
                    }
                    if (result.count("thresholds")) {
                        movi_options.set_thresholds(true);
                    }
#if TALLY_MODE
                    if (result.count("tally") >= 1) {
                        movi_options.set_tally_checkpoints(static_cast<uint32_t>(result["tally"].as<uint32_t>()));
                    }
#endif
#if USE_THRESHOLDS
                    // In these modes, thresholds are always stored
                    movi_options.set_thresholds(true);
#endif
                } else {
                    const std::string message = "Please include one index directory and one fasta file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "build-SA") {
                if (result.count("index") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    if (result.count("sample-rate") >= 1) {
                        movi_options.set_SA_sample_rate(static_cast<uint64_t>(result["sample-rate"].as<uint64_t>()));
                    }
                } else {
                    const std::string message = "Please specify the index directory file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "color") { 
                if (result.count("index") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    if (result.count("color-vectors") == 1) {
                        movi_options.set_doc_sets_vector_of_vectors(true);
                    }
                    if (result.count("full")) {
                        movi_options.set_full_color(true);
                    }
                    if (result.count("compress")) {
                        movi_options.set_compressed(true);
                    }
                    if (result.count("threads") == 1) {
                        if (!supress_messages) {
                            std::cerr << "threads: " << result["threads"].as<int>() << "\n";
                        }
                        movi_options.set_threads(static_cast<size_t>(result["threads"].as<int>()));
                    }
                } else {
                    const std::string message = "Please include one index directory and one fasta file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "query") {
                try {
                    if (result.count("index") == 1 and result.count("read") == 1) {
                        movi_options.set_index_dir(result["index"].as<std::string>());
                        movi_options.set_read_file(result["read"].as<std::string>());
                        if (result.count("k") >= 1) { movi_options.set_k(static_cast<uint32_t>(result["k"].as<uint32_t>())); }
                        if (result.count("ftab-k") >= 1) { movi_options.set_ftab_k(static_cast<uint32_t>(result["ftab-k"].as<uint32_t>())); }
                        if (result.count("bin-width") >= 1) { movi_options.set_bin_width(static_cast<uint32_t>(result["bin-width"].as<uint32_t>())); }
                        if (result.count("multi-ftab") >= 1) { movi_options.set_multi_ftab(true); }
                        if (result.count("kmer") >= 1) { movi_options.set_kmer(); }
                        if (result.count("kmer-count") >= 1) { movi_options.set_kmer(); movi_options.set_kmer_count(true); }
                        if (result.count("count") >= 1) { movi_options.set_count(); }
                        if (result.count("zml") >= 1) { movi_options.set_zml(); }
                        if (result.count("pml") >= 1) { movi_options.set_pml(); }
                        if (result.count("sa-entries") >= 1) { movi_options.set_get_sa_entries(true); }
                        if (result.count("rpml") >= 1) { movi_options.set_random_repositioning(true); }
                        if (result.count("classify") >= 1) { movi_options.set_classify(true); }
                        if (result.count("filter") >= 1) { movi_options.set_filter(true); }
                        if (result.count("multi-classify") >= 1) { movi_options.set_multi_classify(true); }
                        if (result.count("early-stop") >= 1) { movi_options.set_early_stop(true); }
                        if (result.count("report-all") >= 1) { movi_options.set_report_all(true); }
                        if (result.count("report-colors") >= 1) { movi_options.set_report_colors(true); }
                        if (result.count("report-color-ids") >= 1) { movi_options.set_report_color_ids(true); }
                        if (result.count("min-diff-frac") >= 1) { movi_options.set_min_diff_frac(static_cast<float>(result["min-diff-frac"].as<float>())); }
                        if (result.count("min-score-frac") >= 1) { movi_options.set_min_score_frac(static_cast<float>(result["min-score-frac"].as<float>())); }
                        if (result.count("min-len")) { movi_options.set_min_match_len(result["min-len"].as<uint8_t>()); }
                        if (result.count("pvalue-scoring") == 1) { movi_options.set_pvalue_scoring(true); }
                        if (result.count("color-move-rows") == 1) { movi_options.set_color_move_rows(true); }
                        if (result.count("color-vectors") == 1) { movi_options.set_doc_sets_vector_of_vectors(true); }
                        if (result.count("reverse") == 1) { movi_options.set_reverse(true); }
                        if (result.count("out-file")) { movi_options.set_out_file(result["out-file"].as<std::string>()); }
                        if (result.count("multi-classify")) {
                            if (result.count("full")) {
                                movi_options.set_full_color(true);
                            }
                            if (result.count("freq-compress")) {
                                movi_options.set_freq_compressed(true);
                                movi_options.set_tree_compressed(false);
                            }
                            if (result.count("tree-compress")) {
                                movi_options.set_tree_compressed(true);
                                movi_options.set_freq_compressed(false);
                            }

                            if (result.count("out-file") == 1) {
                                movi_options.set_out_file(result["out-file"].as<std::string>());
                            } else {
                                const std::string message = "Please include the out file for multi-classify.";
                                cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                            }
                        }
                        if (result.count("ignore-illegal-chars") == 1) {
                            if (!movi_options.set_ignore_illegal_chars(result["ignore-illegal-chars"].as<int>())) {
                                const std::string message = "ignore-illegal-chars should be either 1 (set illegal chars to \'A\') or 2 (set illegal chars to a random char).";
                                cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                            }
                        }
                        if (movi_options.is_pml() and movi_options.is_count()) {
                            const std::string message = "Please only specify count or pml as the type of queries.";
                            cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                        }
                        if (result.count("mmap") == 1) {
                            movi_options.set_mmap(true);
                        }
                        if (result.count("no-prefetch") == 1) {
                            movi_options.set_prefetch(false);
                        }
                        if (result.count("strands") == 1) {
                            if (!supress_messages) {
                                std::cerr << "Number of strands: " << result["strands"].as<int>() << "\n";
                            }
                            movi_options.set_strands(static_cast<size_t>(result["strands"].as<int>()));
                        }
                        if (result.count("threads") == 1) {
                            if (!supress_messages) {
                                std::cerr << "Number of threads: " << result["threads"].as<int>() << "\n";
                            }
                            movi_options.set_threads(static_cast<size_t>(result["threads"].as<int>()));
                        }
                        if (result.count("stdout")) {
                            // Set global verbose flag
                            movi_options.set_stdout(true);
                        }
                        if (result.count("no-output")) {
                            movi_options.set_no_output(true);
                        }
                    } else {
                        const std::string message = "Please include one index directory and one read file.";
                        cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                    }
                } catch (const cxxopts::exceptions::exception& e) {
                    if (!supress_messages) {
                        std::cerr << "Error parsing command line options: " << e.what() << "\n";
                        std::cerr << options.help(help_groups) << "\n";
                    }
                    return false;
                }
            } else if (command == "rlbwt") {
                if (result.count("bwt-file") == 1) {
                    movi_options.set_bwt_file(result["bwt-file"].as<std::string>());
                } else {
                    const std::string message = "Please specify one bwt file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "view") {
                if (result.count("mls-file") == 1) {
                    movi_options.set_mls_file(result["mls-file"].as<std::string>());
                    if (result.count("small-pml") >= 1) {
                        movi_options.set_small_pml_lens(true);
                    }

                    if (result.count("large-pml") >= 1) {
                        movi_options.set_large_pml_lens(true);
                    }

                    if (result.count("classify") >= 1) {
                        movi_options.set_classify(true);
                        if (result.count("bin-width") >= 1) {
                            movi_options.set_bin_width(static_cast<uint32_t>(result["bin-width"].as<uint32_t>()));
                        }
                        if (result.count("zml") >= 1) {
                            movi_options.set_zml();
                        }
                        if (result.count("pml") >= 1) {
                            movi_options.set_pml();
                        }
                        if (result.count("index") == 1) {
                            movi_options.set_index_dir(result["index"].as<std::string>());
                        } else {
                            const std::string message = "Please specify the index directory file.";
                            cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                        }
                    }
                } else {
                    const std::string message = "Please specify one mls.bin file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "LF") {
                if (result.count("index") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    if (result.count("lf-type")) {
                        if (!movi_options.set_LF_type(result["lf-type"].as<std::string>())) {
                            const std::string message = "The LF type is not defined, please choose from: \"reconstruct\", \"sequential\", or \"random\"";
                            cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                        }
                    }
                } else {
                    const std::string message = "Please specify the index directory file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "inspect") {
                if (result.count("index")) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    if (result.count("output-ids") >= 1) { movi_options.set_output_ids(true); }
                    if (result.count("flat-color-vectors") == 1) { movi_options.set_flat_color_vectors(true); }
                } else {
                    const std::string message = "Please specify the index directory file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "ftab") {
                if (result.count("index") and result.count("ftab-k")) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    movi_options.set_ftab_k(static_cast<uint32_t>(result["ftab-k"].as<uint32_t>()));
                    if (result.count("multi-ftab") >= 1) { movi_options.set_multi_ftab(true); }
                } else {
                    const std::string message = "Please specify the index directory file and the k length for the ftab.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "null") {
                if (result.count("index") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    if (result.count("gen-reads") >= 1) {
                        movi_options.set_generate_null_reads(true);
                        if (result.count("fasta") >= 1) {
                            movi_options.set_ref_file(result["fasta"].as<std::string>());
                        } else {
                            const std::string message = "Please specify the reference fasta file.";
                            cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                        }
                    }
                    if (result.count("zml") >= 1) { movi_options.set_zml(); }
                    if (result.count("pml") >= 1) { movi_options.set_pml(); }
                } else {
                    const std::string message = "Please specify the index directory file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else {
                const std::string message = "Invalid action: \"" + command + "\"";
                cxxopts::throw_or_mimic<cxxopts::exceptions::no_such_option>(message);
            }

            if (result.count("help")) {
                if (!supress_messages) {
                    std::cerr << options.help(help_groups) << std::endl;
                }
                return 0;
            }

        } else {
            const std::string message = "No action specified.";
            help_groups.push_back(" ");
            cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
        }
    } catch (const cxxopts::exceptions::exception& e) {
        if (!supress_messages) {
            std::cerr << "Error parsing command line options: " << e.what() << "\n";
            std::cerr << options.help(help_groups) << "\n";
        }
        return false;
    }
    return true;
}