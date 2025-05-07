#include "movi_parser.hpp"

bool parse_command(int argc, char** argv, MoviOptions& movi_options) {
    // movi_options.print_options();

    // cxxopts::Options options("movi-" + program(), "Please use the following format:");
    cxxopts::Options options("movi", "");

    options.add_options()
        ("command", "Command to execute", cxxopts::value<std::string>())
        ("type", "Which index type should be built or used.", cxxopts::value<std::string>())
        ("h,help", "Print help")
        ("no-header", "Header information in not stored")
        ("d,dbg", "Enable debug mode")
        ("v,verbose", "Enable verbose mode")
        ("l,logs", "Enable logs");

    auto buildOptions = options.add_options("build")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("f,fasta", "Reference file", cxxopts::value<std::string>())
        ("l,list", "List of fasta files, only works with 'movi' binary", cxxopts::value<std::string>())
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

    auto queryOptions = options.add_options("query")
        ("pml", "Compute the pseudo-matching lengths (PMLs)")
        ("rpml", "Compute the pseudo-matching lengths using random repositioning (RPMLs)")
        ("zml", "Compute the Ziv-Merhav cross parsing length (ZMLs)")
        ("count", "Compute the count queries")
        ("kmer", "Search all the kmers")
        ("kmer-count", "Find the count of every kmer")
        ("classify", "Classify the reads")
        ("multi-classify", "Multi-class classification with PMLs")
        ("thres", "Threshold for classification (only consider PMLs above thres)", cxxopts::value<uint8_t>())
        ("scale", "Scale for p-value classification", cxxopts::value<double>())
        ("full", "Use full coloring information to compute pseudo-matching lengths (PMLs)")
        ("compress", "Use compressed document sets for classification")
        ("color-move-rows", "Color the move rows, query is not performed")
        ("bin-width", "The width of the bin used for classification", cxxopts::value<uint32_t>())
        ("reverse", "Use the reverse (not reverse complement) of the reads to perform queries")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("r,read", "fasta/fastq Read file for query", cxxopts::value<std::string>())
        ("n,no-prefetch", "Disable prefetching for query")
        ("k,k-length", "The length of the kmer", cxxopts::value<uint32_t>())
        ("ftab-k", "The length of the ftba kmer", cxxopts::value<uint32_t>())
        ("multi-ftab", "Use ftabs with smaller k values if the largest one fails")
        ("s,strands", "Number of strands for query", cxxopts::value<int>())
        ("t,threads", "Number of threads for query", cxxopts::value<int>())
        ("out_file", "Output file if computing PMLs for classification", cxxopts::value<std::string>())
        ("stdout", "Write the output to stdout")
        ("ignore-illegal-chars", "In the case of illegal characters (i.e., non-ACGT for genomic data), substitute the character with \'A\'(1) or a random character from the alphabet (2).", cxxopts::value<int>());

    auto colorOptions = options.add_options("color")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("compress", "Whether or not we compress doc sets (only keep most frequent few)")
        ("full", "Whether or not to store all document information (or just the sets for each run)")
        ("t,threads", "Number of threads for query", cxxopts::value<int>());

    auto viewOptions = options.add_options("view")
        ("mls-file", "The matching lengths (PML or ZML) file in the binary format", cxxopts::value<std::string>());

    auto rlbwtOptions = options.add_options("rlbwt")
        ("bwt-file", "BWT file", cxxopts::value<std::string>());

    auto LFOptions = options.add_options("LF")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("lf-type", "type of the LF query: \"reconstruct\", \"sequential\", or \"random\"", cxxopts::value<std::string>());

    auto inspectOptions = options.add_options("inspect")
        ("output-ids", "Output the adjusted ids of all the runs to ids.* files, one file per character")
        ("i,index", "Index directory", cxxopts::value<std::string>());

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

        if (result.count("no-header")) {
            // Set global verbose flag
            movi_options.set_no_header(true);
        }

        if (result.count("verbose")) {
            // Set global verbose flag
            movi_options.set_verbose(true);
        }

        if (result.count("logs")) {
            // Set global logs flag
            movi_options.set_logs(true);
        }

        if (result.count("dbg")) {
            // Set global debug flag
            movi_options.set_debug(true);
        }

        if (result.count("command")) {
            std::string command = result["command"].as<std::string>();
            movi_options.set_command(command);
            help_groups.push_back("");
            help_groups.push_back(command);

            if (command == "build") {
                if (result.count("index") == 1 and result.count("fasta") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    movi_options.set_ref_file(result["fasta"].as<std::string>());
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
            } else if (command == "color") { 
                if (result.count("index") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    if (result.count("full")) {
                        movi_options.set_full_color(true);
                    }
                    if (result.count("compress")) {
                        movi_options.set_compress(true);
                    }
                    if (result.count("threads") == 1) {
                        std::cerr << "threads: " << result["threads"].as<int>() << "\n";
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
                        if (result.count("rpml") >= 1) { movi_options.set_random_repositioning(true); }
                        if (result.count("classify") >= 1) { movi_options.set_classify(true); }
                        if (result.count("multi-classify") >= 1) { movi_options.set_multi_classify(true); }
                        if (result.count("thres")) { movi_options.set_thres(result["thres"].as<uint8_t>()); }
                        if (result.count("scale")) { movi_options.set_scale(result["scale"].as<double>()); }
                        if (result.count("color-move-rows") == 1) { movi_options.set_color_move_rows(true); }
                        if (result.count("reverse") == 1) { movi_options.set_reverse(true); }
                        if (result.count("pml") || result.count("zml")) {
                            if (result.count("full")) {
                                movi_options.set_full_color(true);
                            }
                            if (result.count("compress")) {
                                movi_options.set_compress(true);
                            }
                            if (result.count("out_file")) {
                                movi_options.set_out_file(result["out_file"].as<std::string>());
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
                        if (result.count("no-prefetch") == 1) {
                            movi_options.set_prefetch(false);
                        }
                        if (result.count("strands") == 1) {
                            std::cerr << "strands: " << result["strands"].as<int>() << "\n";
                            movi_options.set_strands(static_cast<size_t>(result["strands"].as<int>()));
                        }
                        if (result.count("threads") == 1) {
                            std::cerr << "threads: " << result["threads"].as<int>() << "\n";
                            movi_options.set_threads(static_cast<size_t>(result["threads"].as<int>()));
                        }
                        if (result.count("stdout")) {
                            // Set global verbose flag
                            movi_options.set_stdout(true);
                        }
                    } else {
                        const std::string message = "Please include one index directory and one read file.";
                        cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                    }
                } catch (const cxxopts::exceptions::exception& e) {
                    std::cerr << "Error parsing command line options: " << e.what() << "\n";
                    std::cerr << options.help(help_groups) << "\n";
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
                    std::cerr << options.help(help_groups) << std::endl;
                    return 0;
                }

        } else {
            const std::string message = "No action specified.";
            cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
        }
    } catch (const cxxopts::exceptions::exception& e) {
        std::cerr << "Error parsing command line options: " << e.what() << "\n";
        std::cerr << options.help(help_groups) << "\n";
        return false;
    }
    return true;
}