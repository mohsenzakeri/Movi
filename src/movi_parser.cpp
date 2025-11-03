#include "movi_parser.hpp"

void add_all_groups(std::vector<std::string>& all_groups, std::string command) {
    if (command == "query" or command == "") {
        all_groups.push_back("query (advanced)");
        all_groups.push_back("query (color)");
    }

    if (command == "build" or command == "") {
        all_groups.push_back("build (advanced)");
        all_groups.push_back("build (color)");
    }

    all_groups.push_back("developer");
}

bool handle_help(int argc, char** argv, cxxopts::Options& options, std::vector<std::string>& help_groups,
                 std::vector<std::string>& main_actions, std::vector<std::string>& all_actions) {

    bool help = false;
    bool help_all = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            help = true;
            break;
        } else if (arg == "--help-all") {
            help_all = true;
            break;
        }
    }

    std::string command = get_action(argc, argv);

    // If the command is valid but not in the main actions, add it to the main actions for the help message
    if (command != "" and
        std::find(main_actions.begin(), main_actions.end(), command) == main_actions.end() and
        std::find(all_actions.begin(), all_actions.end(), command) != all_actions.end()) {
        // Add the command to the main actions
        main_actions.push_back(command);
    }

    // Check if the command is a valid action
    if (command != "" and std::find(all_actions.begin(), all_actions.end(), command) != all_actions.end()) {
        // Print the hlep message on for this action
        help_groups.push_back("");
        help_groups.push_back(command);
        if (help_all) {
            add_all_groups(help_groups, command);
        }
        options.set_actions(main_actions);
    } else {
        // When the action is missing or is not valid
        if (help_all) {
            help_groups = all_actions;
            add_all_groups(help_groups, "");
            options.set_actions(all_actions);
        } else {
            help_groups = main_actions;
            options.set_actions(main_actions);
        }
    }

    if (help or help_all) {
        return true;
    }

    return false;
}

bool parse_command(int argc, char** argv, MoviOptions& movi_options, bool supress_messages) {

    // cxxopts::Options options("movi-" + program(), "Please use the following format:");
    cxxopts::Options options("movi", "Movi: A tool for indexing and querying pangenomes\n");

    std::vector<std::string> all_actions;
    std::vector<std::string> main_actions;

    all_actions.push_back("");
    main_actions.push_back("");
    options.add_options()
        ("command", "Command to execute", cxxopts::value<std::string>())
        ("type", "The type of the index: regular-thresholds (default), regular, blocked-thresholds, blocked, sampled-thresholds, sampled", cxxopts::value<std::string>())
        ("h,help", "Print help")
        ("help-all", "Print help with all options (including advanced)")
        ("verbose", "Enable verbose mode");

    options.add_options("developer")
        ("validate-flags", "Validate the command line flags are valid")
        ("no-header", "Header information is not stored")
        ("legacy-header", "Use legacy single-byte header format")
        ("default-blocks", "The block size should not be read from the index (for blocked indexes)")
        ("d,debug", "Enable debug messages")
        ("logs", "Enable logging");

    all_actions.push_back("build");
    main_actions.push_back("build");
    auto buildOptions = options.add_options("build")
        ("i,index", "Index directory [REQUIRED]", cxxopts::value<std::string>())
        ("f,fasta", "Reference file [REQUIRED unless -l is passed]", cxxopts::value<std::string>())
        ("l,list", "List of fasta files, only works with 'movi' binary [REQUIRED unless -f is passed]", cxxopts::value<std::string>())
        ("separators", "Use separators in between fasta entries")
        ("checkpoint", "Create checkpoint for id field of sampled and sampled-thresholds indexes every n move rows", cxxopts::value<uint32_t>())
        ("verify", "Verify if all the LF-move operations are correct");

    auto buildAdvancedOptions = options.add_options("build (advanced)")
        ("thresholds", "Store the threshold values by splitting the runs at threshold boundaries (default)")
        ("output-ids", "Output the adjusted ids of all the runs to ids.* files, one file per character")
        ("max-run-length", "The maximum length of the run", cxxopts::value<uint64_t>())
        ("ftab-k", "The length of the ftab kmer", cxxopts::value<uint32_t>())
        ("multi-ftab", "Use ftabs with smaller k values if the largest one fails")
        ("preprocessed", "The BWT is preprocessed into heads and lens files")
        ("keep", "Keep the extra files after the build step")
        ("skip-prepare", "Skip the prepare_ref step")
        ("skip-pfp", "Skip the pfp_thresholds step")
        ("skip-rlbwt", "Skip the rlbwt step -- constant and split index")
        ("skip-r-permute", "Skip the r-permute step -- constant and split index");

    auto buildColorOptions = options.add_options("build (color)")
        ("color", "Add colors to the index for multi-class classification")
        ("color-vectors", "Build a vector of vectors for colors (builds \"ref.fa.doc_sets.bin\")");

    all_actions.push_back("inspect");
    main_actions.push_back("inspect");
    auto inspectOptions = options.add_options("inspect")
        ("i,index", "Index directory [REQUIRED]", cxxopts::value<std::string>())
        ("output-ids", "Output the adjusted ids of all the runs to ids.* files, one file per character")
        ("flat-color-vectors", "Flat and serialize the colors vectors");

    all_actions.push_back("query");
    main_actions.push_back("query");
    auto queryOptions = options.add_options("query")
        ("i,index", "Index directory [REQUIRED]", cxxopts::value<std::string>())
        ("r,read", "fasta/fastq Read file for query [REQUIRED]", cxxopts::value<std::string>())
        ("o,out-file", "Output file prefix if computing PMLs for classification", cxxopts::value<std::string>())
        ("t,threads", "Number of threads for query", cxxopts::value<int>())
        ("pml", "Compute the pseudo-matching lengths (default)")
        ("zml", "Compute the Ziv-Merhav cross parsing length)")
        ("count", "Compute the count queries")
        ("classify", "Enable binary classification of the reads")
        ("filter", "Filter the reads based on the matching lengths, output the filtered reads to stdout")
        ("v,invert", "Output the not found reads during filtering")
        ("stdout", "Write the output to stdout, writes the matching lengths by default, or the report of classification if --classify is passed");

    // Advanced query options (hidden by default, shown with --help-all)
    auto queryAdvancedOptions = options.add_options("query (advanced)")
        ("mem", "Compute the maximal exact matches (MEMs)")
        ("l,min-mem-length", "The minimum length of the MEMs", cxxopts::value<uint32_t>())
        ("rpml", "Compute the pseudo-matching lengths using random repositioning (RPMLs)")
        ("kmer", "Search all the kmers")
        ("kmer-count", "Find the count of every kmer")
        ("k,k-length", "The length of the kmer", cxxopts::value<uint32_t>())
        ("bin-width", "The width of the bin used for binary classification", cxxopts::value<uint32_t>())
        ("ignore-illegal-chars", "In the case of illegal characters (i.e., non-ACGT for genomic data), substitute the character with \'A\'(1) or a random character from the alphabet (2).", cxxopts::value<int>())
        ("sa-entries", "Find the SA entries for each read (requires build-SA step before executing this command)")
        ("reverse", "Use the reverse (not reverse complement) of the reads to perform queries")
        ("s,strands", "Number of strands for query", cxxopts::value<int>())
        ("mmap", "Use memory mapping to read the index")
        ("n,no-prefetch", "Disable prefetching for query")
        ("ftab-k", "The length of the ftba kmer", cxxopts::value<uint32_t>())
        ("multi-ftab", "Use ftabs with smaller k values if the largest one fails")
        ("no-output", "Do not write any output, ignores other options about the output");

    auto moviColorOptions = options.add_options("query (color)")
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
        ("color-vectors", "Use vector of vectors for colors (requires \"ref.fa.doc_sets.bin\")");

    all_actions.push_back("view");
    main_actions.push_back("view");
    auto viewOptions = options.add_options("view")
        ("bpf", "The base profile format (BPF) file to view [REQUIRED]", cxxopts::value<std::string>())
        ("small-bpf", "Read the file with PMLs stored as uint16_t (default: uint32_t).")
        ("large-bpf", "Read the file with PMLs stored as uint64_t (default: uint32_t).");

    all_actions.push_back("build-SA");
    auto buildSAOptions = options.add_options("build-SA")
        ("i,index", "Index directory [REQUIRED]", cxxopts::value<std::string>())
        ("sample-rate", "The sample rate for storing the sampled SA (default: 100)", cxxopts::value<uint64_t>());

    all_actions.push_back("color");
    auto colorOptions = options.add_options("color")
        ("i,index", "Index directory [REQUIRED]", cxxopts::value<std::string>())
        ("color-vectors", "Build a vector of vectors for colors (builds \"ref.fa.doc_sets.bin\")")
        ("compress", "Whether or not we compress doc sets (only keep most frequent few)")
        ("full", "Whether or not to store all document information (or just the sets for each run)")
        ("t,threads", "Number of threads for query", cxxopts::value<int>());

    all_actions.push_back("rlbwt");
    auto rlbwtOptions = options.add_options("rlbwt")
        ("bwt-file", "BWT file", cxxopts::value<std::string>());

    all_actions.push_back("color-move-rows");
    auto colorMoveRowsOptions = options.add_options("color-move-rows")
        ("i,index", "Index directory [REQUIRED]", cxxopts::value<std::string>());

    all_actions.push_back("LF");
    auto LFOptions = options.add_options("LF")
        ("i,index", "Index directory [REQUIRED]", cxxopts::value<std::string>())
        ("lf-type", "type of the LF query: \"reconstruct\", \"sequential\", or \"random\"", cxxopts::value<std::string>());

    all_actions.push_back("ftab");
    auto ftabOptions = options.add_options("ftab")
        ("i,index", "Index directory [REQUIRED]", cxxopts::value<std::string>())
        ("ftab-k", "The length of the ftab kmer", cxxopts::value<uint32_t>())
        ("multi-ftab", "Use ftabs with smaller k values if the largest one fails");

    all_actions.push_back("null");
    auto nullOptions = options.add_options("null")
        ("i,index", "Index directory [REQUIRED]", cxxopts::value<std::string>())
        ("gen-reads", "Generate null reads");

    options.parse_positional({ "command" });
    std::vector<std::string> help_groups;

    std::string command = get_action(argc, argv);

    if (handle_help(argc, argv, options, help_groups, main_actions, all_actions)) {
        INFO_MSG(options.help(help_groups));
        return false;
    }

    try {
        cxxopts::ParseResult result = options.parse(argc, argv);

        if (result.count("validate-flags")) {
            movi_options.set_validate_flags(true);
        }

        if (result.count("no-header")) {
            // For hadnling old indexes with no header
            movi_options.set_no_header(true);
        }

        if (result.count("legacy-header")) {
            // For handling indexes with old header format
            movi_options.set_legacy_header(true);
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

        if (result.count("debug")) {
            movi_options.set_debug(true);
        }

        // Set command-specific flags
        if (result.count("command")) {
            // Get the command from the result (must be the same as what was computed by get_action)
            command = result["command"].as<std::string>();
            movi_options.set_command(command);


            if (command == "build") {
                if (result.count("index") == 1 and result.count("fasta") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    movi_options.set_ref_file(result["fasta"].as<std::string>());
                    if (result.count("separators")) { movi_options.set_use_separators(true); }
                    if (result.count("color")) { movi_options.set_color(true); }
                    if (result.count("color-vectors") == 1) { movi_options.set_doc_sets_vector_of_vectors(true); }
                    if (result.count("ftab-k") >= 1) { movi_options.set_ftab_k(static_cast<uint32_t>(result["ftab-k"].as<uint32_t>())); }
                    if (result.count("multi-ftab") >= 1) { movi_options.set_multi_ftab(true); }
                    if (result.count("output-ids") >= 1) { movi_options.set_output_ids(true); }
                    if (result.count("max-run-length") >= 1) {
                        if (static_cast<uint64_t>(result["max-run-length"].as<uint64_t>()) > MAX_RUN_LENGTH) {
                            WARNING_MSG("The input maximum run length is greater than the maximum allowed, setting it to " + std::to_string(MAX_RUN_LENGTH) + ".");
                            movi_options.set_max_run_length(MAX_RUN_LENGTH);
                        } else {
                            movi_options.set_max_run_length(static_cast<uint64_t>(result["max-run-length"].as<uint64_t>()));
                        }
                    }
                    if (result.count("verify")) {
                        movi_options.set_verify(true);
                    }
                    if (result.count("preprocessed")) {
                        movi_options.set_preprocessed(true);
                    }
                    if (result.count("thresholds")) {
                        movi_options.set_thresholds(true);
                    }
#if TALLY_MODES
                    if (result.count("checkpoint") >= 1) {
                        movi_options.set_tally_checkpoints(static_cast<uint32_t>(result["checkpoint"].as<uint32_t>()));
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
                            INFO_MSG("threads: " + std::to_string(result["threads"].as<int>()));
                        }
                        movi_options.set_threads(static_cast<size_t>(result["threads"].as<int>()));
                    }
                } else {
                    const std::string message = "Please include one index directory and one fasta file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "query") {
                if (result.count("index") == 1 and result.count("read") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    movi_options.set_read_file(result["read"].as<std::string>());
                    if (result.count("out-file")) { movi_options.set_out_file(result["out-file"].as<std::string>()); }
                    if (result.count("k") >= 1) { movi_options.set_k(static_cast<uint32_t>(result["k"].as<uint32_t>())); }
                    if (result.count("min-mem-length") >= 1) { movi_options.set_min_mem_length(static_cast<uint32_t>(result["min-mem-length"].as<uint32_t>())); }
                    if (result.count("ftab-k") >= 1) { movi_options.set_ftab_k(static_cast<uint32_t>(result["ftab-k"].as<uint32_t>())); }
                    if (result.count("bin-width") >= 1) { movi_options.set_bin_width(static_cast<uint32_t>(result["bin-width"].as<uint32_t>())); }
                    if (result.count("multi-ftab") >= 1) { movi_options.set_multi_ftab(true); }
                    if (result.count("kmer") >= 1) { movi_options.set_kmer(); movi_options.set_prefetch(false); }
                    if (result.count("kmer-count") >= 1) { movi_options.set_kmer(); movi_options.set_kmer_count(true); movi_options.set_prefetch(false); }
                    if (result.count("mem") >= 1) { movi_options.set_mem(); }
                    if (result.count("count") >= 1) { movi_options.set_count(); }
                    if (result.count("zml") >= 1) { movi_options.set_zml(); }
                    if (result.count("pml") >= 1) { movi_options.set_pml(); }
                    if (result.count("sa-entries") >= 1) { movi_options.set_get_sa_entries(true); }
                    if (result.count("rpml") >= 1) { movi_options.set_random_repositioning(true); }
                    if (result.count("classify") >= 1) { movi_options.set_classify(true); }
                    if (result.count("filter") >= 1) { movi_options.set_filter(true); }
                    if (result.count("invert") >= 1) { movi_options.set_invert(true); }
                    if (result.count("multi-classify") >= 1) { movi_options.set_multi_classify(true); }
                    if (result.count("early-stop") >= 1) { movi_options.set_early_stop(true); }
                    if (result.count("report-all") >= 1) { movi_options.set_report_all(true); }
                    if (result.count("min-diff-frac") >= 1) { movi_options.set_min_diff_frac(static_cast<float>(result["min-diff-frac"].as<float>())); }
                    if (result.count("min-score-frac") >= 1) { movi_options.set_min_score_frac(static_cast<float>(result["min-score-frac"].as<float>())); }
                    if (result.count("min-len")) { movi_options.set_min_match_len(result["min-len"].as<uint8_t>()); }
                    if (result.count("pvalue-scoring") == 1) { movi_options.set_pvalue_scoring(true); }
                    if (result.count("color-move-rows") == 1) { movi_options.set_color_move_rows(true); }
                    if (result.count("color-vectors") == 1) { movi_options.set_doc_sets_vector_of_vectors(true); }
                    if (result.count("reverse") == 1) { movi_options.set_reverse(true); }
                    if (result.count("report-colors") >= 1) {
                        movi_options.set_report_colors(true);
                        // Needs to be in the multi-classify mode to report colors
                        movi_options.set_multi_classify(true);
                    }
                    if (result.count("report-color-ids") >= 1) {
                        movi_options.set_report_color_ids(true);
                        // Needs to be in the multi-classify mode to report color ids
                        movi_options.set_multi_classify(true);
                    }
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
                            INFO_MSG("Number of strands: " + std::to_string(result["strands"].as<int>()));
                        }
                        movi_options.set_strands(static_cast<size_t>(result["strands"].as<int>()));
                    }
                    if (result.count("threads") == 1) {
                        if (!supress_messages) {
                            INFO_MSG("Number of threads: " + std::to_string(result["threads"].as<int>()));
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
            } else if (command == "rlbwt") {
                if (result.count("bwt-file") == 1) {
                    movi_options.set_bwt_file(result["bwt-file"].as<std::string>());
                } else {
                    const std::string message = "Please specify one bwt file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "view") {
                if (result.count("bpf") == 1) {
                    movi_options.set_mls_file(result["bpf"].as<std::string>());

                    if (result.count("no-header") >= 1) {
                        // To handle legacy BPF files (mls files, e.g., pmls.bin)
                        // Might use small-bfp or large-bfp options to read the file
                        movi_options.set_no_header(true);
                    }

                    if (result.count("small-bpf") >= 1) {
                        movi_options.set_small_pml_lens(true);
                    } else if (result.count("large-bpf") >= 1) {
                        movi_options.set_small_pml_lens(false);
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
                    const std::string message = "Please specify a base profile format (BPF) file to view.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "color-move-rows") {
                if (result.count("index") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                } else {
                    const std::string message = "Please specify the index directory file.";
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
                cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
            }

        } else {
            const std::string message = "The action is missing.";
            cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
        }
    } catch (const cxxopts::exceptions::exception& e) {
        if (!supress_messages) {
            std::cerr << ERROR_MSG("Error parsing command line options: " + std::string(e.what())) << std::endl;
            INFO_MSG(options.help(help_groups));
        }

        return false;
    }
    return true;
}