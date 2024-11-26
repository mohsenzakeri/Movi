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
#include "cxxopts.hpp"

#include "move_structure.hpp"
#include "move_query.hpp"
#include "read_processor.hpp"
#include "movi_options.hpp"

// STEP 1: declare the type of file handler and the read() function
// KSEQ_INIT(gzFile, gzread)
std::string program() {
#if MODE == 0
    return "default";
#endif
#if MODE == 1
    return "constant";
#endif
#if MODE == 3
    return "compact";
#endif
#if MODE == 4
    return "split";
#endif
#if MODE == 5
    return "tally";
#endif
#if MODE == 6
    return "compact-thresholds";
#endif
#if MODE == 7
    return "tally-thresholds";
#endif
}

kseq_t* open_kseq(gzFile& fp, std::string file_address) {
    std::cerr << "file_address: " << file_address << "\n";
    kseq_t *seq;
    // Solution for handling the stdin input: https://biowize.wordpress.com/2013/03/05/using-kseq-h-with-stdin/
    if (file_address == "-") {
        FILE *instream = stdin;
        fp = gzdopen(fileno(instream), "r"); // STEP 2: open the file handler
    } else {
        fp = gzopen(file_address.c_str(), "r"); // STEP 2: open the file handler
    }
    seq = kseq_init(fp); // STEP 3: initialize seq
    return seq;
}

void close_kseq(kseq_t *seq, gzFile& fp) {
    kseq_destroy(seq); // STEP 5: destroy seq
    std::cerr << "kseq destroyed!\n";
    gzclose(fp); // STEP 6: close the file handler
    std::cerr << "fp file closed!\n";
}

bool parse_command(int argc, char** argv, MoviOptions& movi_options) {
    // movi_options.print_options();

    cxxopts::Options options("movi-" + program(), "Please use the following format.");

    options.add_options()
        ("command", "Command to execute", cxxopts::value<std::string>())
        ("h,help", "Print help")
        ("d,dbg", "Enable debug mode")
        ("v,verbose", "Enable verbose mode")
        ("l,logs", "Enable logs");

    auto buildOptions = options.add_options("build")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("f,fasta", "Reference file", cxxopts::value<std::string>())
        ("thresholds", "Store the threshold values in the compact mode by splitting the runs at threshold boundaries")
        ("preprocessed", "The BWT is preprocessed into heads and lens files")
        ("verify", "Verify if all the LF_move operations are correct")
        ("output-ids", "Output the adjusted ids of all the runs to ids.* files, one file per character")
        ("ftab-k", "The length of the ftab kmer", cxxopts::value<uint32_t>())
        ("tally", "Sample id at every tally runs", cxxopts::value<uint32_t>())
        ("multi-ftab", "Use ftabs with smaller k values if the largest one fails");

    auto queryOptions = options.add_options("query")
        ("pml", "Compute the pseudo-matching lengths (PMLs)")
        ("zml", "Compute the Ziv-Merhav cross parsing length (ZMLs)")
        ("count", "Compute the count queries")
        ("kmer", "Search all the kmers")
        ("kmer-count", "Find the count of every kmer")
        ("reverse", "Use the reverse (not reverse complement) of the reads to perform queries")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("r,read", "fasta/fastq Read file for query", cxxopts::value<std::string>())
        ("n,no-prefetch", "Disable prefetching for query")
        ("k,k-length", "The length of the kmer", cxxopts::value<uint32_t>())
        ("ftab-k", "The length of the ftba kmer", cxxopts::value<uint32_t>())
        ("multi-ftab", "Use ftabs with smaller k values if the largest one fails")
        ("s,strands", "Number of strands for query", cxxopts::value<int>())
        ("stdout", "Write the output to stdout")
        ("ignore-illegal-chars", "In the case of illegal characters (i.e., non-ACGT for genomic data), substitute the character with \'A\'(1) or a random character from the alphabet (2).", cxxopts::value<int>());

    auto viewOptions = options.add_options("view")
        ("mls-file", "The matching lengths (PML or ZML) file in the binary format", cxxopts::value<std::string>());

    auto rlbwtOptions = options.add_options("rlbwt")
        ("bwt-file", "BWT file", cxxopts::value<std::string>());

    auto LFOptions = options.add_options("LF")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("type", "type of the LF query: \"reconstruct\", \"sequential\", or \"random\"", cxxopts::value<std::string>());

    auto statsOptions = options.add_options("stats")
        ("output-ids", "Output the adjusted ids of all the runs to ids.* files, one file per character")
        ("i,index", "Index directory", cxxopts::value<std::string>());

    auto ftabOptions = options.add_options("ftab")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("ftab-k", "The length of the ftab kmer", cxxopts::value<uint32_t>())
        ("multi-ftab", "Use ftabs with smaller k values if the largest one fails");

    options.parse_positional({ "command" });

    try {
        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cerr << options.help() << std::endl;
            return 0;
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
#if MODE == 5 or MODE == 7
                    if (result.count("tally") >= 1) {
                        movi_options.set_tally_checkpoints(static_cast<uint32_t>(result["tally"].as<uint32_t>()));
                    }
#endif
#if MODE == 0 or MODE == 1 or MODE == 4 or MODE == 6 or MODE == 7
                    // In these modes, thresholds are always stored
                    movi_options.set_thresholds(true);
#endif
                } else {
                    const std::string message = "Please include one index directory and one fasta file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "query") {
                if (result.count("index") == 1 and result.count("read") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    movi_options.set_read_file(result["read"].as<std::string>());
                    if (result.count("k") >= 1) { movi_options.set_k(static_cast<uint32_t>(result["k"].as<uint32_t>())); }
                    if (result.count("ftab-k") >= 1) { movi_options.set_ftab_k(static_cast<uint32_t>(result["ftab-k"].as<uint32_t>())); }
                    if (result.count("multi-ftab") >= 1) { movi_options.set_multi_ftab(true); }
                    if (result.count("kmer") >= 1) { movi_options.set_kmer(); }
                    if (result.count("kmer-count") >= 1) { movi_options.set_kmer(); movi_options.set_kmer_count(true); }
                    if (result.count("count") >= 1) { movi_options.set_count(); }
                    if (result.count("zml") >= 1) { movi_options.set_zml(); }
                    if (result.count("pml") >= 1) { movi_options.set_pml(); }
                    if (result.count("reverse") == 1) { movi_options.set_reverse(true); }
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
                    if (result.count("stdout")) {
                        // Set global verbose flag
                        movi_options.set_stdout(true);
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
                if (result.count("mls-file") == 1) {
                    movi_options.set_mls_file(result["mls-file"].as<std::string>());
                } else {
                    const std::string message = "Please specify one pml file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "LF") {
                if (result.count("index") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    if (result.count("type")) {
                        if (!movi_options.set_LF_type(result["type"].as<std::string>())) {
                            const std::string message = "The LF type is not defined, please choose from: \"reconstruct\", \"sequential\", or \"random\"";
                            cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                        }
                    }
                } else {
                    const std::string message = "Please specify the index directory file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "stats") {
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
            } else {
                const std::string message = "Invalid command: \"" + command + "\"";
                cxxopts::throw_or_mimic<cxxopts::exceptions::no_such_option>(message);
            }
        } else {
            const std::string message = "No command specified.";
            cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
        }
    } catch (const cxxopts::exceptions::exception& e) {
        std::cerr << "Error parsing command line options: " << e.what() << "\n";
        std::cerr << options.help() << "\n";
        return false;
    }
    return true;
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

void query(MoveStructure& mv_, MoviOptions& movi_options) {
    if (movi_options.get_ftab_k() != 0) {
        mv_.read_ftab();
        std::cerr<<"Ftab was read!\n";
    }

    if (!movi_options.no_prefetch()) {
        ReadProcessor rp(movi_options.get_read_file(), mv_, movi_options.get_strands(), movi_options.is_verbose(), movi_options.is_reverse());
        if (movi_options.is_pml() or movi_options.is_zml() or movi_options.is_count()) {
#if MODE == 5 or MODE == 7
            rp.process_latency_hiding_tally();
#endif
#if MODE == 0 or MODE == 1 or MODE == 4 or MODE == 3 or MODE == 6
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
        std::string index_type = mv_.index_type();
        if (movi_options.is_logs()) {
            costs_file = std::ofstream(movi_options.get_read_file() + "." + index_type + ".costs");
            scans_file = std::ofstream(movi_options.get_read_file() + "." + index_type + ".scans");
            fastforwards_file = std::ofstream(movi_options.get_read_file() + "." + index_type + ".fastforwards");
        }
        uint64_t total_ff_count = 0;

        std::ofstream mls_file;
        std::ofstream count_file;
        if (!movi_options.is_stdout()) {
            if (movi_options.is_pml())
                mls_file = std::ofstream(movi_options.get_read_file() + "." + index_type + ".pml.bin", std::ios::out | std::ios::binary);
            else if (movi_options.is_zml())
                mls_file = std::ofstream(movi_options.get_read_file() + "." + index_type + ".zml.bin", std::ios::out | std::ios::binary);
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
                if (movi_options.is_pml())
                    total_ff_count += mv_.query_pml(mq, random_jump);
                else if (movi_options.is_zml()) {
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
    }
}

int main(int argc, char** argv) {
    MoviOptions movi_options;
    if (!parse_command(argc, argv, movi_options)) {
        return 0;
    }
    std::string command = movi_options.get_command();
    if (command == "build") {
        MoveStructure mv_(&movi_options, MODE == 1 or MODE == 4, MODE == 1);
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
    } else if (command == "stats") {
        MoveStructure mv_(&movi_options);
        mv_.deserialize();
        mv_.print_stats();
        // mv_.compute_run_lcs();
        // mv_.analyze_rows();
    } else if (command == "ftab") {
        MoveStructure mv_(&movi_options);
        mv_.deserialize();
        build_ftab(mv_, movi_options);
    }
}
