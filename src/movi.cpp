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
}

kseq_t* open_kseq(gzFile& fp, std::string file_address) {
    kseq_t *seq;
    fp = gzopen(file_address.c_str(), "r"); // STEP 2: open the file handler
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
        ("v,verbose", "Enable verbose mode")
        ("l,logs", "Enable logs");

    auto buildOptions = options.add_options("build")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("f,fasta", "Reference file", cxxopts::value<std::string>());

    auto queryOptions = options.add_options("query")
        ("pml", "Compute the pseudo-matching lengths (PMLs)")
        ("count", "Compute the count queries")
        ("reverse", "Use the reverse (not reverse complement) of the reads to perform queries")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("r,read", "fasta/fastq Read file for query", cxxopts::value<std::string>())
        ("n,no-prefetch", "Disable prefetching for query")
        ("s,strands", "Number of strands for query", cxxopts::value<int>());

    auto viewOptions = options.add_options("view")
        ("pml-file", "PML file in the binary format", cxxopts::value<std::string>());

    auto rlbwtOptions = options.add_options("rlbwt")
        ("bwt-file", "BWT file", cxxopts::value<std::string>());

    auto LFOptions = options.add_options("LF")
        ("i,index", "Index directory", cxxopts::value<std::string>())
        ("type", "type of the LF query: \"reconstruct\", \"sequential\", or \"random\"", cxxopts::value<std::string>());

    auto statsOptions = options.add_options("stats")
        ("i,index", "Index directory", cxxopts::value<std::string>());

    options.parse_positional({ "command" });

    try {
        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help() << std::endl;
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

        if (result.count("command")) {
            std::string command = result["command"].as<std::string>();
            movi_options.set_command(command);

            if (command == "build") {
                if (result.count("index") == 1 and result.count("fasta") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    movi_options.set_ref_file(result["fasta"].as<std::string>());
                } else {
                    const std::string message = "Please include one index directory and one fasta file.";
                    cxxopts::throw_or_mimic<cxxopts::exceptions::invalid_option_format>(message);
                }
            } else if (command == "query") {
                if (result.count("index") == 1 and result.count("read") == 1) {
                    movi_options.set_index_dir(result["index"].as<std::string>());
                    movi_options.set_read_file(result["read"].as<std::string>());
                    if (result.count("count") >= 1) { movi_options.set_count(true); }
                    if (result.count("pml") >= 1) { movi_options.set_pml(true); }
                    if (result.count("reverse") == 1) { movi_options.set_reverse(true); }
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
                if (result.count("pml-file") == 1) {
                    movi_options.set_pml_file(result["pml-file"].as<std::string>());
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
                } else {
                    const std::string message = "Please specify the index directory file.";
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
        std::cout << options.help() << "\n";
        return false;
    }
    return true;
}

void query(MoveStructure& mv_, MoviOptions& movi_options) {
    if (!movi_options.no_prefetch()) {
        ReadProcessor rp(movi_options.get_read_file(), mv_, movi_options.get_strands(), movi_options.is_pml(), movi_options.is_reverse());
        if (movi_options.is_pml()) {
            rp.process_latency_hiding(mv_);
        } else if (movi_options.is_count()) {
            rp.backward_search_latency_hiding(mv_);
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

        std::ofstream pmls_file;
        std::ofstream count_file;
        if (movi_options.is_pml())
            pmls_file = std::ofstream(movi_options.get_read_file() + "." + index_type + ".mpml.bin", std::ios::out | std::ios::binary);
        else if (movi_options.is_count())
            count_file = std::ofstream(movi_options.get_read_file() + "." + index_type + ".matches");

        uint64_t read_processed = 0;
        while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
            if (read_processed % 1000 == 0)
                std::cerr << read_processed << "\r";
            read_processed += 1 ;
            std::string query_seq = seq->seq.s;
            MoveQuery mq;
            if (movi_options.is_pml()) {
                if (movi_options.is_reverse())
                    std::reverse(query_seq.begin(), query_seq.end());
                mq = MoveQuery(query_seq);
                bool random_jump = false;
                total_ff_count += mv_.query_pml(mq, random_jump);
                uint16_t st_length = seq->name.m;
                pmls_file.write(reinterpret_cast<char*>(&st_length), sizeof(st_length));
                pmls_file.write(reinterpret_cast<char*>(&seq->name.s[0]), st_length);
                auto& pml_lens = mq.get_pml_lens();
                uint64_t mq_pml_lens_size = pml_lens.size();
                pmls_file.write(reinterpret_cast<char*>(&mq_pml_lens_size), sizeof(mq_pml_lens_size));
                pmls_file.write(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));
            } else if (movi_options.is_count()) {
                // std::string R = std::string(seq->seq.s);
                int32_t pos_on_r = query_seq.length() - 1;
                uint64_t match_count = mv_.backward_search(query_seq, pos_on_r);
                count_file << seq->name.s << "\t";
                if (pos_on_r != 0) pos_on_r += 1;
                count_file << query_seq.length() - pos_on_r << "/" << query_seq.length() << "\t" << match_count << "\n";                
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
        
        if (movi_options.is_pml()) {
            std::cerr << "all fast forward counts: " << total_ff_count << "\n";
            pmls_file.close();
            std::cerr << "pmls file closed.\n";
        } else if (movi_options.is_count()) {
            count_file.close();
            std::cerr << "count file is closed.\n";
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
    std::ifstream pmls_file(movi_options.get_pml_file(), std::ios::in | std::ios::binary);
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
}

int main(int argc, char** argv) {
    // std::ios_base::sync_with_stdio(false);
    MoviOptions movi_options;
    if (!parse_command(argc, argv, movi_options)) {
        return 0;
    }
    std::string command = movi_options.get_command();
    if (command == "build") {
        MoveStructure mv_(movi_options.get_ref_file(), false, movi_options.is_verbose(), movi_options.is_logs(), false, MODE == 1);
        std::cerr << "The move structure is successfully stored at " << movi_options.get_index_dir() << "\n";
        mv_.serialize(movi_options.get_index_dir());
        std::cerr << "The move structure is read from the file successfully.\n";
    } else if (command == "query") {
        MoveStructure mv_(movi_options.is_verbose(), movi_options.is_logs());
        auto begin = std::chrono::system_clock::now();
        mv_.deserialize(movi_options.get_index_dir());
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::printf("Time measured for loading the index: %.3f seconds.\n", elapsed.count() * 1e-9);
        begin = std::chrono::system_clock::now();
        query(mv_, movi_options);
        end = std::chrono::system_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::printf("Time measured for processing the reads: %.3f seconds.\n", elapsed.count() * 1e-9);
    } else if (command == "view") {
        view(movi_options);
    } else if (command == "rlbwt") {
        std::cerr << "The run and len files are being built.\n";
        MoveStructure mv_(movi_options.is_verbose(), movi_options.is_logs());
        mv_.build_rlbwt(movi_options.get_bwt_file());
    } else if (command == "LF") {
        MoveStructure mv_(movi_options.is_verbose(), movi_options.is_logs());
        mv_.deserialize(movi_options.get_index_dir());
        std::cerr << "The move structure is read from the file successfully.\n";
        if (movi_options.get_LF_type() == "sequential")
            mv_.sequential_lf();
        else if (movi_options.get_LF_type() == "random")
            mv_.random_lf();
        else if (movi_options.get_LF_type() == "reconstruct")
            mv_.reconstruct_lf();
    } else if (command == "stats") {
        MoveStructure mv_(movi_options.is_verbose(), movi_options.is_logs());
        mv_.deserialize(movi_options.get_index_dir());
        mv_.print_stats();
    }
}

/*int main(int argc, char** argv) {
    // std::ios_base::sync_with_stdio(false);
    std::string command = (argc >= 2) ? argv[1] : "";
    std::cerr << "command: " << command << "\n";
    if (command == "build") {
        std::cerr << "The move structure is being built.\n";

        bool onebit = std::string(argv[2]) == "onebit" ? true : false;
        uint16_t splitting = std::string(argv[2]) == "split" ? 5 : 0;
        bool constant = std::string(argv[2]) == "constant" ? true : false;
        if (constant or onebit) {
            splitting = 5;
        }
        std::cerr << "splitting: " << splitting << "\n";
        std::cerr << "onebit: " << onebit << "\n";
        std::cerr << "constant: " << constant << "\n";

        bool verbose = (argc > 5 and std::string(argv[5]) == "verbose");
        bool logs = (argc > 5 and std::string(argv[5]) == "logs");

        MoveStructure mv_(argv[3], onebit, verbose, logs, splitting, constant);
        std::cerr << "The move structure is successfully built!\n";

        // mv_.reconstruct();
        // std::cerr << "The original string is reconstructed.\n";
        // std::cerr << "The original string is:\n" << mv_.R() << "\n";

        mv_.serialize(argv[4]);
        std::cerr << "The move structure is successfully stored at " << argv[4] << "\n";
        // if (logs) {
        //     std::ofstream rl_file(static_cast<std::string>(argv[4]) + "/run_lengths");
        //     for (auto& run_length : mv_.run_lengths) {
        //         rl_file <<run_length.first << "\t" << run_length.second << "\n";
        //     }
        //     rl_file.close();
        // }
    } else if (command == "query-pf") {
        bool verbose = (argc > 5 and std::string(argv[5]) == "verbose");
        bool logs = (argc > 5 and std::string(argv[5]) == "logs");
        MoveStructure mv_(verbose, logs);

        auto begin = std::chrono::system_clock::now();
        mv_.deserialize(argv[2]);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::printf("Time measured for loading the index: %.3f seconds.\n", elapsed.count() * 1e-9);
        std::cerr << "The move structure is read from the file successfully.\n";

        ReadProcessor rp(argv[3], mv_, atoi(argv[4]), true);
        rp.process_latency_hiding(mv_);
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
        // std::cerr << "The original string is: " << mv_.reconstruct() << "\n";
        // std::string query = argv[3];

        // Fasta/q reader from http://lh3lh3.users.sourceforge.net/parsefastq.shtml
        gzFile fp;
        int l;
        // kseq_t *seq;
        // fp = gzopen(argv[3], "r"); // STEP 2: open the file handler
        // seq = kseq_init(fp); // STEP 3: initialize seq
        kseq_t* seq = open_kseq(fp, argv[3]);
        // std::ofstream pmls_file(static_cast<std::string>(argv[3]) + ".mpml");
        std::ofstream costs_file;
        std::ofstream scans_file;
        std::ofstream fastforwards_file;
        std::string index_type = mv_.index_type();
        if (logs) {
            costs_file = std::ofstream(static_cast<std::string>(argv[3]) + "." + index_type + ".costs");
            scans_file = std::ofstream(static_cast<std::string>(argv[3]) + "." + index_type + ".scans");
            fastforwards_file = std::ofstream(static_cast<std::string>(argv[3]) + "." + index_type + ".fastforwards");
        }
        std::ofstream pmls_file(static_cast<std::string>(argv[3]) + "." + index_type + ".mpml.bin", std::ios::out | std::ios::binary);
        uint64_t total_ff_count = 0;
        // uint64_t query_pml_tot_time = 0;
        // uint64_t iteration_tot_time = 0;
        uint64_t read_processed = 0;
        while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
            if (read_processed % 1000 == 0)
                std::cerr << read_processed << "\r";
            read_processed += 1 ;
            // printf("name: %s\n", seq->name.s);
            // if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            // printf("seq: %s\n", seq->seq.s);
            // if (seq->qual.l) printf("qual: %s\n", seq->qual.s); 

            std::string query_seq = seq->seq.s;
            // reverse for the null reads
            if (reverse_query)
                std::reverse(query_seq.begin(), query_seq.end());
            // auto t1 = std::chrono::high_resolution_clock::now();    
            MoveQuery mq(query_seq);
            bool random_jump = false;
            // std::cerr << seq->name.s << "\n";
            total_ff_count += mv_.query_pml(mq, random_jump);
            // auto t2 = std::chrono::high_resolution_clock::now();
            // query_pml_tot_time += static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

            // pmls_file << ">" << seq->name.s << "\n";
            // for (int64_t i = mq.pml_lens.size() - 1; i >= 0; i--) {
            //     pmls_file << mq.pml_lens[i] << " ";
            // }
            // pmls_file << "\n";
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
        std::cerr << "all fast forward counts: " << total_ff_count << "\n";
        // std::cerr << "query_pml_tot_time: " << query_pml_tot_time << "\n";
        // std::cerr << "iteration_tot_time: " << iteration_tot_time << "\n";
        if (logs) {
            costs_file.close();
            scans_file.close();
            fastforwards_file.close();
        }
        pmls_file.close();
        std::cerr << "pmls file closed!\n";
        // printf("return value: %d\n", l);
        close_kseq(seq, fp);
        // kseq_destroy(seq); // STEP 5: destroy seq
        // std::cerr << "kseq destroyed!\n";
        // gzclose(fp); // STEP 6: close the file handler
        // std::cerr << "fp file closed!\n";

        if (logs) {
            std::ofstream jumps_file(static_cast<std::string>(argv[3]) + "." + index_type + ".jumps");
            for (auto& jump : mv_.jumps) {
                jumps_file <<jump.first << "\t" << jump.second << "\n";
            }
            jumps_file.close();
        }
    } else if (command == "match" or command == "match-onebit") {
        bool verbose = (argc > 4 and std::string(argv[4]) == "verbose");
        bool logs = (argc > 4 and std::string(argv[4]) == "logs");
        bool reverse_query = (argc > 4 and std::string(argv[4]) == "reverse");
        std::cerr << verbose << " " << logs << "\n";
        MoveStructure mv_(verbose, logs); 
        if (command == "match-onebit")
            mv_.set_onebit();
        auto begin = std::chrono::system_clock::now();
        mv_.deserialize(argv[2]);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::printf("Time measured for loading the index: %.3f seconds.\n", elapsed.count() * 1e-9);
        std::cerr << "The move structure is read from the file successfully.\n";

        // Fasta/q reader from http://lh3lh3.users.sourceforge.net/parsefastq.shtml
        gzFile fp;
        int l;
        // kseq_t *seq;
        // fp = gzopen(argv[3], "r"); // STEP 2: open the file handler
        // seq = kseq_init(fp); // STEP 3: initialize seq
        kseq_t* seq = open_kseq(fp, argv[3]);
        std::string index_type = mv_.index_type();
        std::ofstream matches_file(static_cast<std::string>(argv[3]) + "." + index_type + ".matches.bin", std::ios::out | std::ios::binary);
        uint64_t total_ff_count = 0;
        while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence

            std::string query_seq = seq->seq.s;
            // reverse for the null reads
            if (reverse_query)
                std::reverse(query_seq.begin(), query_seq.end());
            MoveQuery mq(query_seq);
            total_ff_count += mv_.exact_matches(mq);
            uint16_t st_length = seq->name.m;
            matches_file.write(reinterpret_cast<char*>(&st_length), sizeof(st_length));
            matches_file.write(reinterpret_cast<char*>(&seq->name.s[0]), st_length);
            auto& pml_lens = mq.get_pml_lens();
            uint64_t mq_pml_lens_size = pml_lens.size();
            matches_file.write(reinterpret_cast<char*>(&mq_pml_lens_size), sizeof(mq_pml_lens_size));
            matches_file.write(reinterpret_cast<char*>(&pml_lens[0]), mq_pml_lens_size * sizeof(pml_lens[0]));
        }
        matches_file.close();
        std::cerr << "pmls file closed!\n";
        // printf("return value: %d\n", l);
        close_kseq(seq, fp);
        // kseq_destroy(seq); // STEP 5: destroy seq
        // std::cerr << "kseq destroyed!\n";
        // gzclose(fp); // STEP 6: close the file handler
        // std::cerr << "fp file closed!\n";
    } else if (command == "count" or command == "count-onebit") {
        bool verbose = (argc > 4 and std::string(argv[4]) == "verbose");
        bool logs = (argc > 4 and std::string(argv[4]) == "logs");
        std::cerr << verbose << " " << logs << "\n";
        MoveStructure mv_(verbose, logs);
        if (command == "count-onebit")
            mv_.set_onebit();
        auto begin = std::chrono::system_clock::now();
        mv_.deserialize(argv[2]);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::printf("Time measured for loading the index: %.3f seconds.\n", elapsed.count() * 1e-9);
        std::cerr << "The move structure is read from the file successfully.\n";

        // Fasta/q reader from http://lh3lh3.users.sourceforge.net/parsefastq.shtml
        gzFile fp;
        int l;
        // kseq_t *seq;
        // fp = gzopen(argv[3], "r"); // STEP 2: open the file handler
        // seq = kseq_init(fp); // STEP 3: initialize seq
        kseq_t* seq = open_kseq(fp, argv[3]);
        std::string index_type = mv_.index_type();
        std::ofstream count_file(static_cast<std::string>(argv[3]) + "." + index_type + ".matches");
        uint64_t total_ff_count = 0;
        uint64_t read_processed = 0;
        while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
            if (read_processed % 1000 == 0)
                std::cerr << read_processed << "\r";
            std::string R = std::string(seq->seq.s);
            int32_t pos_on_r = R.length() - 1;
            uint64_t match_count = mv_.backward_search(R, pos_on_r);
            // count_file << seq->name.s << "\t" << (pos_on_r == 0 ? "Found\t" : "Not-Found\t");
            count_file << seq->name.s << "\t";
            if (pos_on_r != 0) pos_on_r += 1;
            count_file << R.length() - pos_on_r << "/" << R.length() << "\t" << match_count << "\n";
            read_processed += 1;
        }
        count_file.close();
        std::cerr << "output file closed!\n";
        close_kseq(seq, fp);
        // kseq_destroy(seq); // STEP 5: destroy seq
        // std::cerr << "kseq destroyed!\n";
        // gzclose(fp); // STEP 6: close the file handler
        // std::cerr << "fp file closed!\n";
    } else if (command == "count-pf") {
        bool verbose = (argc > 5 and std::string(argv[5]) == "verbose");
        bool logs = (argc > 5 and std::string(argv[5]) == "logs");
        MoveStructure mv_(verbose, logs);

        auto begin = std::chrono::system_clock::now();
        mv_.deserialize(argv[2]);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::printf("Time measured for loading the index: %.3f seconds.\n", elapsed.count() * 1e-9);
        std::cerr << "The move structure is read from the file successfully.\n";

        ReadProcessor rp(argv[3], mv_, atoi(argv[4]), false);
        rp.backward_search_latency_hiding(mv_);
    } else if (command == "rlbwt") {
        std::cerr << "The run and len files are being built.\n";
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
            std::string index_type = mv_.index_type();
            std::ofstream ff_counts_file(static_cast<std::string>(argv[3]) + "." + index_type + ".fastforwards");
            for (auto& ff_count : mv_.ff_counts) {
                ff_counts_file <<ff_count.first << "\t" << ff_count.second << "\n";
            }
            ff_counts_file.close();
        }
    } else {
        std::cerr << "clean fasta files:\t\t./prepare_ref <fasta list file> <output fasta> list\n\n";
        std::cerr << ">>>>>> You can run the following if the .thr_pos and .bwt files are provided for <output fasta>.\n";
        std::cerr << "build default index:\t\t./movi-default build default <output fasta> <index dir>\n\n";
        std::cerr << "query (compute PMLs):\t\t./movi-default query <index dir> <reads file>\n";
        std::cerr << "view the mpml.bin file:\t\t./movi-default view <mpml.bin file> | less\n";
    }
}*/
