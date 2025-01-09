#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <filesystem>
#include <sys/stat.h>
#include <regex>

std::string construct_command(const std::string& binary, const std::vector<std::string>& args) {
    std::ostringstream cmd;
    cmd << binary << " build";
    for (const auto& arg : args) {
        cmd << " " << arg;
    }
    return cmd.str();
}


// Function to parse the index type and filter the remaining arguments
bool parse_arguments(int argc, char* argv[],
                    std::string& index_type,
                    std::string& fasta_file,
                    std::string& index_path,
                    std::vector<std::string>& all_args) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--index-type") {
            if (i + 1 < argc) {
                index_type = argv[++i]; // Get the next argument as the index type
            } else {
                throw std::runtime_error("Error: --index-type flag requires a value.\n");
            }
        } else if (std::string(argv[i]) == "-i" || std::string(argv[i]) == "--index") {
            if (i + 1 < argc) {
                all_args.push_back(argv[i]);
                all_args.push_back(argv[i + 1]);
                index_path = argv[++i]; // Get the next argument as the index path
            } else {
                throw std::runtime_error("Error: " + std::string(argv[i]) + " flag requires a value.");
            }
        } else if (std::string(argv[i]) == "-f" || std::string(argv[i]) == "--fasta") {
            if (i + 1 < argc) {
                fasta_file = argv[++i]; // Get the next argument as the fasta file
            } else {
                throw std::runtime_error("Error: " + std::string(argv[i]) + " flag requires a value.\n");
            }
        } else {
            all_args.push_back(arg); // Collect other arguments
        }
    }

    if (index_type.empty()) {
        // std::cerr << "Error: --index-type flag is required.\n";
        std::cerr << "Default index is selected: default.\n";
        index_type = "default";
        return true;
    }

    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: movi-builder --index <index directory> --fasta <fasta file> --index-type <type> [other options]\n";
        return 0;
    }

    std::string index_type;
    std::string index_path;
    std::string fasta_file;
    std::vector<std::string> all_args;

    std::string binary_dir = std::filesystem::path(argv[0]).parent_path().string();

    try {
        // Parse the index type and other arguments
        if (!parse_arguments(argc, argv, index_type, fasta_file, index_path, all_args)) {
            return 0;
        }

        // Map index type to binary
        std::unordered_map<std::string, std::string> index_type_to_binary = {
            {"default", "movi-default"},
            {"constant", "movi-constant"},
            {"split", "movi-split"},
            {"compact", "movi-compact"},
            {"compact-thresholds", "movi-compact-thresholds"},
            {"blocked", "movi-blocked"},
            {"blocked-thresholds", "movi-blocked-thresholds"},
            {"movi-tally", "movi-tally"},
            {"movi-tally-thresholds", "movi-tally-thresholds"}
        };

        if (index_type_to_binary.find(index_type) == index_type_to_binary.end()) {
            throw std::runtime_error("Error: Unrecognized index_type '" + index_type + "'");
        }
        std::string binary = index_type_to_binary[index_type];

        // Execute the prepare_ref step
        std::string preprocess_command = binary_dir + "/prepare_ref " + fasta_file + " " + index_path + "/ref.fa ";
        std::cout << "Executing: " << preprocess_command << "\n";
        mkdir(index_path.c_str(),0777);
        if (!std::system(preprocess_command.c_str())) {
            std::runtime_error("Exiting in the preppare fasta step.");
        }
        std::string clean_fasta = index_path + "/ref.fa";
        all_args.push_back("--fasta");
        all_args.push_back(clean_fasta);

        // Execute the pfp_thresholds step
        std::string pfp_command = binary_dir + "/pfp-thresholds-prefix/src/pfp-thresholds-build/pfp_thresholds -f " + clean_fasta;
        std::cout << "Executing: " << pfp_command << "\n";
        if (!std::system(pfp_command.c_str())) {
            std::runtime_error("Exiting in the pfp_thresholds step.");
        }

        // Execute the binary with the appropriate options
        std::string command = construct_command(binary_dir + "/" + binary, all_args);
        std::cout << "Executing: " << command << "\n";
        if (!std::system(command.c_str())) {
            std::runtime_error("Exiting in the movi build step.");
        }

        // Remove extra files
        std::filesystem::remove(clean_fasta);
        std::filesystem::remove(clean_fasta + ".occ");
        std::filesystem::remove(clean_fasta + ".parse");
        std::filesystem::remove(clean_fasta + ".thr");
        std::filesystem::remove(clean_fasta + ".thr_pos");
        std::filesystem::remove(clean_fasta + ".thr.log");
        std::filesystem::remove(clean_fasta + ".bwt");
        std::filesystem::remove(clean_fasta + ".dict");
        std::filesystem::remove(clean_fasta + ".esa");
        std::filesystem::remove(clean_fasta + ".ssa");
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 0;
    }
}