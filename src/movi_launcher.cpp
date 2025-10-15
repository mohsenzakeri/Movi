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

#include "utils.hpp"

#define PFP_THRESHOLDS_PATH "/external_repos/pfp-thresholds-build/pfp_thresholds"
#define R_PERMUTE_PATH "/external_repos/r-permute-build/test/src/"
#define DEFAULT_INDEX_TYPE "regular-thresholds"
#define DEFAULT_INDEX_TYPE_CODE 6

/* * * * * * * * * * * * *
arguments specific to the launcher:
--keep: keep the extra files after the build step
--skip-prepare: skip the prepare_ref step
--skip-pfp: skip the pfp_thresholds step
--skip-rlbwt: skip the rlbwt step
--skip-r-permute: skip the r-permute step
--type: the index type to build
--list or -l: the list of fasta files
--non-preprocessed: undo the default behavior of writing the BWT in the run length compressed format (ignored if --preprocessed is set)
 * * * * * * * * * * * * */

std::string binary_dir; // set from argv[0]

struct Args {
    std::string action;
    std::string index_path;
    std::string index_type = DEFAULT_INDEX_TYPE;
    std::string fasta_file;
    bool verbose = false;
    bool index_provided = true;
    bool fasta_provided = false;
    bool fasta_list_provided = false;
    bool remove_extra_files = true;
    bool skip_prepare = false;
    bool skip_pfp = false;
    bool skip_rlbwt = false;
    bool skip_r_permute = false;
    // By default, write the BWT in the run length compressed format
    bool preprocessed = true;
    bool preprocessed_flag = false;
    bool non_preprocessed_flag = false;
};

// Index type mapping
const std::unordered_map<std::string, char> type_to_char = {
    {"large", 0},
    {"constant", 1},
    {"split", 4},
    {"regular", 3},
    {"regular-thresholds", 6},
    {"blocked", 2},
    {"blocked-thresholds", 8},
    {"tally", 5},
    {"tally-thresholds", 7}
};

const std::unordered_map<char, std::string> char_to_type = {
    {0, "large"},
    {1, "constant"},
    {4, "split"},
    {3, "regular"},
    {6, "regular-thresholds"},
    {2, "blocked"},
    {8, "blocked-thresholds"},
    {5, "tally"},
    {7, "tally-thresholds"}
};

// Function declarations
std::string construct_movi_command(const std::string& binary, const std::vector<std::string>& all_args);
void execute_command_line(const std::string& command, const std::string& error_msg, const Args& args);
void throw_missing_argument(const std::string& argument);
void print_error(const std::string& error_msg);
void print_warning(const std::string& warning_msg);
void parse_build_arguments(int argc, char* argv[], Args& args, std::vector<std::string>& all_args);
void parse_arguments(int argc, char* argv[], Args& args, std::vector<std::string>& all_args);
void handle_build(const Args& args, const std::vector<std::string>& all_args);
void handle(const Args& args, const std::vector<std::string>& all_args);
char get_index_type(const std::string& index_dir);
void cleanup_files(const std::string& clean_fasta, const std::string& index_type);

int main(int argc, char* argv[]) {
    binary_dir = std::filesystem::path(argv[0]).parent_path().string();

    Args args;
    std::vector<std::string> all_args;

    try {
        if (argc < 2) {
            throw std::runtime_error("Error parsing command line options: The action is missing.");
        } else {
            args.action = argv[1];
        }

        if (args.action == "build") {
            parse_build_arguments(argc, argv, args, all_args);
            handle_build(args, all_args);
        } else {
            parse_arguments(argc, argv, args, all_args);
            handle(args, all_args);
        }
    } catch (const std::exception& e) {
        if (std::filesystem::exists(binary_dir + "/bin/movi-" + args.index_type)) {
            print_error(std::string(e.what()));
            std::string help_command = binary_dir + "/bin/movi-" + args.index_type + " " + args.action + " -h";
            execute_command_line(help_command, "Failed in movi " + args.action + " step", args);
        } else {
            print_error("Error parsing command line options: The binary '" + binary_dir + "/bin/movi-" + args.index_type + "' does not exist.");
        }
        return 1;
    }

    return 0;
}

std::string construct_movi_command(const std::string& binary,
                              const std::vector<std::string>& all_args) {
    // Checke if the binary exists
    if (!std::filesystem::exists(binary_dir + "/" + binary)) {
        throw std::runtime_error("Error parsing command line options: The binary '" + binary_dir + "/" + binary + "' does not exist.");
    }

    std::ostringstream cmd;
    cmd << binary_dir << "/" << binary;
    for (const auto& arg : all_args) {
        cmd << " " << arg;
    }
    return cmd.str();
}

void execute_command_line(const std::string& command, const std::string& error_msg, const Args& args) {
    if (args.verbose) {
        INFO_MSG("Executing: " + command);
    }
    if (std::system(command.c_str()) != 0 and args.action != "view") {
        print_error("Error executing the command: " + error_msg);
        exit(1);
    }
}

void throw_missing_argument(const std::string& argument) {
    throw std::runtime_error("Error parsing command line options: Option '" + argument + "' is missing an argument.");
}

void print_error(const std::string& error_msg) {
    std::cerr << ERROR_MSG(error_msg);
}

void print_warning(const std::string& warning_msg) {
    WARNING_MSG(warning_msg);
}

// Implementation of build-specific functionality
void handle_build(const Args& args, const std::vector<std::string>& all_args) {
    // Create index directory
    mkdir(args.index_path.c_str(), 0777);

    // Prepare reference
    std::string clean_fasta = args.index_path + "/ref.fa";
    if (!args.skip_prepare) {
        std::string prepare_ref_command = binary_dir + "/bin/movi-prepare-ref "
                                        + args.fasta_file + " " + clean_fasta
                                        + (args.fasta_list_provided ? " list" : "");
        execute_command_line(prepare_ref_command, "Failed in prepare fasta step", args);
    }

    // Execute pfp_thresholds
    if (!args.skip_pfp) {
        std::string pfp_command = binary_dir + PFP_THRESHOLDS_PATH + " -P ";
        if (args.preprocessed) {
            pfp_command += " -r ";
        }
        pfp_command += " -f " + clean_fasta;
        execute_command_line(pfp_command, "Failed in pfp_thresholds step", args);
    }

    // Handle special cases for constant/split indexes
    if (args.index_type == "constant" || args.index_type == "split") {
        std::string binary = "bin/movi-" + args.index_type;
        if (!args.skip_rlbwt) {
            std::string rlbwt_command = binary_dir + "/" + binary + " rlbwt --bwt-file " + clean_fasta + ".bwt";
            execute_command_line(rlbwt_command, "Failed in rlbwt step", args);
        }
        if (!args.skip_r_permute) {
            std::string build_constructor_command = binary_dir + R_PERMUTE_PATH + "/build_constructor " + clean_fasta;
            execute_command_line(build_constructor_command, "Failed in r-permute build_constructor step", args);
            std::string run_constructor_command = binary_dir + R_PERMUTE_PATH + "/run_constructor " + clean_fasta + " -d 5";
            execute_command_line(run_constructor_command, "Failed in r-permute run_constructor step", args);
        }
    }

    // Build the index
    std::vector<std::string> build_args = all_args;
    build_args.push_back("--fasta");
    build_args.push_back(clean_fasta);
    
    std::string binary = "bin/movi-" + args.index_type;
    std::string build_command = construct_movi_command(binary, build_args);
    execute_command_line(build_command, "Failed in movi build step", args);

    // Cleanup if requested
    if (args.remove_extra_files) {
        cleanup_files(clean_fasta, args.index_type);
    }
}

void handle(const Args& args, const std::vector<std::string>& all_args) {
    char index_type_char = args.index_provided ? get_index_type(args.index_path) : type_to_char.at(args.index_type);

    if (char_to_type.find(index_type_char) == char_to_type.end()) {
        throw std::runtime_error("Error parsing command line options: Unrecognized index_type '" + std::to_string(index_type_char) + "'");
    }

    std::string binary = "bin/movi-" + char_to_type.at(index_type_char);
    std::string command = construct_movi_command(binary, all_args);
    execute_command_line(command, "Failed in movi " + args.action + " step", args);
}

// Function to parse the index type and filter the remaining arguments
void parse_build_arguments(int argc, char* argv[],
                     Args& script_args,
                     std::vector<std::string>& all_args) {

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--verbose") {
            script_args.verbose = true;
            all_args.push_back(arg);
        } else if (arg == "--keep") {
            script_args.remove_extra_files = false;
        } else if (arg == "--skip-prepare") {
            script_args.skip_prepare = true;
        } else if (arg == "--skip-pfp") {
            script_args.skip_pfp = true;
        } else if (arg == "--skip-rlbwt") {
            script_args.skip_rlbwt = true;
        } else if (arg == "--skip-r-permute") {
            script_args.skip_r_permute = true;
        } else if (arg == "--non-preprocessed") {
            script_args.non_preprocessed_flag = true;
        } else if (arg == "--preprocessed") {
            script_args.preprocessed_flag = true;
        } else if (arg == "--type") {
            if (i + 1 < argc) {
                script_args.index_type = argv[++i]; // Get the next argument as the index type
            } else {
                throw_missing_argument(arg);
            }
        } else if (arg == "-i" || arg == "--index") {
            if (i + 1 < argc) {
                all_args.push_back(argv[i]);
                all_args.push_back(argv[i + 1]);
                script_args.index_path = argv[++i]; // Get the next argument as the index path
            } else {
                throw_missing_argument(arg);
            }
        } else if (arg == "-f" || arg == "--fasta") {
            if (i + 1 < argc) {
                script_args.fasta_file = argv[++i]; // Get the next argument as the fasta file
                script_args.fasta_provided = true;
            } else {
                throw_missing_argument(arg);
            }
        } else if (arg == "-l" || arg == "--list") {
            if (i + 1 < argc) {
                script_args.fasta_file = argv[++i]; // Get the next argument as the fasta file (list of fastas)
                script_args.fasta_list_provided = true;
            } else {
                throw_missing_argument(arg);
            }
        } else {
            all_args.push_back(arg); // Collect other arguments
        }
    }

    if (script_args.index_type.empty()) {
        script_args.index_type = DEFAULT_INDEX_TYPE;
        print_warning("Default index is selected (" + std::string(DEFAULT_INDEX_TYPE) + ").\n");
    }

    if (script_args.fasta_file.empty()) {
        throw std::runtime_error("The fasta file or list of fasta files should be provided");
    }

    if (script_args.index_path.empty()) {
        throw std::runtime_error("The index directory should be provided");
    }

    // By default, write the BWT in the run length compressed format
    // Only write the uncompressed BWT if the --non-preprocessed flag is set and the --preprocessed flag is not set
    script_args.preprocessed = script_args.preprocessed_flag or !script_args.non_preprocessed_flag;
    if (script_args.preprocessed) {
        all_args.push_back("--preprocessed");
        script_args.skip_rlbwt = true;
    }
}

void parse_arguments(int argc, char* argv[],
                     Args& script_args,
                     std::vector<std::string>& all_args) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--verbose") {
            script_args.verbose = true;
            all_args.push_back(arg);
        } else if (arg == "-i" || arg == "--index") {
            if (i + 1 < argc) {
                all_args.push_back(argv[i]);
                all_args.push_back(argv[i + 1]);
                script_args.index_path = argv[++i]; // Get the next argument as the index type
            } else {
                throw_missing_argument(arg);
            }
        } else {
            all_args.push_back(arg); // Collect other arguments
        }
    }

    if (script_args.index_path.empty()) {
        script_args.index_provided = false;
        if (script_args.index_type.empty()) {
            script_args.index_type = DEFAULT_INDEX_TYPE;
        }
    }
}

char get_index_type(const std::string& index_dir) {
    std::vector<std::string> possible_filenames = {
        index_dir + "/index.movi",
        index_dir + "/movi_index.bin"
    };

    for (const auto& fname : possible_filenames) {
        std::ifstream fin(fname, std::ios::binary);
        if (fin) {
            // Try to read as modern header first
            MoviHeader header;
            fin.read(reinterpret_cast<char*>(&header), sizeof(MoviHeader));

            // Check if it's a modern header by verifying magic number
            if (header.magic == MOVI_MAGIC) {
                return static_cast<char>(header.type);
            }

            // If not a modern header, try as legacy header
            fin.seekg(0, std::ios::beg);
            char index_type;
            fin.read(reinterpret_cast<char*>(&index_type), sizeof(index_type));
            return index_type;
        }
    }
    throw std::runtime_error("Error parsing command line options: Unable to open the index file at: " + index_dir);
}

// Clean up temporary files
void cleanup_files(const std::string& clean_fasta, const std::string& index_type) {
    std::vector<std::string> files_to_remove = {
        clean_fasta,
        clean_fasta + ".occ",
        clean_fasta + ".parse", 
        clean_fasta + ".thr",
        clean_fasta + ".thr_pos",
        clean_fasta + ".thr.log",
        clean_fasta + ".bwt",
        clean_fasta + ".dict",
        clean_fasta + ".esa",
        clean_fasta + ".ssa"
    };

    if (index_type == "constant" || index_type == "split") {
        files_to_remove.insert(files_to_remove.end(), {
            clean_fasta + ".d_col",
            clean_fasta + ".d_construct", 
            clean_fasta + ".bwt.heads",
            clean_fasta + ".bwt.len"
        });
    }

    for (const auto& file : files_to_remove) {
        std::filesystem::remove(file);
    }
}
