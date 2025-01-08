#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <stdexcept>

char get_index_type(const std::string& index_dir) {
    // Open the index file and read the index type
    std::string fname = index_dir + "/movi_index.bin";
    std::ifstream fin(fname, std::ios::binary);
    if (!fin) {
        throw std::runtime_error("Error: Unable to open index file " + fname);
    }
    fin.seekg(0, std::ios::beg);
    char index_type;
    fin.read(reinterpret_cast<char*>(&index_type), sizeof(index_type));
    return index_type;
}

std::string construct_command(const std::string& binary, const std::vector<std::string>& args) {
    std::ostringstream cmd;
    cmd << "./" << binary;
    for (const auto& arg : args) {
        cmd << " " << arg;
    }
    return cmd.str();
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " -i <index_path> [other_args...]\n";
        return 1;
    }

    // Parse arguments
    std::string index_path;
    std::vector<std::string> all_args;
    for (int i = 1; i < argc; ++i) {
        if ( (std::string(argv[i]) == "-i" || std::string(argv[i]) == "--index") && i + 1 < argc) {
            all_args.push_back(argv[i]);
            all_args.push_back(argv[i + 1]);
            index_path = argv[++i];
        } else {
            all_args.push_back(argv[i]);
        }
    }

    if (index_path.empty()) {
        std::cerr << "Error: Index file path (-i) is required.\n";
        return 1;
    }

    try {
        // Step 1: Extract index_type
        char index_type = get_index_type(index_path);

        // Step 2: Map index_type to binary
        std::unordered_map<char, std::string> index_type_to_binary = {
            {0, "movi-default"},
            {1, "movi-constant"},
            {4, "movi-split"},
            {3, "movi-compact"},
            {6, "movi-compact-thresholds"},
            {2, "movi-blocked"},
            {8, "movi-blocked-thresholds"},
            {5, "movi-tally"},
            {7, "movi-tally-thresholds"}
        };

        if (index_type_to_binary.find(index_type) == index_type_to_binary.end()) {
            throw std::runtime_error("Error: Unrecognized index_type '" + std::to_string(index_type) + "'");
        }

        std::string binary = index_type_to_binary[index_type];

        // Step 3: Construct the command
        std::string command = construct_command(binary, all_args);
        std::cerr << command << "\n";

        // Step 4: Execute the command
        int ret_code = std::system(command.c_str());
        return ret_code;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}