#include "test_utils.hpp"

using namespace Catch::Matchers;

// Helper function to test index creation
void test_index_creation(const std::string& index_type, const std::string& index_suffix, size_t expected_size) {
    ensure_test_data_dir();
    std::string fasta_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/ref.fasta";
    REQUIRE(std::filesystem::exists(fasta_file));
    
    const std::string index_dir = TEST_DATA_DIR + "/index_" + index_suffix;
    std::string cmd = std::string(BINARY_DIR) + "/movi build --type " + index_type + " --index " + index_dir +
              " --fasta " + fasta_file + " --keep --verify > /dev/null 2>&1";

    int exit_code = system(cmd.c_str());

    // Check that command executed successfully
    REQUIRE(exit_code == 0);

    // Check that the index file exists
    const std::string index_file = index_dir + "/index.movi";
    REQUIRE(std::filesystem::exists(index_file));
    
    // Check that the index file has the expected size
    auto file_size = std::filesystem::file_size(index_file);
    REQUIRE(file_size == expected_size);
}

TEST_CASE("MoveStructure - index building", "[move_structure_index_building]") {
    SECTION("Regular index creation") {
        test_index_creation("regular", "regular", 871482);
    }

    SECTION("Regular-thresholds index creation") {
        test_index_creation("regular-thresholds", "regular_thresholds", 948122);
    }

    SECTION("Tally index creation") {
        test_index_creation("tally", "tally", 437009);
    }

    SECTION("Tally-thresholds index creation") {
        test_index_creation("tally-thresholds", "tally_thresholds", 475329);
    }

    SECTION("Blocked index creation") {
        test_index_creation("blocked", "blocked", 654256);
    }

    SECTION("Blocked-thresholds index creation") {
        test_index_creation("blocked-thresholds", "blocked_thresholds", 711736);
    }

    SECTION("Large index creation") {
        test_index_creation("large", "large", 1305998);
    }
}