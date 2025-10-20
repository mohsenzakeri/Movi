#include "test_utils.hpp"

using namespace Catch::Matchers;

// Helper function to test index creation
void test_index_creation(const std::string& index_type, const std::string& index_suffix, size_t expected_size, std::string flags = "") {
    ensure_test_data_dir();
    std::string fasta_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/ref.fasta";
    REQUIRE(std::filesystem::exists(fasta_file));
    
    const std::string index_dir = TEST_DATA_DIR + "/index_" + index_suffix;
    std::string cmd = std::string(BINARY_DIR) + "/movi build --type " + index_type + " --index " + index_dir +
                                 " --fasta " + fasta_file + " " + flags + " --verify > /dev/null 2>&1";
    int exit_code = system(cmd.c_str());

    // Check that command executed successfully
    REQUIRE(exit_code == 0);

    // Check that the index file exists
    const std::string index_file = index_dir + "/index.movi";
    REQUIRE(std::filesystem::exists(index_file));
    
    // Check that the index file has the expected size
    auto file_size = std::filesystem::file_size(index_file);
    REQUIRE(file_size == expected_size);

    // Remove the index file
    std::filesystem::remove_all(index_file);
}

TEST_CASE("MoveStructure - index building", "[move_structure_index_building]") {
    SECTION("Regular index creation") {
        test_index_creation("regular", "regular", 871479);
    }

    SECTION("Regular-thresholds index creation") {
        test_index_creation("regular-thresholds", "regular_thresholds", 948119);
    }

    SECTION("Sampled index creation") {
        test_index_creation("sampled", "sampled", 437006);
    }

    SECTION("Sampled-thresholds index creation") {
        test_index_creation("sampled-thresholds", "sampled_thresholds", 475326);
    }

    SECTION("Blocked index creation") {
        test_index_creation("blocked", "blocked", 654253);
    }

    SECTION("Blocked-thresholds index creation") {
        test_index_creation("blocked-thresholds", "blocked_thresholds", 711733);
    }

    SECTION("Large index creation") {
        test_index_creation("large", "large", 1305995);
    }

    SECTION("Sampled-thresholds index creation (preprocessed)") {
        test_index_creation("sampled-thresholds", "sampled_thresholds", 475326, "--preprocessed");
    }

    SECTION("Sampled-thresholds index creation (non-preprocessed)") {
        test_index_creation("large", "large", 1305995, "--non-preprocessed");
    }

    SECTION("Sampled-thresholds index creation (non-preprocessed and preprocessed)") {
        test_index_creation("sampled-thresholds", "sampled_thresholds", 475326, "--non-preprocessed --preprocessed");
    }
}

TEST_CASE("MoveStructure - index building with separators", "[move_structure_index_building_with_separators]") {
    SECTION("Regular index creation") {
        test_index_creation("regular", "regular", 871496, "--separators");
    }

    SECTION("Regular-thresholds index creation") {
        test_index_creation("regular-thresholds", "regular_thresholds", 948232, "--separators");
    }

    SECTION("Sampled index creation") {
        test_index_creation("sampled", "sampled", 464203, "--separators");
    }

    SECTION("Sampled-thresholds index creation") {
        test_index_creation("sampled-thresholds", "sampled_thresholds", 505009, "--separators");
    }

    SECTION("Blocked index creation") {
        test_index_creation("blocked", "blocked", 654280, "--separators");
    }

    SECTION("Blocked-thresholds index creation") {
        test_index_creation("blocked-thresholds", "blocked_thresholds", 711854, "--separators");
    }
}