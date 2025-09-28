#include "test_utils.hpp"

using namespace Catch::Matchers;

// Helper function to test PML computation
void test_pml_computation(const std::string& index_type, const std::string& index_suffix, 
                         const std::string& output_suffix, const std::string& extra_query_args = "") {
    ensure_test_data_dir();
    std::string fasta_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/ref.fasta";
    const std::string index_dir = TEST_DATA_DIR + "/index_" + index_suffix;
    
    // Build index
    std::string cmd = std::string(BINARY_DIR) + "/movi build --type " + index_type + " --index " + index_dir +
              " --fasta " + fasta_file + " --verify > /dev/null 2>&1";
    int exit_code = system(cmd.c_str());
    REQUIRE(exit_code == 0);

    // Run PML query
    std::string reads_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/reads.fasta";
    std::string query_output = TEST_DATA_DIR + "/" + output_suffix + ".pmls";
    std::string query_cmd = std::string(BINARY_DIR) + "/movi query --index " + index_dir +
              " --read " + reads_file + " --pml " + extra_query_args + " --stdout > " + query_output + " 2>/dev/null";
    int query_exit_code = system(query_cmd.c_str());
    REQUIRE(query_exit_code == 0);

    // Sort output using C locale for consistent sorting across environments
    std::string query_sort_cmd = "LC_ALL=C sort " + query_output + " > " + query_output + ".sorted";
    int sort_exit_code = system(query_sort_cmd.c_str());
    REQUIRE(sort_exit_code == 0);

    // Compare with expected results
    std::string pmls_source_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/reads.fasta.pmls.sorted";
    REQUIRE(std::filesystem::exists(query_output + ".sorted"));
    REQUIRE(std::filesystem::exists(pmls_source_file));
    
    std::string compare_cmd = "diff " + query_output + ".sorted " + pmls_source_file + " > /dev/null 2>&1";
    int compare_exit_code = system(compare_cmd.c_str());
    REQUIRE(compare_exit_code == 0);

    // Remove the outputs
    std::filesystem::remove_all(index_dir + "/index.movi");
    std::filesystem::remove_all(query_output + ".sorted ");
    std::filesystem::remove_all(query_output);
}

TEST_CASE("MoveStructure - PML query", "[move_structure_pml]") {
    SECTION("Regular thresholds PML computation") {
        test_pml_computation("regular-thresholds", "regular_thresholds", "regular_thresholds", "--no-prefetch -t1");
    }

    SECTION("Blocked-thresholds PML computation") {
        test_pml_computation("blocked-thresholds", "blocked_thresholds", "blocked_thresholds", "--no-prefetch -t1");
    }

    SECTION("Tally-thresholds PML computation") {
        test_pml_computation("tally-thresholds", "tally_thresholds", "tally_thresholds", "--no-prefetch -t1");
    }

    SECTION("Large PML computation") {
        test_pml_computation("large", "large", "large", "--no-prefetch -t1");
    }

    SECTION("Tally-thresholds PML computation (threaded)") {
        test_pml_computation("tally-thresholds", "tally_thresholds", "tally_thresholds", "--no-prefetch -t16");
    }
}

TEST_CASE("ReadProcessor - PML query", "[read_processor]") {
    SECTION("Regular thresholds PML computation with 16 strands and 1 thread") {
        test_pml_computation("regular-thresholds", "regular_thresholds", "regular_thresholds", "-s16 -t1");
    }

    SECTION("Regular thresholds PML computation with 16 strands and 16 thread") {
        test_pml_computation("regular-thresholds", "regular_thresholds", "regular_thresholds", "-s16 -t16");
    }
}

// Cleanup function
// TEST_CASE("Cleanup test data", "[cleanup]") {
//     SECTION("Remove test data directory") {
//         if (std::filesystem::exists(TEST_DATA_DIR)) {
//             std::filesystem::remove_all(TEST_DATA_DIR);
//         }
//         REQUIRE_FALSE(std::filesystem::exists(TEST_DATA_DIR));
//     }
// }