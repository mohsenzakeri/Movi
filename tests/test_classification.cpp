#include "test_utils.hpp"

using namespace Catch::Matchers;

// Helper function to test PML filtering
void test_pml_filtering(const std::string& index_type, const std::string& index_suffix,
                        const std::string& output_suffix, const std::string& extra_query_args = "") {
    ensure_test_data_dir();
    std::string fasta_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/ref.fasta";
    const std::string index_dir = TEST_DATA_DIR + "/index_" + index_suffix;

    // Build index
    std::string cmd = std::string(BINARY_DIR) + "/movi build --type " + index_type + " --index " + index_dir +
                      " --fasta " + fasta_file + " --verify > /dev/null 2>&1";
    int exit_code = system(cmd.c_str());
    REQUIRE(exit_code == 0);

    // Run PML query with classification
    std::string reads_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/classification.fasta";
    std::string reads_filtered = TEST_DATA_DIR + "/" + output_suffix + ".pmls.filtered";
    std::string query_cmd = std::string(BINARY_DIR) + "/movi query --index " + index_dir +
                            " --read " + reads_file + " --pml --filter " + extra_query_args + " --stdout > " + reads_filtered + " 2>/dev/null";
    int query_exit_code = system(query_cmd.c_str());
    REQUIRE(query_exit_code == 0);

    // Sort output using C locale for consistent sorting across environments
    std::string query_sort_cmd = "LC_ALL=C sort " + reads_filtered + " > " + reads_filtered + ".sorted";
    int sort_exit_code = system(query_sort_cmd.c_str());
    REQUIRE(sort_exit_code == 0);

    // Compare with expected results
    std::string filtered_source_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/classification.fasta.pmls.filtered.sorted";
    REQUIRE(std::filesystem::exists(reads_filtered + ".sorted"));
    REQUIRE(std::filesystem::exists(filtered_source_file));

    std::string compare_cmd = "diff " + reads_filtered + ".sorted " + filtered_source_file + " > /dev/null 2>&1";
    int compare_exit_code = system(compare_cmd.c_str());

    if (compare_exit_code != 0) {
        std::cerr << "----------------------------------------\n";
        std::cerr << query_cmd << "\n";
        std::cerr << compare_cmd << "\n";
        std::cerr << "----------------------------------------\n";
    }

    REQUIRE(compare_exit_code == 0);

    // Remove the outputs
    std::filesystem::remove_all(index_dir + "/index.movi");
    std::filesystem::remove_all(reads_filtered + ".sorted ");
    std::filesystem::remove_all(reads_filtered);
}

// Helper function to test PML filtering (invert)
void test_pml_filtering_invert(const std::string& index_type, const std::string& index_suffix,
                               const std::string& output_suffix, const std::string& extra_query_args = "") {
    ensure_test_data_dir();
    std::string fasta_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/ref.fasta";
    const std::string index_dir = TEST_DATA_DIR + "/index_" + index_suffix;

    // Build index
    std::string cmd = std::string(BINARY_DIR) + "/movi build --type " + index_type + " --index " + index_dir +
                      " --fasta " + fasta_file + " --verify > /dev/null 2>&1";
    int exit_code = system(cmd.c_str());
    REQUIRE(exit_code == 0);

    // Run PML query with classification
    std::string reads_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/sample.fasta";
    std::string reads_filtered = TEST_DATA_DIR + "/" + output_suffix + ".pmls.filtered_notfound";
    std::string query_cmd = std::string(BINARY_DIR) + "/movi query --index " + index_dir +
                            " --read " + reads_file + " --pml --filter --invert " + extra_query_args + " --stdout > " + reads_filtered + " 2>/dev/null";
    int query_exit_code = system(query_cmd.c_str());
    REQUIRE(query_exit_code == 0);

    // Sort output using C locale for consistent sorting across environments
    std::string query_sort_cmd = "LC_ALL=C sort " + reads_filtered + " > " + reads_filtered + ".sorted";
    int sort_exit_code = system(query_sort_cmd.c_str());
    REQUIRE(sort_exit_code == 0);

    // Compare with expected results
    std::string filtered_source_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/sample.fasta.pmls.filtered_notfound.sorted";
    REQUIRE(std::filesystem::exists(reads_filtered + ".sorted"));
    REQUIRE(std::filesystem::exists(filtered_source_file));

    std::string compare_cmd = "diff " + reads_filtered + ".sorted " + filtered_source_file + " > /dev/null 2>&1";
    int compare_exit_code = system(compare_cmd.c_str());

    if (compare_exit_code != 0) {
        std::cerr << "----------------------------------------\n";
        std::cerr << query_cmd << "\n";
        std::cerr << compare_cmd << "\n";
        std::cerr << "----------------------------------------\n";
    }

    REQUIRE(compare_exit_code == 0);

    // Remove the outputs
    std::filesystem::remove_all(index_dir + "/index.movi");
    std::filesystem::remove_all(reads_filtered + ".sorted ");
    std::filesystem::remove_all(reads_filtered);
}


// Helper function to test PML classification
void test_classification_computation(const std::string& index_type, const std::string& index_suffix,
                         const std::string& output_suffix, const std::string& extra_query_args = "") {
    ensure_test_data_dir();
    std::string fasta_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/ref.fasta";
    const std::string index_dir = TEST_DATA_DIR + "/index_" + index_suffix;
    
    // Build index
    std::string cmd = std::string(BINARY_DIR) + "/movi build --type " + index_type + " --index " + index_dir +
              " --fasta " + fasta_file + " --keep --verify > /dev/null 2>&1";
    int exit_code = system(cmd.c_str());
    REQUIRE(exit_code == 0);

    // Run PML query with classification
    std::string reads_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/classification.fasta";
    std::string report_output = TEST_DATA_DIR + "/" + output_suffix + ".pmls.report";
    std::string query_cmd = std::string(BINARY_DIR) + "/movi query --index " + index_dir +
              " --read " + reads_file + " --pml --classify " + extra_query_args + " --stdout > " + report_output + " 2>/dev/null";
    int query_exit_code = system(query_cmd.c_str());
    REQUIRE(query_exit_code == 0);

    // Sort output using C locale for consistent sorting across environments
    std::string query_sort_cmd = "LC_ALL=C sort " + report_output + " > " + report_output + ".sorted";
    int sort_exit_code = system(query_sort_cmd.c_str());
    REQUIRE(sort_exit_code == 0);

    REQUIRE(std::filesystem::exists(report_output + ".sorted"));

    // Compare with expected results
    std::string summary_report_source_file = std::string(BINARY_DIR) + "/../" + TEST_DATA_DIR + "/classification.fasta.pmls.report.summary.sorted";
    REQUIRE(std::filesystem::exists(summary_report_source_file));

    std::string summary_cmd = "cat " + report_output + ".sorted | awk '{print $1,$2}' > " + report_output + ".summary.sorted";
    int summary_exit_code = system(summary_cmd.c_str());
    REQUIRE(summary_exit_code == 0);
    
    std::string compare_cmd = "diff " + report_output + ".summary.sorted " + summary_report_source_file + " > /dev/null 2>&1";
    int compare_exit_code = system(compare_cmd.c_str());

    if (compare_exit_code != 0) {
        std::cerr << "----------------------------------------\n";
        std::cerr << query_cmd << "\n";
        std::cerr << compare_cmd << "\n";
        std::cerr << "----------------------------------------\n";
    }

    REQUIRE(compare_exit_code == 0);

    // Remove the outputs
    std::filesystem::remove_all(index_dir + "/index.movi");
    std::filesystem::remove_all(report_output + ".sorted ");
    std::filesystem::remove_all(report_output + ".summary.sorted");
    std::filesystem::remove_all(report_output);
}

TEST_CASE("MoveStructure - PML query with classification", "[move_structure_pml_classification]") {
    SECTION("Regular thresholds PML computation") {
        test_classification_computation("regular-thresholds", "regular_thresholds", "regular_thresholds", "--no-prefetch -t1");
    }

    SECTION("Blocked-thresholds PML classification") {
        test_classification_computation("blocked-thresholds", "blocked_thresholds", "blocked_thresholds", "--no-prefetch -t1");
    }

    SECTION("Tally-thresholds PML classification") {
        test_classification_computation("tally-thresholds", "tally_thresholds", "tally_thresholds", "--no-prefetch -t1");
    }

    SECTION("Large PML classification") {
        test_classification_computation("large", "large", "large", "--no-prefetch -t1");
    }

    SECTION("Tally-thresholds PML classification (threaded)") {
        test_classification_computation("tally-thresholds", "tally_thresholds", "tally_thresholds", "--no-prefetch -t16");
    }

    SECTION("Regular thresholds PML classification with 16 strands and 1 thread") {
        test_classification_computation("regular-thresholds", "regular_thresholds", "regular_thresholds", "-s16 -t1");
    }

    SECTION("Regular thresholds PML classification with 16 strands and 16 thread") {
        test_classification_computation("regular-thresholds", "regular_thresholds", "regular_thresholds", "-s16 -t16");
    }
}

TEST_CASE("ReadProcessor - PML filtering", "[read_processor]") {
    SECTION("Regular thresholds PML filtering with 16 strands and 2 thread") {
        test_pml_filtering("regular-thresholds", "regular_thresholds", "regular_thresholds", " -s16 -t2");
    }

    SECTION("Tally thresholds PML filtering with no-prefetch and 16 thread") {
        test_pml_filtering("tally-thresholds", "tally_thresholds", "tally_thresholds", " --no-prefetch -t1");
    }

    SECTION("Large PML invert filtering with 16 strands and 1 thread") {
        test_pml_filtering_invert("large", "large", "large", " -s16 -t1");
    }

    // SECTION("Blocked thresholds PML invert filtering with no-prefetch and 1 thread") {
    //     test_pml_filtering_invert("blocked-thresholds", "blocked_thresholds", "blocked_thresholds", " --no-prefetch -t1");
    // }

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