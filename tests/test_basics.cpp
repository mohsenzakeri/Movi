#include "test_utils.hpp"

using namespace Catch::Matchers;

// Test suite for MoviOptions class
TEST_CASE("MoviOptions", "[movi_options]") {
    MoviOptions options;
    
    SECTION("Default values") {
        REQUIRE(options.get_k() == 31);
        REQUIRE(options.get_ftab_k() == 0);
        REQUIRE(options.is_pml());
        REQUIRE_FALSE(options.is_zml());
        REQUIRE_FALSE(options.is_count());
        REQUIRE_FALSE(options.is_kmer());
        REQUIRE_FALSE(options.is_verbose());
        REQUIRE_FALSE(options.is_debug());
        REQUIRE_FALSE(options.is_stdout());
        REQUIRE_FALSE(options.is_no_output());
        REQUIRE_FALSE(options.is_small_pml_lens());
        REQUIRE_FALSE(options.is_large_pml_lens());
        REQUIRE_FALSE(options.is_output_ids());
        REQUIRE_FALSE(options.is_logs());
        REQUIRE_FALSE(options.is_random_repositioning());
        REQUIRE_FALSE(options.is_get_sa_entries());
        REQUIRE_FALSE(options.is_pvalue_scoring());
        REQUIRE_FALSE(options.is_multi_classify());
        REQUIRE_FALSE(options.is_generate_null_reads());
        REQUIRE_FALSE(options.is_filter());
        REQUIRE_FALSE(options.is_early_stop());
        REQUIRE_FALSE(options.is_report_colors());
        REQUIRE_FALSE(options.is_report_color_ids());
        REQUIRE_FALSE(options.is_report_all());
        REQUIRE(options.ignore_illegal_chars_status() == 0);
        REQUIRE(options.get_strands() == 16);
        REQUIRE_FALSE(options.is_full_color());
        REQUIRE_FALSE(options.is_compressed());
        REQUIRE_FALSE(options.is_freq_compressed());
        REQUIRE_FALSE(options.is_tree_compressed());
        REQUIRE_FALSE(options.is_color_move_rows());
        REQUIRE_FALSE(options.is_flat_color_vectors());
        REQUIRE_FALSE(options.is_color());
        REQUIRE_FALSE(options.is_doc_sets_vector_of_vectors());
        REQUIRE(options.get_min_match_len() == 1);
        REQUIRE_FALSE(options.is_pvalue_scoring());
        REQUIRE(options.get_threads() == 1);
        REQUIRE(options.get_tally_checkpoints() == 20);
        REQUIRE(options.get_SA_sample_rate() == 100);
        REQUIRE_FALSE(options.is_multi_ftab());
        REQUIRE_FALSE(options.is_verify());
        REQUIRE_FALSE(options.is_classify());
        REQUIRE_FALSE(options.is_filter());
        REQUIRE_FALSE(options.is_multi_classify());
        REQUIRE_FALSE(options.is_early_stop());
        REQUIRE_FALSE(options.is_report_colors());
    }
    
    SECTION("Setters and getters") {
        options.set_command("build");
        REQUIRE(options.get_command() == "build");
        
        options.set_ref_file("test.fa");
        REQUIRE(options.get_ref_file() == "test.fa");
        
        options.set_read_file("reads.fq");
        REQUIRE(options.get_read_file() == "reads.fq");
        
        options.set_index_dir("index");
        REQUIRE(options.get_index_dir() == "index");
        
        options.set_k(31);
        REQUIRE(options.get_k() == 31);
        
        options.set_ftab_k(5);
        REQUIRE(options.get_ftab_k() == 5);
        
        options.set_threads(4);
        REQUIRE(options.get_threads() == 4);

        options.set_verbose(true);
        REQUIRE(options.is_verbose());
        
        options.set_debug(true);
        REQUIRE(options.is_debug());
        
        options.set_stdout(true);
        REQUIRE(options.is_stdout());
        
        options.set_no_output(true);
        REQUIRE(options.is_no_output());

    }
    
    SECTION("Query type flags") {
        options.set_pml();
        REQUIRE(options.is_pml());
        REQUIRE_FALSE(options.is_zml());
        REQUIRE_FALSE(options.is_count());
        REQUIRE_FALSE(options.is_kmer());
        
        options.set_zml();
        REQUIRE(options.is_zml());
        REQUIRE_FALSE(options.is_pml());
        REQUIRE_FALSE(options.is_count());
        REQUIRE_FALSE(options.is_kmer());
        
        options.set_count();
        REQUIRE(options.is_count());
        REQUIRE_FALSE(options.is_pml());
        REQUIRE_FALSE(options.is_zml());
        REQUIRE_FALSE(options.is_kmer());
        
        options.set_kmer();
        REQUIRE(options.is_kmer());
        REQUIRE_FALSE(options.is_pml());
        REQUIRE_FALSE(options.is_zml());
        REQUIRE_FALSE(options.is_count());
    }
}

// Test suite for command line parsing
TEST_CASE("Command line parsing", "[parser]") {
    SECTION("Build command parsing") {
        const char* argv[] = {"movi", "build", "-i", "test_index", "-f", "test.fa"};
        int argc = 6;
        
        MoviOptions options;
        bool result = parse_command(argc, const_cast<char**>(argv), options, true);
        
        REQUIRE(result == true);
        REQUIRE(options.get_command() == "build");
        REQUIRE(options.get_index_dir() == "test_index");
        REQUIRE(options.get_ref_file() == "test.fa");
    }
    
    SECTION("Query command parsing") {
        const char* argv[] = {"movi", "query", "-i", "test_index", "-r", "reads.fq", "--pml"};
        int argc = 7;
        
        MoviOptions options;
        bool result = parse_command(argc, const_cast<char**>(argv), options, true);
        
        REQUIRE(result == true);
        REQUIRE(options.get_command() == "query");
        REQUIRE(options.get_index_dir() == "test_index");
        REQUIRE(options.get_read_file() == "reads.fq");
        REQUIRE_FALSE(options.is_verbose());
        REQUIRE(options.is_pml());
        REQUIRE_FALSE(options.is_zml());
        REQUIRE_FALSE(options.is_count());
        REQUIRE_FALSE(options.is_kmer());
        REQUIRE_FALSE(options.is_kmer_count());
        REQUIRE_FALSE(options.is_reverse());
        REQUIRE_FALSE(options.is_multi_ftab());
        REQUIRE_FALSE(options.is_classify());
    }
    
    SECTION("Query with multiple options") {
        const char* argv[] = {"movi", "query", "-i", "test_index", "-r", "reads.fq", 
                             "--zml", "--verbose", "--threads", "4"};
        int argc = 10;
        
        MoviOptions options;
        bool result = parse_command(argc, const_cast<char**>(argv), options, true);
        
        REQUIRE(result == true);
        REQUIRE(options.get_command() == "query");
        REQUIRE(options.is_verbose());
        REQUIRE(options.get_threads() == 4);
        REQUIRE(options.is_zml());
        REQUIRE_FALSE(options.is_pml());
        REQUIRE_FALSE(options.is_count());
        REQUIRE_FALSE(options.is_kmer());
        REQUIRE_FALSE(options.is_kmer_count());
        REQUIRE_FALSE(options.is_reverse());
        REQUIRE_FALSE(options.is_multi_ftab());
        REQUIRE_FALSE(options.is_classify());
    }
    
    SECTION("Invalid command parsing") {
        const char* argv[] = {"movi", "invalid_command"};
        int argc = 2;
        
        MoviOptions options;
        bool result = parse_command(argc, const_cast<char**>(argv), options, true);
        
        REQUIRE(result == false);
    }
    
    SECTION("Missing required arguments") {
        const char* argv[] = {"movi", "build", "-i", "test_index"};
        int argc = 4;
        
        MoviOptions options;
        bool result = parse_command(argc, const_cast<char**>(argv), options, true);
        
        REQUIRE(result == false);
    }
}

// Test suite for MoveInterval
TEST_CASE("MoveInterval", "[move_interval]") {
    SECTION("Default constructor") {
        MoveInterval interval;
        // Default constructor doesn't initialize values, so we just test that it exists
        // and can be assigned to
        interval.run_start = 0;
        interval.offset_start = 0;
        interval.run_end = 1;
        interval.offset_end = 0;
        REQUIRE(interval.run_start == 0);
        REQUIRE(interval.offset_start == 0);
        REQUIRE(interval.run_end == 1);
        REQUIRE(interval.offset_end == 0);
    }
    
    SECTION("Parameterized constructor") {
        MoveInterval interval(10, 5, 20, 15);
        REQUIRE(interval.run_start == 10);
        REQUIRE(interval.offset_start == 5);
        REQUIRE(interval.run_end == 20);
        REQUIRE(interval.offset_end == 15);
    }
    
    SECTION("Assignment operator") {
        MoveInterval interval1(10, 5, 20, 15);
        MoveInterval interval2;
        
        interval2 = interval1;
        
        REQUIRE(interval2.run_start == 10);
        REQUIRE(interval2.offset_start == 5);
        REQUIRE(interval2.run_end == 20);
        REQUIRE(interval2.offset_end == 15);
    }
    
    SECTION("Make empty") {
        MoveInterval interval(10, 5, 20, 15);
        interval.make_empty();
        
        REQUIRE(interval.run_start == 1);
        REQUIRE(interval.offset_start == 0);
        REQUIRE(interval.run_end == 0);
        REQUIRE(interval.offset_end == 0);
        REQUIRE(interval.is_empty());
    }
}

// Test suite for DocSet (Colors)
TEST_CASE("DocSet", "[docset]") {
    SECTION("Default constructor") {
        DocSet docset;
        REQUIRE(docset.docs.empty());
    }
    
    SECTION("Constructor with vector") {
        std::vector<uint16_t> docs = {1, 5, 10};
        DocSet docset(docs);
        
        REQUIRE(docset.docs.size() == 3);
        REQUIRE(docset.docs[0] == 1);
        REQUIRE(docset.docs[1] == 5);
        REQUIRE(docset.docs[2] == 10);
    }
    
    SECTION("Move constructor") {
        std::vector<uint16_t> docs = {1, 5, 10};
        DocSet docset(std::move(docs));
        
        REQUIRE(docset.docs.size() == 3);
        REQUIRE(docset.docs[0] == 1);
        REQUIRE(docset.docs[1] == 5);
        REQUIRE(docset.docs[2] == 10);
    }
    
    SECTION("Equality operator") {
        std::vector<uint16_t> docs1 = {1, 5, 10};
        std::vector<uint16_t> docs2 = {1, 5, 10};
        std::vector<uint16_t> docs3 = {1, 5, 11};
        
        DocSet docset1(docs1);
        DocSet docset2(docs2);
        DocSet docset3(docs3);
        
        REQUIRE(docset1 == docset2);
        REQUIRE_FALSE(docset1 == docset3);
    }

    SECTION("DocSet with many documents") {
        std::vector<uint16_t> docs;
        for (int i = 0; i < 1000; ++i) {
            docs.push_back(i);
        }
        
        DocSet docset(docs);
        REQUIRE(docset.docs.size() == 1000);
        REQUIRE(docset.docs[0] == 0);
        REQUIRE(docset.docs[999] == 999);
    }
}

// Test suite for MoveQuery advanced functionality
TEST_CASE("MoveQuery", "[move_query]") {
    SECTION("Large matching lengths") {
        MoveQuery query("ATCGATCGATCG");
        
        // Test with large values
        query.add_ml(65535, false);
        query.add_ml(65536, false);
        
        auto& lengths = query.get_matching_lengths();
        REQUIRE(lengths.size() == 2);
        REQUIRE(lengths[0] == 65535);
        REQUIRE(lengths[1] == 65535);
    }
    
    SECTION("Matching lengths with stdout output") {
        MoveQuery query("ATCG");
        
        query.add_ml(123, true);
        query.add_ml(456, true);
        
        auto& lengths = query.get_matching_lengths();
        REQUIRE(lengths.size() == 2);
        REQUIRE(lengths[0] == 123);
        REQUIRE(lengths[1] == 456);
        
        // Check that stdout string contains the lengths
        std::string output = query.get_matching_lengths_string();
        REQUIRE_THAT(output, Catch::Matchers::ContainsSubstring("321")); // 123 reversed
        REQUIRE_THAT(output, Catch::Matchers::ContainsSubstring("654")); // 456 reversed
    }
    
    SECTION("Multiple kmer additions") {
        MoveQuery query("ATCGATCG");
        
        query.add_kmer(0, 10);
        query.add_kmer(1, 20);
        query.add_kmer(2, 30);
        
        REQUIRE(query.found_kmer_count == 60);
        
        std::string output = query.get_matching_lengths_string();
        REQUIRE_THAT(output, Catch::Matchers::ContainsSubstring("0:10"));
        REQUIRE_THAT(output, Catch::Matchers::ContainsSubstring("1:20"));
        REQUIRE_THAT(output, Catch::Matchers::ContainsSubstring("2:30"));
    }

    SECTION("Adding colors") {
        MoveQuery query("ATCG");
        
        query.add_color(1);
        query.add_color(5);
        query.add_color(10);
        
        auto& colors = query.get_matching_colors();
        REQUIRE(colors.size() == 3);
        REQUIRE(colors[0] == 1);
        REQUIRE(colors[1] == 5);
        REQUIRE(colors[2] == 10);
    }

    SECTION("Adding SA entries") {
        MoveQuery query("ATCG");
        
        query.add_sa_entries(100);
        query.add_sa_entries(200);
        
        auto& sa_entries = query.get_sa_entries();
        REQUIRE(sa_entries.size() == 2);
        REQUIRE(sa_entries[0] == 100);
        REQUIRE(sa_entries[1] == 200);
    }

    SECTION("Performance metrics accumulation") {
        MoveQuery query("ATCG");
        
        // Add multiple metrics
        for (int i = 0; i < 100; ++i) {
            query.add_scan(i);
            query.add_fastforward(i * 2);
            query.add_cost(std::chrono::nanoseconds(i * 1000));
        }
        
        auto& scans = query.get_scans();
        auto& fastforwards = query.get_fastforwards();
        auto& costs = query.get_costs();
        
        REQUIRE(scans.size() == 100);
        REQUIRE(fastforwards.size() == 100);
        REQUIRE(costs.size() == 100);
        
        // Check some specific values
        REQUIRE(scans[0] == 0);
        REQUIRE(scans[99] == 99);
        REQUIRE(fastforwards[0] == 0);
        REQUIRE(fastforwards[99] == 198);
        REQUIRE(costs[0].count() == 0);
        REQUIRE(costs[99].count() == 99000);
    }

    SECTION("MoveQuery with many matching lengths") {
        MoveQuery query("ATCGATCGATCG");
        
        // Add many matching lengths
        for (int i = 0; i < 10000; ++i) {
            query.add_ml(i % 100, false);
        }
        
        auto& lengths = query.get_matching_lengths();
        REQUIRE(lengths.size() == 10000);
        
        // Verify some specific values
        REQUIRE(lengths[0] == 0);
        REQUIRE(lengths[100] == 0);
        REQUIRE(lengths[101] == 1);
    }
}