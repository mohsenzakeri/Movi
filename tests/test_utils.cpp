#include "test_utils.hpp"

// Helper function to create test data directory
void ensure_test_data_dir() {
    if (!std::filesystem::exists(TEST_DATA_DIR)) {
        std::filesystem::create_directory(TEST_DATA_DIR);
    }
}