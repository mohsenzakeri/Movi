#pragma once

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/matchers/catch_matchers.hpp>

#include "movi_options.hpp"
#include "move_query.hpp"
#include "move_structure.hpp"
#include "movi_parser.hpp"
#include "utils.hpp"
#include "batch_loader.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// Test data directory
const std::string TEST_DATA_DIR = "tests_data";

// Helper function to create test data directory
void ensure_test_data_dir();