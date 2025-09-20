# Movi Testing Framework

This document describes the comprehensive testing framework for Movi.

## Overview

The testing framework uses [Catch2](https://github.com/catchorg/Catch2) as the testing library and provides comprehensive coverage of the Movi codebase including:


## Running Tests

### Quick Start

```bash
# Run all tests with the provided script
./run_tests.sh
```

### Manual Testing

```bash
# Create build directory
mkdir -p build_test
cd build_test

# Configure and build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j

# Run all tests using the test runner script
../run_tests.sh

# Or run individual test executables
./basics-tests
./build-tests
./pml-tests

```

### Running Specific Tests

```bash
# Run individual test suites
./basics-tests      # MoviOptions, Command line parsing, MoveInterval, DocSet, MoveQuery
./build-tests       # MoveStructure index building tests
./pml-tests         # PML query tests and cleanup

# Run with specific tags
./basics-tests "[movi_options]"
./basics-tests "[parser]"
./basics-tests "[move_interval]"
./basics-tests "[docset]"
./basics-tests "[move_query]"
./build-tests "[move_structure_index_building]"
./pml-tests "[move_structure_pml]"
./pml-tests "[read_processor]"
```

### Test Data

The testing framework relies on the following files found in the `tests_data/` directory:

- Reference sequence: ref.fasta
- Read sequence: reads.fasta
- Expected PML values: reads.fasta.pmls.sorted

## Adding New Tests

### 1. Create Test File

```cpp
#include "test_utils.hpp"

using namespace Catch::Matchers;

TEST_CASE("Your Test Description", "[tag]") {
    SECTION("Specific test case") {
        // Your test code here
        REQUIRE(condition);
    }
}
```

### 2. Update CMakeLists.txt

Add your test file to `tests/CMakeLists.txt`:

```cmake
add_executable(your-tests your_test_file.cpp test_utils.cpp)
target_link_libraries(your-tests PRIVATE movi_lib_test Catch2::Catch2WithMain)
target_include_directories(your-tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(your-tests PRIVATE -Wno-deprecated-declarations)
target_compile_definitions(your-tests PRIVATE BINARY_DIR="${CMAKE_BINARY_DIR}")
```

### 3. Update Test Runner

Add your new test executable to `run_tests.sh`:

```bash
if ! run_test "Your Tests" "your-tests"; then
    ((FAILED_TESTS++))
fi
((TOTAL_TESTS++))
```