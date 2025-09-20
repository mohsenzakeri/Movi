#!/bin/bash

# Movi Test Runner Script
# This script builds and runs all tests for the Movi project

set -e  # Exit on any error

# Parse command line arguments
JOBS=1  # Default to 1 job

while [[ $# -gt 0 ]]; do
    case $1 in
        -j|--jobs)
            JOBS="$2"
            shift 2
            ;;
        -j*)
            JOBS="${1#-j}"
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  -j, --jobs N     Number of parallel jobs for compilation (default: 1)"
            echo "  -h, --help       Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# Validate jobs parameter
if ! [[ "$JOBS" =~ ^[0-9]+$ ]] || [ "$JOBS" -lt 1 ]; then
    echo "Error: Jobs must be a positive integer, got: $JOBS"
    exit 1
fi

echo "=========================================="
echo "Movi Test Runner"
echo "=========================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if we're in the right directory
if [ ! -f "CMakeLists.txt" ]; then
    print_error "CMakeLists.txt not found. Please run this script from the Movi project root directory."
    exit 1
fi

# Create build directory if it doesn't exist
BUILD_DIR="build_test"
if [ ! -d "$BUILD_DIR" ]; then
    print_status "Creating build directory: $BUILD_DIR"
    mkdir -p "$BUILD_DIR"
fi

# Change to build directory
cd "$BUILD_DIR"

# Configure CMake
print_status "Configuring CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_STANDARD=20 -DBUILD_TEST=1

# Build the project
print_status "Building project with $JOBS parallel job(s)..."
make -j$JOBS

# Check if test executables were built
TEST_EXECUTABLES=("basics-tests" "build-tests" "pml-tests")
MISSING_TESTS=()

for test_exe in "${TEST_EXECUTABLES[@]}"; do
    if [ ! -f "$test_exe" ]; then
        MISSING_TESTS+=("$test_exe")
    fi
done

if [ ${#MISSING_TESTS[@]} -gt 0 ]; then
    print_warning "Some test executables were not built:"
    for test in "${MISSING_TESTS[@]}"; do
        echo "  - $test"
    done
fi

# Run tests
print_status "Running tests..."

# Function to run a test and report results
run_test() {
    local test_name="$1"
    local test_exe="$2"
    
    if [ -f "$test_exe" ]; then
        print_status "Running $test_name..."
        echo "----------------------------------------"
        
        # Run tests with verbose output showing individual test cases
        if ./"$test_exe" --verbosity high --reporter compact; then
            echo "----------------------------------------"
            print_success "$test_name passed"
            echo ""
            return 0
        else
            echo "----------------------------------------"
            print_error "$test_name failed"
            echo ""
            echo "Detailed failure information:"
            echo "Running $test_exe with full output..."
            ./"$test_exe" --verbosity high --reporter console
            echo ""
            return 1
        fi
    else
        print_warning "$test_name executable not found, skipping"
        return 1
    fi
}

# Track test results
FAILED_TESTS=0
TOTAL_TESTS=0

# Temporarily disable exit on error for test execution
set +e

# Run individual test suites
if ! run_test "Basics Tests" "basics-tests"; then
    ((FAILED_TESTS++))
fi
((TOTAL_TESTS++))

if ! run_test "Build Tests" "build-tests"; then
    ((FAILED_TESTS++))
fi
((TOTAL_TESTS++))

if ! run_test "PML Tests" "pml-tests"; then
    ((FAILED_TESTS++))
fi
((TOTAL_TESTS++))

# Re-enable exit on error
set -e

# Print test summary
echo ""
echo "========================================"
echo "TEST SUMMARY"
echo "========================================"
echo "Total test suites: $TOTAL_TESTS"
echo "Passed test suites: $((TOTAL_TESTS - FAILED_TESTS))"
echo "Failed test suites: $FAILED_TESTS"

if [ $FAILED_TESTS -eq 0 ]; then
    print_success "All test suites passed successfully!"
    echo ""
    exit 0
else
    print_error "$FAILED_TESTS out of $TOTAL_TESTS test suites failed"
    echo ""
    exit 1
fi