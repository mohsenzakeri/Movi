#ifndef COMMONS_HPP
#define COMMONS_HPP

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdint>

#define LARGE_INDEX MODE == 0
#define CONSTANT_INDEX MODE == 1
#define SPLIT_INDEX MODE == 4
#define REGULAR_INDEX MODE == 3
#define REGULAR_THRESHOLDS_INDEX MODE == 6
#define BLOCKED_INDEX MODE == 2
#define BLOCKED_THRESHOLDS_INDEX MODE == 8
// Note: External interface uses "sampled" naming, internal code still uses in many places "tally"
#define TALLY_INDEX MODE == 5
#define TALLY_THRESHOLDS_INDEX MODE == 7

#define NO_SAMPPLED_ID MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4 or MODE == 8 or MODE == 3 or MODE == 6
#define USE_THRESHOLDS MODE == 0 or MODE == 1 or MODE == 4 or MODE == 6 or MODE == 7 or MODE == 8
#define NO_THRESHOLDS MODE == 2 or MODE == 3 or MODE == 5
#define THRESHOLDS_WITHOUT_NEXTS MODE == 0 or MODE == 4 or MODE == 8 or MODE == 7 or MODE == 6
#define USE_NEXT_POINTERS MODE == 1
#define SPLIT_THRESHOLDS_FALSE MODE == 0 or MODE == 1 or MODE == 4
#define SPLIT_THRESHOLDS_TRUE MODE == 8 or MODE == 7 or MODE == 6
#define SPLIT_MAX_RUN MODE == 3 or MODE == 6 or MODE == 2 or MODE == 8 or MODE == 5 or MODE == 7
#define SPLIT_ARRAY MODE == 1 or MODE == 4
#define NO_EXTRA_TABLE MODE == 0 or MODE == 1 or MODE == 4 or MODE == 3 or MODE == 6
#define REGULAR_MODES MODE == 3 or MODE == 6
#define BLOCKED_MODES MODE == 2 or MODE == 8
#define TALLY_MODES MODE == 5 or MODE == 7
#define MOVI1_STYLE MODE == 0 or MODE == 1 or MODE == 4

// Feature compatibility macros
#define SUPPORTS_SEPARATORS (MODE == 2 or MODE == 3 or MODE == 5 or MODE == 6 or MODE == 7 or MODE == 8)

#define END_CHARACTER 0
#define THRBYTES 5
#define MIN_MATCHING_LENGTH 3
#define NULL_READ_CHUNK 150
#define NUM_NULL_READS 800 // 150,000 = 150 bp * 1000 reads
#define NULL_READ_BOUND 1000

#define UNCLASSIFIED_THRESHOLD 0.4

// Utility function to get current timestamp
inline std::string get_timestamp() {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto tm = *std::localtime(&time_t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%H:%M:%S");
    return oss.str();
}

#define ERROR_MSG(msg) (std::string("\033[31m\n\n[ERROR] ") + get_timestamp() + " - " + msg + "\n\033[0m")
#define WARNING_MSG(msg) (std::cerr << "\033[33m\n\n[WARN] " << get_timestamp() << " - " << msg << "\n\033[0m")
#define DEBUG_MSG(msg) (std::cerr << "\033[36m[DEBUG] " << get_timestamp() << " - " << msg << "\n\033[0m")
#define INFO_MSG(msg) (std::cerr << "\033[37m[INFO] " << get_timestamp() << " - " << msg << "\n\033[0m")
#define PROGRESS_MSG(msg) (std::cerr << "\033[95m\n[PROGRESS] " << get_timestamp() << " - " << msg << "\n\033[0m")
#define QUERY_PROGRESS_MSG(msg) (std::cerr << "\033[95m[PROGRESS] " << get_timestamp() << " - " << msg << "\r\033[0m")
#define SUCCESS_MSG(msg) (std::cerr << "\033[92m\n[SUCCESS] " << get_timestamp() << " - " << msg << "\n\n\033[0m")

// Function for timing messages: formats elapsed time and passes to PROGRESS_MSG
// Usage: TIMING_MSG(elapsed, "loading the index")
// where elapsed is std::chrono::nanoseconds
inline void TIMING_MSG(const std::chrono::nanoseconds& elapsed_ns, const std::string& description) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3);
    oss << "Time measured for " << description << ": " << (elapsed_ns.count() * 1e-9) << " seconds.\n";
    INFO_MSG(oss.str());
}

// In case we want to call the timing message with a uint64_t instead of a std::chrono::nanoseconds
inline void TIMING_MSG(const uint64_t& elapsed_ns_int, const std::string& description) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3);
    oss << "Time measured for " << description << ": " << (elapsed_ns_int * 1e-9) << " seconds.\n";
    INFO_MSG(oss.str());
}

// Function to format numbers with commas (e.g., 1334 -> "1,334")
inline std::string format_number_with_commas(uint64_t number) {
    std::string num_str = std::to_string(number);
    std::string result;

    int count = 0;
    for (int i = num_str.length() - 1; i >= 0; i--) {
        if (count > 0 && count % 3 == 0) {
            result = "," + result;
        }
        result = num_str[i] + result;
        count++;
    }

    return result;
}

#define SEPARATOR '%'
#define SEPARATOR_INDEX 0

// Progress bar utility function
void print_progress_bar(uint64_t current, uint64_t total, const std::string& operation, uint64_t step = 0, uint64_t total_steps = 0);

std::string get_action(int argc, char* argv[]);

#endif