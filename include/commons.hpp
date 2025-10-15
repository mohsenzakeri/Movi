#ifndef COMMONS_HPP
#define COMMONS_HPP

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

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

#endif