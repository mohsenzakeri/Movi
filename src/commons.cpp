#include "commons.hpp"

void print_progress_bar(uint64_t current, uint64_t total, const std::string& operation, uint64_t step, uint64_t total_steps) {
    const int bar_width = 50;
    float progress = (float)current / total;
    int pos = bar_width * progress;
    std::cerr << "\033[95m";
    if (current == 0) {
        if (step != 0) {
            std::cerr << "\n[" << step << "/" << total_steps << "] " << get_timestamp() << " - "<< operation << "..\n";
        } else {
            std::cerr << "\n[PROGRESS] " << get_timestamp() << " - " << operation << "..\n";
        }
    }
    std::cerr << "\r" << "[";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cerr << "█";
        else if (i == pos) std::cerr << "▌";
        else std::cerr << " ";
    }
    std::cerr << "] " << int(progress * 100.0) << "% (" << format_number_with_commas(current) << "/" << format_number_with_commas(total) << ")";
    std::cerr.flush();
}