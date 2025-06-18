#ifndef MEM_FINDER_HPP
#define MEM_FINDER_HPP

struct MEMStatistics {
    uint64_t total_steps() {
        return ftab_backward_extensions + ftab_forward_extensions + ftab_bidirectional_extensions + 
                backward_extensions + forward_extensions + bidirectional_left_extensions + bidirectional_right_extensions;
    }

    uint64_t total_bidirectional_extensions() {
        return bidirectional_left_extensions + bidirectional_right_extensions;
    }
    
    void print() {
        std::cout << "\n\nNumber of MEMs found in the index:\t" << mems << "\n";
        std::cout << "Total length of patterns:             \t" << total_length << "\n";
        std::cout << "Sum of the counts of found MEMs:      \t" << total_counts << "\n\n\n";

        std::cout << "\n- - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - -\n";
        std::cout << "ftab_skips:               \t" << ftab_skips << "\n";
        std::cout << "ftab_backward_fail:       \t" << ftab_backward_fail << "\n";
        std::cout << "ftab_forward_fail:        \t" << ftab_forward_fail << "\n";
        std::cout << "ftab_bidirectional_fail:  \t" << ftab_bidirectional_fail << "\n";
        std::cout << "bidirectional_count_scans:\t" << bidirectional_count_scans << "\t("
                    << std::setprecision(4) << static_cast<double>(bidirectional_count_scans) / static_cast<double>(total_bidirectional_extensions()) << " avg.)\n";
        std::cout << "bml_skips:                \t" << bml_skips << "\n";
        std::cout << "bml_chars_skipped:        \t" << bml_chars_skipped << "\n\n\n";
        
        std::cout << "\n- - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - -\n";
        std::cout << "total_steps:                   \t" << total_steps() << "\n";
        std::cout << "ftab_backward_extensions:      \t" << ftab_backward_extensions << "\t"
                    << std::setprecision(4) << 100 * static_cast<double>(ftab_backward_extensions) / static_cast<double>(total_steps()) << "%\n";
        std::cout << "ftab_forward_extensions:       \t" << ftab_forward_extensions << "\t"
                    << std::setprecision(4) << 100 * static_cast<double>(ftab_forward_extensions) / static_cast<double>(total_steps()) << "%\n";
        std::cout << "ftab_bidirectional_extensions: \t" << ftab_bidirectional_extensions << "\t"
                    << std::setprecision(4) << 100 * static_cast<double>(ftab_bidirectional_extensions) / static_cast<double>(total_steps()) << "%\n";
        std::cout << "backward_extensions:           \t" << backward_extensions << "\t"
                    << std::setprecision(4) << 100 * static_cast<double>(backward_extensions) / static_cast<double>(total_steps()) << "%\n";
        std::cout << "forward_extensions:            \t" << forward_extensions << "\t"
                    << std::setprecision(4) << 100 * static_cast<double>(forward_extensions) / static_cast<double>(total_steps()) << "%\n";
        std::cout << "bidirectional_left_extensions: \t" << bidirectional_left_extensions << "\t"
                    << std::setprecision(4) << 100 * static_cast<double>(bidirectional_left_extensions) / static_cast<double>(total_steps()) << "%\n";
        std::cout << "bidirectional_right_extensions:\t" << bidirectional_right_extensions << "\t"
                    << std::setprecision(4) << 100 * static_cast<double>(bidirectional_right_extensions) / static_cast<double>(total_steps()) << "%\n";
        std::cout << "- - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - -\n\n";
    }
    uint64_t mems = 0;
    uint64_t total_length = 0;
    uint64_t total_counts = 0;
    uint64_t ftab_skips = 0;
    uint64_t ftab_backward_fail = 0;
    uint64_t ftab_forward_fail = 0;
    uint64_t ftab_bidirectional_fail = 0;
    uint64_t bml_skips = 0;
    uint64_t bml_chars_skipped = 0;
    uint64_t ftab_backward_extensions = 0;
    uint64_t ftab_forward_extensions = 0;
    uint64_t ftab_bidirectional_extensions = 0;
    uint64_t backward_extensions = 0;
    uint64_t forward_extensions = 0;
    uint64_t bidirectional_left_extensions = 0;
    uint64_t bidirectional_right_extensions = 0;
    uint64_t bidirectional_count_scans = 0;
};

#endif
