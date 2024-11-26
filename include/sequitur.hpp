
#ifndef __SEQUITUR__
#define __SEQUITUR__

struct KmerStatistics {
    uint64_t total_kmers() {
        uint64_t negative_kmers = look_ahead_skipped + initialize_skipped + backward_search_failed + backward_search_empty +
                                  right_extension_failed;
        return positive_kmers + negative_kmers;
    }

    void print(bool count_query = false) {
        std::cout << "\n\nNumber of kmers found in the index:\t" << positive_kmers << "\n";
        if (count_query)
            std::cout << "Sum of the counts of found k-mers:\t" << total_counts << "\n\n\n";

        std::cout << "\n- - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - -\n";
        std::cout << "total_kmers:\t\t" <<  total_kmers() << "\n";
        std::cout << "positive_skipped:\t" << positive_skipped << "\t"
                    << std::setprecision(4) << 100 * static_cast<double>(positive_skipped) / static_cast<double>(total_kmers()) << "%\n";
        std::cout << "backward_search_failed:\t" << backward_search_failed << "\t"
                    << std::setprecision(4) << 100 * static_cast<double>(backward_search_failed) / static_cast<double>(total_kmers()) << "%\n";
        std::cout << "look_ahead_skipped:\t" << look_ahead_skipped << "\t"
                    << std::setprecision(4) << 100 * static_cast<double>(look_ahead_skipped) / static_cast<double>(total_kmers()) << "%\n";
        std::cout << "initialize_skipped:\t" << initialize_skipped << "\t"
                    << std::setprecision(2) << 100 * static_cast<double>(initialize_skipped) / static_cast<double>(total_kmers()) << "%\n";
        std::cout << "backward_search_empty:\t" << backward_search_empty << "\t"
                    << std::setprecision(2) << 100 * static_cast<double>(backward_search_empty) / static_cast<double>(total_kmers()) << "%\n";
        std::cout << "right_extension_failed:\t" << right_extension_failed << "\t"
                    << std::setprecision(2) << 100 * static_cast<double>(right_extension_failed) / static_cast<double>(total_kmers()) << "%\n";

        std::cout << "- - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - -\n\n";
    }
    uint64_t total_counts = 0;
    uint64_t positive_kmers = 0;
    uint64_t positive_skipped = 0;
    uint64_t look_ahead_skipped = 0;
    uint64_t initialize_skipped = 0;
    uint64_t backward_search_failed = 0;
    uint64_t backward_search_empty = 0;
    uint64_t right_extension_failed = 0;
};

#endif