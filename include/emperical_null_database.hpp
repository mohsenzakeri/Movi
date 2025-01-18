#ifndef EMPERICAL_NULL_DATABASE_HPP
#define EMPERICAL_NULL_DATABASE_HPP

#include "move_structure.hpp"
#include "movi_options.hpp"
#include "utils.hpp"

// Borrowed from spumoni written by Omar Ahmed: https://github.com/oma219/spumoni/tree/main
class EmpNullDatabase {
private:
    std::vector<size_t> pml_stats;
    uint32_t max_stat_width;
    size_t num_values;
    std::vector<size_t> null_stats;
    double mean_null_stat;
    size_t percentile_value;

public:
    size_t get_percentile_value();

    void generate_stats(MoviOptions& movi_options, MoveStructure& mv_, std::string pattern_file);

    void compute_stats();

    void serialize(MoviOptions& movi_options);

    void deserialize(MoviOptions& movi_options);
};

#endif