/*
 * This file is based on:
 * https://github.com/oma219/spumoni/blob/main/include/emp_null_database.hpp
 *
 * Original author: Omar Ahmed
 * License: GPL-3.0
 * Minor modifications made by Mohsen Zakeri
 */

#ifndef EMPERICAL_NULL_DATABASE_HPP
#define EMPERICAL_NULL_DATABASE_HPP

#include "move_structure.hpp"
#include "movi_options.hpp"
#include "utils.hpp"

class EmpNullDatabase {
private:
    std::vector<size_t> ml_stats;
    uint32_t max_stat_width;
    size_t num_values;
    std::vector<size_t> null_stats;
    double mean_null_stat;
    size_t percentile_value;

public:
    size_t get_percentile_value();

    void generate_stats(MoviOptions& movi_options, MoveStructure& mv_, std::string pattern_file);

    void compute_stats(MoviOptions& movi_options);

    void serialize(MoviOptions& movi_options);

    void deserialize(MoviOptions& movi_options);
};

#endif