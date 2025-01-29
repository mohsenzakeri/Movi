#ifndef CLASSIFIER_HPP
#define CLASSIFIER_HPP

#include <iostream>

#include "utils.hpp"
#include "movi_options.hpp"
#include "move_structure.hpp"
#include "emperical_null_database.hpp"

class Classifier {
public:
    void generate_null_statistics(MoveStructure& mv_, MoviOptions& movi_options);
    size_t initialize_report_file(MoviOptions& movi_options);
    void close_report_file();
    void classify(std::string read_name, std::vector<uint16_t>& matching_lens, MoviOptions& movi_options);

private:
    MoviOptions* movi_options_;
    std::ofstream report_file;
    uint16_t max_value_thr;
};

#endif