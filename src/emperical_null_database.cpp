#include "emperical_null_database.hpp"

size_t EmpNullDatabase::get_percentile_value() {
    return percentile_value;
}

void EmpNullDatabase::generate_stats(MoviOptions& movi_options, MoveStructure& mv_, std::string pattern_file) {
    bool random_jump = (USE_THRESHOLDS ? false : true);
    std::cerr << "random_jump: " << random_jump << "\n";
    gzFile fp;
    int l;
    kseq_t* seq = open_kseq(fp, pattern_file);
    while ((l = kseq_read(seq)) >= 0) { 
        std::string query_seq = seq->seq.s;
        std::reverse(query_seq.begin(), query_seq.end());
        MoveQuery mq = MoveQuery(query_seq);

        bool random_jump = false;
        mv_.query_pml(mq, random_jump);
        auto& pml_lens = mq.get_matching_lengths();
        pml_stats.insert(pml_stats.end(), pml_lens.begin(), pml_lens.end());
    }
    close_kseq(seq, fp);
}

void EmpNullDatabase::compute_stats() {
    // Determine vector size needed
    auto max_null_stat = std::max_element(pml_stats.begin(), pml_stats.end()); 
    max_stat_width = std::max(static_cast<int>(std::ceil(std::log2(*max_null_stat))), 1);

    std::cerr << "Maximum null statistic: " << *max_null_stat << "\n";
    std::cerr << "Number of bits used per null statistic: " << max_stat_width << "\n";

    // Initialize attributes
    num_values = pml_stats.size();
    null_stats.resize(num_values, 0);

    // Compute mean
    double sum_values = 0.0;
    for (size_t i = 0; i < pml_stats.size(); i++) {
        sum_values += pml_stats[i];
        null_stats[i] = pml_stats[i];
    }
    mean_null_stat = sum_values/pml_stats.size();

    // Find largest common value
    std::sort(pml_stats.begin(), pml_stats.end());

    size_t largest_val = 0, curr_val = pml_stats[0];
    size_t num_occs = 0;
    for (auto x: pml_stats) {
        if (x == curr_val)
            num_occs++;
        else {
            if (num_occs >= 5)
                largest_val = curr_val;
            curr_val = x;
            num_occs = 1;
        }
    }
    if (num_occs >= 5)  
        largest_val = curr_val;
    percentile_value = largest_val;

    std::cerr << "Largest common null statistic: " << largest_val << "\n";
    std::cerr << "Mean null statistic: " << mean_null_stat << "\n";
    std::cerr << "Percentile value: " << percentile_value << "\n";
}

void EmpNullDatabase::serialize(MoviOptions& movi_options) {
    std::string output_nulldb_name = movi_options.get_index_dir() + "/movi.pmlnulldb";
    std::ofstream output_nulldb(output_nulldb_name, std::ios::out | std::ios::binary);

    output_nulldb.write(reinterpret_cast<char*>(&num_values), sizeof(num_values));
    output_nulldb.write(reinterpret_cast<char*>(&mean_null_stat), sizeof(mean_null_stat));
    output_nulldb.write(reinterpret_cast<char*>(&percentile_value), sizeof(percentile_value));
    output_nulldb.write(reinterpret_cast<char*>(&null_stats[0]), num_values * sizeof(null_stats[0]));

    output_nulldb.close();
}

void EmpNullDatabase::deserialize(MoviOptions& movi_options) {
    std::string input_nulldb_name = movi_options.get_index_dir() + "/movi.pmlnulldb";
    std::ifstream input_nulldb(input_nulldb_name, std::ios::in | std::ios::binary);

    input_nulldb.read(reinterpret_cast<char*>(&num_values), sizeof(num_values));
    input_nulldb.read(reinterpret_cast<char*>(&mean_null_stat), sizeof(mean_null_stat));
    input_nulldb.read(reinterpret_cast<char*>(&percentile_value), sizeof(percentile_value));
    null_stats.resize(num_values);
    input_nulldb.read(reinterpret_cast<char*>(&null_stats[0]), num_values * sizeof(null_stats[0]));

    input_nulldb.close();
}