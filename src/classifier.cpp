#include "classifier.hpp"

void Classifier::generate_null_statistics(MoveStructure& mv_, MoviOptions& movi_options) {
    std::string pattern_file = movi_options.get_index_dir() + "/null_reads.fasta";
    if (movi_options.is_generate_null_reads()) {
        parse_null_reads(movi_options.get_ref_file().c_str(), pattern_file.c_str());
    }

    EmpNullDatabase nulldb;
    nulldb.generate_stats(movi_options, mv_, pattern_file);
    nulldb.compute_stats();
    nulldb.serialize(movi_options);
}

size_t Classifier::initialize_report_file(MoviOptions& movi_options) {

    std::string index_type = program();

    // load the threshold from the null database
    EmpNullDatabase null_db;
    null_db.deserialize(movi_options);
    max_value_thr = std::max(null_db.get_percentile_value(), static_cast<size_t>(MIN_MATCHING_LENGTH)) + 1;

    // initial the report file
    if (!movi_options.is_stdout()) {
        std::string report_file_name = "";
        if (movi_options.get_read_file() != "") {
            report_file_name = movi_options.get_read_file() + "." + index_type + "." + query_type(movi_options) + ".report";;
        } else {
            report_file_name = movi_options.get_mls_file() + ".report";
        }
        std::cerr << "Report file name: " << report_file_name << "\n";
        report_file = std::ofstream(report_file_name);
    }

    std::ostream& out = movi_options.is_stdout() ? std::cout : report_file;
    out.precision(4);
    out << std::setw(30) << std::left << "read id:"
        << std::setw(15) << std::left << "status:"
        << std::setw(19) << std::left << "avg max-value (thr="
        << std::setw(2) << std::left << max_value_thr
        << std::setw(5) << std::left << "):"
        << std::setw(12) << std::left << "above thr:"
        << std::setw(12) << std::left << "below thr:" << std::endl;

    return max_value_thr;
}

void Classifier::close_report_file() {
    report_file.close();
}

// Borrowed from spumoni written by Omar Ahmed: https://github.com/oma219/spumoni/tree/main
void Classifier::classify(std::string read_name, std::vector<uint16_t>& matching_lens, MoviOptions& movi_options) {

    std::vector<size_t> bins_max_value;
    std::string status = "";
    size_t sum_max_bin_values = 0.0;
    size_t bins_above = 0, bins_below = 0;
    size_t start_pos = 0, end_pos = 0;
    size_t bin_width = movi_options.get_bin_width();

    while (start_pos < matching_lens.size()) {
        end_pos = (start_pos + bin_width < matching_lens.size()) ? start_pos + bin_width : matching_lens.size();

        // avoids small regions at the end of read
        if (matching_lens.size() - end_pos < bin_width)
            end_pos = matching_lens.size();

        // grab maximum value in this region and update variables
        auto max_val = *std::max_element(matching_lens.begin()+start_pos, matching_lens.begin()+end_pos);
        if (max_val >= max_value_thr)
            bins_above++;
        else
            bins_below++;
        bins_max_value.push_back(max_val);
        start_pos += (end_pos - start_pos);
    }
    std::for_each(bins_max_value.begin(), bins_max_value.end(), [&] (double n) {sum_max_bin_values += n;});
    bool read_found = (bins_above/(bins_above+bins_below+0.0) > 0.50);
    status = (read_found) ? "FOUND" : "NOT_PRESENT";

    #pragma omp critical
    {

        std::ostream& out = movi_options.is_stdout() ? std::cout : report_file;
        out.precision(3);
        out << std::setw(30) << std::left << read_name
                    << std::setw(15) << std::left << status
                    << std::setw(26) << std::left <<  (sum_max_bin_values+0.0)/bins_max_value.size()
                    << std::setw(12) << std::left <<  bins_above
                    << std::setw(12) << std::left <<  bins_below
                    << "\n";
    }
}