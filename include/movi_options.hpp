#ifndef MOVI_OPTIONS_HPP
#define MOVI_OPTIONS_HPP

class MoviOptions {
    public:
        MoviOptions() {
            command = "";
            ref_file = "";
            read_file = "";
            mls_file = "";
            index_dir = "";
            LF_type = "reconstruct";
            pml_query = true;
        }
        bool is_no_header() { return no_header; }
        bool is_preprocessed() { return preprocessed; }
        bool is_split() { return split; }
        bool is_thresholds() { return thresholds; }
        bool is_mmap() { return mmap; }
        bool no_prefetch() { return !prefetch; }
        bool is_no_output() { return no_output; }
        bool is_output_ids() { return output_ids; }
        bool is_stdout() { return write_stdout; }
        bool is_verbose() { return verbose; }
        bool is_logs() { return logs; }
        bool is_debug() { return debug; }
        bool if_verify() { return verify; }
        bool is_random_repositioning() { return random_repositioning; }
        bool is_get_sa_entries() { return get_sa_entries; }
        bool is_pml() { return pml_query; }
        bool is_zml() { return zml_query; }
        bool is_count() { return count_query; }
        bool is_kmer() { return kmer_query; }
        bool is_kmer_count() { return kmer_count; }
        bool is_reverse() { return reverse; }
        bool is_multi_ftab() { return multi_ftab; }
        bool is_classify() { return classify; }
        bool is_filter() { return filter; }
        bool is_early_stop() { return early_stop; }
        bool is_report_all() { return report_all; }
        bool is_multi_classify() { return multi_classify; }
        bool is_generate_null_reads() { return generate_null_reads; }
        size_t get_bin_width() { return bin_width; }
        int ignore_illegal_chars_status() { return ilc; }
        size_t get_strands() { return strands; }
        bool is_full_color() { return full_color; }
        bool is_compressed() { return compress; }
        bool is_color_move_rows() { return color_move_rows; }
        int get_thres() { return cl_thres; }
        size_t get_threads() { return threads; }
        uint32_t get_k () { return k; }
        uint32_t get_ftab_k () { return ftab_k; }
        uint32_t get_tally_checkpoints () { return tally_checkpoints; }
        uint64_t get_SA_sample_rate() { return SA_sample_rate; }
        std::string get_command() { return command; }
        std::string get_LF_type() { return LF_type; }

        std::string get_ref_file() { return ref_file; }
        std::string get_bwt_file() { return bwt_file; }
        std::string get_read_file() { return read_file; }
        std::string get_mls_file() { return mls_file; }
        std::string get_index_dir() { return index_dir; }
        std::string get_out_file() { return out_file; }

        void set_no_header(bool no_header_) { no_header = no_header_; }
        void set_preprocessed(bool preprocessed_) { preprocessed = preprocessed_; }
        void set_thresholds(bool thresholds_) { thresholds = thresholds_; }
        void set_no_output(bool no_output_) { no_output = no_output_; }
        void set_output_ids(bool output_ids_) { output_ids = output_ids_; }
        void set_stdout(bool write_stdout_) { write_stdout = write_stdout_; }
        void set_verbose(bool verbose_) { verbose = verbose_; }
        void set_logs(bool logs_) { logs = logs_; }
        void set_debug(bool debug_) { debug = debug_; }
        void set_verify(bool verify_) { verify = verify_; }
        void set_random_repositioning(bool random_repositioning_) { random_repositioning = random_repositioning_; }
        void set_get_sa_entries(bool get_sa_entries_) { get_sa_entries = get_sa_entries_; }
        void set_pml()   { pml_query = true; count_query = false; kmer_query = false; zml_query = false; }
        void set_zml()   { zml_query = true; pml_query = false; count_query = false; kmer_query = false; }
        void set_count() { count_query = true; pml_query = false; kmer_query = false; zml_query = false; }
        void set_kmer()  { kmer_query = true; pml_query = false; count_query = false; zml_query = false; }
        void set_kmer_count(bool kmer_count_) { kmer_count = kmer_count_; }
        void set_k(uint32_t k_) { k = k_; }
        void set_ftab_k(uint32_t ftab_k_) { ftab_k = ftab_k_; }
        void set_tally_checkpoints(uint32_t tally_checkpoints_) { tally_checkpoints = tally_checkpoints_; }
        void set_SA_sample_rate(uint64_t SA_sample_rate_) { SA_sample_rate = SA_sample_rate_; }
        void set_multi_ftab(bool multi_ftab_) { multi_ftab = multi_ftab_; }
        void set_reverse(bool reverse_) { reverse = reverse_; }
        void set_classify(bool classify_) { classify = classify_; }
        void set_filter(bool filter_) { filter = filter_; }
        void set_multi_classify(bool multi_classify_) { multi_classify = multi_classify_; }
        void set_early_stop(bool val) { early_stop = val; }
        void set_report_all(bool report_all_) { report_all = report_all_; }
        void set_generate_null_reads(bool generate_null_reads_) { generate_null_reads = generate_null_reads_; }
        void set_bin_width(size_t bin_width_) { bin_width = bin_width_; }
        bool set_ignore_illegal_chars(int ilc_) {
            if (ilc_ > 2 or ilc_ < 1)
                return false;
            ilc = ilc_;
	        if (ilc == 2) // A random character is generated for each illegal character
                std::srand(time(0));
            return true;
        }
        void set_mmap(bool mmap_) { mmap = mmap_; }
        void set_prefetch(bool prefetch_) { prefetch = prefetch_; }
        void set_threads(size_t threads_) { threads = threads_; }
        void set_strands(size_t strands_) { strands = strands_; }
        void set_command(std::string command_) { command = command_; }
        bool set_LF_type(std::string type) {
            if (type == "reconstruct" or type == "sequential" or type == "random")
                LF_type = type;
            else
                return false;
            return true;
        }
        void set_full_color(bool val) { full_color = val; }
        void set_compress(bool val) { compress = val; }
        void set_color_move_rows(bool val) { color_move_rows = val; }
        void set_thres(uint8_t val) { cl_thres = val; }
        void set_out_file(std::string out_file_) { out_file = out_file_; }

        void set_ref_file(std::string file_address) { ref_file = file_address; }
        void set_bwt_file(std::string file_address) { bwt_file = file_address; }
        void set_read_file(std::string file_address) { read_file = file_address; }
        void set_mls_file(std::string file_address) { mls_file = file_address; }
        void set_index_dir(std::string dir) { index_dir = dir; }

        void print_options() {
            std::cerr << "command:\t" << command << "\n";
            std::cerr << "ref_file:\t" << ref_file << "\n";
            std::cerr << "bwt_file:\t" << bwt_file << "\n";
            std::cerr << "read_file:\t" << read_file << "\n";
            std::cerr << "mls_file:\t" << mls_file << "\n";
            std::cerr << "index_dir:\t" << index_dir << "\n";
            std::cerr << "LF_type:\t" << LF_type << "\n";
            std::cerr << "no_header:\t" << no_header << "\n";
            std::cerr << "preprocessed:\t" << preprocessed << "\n";
            std::cerr << "ilc:\t" << ilc << "\n";
            std::cerr << "split:\t" << split << "\n";
            std::cerr << "thresholds:\t" << thresholds << "\n";
            std::cerr << "random_repositioning:\t" << random_repositioning << "\n";
            std::cerr << "get_sa_entries:\t" << get_sa_entries << "\n";
            std::cerr << "pml_query:\t" << pml_query << "\n";
            std::cerr << "zml_query:\t" << zml_query << "\n";
            std::cerr << "count_query:\t" << count_query << "\n";
            std::cerr << "kmer_query:\t" << kmer_query << "\n";
            std::cerr << "kmer_count:\t" << kmer_count << "\n";
            std::cerr << "reverse:\t" << reverse << "\n";
            std::cerr << "mmap:\t" << mmap << "\n";
            std::cerr << "prefetch:\t" << prefetch << "\n";
            std::cerr << "strands:\t" << strands << "\n";
            std::cerr << "threads:\t" << threads << "\n";
            std::cerr << "k:\t" << k << "\n";
            std::cerr << "ftab_k:\t" << ftab_k << "\n";
            std::cerr << "tally_checkpoints:\t" << tally_checkpoints << "\n";
            std::cerr << "SA_sample_rate:\t" << SA_sample_rate << "\n";
            std::cerr << "multi_ftab:\t" << multi_ftab << "\n";
            std::cerr << "verify:\t" << verify << "\n";
            std::cerr << "classify:\t" << classify << "\n";
            std::cerr << "filter:\t" << filter << "\n";
            std::cerr << "generate_null_reads:\t" << generate_null_reads << "\n";
            std::cerr << "bin_width:\t" << bin_width << "\n";
            std::cerr << "no_output:\t" << no_output << "\n";
            std::cerr << "output_ids:\t" << output_ids << "\n";
            std::cerr << "stdout:\t" << write_stdout << "\n";
            std::cerr << "debug:\t" << debug << "\n";
            std::cerr << "verbose:\t" << verbose << "\n";
            std::cerr << "logs:\t" << logs << "\n";
        }
    private:
        std::string command;
        std::string ref_file;
        std::string bwt_file;
        std::string read_file;
        std::string mls_file;
        std::string index_dir;
        std::string LF_type;
        std::string out_file;
        bool no_header = false;
        bool preprocessed = false;
        int ilc = 0;
        bool split = false;
        bool thresholds = false;
        bool random_repositioning = false;
        bool get_sa_entries = false;
        bool pml_query = true;
        bool zml_query = false;
        bool count_query = false;
        bool kmer_query = false;
        bool kmer_count = false;
        bool reverse = false;
        bool mmap = false;
        bool prefetch = true;
        size_t strands = 16;
        size_t threads = 1;
        uint32_t k = 31;
        uint32_t ftab_k = 0;
        uint32_t tally_checkpoints = 20;
        uint64_t SA_sample_rate = 100;
        bool multi_ftab = false;
        bool verify = false;
        bool classify = false;
        bool filter = false;
        bool multi_classify = false;
        bool early_stop = false;
        bool report_all = false;
        bool generate_null_reads = false;
        size_t bin_width = 150;

        bool no_output = false;
        bool output_ids = false;
        bool write_stdout = false;
        bool debug = false;
        bool verbose = false;
        bool logs = false;
        bool full_color = false;
        bool compress = false;
        bool color_move_rows = false;
        // Classification threshold (only consider PMLs with len >= cl_thres)
        uint8_t cl_thres = 5;
};

#endif
