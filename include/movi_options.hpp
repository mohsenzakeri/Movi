#ifndef __MOVI_OPTIONS__
#define __MOVI_OPTIONS__

class MoviOptions {
    public:
        MoviOptions() {
            command = "";
            ref_file = "";
            read_file = "";
            pml_file = "";
            index_dir = "";
            LF_type = "reconstruct";
            pml_query = true;
        }
        bool is_split() { return split; }
        bool no_prefetch() { return !prefetch; }
        bool is_stdout() { return write_stdout; }
        bool is_verbose() { return verbose; }
        bool is_logs() { return logs; }
        bool is_pml() { return pml_query; }
        bool is_count() { return count_query; }
        bool is_reverse() { return reverse; }
        size_t get_strands() { return strands; }
        std::string get_command() { return command; }
        std::string get_LF_type() { return LF_type; }

        std::string get_ref_file() { return ref_file; }
        std::string get_bwt_file() { return bwt_file; }
        std::string get_read_file() { return read_file; }
        std::string get_pml_file() { return pml_file; }
        std::string get_index_dir() { return index_dir; }

        void set_stdout(bool write_stdout_) { write_stdout = write_stdout_; }
        void set_verbose(bool verbose_) { verbose = verbose_; }
        void set_logs(bool logs_) { logs = logs_; }
        void set_pml(bool pml_) { pml_query = pml_; }
        void set_count(bool count_) { count_query = count_; pml_query = false;}
        void set_reverse(bool reverse_) { reverse = reverse_; }
        void set_prefetch(bool prefetch_) { prefetch = prefetch_; }
        void set_strands(size_t strands_) { strands = strands_; }
        void set_command(std::string command_) { command = command_; }
        bool set_LF_type(std::string type) {
            if (type == "reconstruct" or type == "sequential" or type == "random")
                LF_type = type;
            else
                return false;
            return true;
        }

        void set_ref_file(std::string file_address) { ref_file = file_address; }
        void set_bwt_file(std::string file_address) { bwt_file = file_address; }
        void set_read_file(std::string file_address) { read_file = file_address; }
        void set_pml_file(std::string file_address) { pml_file = file_address; }
        void set_index_dir(std::string dir) { index_dir = dir; }

        void print_options() {
            std::cerr << "split:\t" << split << "\n";
            std::cerr << "ref_file:\t" << ref_file << "\n";
            std::cerr << "bwt_file:\t" << bwt_file << "\n";
            std::cerr << "read_file:\t" << read_file << "\n";
            std::cerr << "pml_file:\t" << pml_file << "\n";
            std::cerr << "index_dir:\t" << index_dir << "\n";
            std::cerr << "LF_type:\t" << LF_type << "\n";
            std::cerr << "pml_query:\t" << pml_query << "\n";
            std::cerr << "count_query:\t" << count_query << "\n";
            std::cerr << "reverse:\t" << reverse << "\n";
            std::cerr << "prefetch:\t" << prefetch << "\n";
            std::cerr << "strands:\t" << strands << "\n";
            std::cerr << "stdout:\t" << write_stdout << "\n";
            std::cerr << "verbose:\t" << verbose << "\n";
            std::cerr << "logs:\t" << logs << "\n";
        }
    private:
        std::string command;
        std::string ref_file;
        std::string bwt_file;
        std::string read_file;
        std::string pml_file;
        std::string index_dir;
        std::string LF_type;
        bool split = false;
        bool pml_query = false;
        bool count_query = false;
        bool reverse = false;
        bool prefetch = true;
        size_t strands = 16;
        bool write_stdout = false;
        bool verbose = false;
        bool logs = false;
};

#endif