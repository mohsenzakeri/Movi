#ifndef MOVE_QUERY_HPP
#define MOVE_QUERY_HPP

class MoveQuery {
    public:
        MoveQuery() {}
        MoveQuery(std::string q) : query_string(q) {}

        std::string& query() { return query_string; }
        uint64_t length() { return query_string.length(); }

        void add_color(uint32_t color_id) {
            // Check if the color id is ever larger than 2^32?

            matching_colors.push_back(color_id);
        }

        void add_ml(uint64_t len_, bool to_stdout) {
            uint32_t len = static_cast<uint32_t>(len_);
            if (len_ > std::numeric_limits<uint32_t>::max()) {
                len = std::numeric_limits<uint32_t>::max();
            }

            matching_lens.push_back(len);
            if (to_stdout) {
                std::string len_str = std::to_string(len);
                std::reverse(len_str.begin(), len_str.end());
                matching_lengths_string += " " + len_str;
            }
        }
        void add_sa_entries(uint64_t sa_entry) {
            sa_entries.push_back(sa_entry);
        }

        void add_cost(std::chrono::nanoseconds cost) { costs.push_back(cost); }
        void add_scan(uint64_t scan) { scans.push_back(scan); }
        void add_fastforward(uint64_t fastforward) { fastforwards.push_back(fastforward); }
        std::vector<uint32_t>& get_matching_lengths() { return matching_lens; }
        std::vector<uint32_t>& get_matching_colors() { return matching_colors; }
        std::string& get_matching_lengths_string() { return matching_lengths_string; }
        std::vector<uint64_t>& get_sa_entries() { return sa_entries; }
        std::vector<uint16_t>& get_scans() { return scans; }
        std::vector<uint16_t>& get_fastforwards() { return fastforwards; }
        std::vector<std::chrono::nanoseconds>& get_costs() { return costs; }

        friend std::ostream& operator<<(std::ostream& output, const MoveQuery& dt) {
            // output << "The matching statistics are:\n";
            for (int64_t i = dt.matching_lens.size() - 1; i >= 0; i--) {
                auto pml = dt.matching_lens[i];
                output << pml << " ";
            }
            return output;
        }
    private:
        std::string query_string;
        std::string matching_lengths_string = "";
        std::vector<uint64_t> sa_entries;
        std::vector<uint32_t> matching_lens;
        std::vector<uint32_t> matching_colors;
        std::vector<uint16_t> scans;
        std::vector<uint16_t> fastforwards;
        std::vector<std::chrono::nanoseconds> costs;
};

#endif