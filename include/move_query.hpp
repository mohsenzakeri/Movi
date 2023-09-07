#ifndef __MOVE_QUERY__
#define __MOVE_QUERY__

class MoveQuery {
    public:
        MoveQuery() {}
        MoveQuery(std::string q) : query_string(q) {}

        std::string query() { return query_string; }
        uint64_t length() { return query_string.length(); }
        void add_pml(uint64_t len) { pml_lens.push_back(len); }
        void add_cost(std::chrono::nanoseconds cost) { costs.push_back(cost); }
        void add_scan(uint64_t scan) { scans.push_back(scan); }
        void add_fastforward(uint64_t fastforward) { fastforwards.push_back(fastforward); }
        std::vector<uint64_t>& get_pml_lens() { return pml_lens; }
        std::vector<uint64_t>& get_scans() { return scans; }
        std::vector<uint64_t>& get_fastforwards() { return fastforwards; }
        std::vector<std::chrono::nanoseconds>& get_costs() { return costs; }

        friend std::ostream& operator<<(std::ostream& output, const MoveQuery& dt) {
            // output << "The matching statistics are:\n";
            for (int64_t i = dt.pml_lens.size() - 1; i >= 0; i--) {
                auto pml = dt.pml_lens[i];
                output << pml << " ";
            }
            return output;
        }
    private:
        std::string query_string;
        std::vector<uint64_t> pml_lens;
        std::vector<uint64_t> scans;
        std::vector<uint64_t> fastforwards;
        std::vector<std::chrono::nanoseconds> costs;
};

#endif