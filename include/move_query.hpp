#ifndef __MOVE_QUERY__
#define __MOVE_QUERY__

class MoveQuery {
    public:
        MoveQuery() {}
        MoveQuery(std::string q) : query_string(q) {}

        std::string query() { return query_string; }
        uint64_t length() { return query_string.length(); }
        void add_ms(uint64_t len) { ms_lens.push_back(len); }

        friend std::ostream& operator<<(std::ostream& output, const MoveQuery& dt) {
            // output << "The matching statistics are:\n";
            for (int64_t i = dt.ms_lens.size() - 1; i >= 0; i--) {
                auto ms = dt.ms_lens[i];
                output << ms << " ";
            }
            return output;
        }
    // private:
        std::string query_string;
        std::vector<uint64_t> ms_lens;
        std::vector<std::chrono::nanoseconds> costs;
};

#endif