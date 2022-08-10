#ifndef __MOVE_QUERY__
#define __MOVE_QUERY__

class MoveQuery {
    public:
        MoveQuery() {}
        MoveQuery(std::string q) : query_string(q) {}

        std::string query() { return query_string; }
        uint32_t length() { return query_string.length(); }
        void add_ms(uint32_t len) { ms_lens.push_back(len); }

        friend std::ostream& operator<<(std::ostream& output, const MoveQuery& dt) {
            output << "The matching statistics are:\n";
            for (auto& ms: dt.ms_lens) {
                output << ms << " ";
            }
            return output;
        }
    private:
        std::string query_string;
        std::vector<uint32_t> ms_lens;
};

#endif