#ifndef __MOVE_STRUCTURE__
#define __MOVE_STRUCTURE__

#include <sdsl/bit_vectors.hpp>

struct move_row{
    move_row(uint32_t p, uint32_t n, 
                 uint32_t pp, uint32_t id, char c) {
        this->p = p;
        this->n = n;
        this->pp = pp;
        this->id = id;
        this->c = c;
    }
    uint32_t p;
    uint32_t n;
    uint32_t pp;
    uint32_t id;
    char c;
};

class MoveStructure {
    public:
        MoveStructure() { }
        MoveStructure(char* input_file);
        uint32_t LF(uint32_t row_number);
        void build(std::ifstream &bwt_file);
    private:
        std::string bwt_string;
        uint32_t length;

        std::vector<unsigned char> alphabet;
        std::vector<uint32_t> counts;
        std::map<unsigned char, uint32_t> alphamap;

        std::vector<sdsl::bit_vector*> occs;
        std::vector<sdsl::rank_support_v<>*> occs_rank;

        std::vector<move_row*> rlbwt;
};

#endif