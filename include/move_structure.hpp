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
        void build();
    private:
        char* bwt;
        std::string bwt_string;
        uint64_t length;

        std::vector<move_row*> rlbwt;
};

#endif