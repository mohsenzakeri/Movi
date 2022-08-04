#ifndef __MOVE_STRUCTURE__
#define __MOVE_STRUCTURE__

#include <fstream>

#include <sdsl/bit_vectors.hpp>

#include "move_query.hpp"

struct move_row{
    move_row () { p = 0; n = 0; pp = 0; id = 0; c = 0; }
    move_row(uint32_t p, uint32_t n, 
                 uint32_t pp, uint32_t id, char c) {
        this->p = p;
        this->n = n;
        this->pp = pp;
        this->id = id;
        this->c = c;
    }
    void init(uint32_t p, uint32_t n, 
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

        void build(std::ifstream &bwt_file);        
        void query_ms(MoveQuery& mq);

        uint32_t LF(uint32_t row_number);
        uint32_t fast_forward(uint32_t pointer, uint32_t index);
        uint32_t jump_up(uint32_t idx, char c);
        uint32_t jump_down(uint32_t idx, char c);
    private:
        std::string bwt_string;
        uint32_t length;
        uint32_t r;

        std::vector<unsigned char> alphabet;
        std::vector<uint32_t> counts;
        std::map<unsigned char, uint32_t> alphamap;

        std::vector<sdsl::bit_vector*> occs;
        std::vector<sdsl::rank_support_v<>*> occs_rank;

        std::vector<move_row> rlbwt;
};

#endif