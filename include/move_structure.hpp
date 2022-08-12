#ifndef __MOVE_STRUCTURE__
#define __MOVE_STRUCTURE__

#include <fstream>

#include <sdsl/bit_vectors.hpp>

#include "move_query.hpp"

#define END_CHARACTER '$'

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
        MoveStructure(bool verbose_ = false) { verbose = verbose_; }
        MoveStructure(char* input_file, bool verbose_ = false);

        void build(std::ifstream &bwt_file);        
        void query_ms(MoveQuery& mq, bool random);
        std::string reconstruct();

        uint32_t LF(uint32_t row_number);
        uint32_t fast_forward(uint32_t pointer, uint32_t index);

        uint32_t naive_lcp(uint32_t row1, uint32_t row2);
        uint32_t naive_sa(uint32_t bwt_row);

        uint32_t jump_up(uint32_t idx, char c);
        uint32_t jump_down(uint32_t idx, char c);
        bool jump_randomly(uint32_t& idx, char r_char);
        bool jump_naive_lcp(uint32_t& idx, uint32_t pointer, char r_char, uint32_t& lcp);

        void seralize(char* output_dir);
        void deseralize(char* index_dir);
    private:
        std::string bwt_string;
        std::string orig_string;
        bool reconstructed;
        uint32_t length;
        uint32_t r;
        uint32_t end_bwt_row;
        bool verbose;

        std::vector<unsigned char> alphabet;
        std::vector<uint32_t> counts;
        std::vector<uint32_t> alphamap;

        std::vector<sdsl::bit_vector*> occs;
        std::vector<sdsl::rank_support_v<>*> occs_rank;
        sdsl::bit_vector bits;
        sdsl::rank_support_v<> rbits;

        std::vector<move_row> rlbwt;
};

#endif