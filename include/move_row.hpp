#ifndef __MOVE_ROW__
#define __MOVE_ROW__

#include <iostream>

class MoveRow{
    public:
        MoveRow () { p = 0; n = 0; pp = 0; id = 0; overflow_bits = 0;}
        MoveRow(uint64_t p_, uint16_t n_, uint64_t pp_, uint64_t id_);
        void init(uint64_t p_, uint16_t n_, uint64_t pp_, uint64_t id_);

        void set_p(uint64_t p_);
        void set_n(uint16_t n_) {n = n_;}
        void set_pp(uint64_t pp_);
        void set_id(uint64_t id_);

        uint64_t get_p();
        uint16_t get_n() {return n;}
        uint64_t get_pp();
        uint64_t get_id();
    private:
        uint32_t p;
        uint16_t n;
        uint32_t pp;
        uint32_t id;
        uint32_t overflow_bits;
};

#endif