#ifndef __MOVE_ROW__
#define __MOVE_ROW__

#include <iostream>
#include <vector>

const uint32_t mask_p = ~((1U << 8) - 1);
const uint32_t mask_pp = ~(((1U << 8) - 1) << 8);
const uint32_t mask_id = ~(((1U << 8) - 1) << 16);

class MoveRow{
    public:
        MoveRow () { p = 0; n = 0; pp = 0; id = 0; overflow_bits = 0;}
        MoveRow(uint64_t p_, uint16_t n_, uint64_t pp_, uint64_t id_);
        void init(uint64_t p_, uint16_t n_, uint64_t pp_, uint64_t id_);
        friend std::ostream& operator<<(std::ostream& os, const MoveRow& mr);

        void set_p(uint64_t p_);
        void set_n(uint16_t n_);
        void set_pp(uint64_t pp_);
        void set_id(uint64_t id_);

        uint64_t get_p() const;
        uint16_t get_n() const;
        uint64_t get_pp() const;
        uint64_t get_id() const;
//    private:
        uint32_t p;
        uint16_t n;
        uint32_t pp;
        uint32_t id;
        uint32_t overflow_bits;
        uint64_t* thresholds;
};

inline uint64_t MoveRow::get_p() const{
    if (overflow_bits != 0) {
        uint32_t a = overflow_bits & (~mask_p);
        uint64_t b = a;
        b = (b << 32);
        uint64_t c = p;
        c = c | b;
        return c;
    } else {
        return static_cast<uint32_t>(p);
    }
}

inline uint16_t MoveRow::get_n() const{
    return n;
}

inline uint64_t MoveRow::get_pp() const{
    if (overflow_bits != 0) {
        uint32_t a = (overflow_bits & (~mask_pp)) >> 8;
        uint64_t b = a;
        b = (b << 32);
        uint64_t c = pp;
        c = c | b;
        return c;
    } else {
        return static_cast<uint32_t>(pp);
    }
}

inline uint64_t MoveRow::get_id() const{
    if (overflow_bits != 0) {
        uint32_t a = (overflow_bits & (~mask_id)) >> 16;
        uint64_t b = a;
        b = (b << 32);
        uint64_t c = id;
        c = c | b;
        return c;
    } else {
        return static_cast<uint32_t>(id);
    }
    return id;
}
#endif