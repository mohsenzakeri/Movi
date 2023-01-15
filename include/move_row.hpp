#ifndef __MOVE_ROW__
#define __MOVE_ROW__

#include <iostream>
#include <vector>

const uint32_t mask_p = ~((1U << 8) - 1);
const uint32_t mask_pp = ~(((1U << 8) - 1) << 8);
const uint32_t mask_id = ~(((1U << 8) - 1) << 16);
const uint32_t mask_c = ~(((1U << 8) - 1) << 24);

class MoveRow{
    public:
        MoveRow () { p = 0; n = 0; pp = 0; id = 0; overflow_bits = 0;}
        MoveRow(uint64_t p_, uint16_t n_, uint64_t pp_, uint64_t id_);
        void init(uint64_t p_, uint16_t n_, uint64_t pp_, uint64_t id_);
        void init(uint64_t p_, uint16_t n_, uint64_t pp_, uint64_t id_, char c_);
        friend std::ostream& operator<<(std::ostream& os, const MoveRow& mr);

        void set_p(uint64_t p_);
        void set_n(uint16_t n_);
        void set_pp(uint64_t pp_);
        void set_id(uint64_t id_);
        void set_c(char c_);

        uint64_t get_p() const;
        uint16_t get_n() const;
        uint64_t get_pp() const;
        uint64_t get_id() const;
        char get_c() const;
//    private:
        uint32_t p; // bwt row of the head before the jump
        uint16_t n; // length of the run
        uint32_t pp; // bwt row of the head after the jump
        uint32_t id; // bwt run after the jump
        uint32_t overflow_bits;

        uint16_t thresholds[3]; // 4*8 = 32 : 1.9 GB ---- 3 * 2 = 6
        uint16_t threshold_1bit; // 4 + 2 + 4 + 4 + 4 = 18 : 534 MB

        // Why the number is close to 2 for 1bit?
        // -- 15.9M vs 51.4M -- tot: 24.6M vs 89.9M
        // --  1.3M vs  4.4M -- tot:  2.2M vs  7.7M
        //
        //
        
        // Store the character bit in the remaining 8 bit of the overflow_bits

        // What's the bug for the big experiment (1bit)

        // can we reduce the size by converting 64 to 16 and then what happens in terms of the query times?
        // -- query time remains similar
        // -- 1.9 GB goes down to 1.3 GB
        //
        //

        // convinced to do the Nate suggestion, we need full bwt offset in order to interpret thresholds? or not? PML vs ML -- new branch?
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

inline char MoveRow::get_c() const{
    uint32_t a = (overflow_bits & (~mask_c)) >> 24;
    char c = static_cast<char>(a);
    return c;
}


#endif