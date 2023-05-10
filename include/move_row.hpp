#ifndef __MOVE_ROW__
#define __MOVE_ROW__

#include <iostream>
#include <vector>
#include <bitset>

// const uint32_t mask_p = ~((1U << 8) - 1);
// const uint32_t mask_pp = ~(((1U << 8) - 1) << 8);
// const uint32_t mask_id = ~(((1U << 8) - 1) << 16);
// const uint32_t mask_c = ~(((1U << 8) - 1) << 24);

const uint16_t mask_id =  static_cast<uint16_t>(~((1U << 8) - 1));
const uint16_t mask_c = static_cast<uint16_t>(~(((1U << 2) - 1) << 8));
const uint16_t mask_overflow_n = static_cast<uint16_t>(~(((1U << 1) - 1) << 10));
const uint16_t mask_overflow_offset = static_cast<uint16_t>(~(((1U << 1) - 1) << 11));


class MoveRow{
    public:
        MoveRow () {n = 0; id = 0; overflow_bits = 0;}
        MoveRow(uint16_t n_, uint16_t offset_, uint64_t id_);
        void init(uint16_t n_, uint16_t offset_, uint64_t id_);
        void init(uint16_t n_, uint16_t offset_, uint64_t id_, char c_);
        friend std::ostream& operator<<(std::ostream& os, const MoveRow& mr);

        // void set_p(uint64_t p_);
        // void set_pp(uint64_t pp_);
        void set_n(uint16_t n_);
        void set_offset(uint16_t offset_);
        void set_id(uint64_t id_);
        void set_c(char c_);

        // uint64_t get_p() const;
        // uint64_t get_pp() const;
        uint16_t get_n() const;
        uint16_t get_offset() const;
        uint64_t get_id() const;
        char get_c() const;

        void set_overflow_n();
        void set_overflow_offset();
        bool is_overflow_n() const;
        bool is_overflow_offset() const;
//    private:
        // offset based: uint32_t p; // bwt row of the head before the jump
        // offset based: uint32_t pp; // bwt row of the head after the jump
        uint16_t offset; // offset of the bwt row head of the current run in the new run after the jump
        uint16_t overflow_bits;
        uint16_t n; // length of the run
        uint16_t thresholds[3];
        uint32_t id; // bwt run after the jump

};

/* inline uint64_t MoveRow::get_p() const{
    /*if (overflow_bits != 0) {
        uint32_t a = overflow_bits & (~mask_p);
        uint64_t b = a;
        b = (b << 32);
        uint64_t c = p;
        c = c | b;
        return c;
    } else {
        return static_cast<uint32_t>(p);
    }
}*/

/*inline uint64_t MoveRow::get_pp() const{
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
}*/

inline uint16_t MoveRow::get_n() const{
    return n;
}

inline uint16_t MoveRow::get_offset() const{
    return offset;
}

inline uint64_t MoveRow::get_id() const{
    if (overflow_bits != 0) {
        uint32_t a = overflow_bits & (~mask_id);
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
    uint32_t a = (overflow_bits & (~mask_c)) >> 8;
    char c = static_cast<char>(a);
    return c;
}

inline bool MoveRow::is_overflow_n() const{
    uint32_t a = (overflow_bits & (~mask_overflow_n)) >> 10;
    bool b = static_cast<bool>(a);
    return !b;
}

inline bool MoveRow::is_overflow_offset() const{
    uint32_t a = (overflow_bits & (~mask_overflow_offset)) >> 11;
    bool b = static_cast<bool>(a);
    return !b;
}


#endif