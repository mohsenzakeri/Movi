#ifndef __MOVE_ROW__
#define __MOVE_ROW__

#include <iostream>
#include <vector>
#include <bitset>

// #ifndef MODE
// #define MODE -1// 0: regular, 1: constant, 2: one-bit
// #endif

const uint16_t mask_thresholds1 = static_cast<uint16_t>(~(((1U << 2) - 1) << 0)); // 00000011
const uint16_t mask_thresholds2 = static_cast<uint16_t>(~(((1U << 2) - 1) << 2)); // 00001100
const uint16_t mask_thresholds3 = static_cast<uint16_t>(~(((1U << 2) - 1) << 4)); // 00110000

const uint16_t mask_id =  static_cast<uint16_t>(~(((1U << 8) - 1) << 0));                   // 00000000 11111111
const uint16_t mask_c = static_cast<uint16_t>(~(((1U << 2) - 1) << 8));                     // 00000011 00000000
const uint16_t mask_overflow_n = static_cast<uint16_t>(~(((1U << 1) - 1) << 10));           // 00000100 00000000
const uint16_t mask_overflow_offset = static_cast<uint16_t>(~(((1U << 1) - 1) << 11));      // 00001000 00000000
const uint16_t mask_overflow_thresholds = static_cast<uint16_t>(~(((1U << 1) - 1) << 12));  // 00010000 00000000

class MoveRow{
    public:
        MoveRow () {n = 0; id = 0; overflow_bits = 0;}
        MoveRow(uint16_t n_, uint16_t offset_, uint64_t id_);
        void init(uint16_t n_, uint16_t offset_, uint64_t id_);
        friend std::ostream& operator<<(std::ostream& os, const MoveRow& mr);

        void set_n(uint16_t n_);
        void set_offset(uint16_t offset_);
        void set_id(uint64_t id_);
        void set_c(char c_, std::vector<uint64_t>& alphamap);

        uint16_t get_n() const;
        uint16_t get_n_ff() const;
        uint16_t get_offset() const;
        uint64_t get_id() const;
        uint8_t get_threshold_status(uint16_t i) const;
        char get_c() const;

        void set_overflow_n();
        void set_overflow_offset();
        void set_overflow_thresholds();
        void set_threshold_status(uint16_t i, uint8_t status);
        bool is_overflow_n() const;
        bool is_overflow_offset() const;
        bool is_overflow_thresholds() const;

#if MODE == 1
        uint16_t get_next_up(uint32_t i) { return next_up[i]; }
        uint16_t get_next_down(uint32_t i) { return next_down[i]; }
        void set_next_up(uint32_t i, uint16_t t) { next_up[i] = t; }
        void set_next_down(uint32_t i, uint16_t t) { next_down[i] = t; }
#endif

        uint16_t get_threshold() { return threshold; }
        void set_threshold(uint16_t t) { threshold = t; }

        uint64_t row_size() {
#if MODE == 0
            return 13;
#endif
#if MODE == 1
            return 25;
#endif
#if MODE == 2
            return 12;
#endif
        }
    private:
        uint16_t offset; // offset of the bwt row head of the current run in the new run after the jump
        uint16_t overflow_bits;
        uint16_t n; // length of the run
        uint32_t id; // bwt run after the jump

        // thresholds for all the rows:
#if MODE == 0 or MODE == 1
        uint8_t thresholds_status;
#endif

#if MODE == 1
        // to store pointers for avoiding scanning
        uint16_t next_up[3];
        uint16_t next_down[3];
#endif

    uint16_t threshold;
};

inline uint8_t MoveRow::get_threshold_status(uint16_t i) const {
    const uint16_t mask_thresholds = static_cast<uint16_t>(~(((1U << 2) - 1) << i*2));
    uint8_t status = static_cast<uint8_t>((thresholds_status & (~mask_thresholds)) >> i*2);
    return status;
}

inline uint16_t MoveRow::get_n() const{
    return n;
}

inline uint16_t MoveRow::get_n_ff() const{
    return n;
}

inline uint16_t MoveRow::get_offset() const{
    return offset;
}

inline uint16_t extract_value(uint16_t source, uint16_t mask, uint16_t offset_bits) {
    uint32_t res = (source & (~mask)) >> offset_bits;
    return res;
}

inline uint64_t MoveRow::get_id() const{
    if (overflow_bits != 0) {
        uint64_t res = static_cast<uint64_t>(extract_value(overflow_bits, mask_id, 0));
        res = res << 32;
        uint64_t c = static_cast<uint64_t>(id);
        c = c | res;
        return c;
    } else {
        return static_cast<uint64_t>(id);
    }
}

inline char MoveRow::get_c() const{
    return static_cast<char>(extract_value(overflow_bits, mask_c, 8));
}

inline bool MoveRow::is_overflow_n() const{
    uint16_t res = extract_value(overflow_bits, mask_overflow_n, 10);
    return !static_cast<bool>(res);
}

inline bool MoveRow::is_overflow_offset() const{
    uint16_t res = extract_value(overflow_bits, mask_overflow_offset, 11);
    return !static_cast<bool>(res);
}

inline bool MoveRow::is_overflow_thresholds() const{
    uint16_t res = extract_value(overflow_bits, mask_overflow_thresholds, 12);
    return !static_cast<bool>(res);
}

#endif
