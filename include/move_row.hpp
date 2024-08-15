#ifndef __MOVE_ROW__
#define __MOVE_ROW__

#include <iostream>
#include <vector>
#include <bitset>

// #ifndef MODE
// #define MODE -1// 0: regular, 1: constant, 2: one-bit
// #endif

#if MODE == 0 or MODE == 1
const uint8_t mask_thresholds1 = static_cast<uint8_t>(~(((1U << 2) - 1) << 0)); // 00000011
const uint8_t mask_thresholds2 = static_cast<uint8_t>(~(((1U << 2) - 1) << 2)); // 00001100
const uint8_t mask_thresholds3 = static_cast<uint8_t>(~(((1U << 2) - 1) << 4)); // 00110000
const uint8_t mask_c = static_cast<uint8_t>(~(((1U << 2) - 1) << 6));           // 11000000
const uint8_t mask_id =  static_cast<uint8_t>(~(((1U << 4) - 1) << 0));                  // 00001111
const uint8_t mask_overflow_n = static_cast<uint8_t>(~(((1U << 1) - 1) << 4));           // 00010000
const uint8_t mask_overflow_offset = static_cast<uint8_t>(~(((1U << 1) - 1) << 5));      // 00100000
const uint8_t mask_overflow_thresholds = static_cast<uint8_t>(~(((1U << 1) - 1) << 6));  // 01000000
#define MAX_RUN_LENGTH 65535 // 2^16-1
#endif
#if MODE == 3
const uint16_t mask_id =  static_cast<uint16_t>(~(((1U << 4) - 1) << 12));                  // 11110000 00000000
const uint16_t mask_offset =  static_cast<uint16_t>(~(((1U << 12) - 1) << 0));              // 00001111 11111111
const uint16_t mask_n =  static_cast<uint16_t>(~(((1U << 12) - 1) << 0));                   // 00001111 11111111
const uint16_t mask_c = static_cast<uint16_t>(~(((1U << 4) - 1) << 12));                    // 11110000 00000000
// const uint16_t mask_overflow_n = static_cast<uint16_t>(~(((1U << 1) - 1) << 4));            // 01000000 00000000
// const uint16_t mask_overflow_offset = static_cast<uint16_t>(~(((1U << 1) - 1) << 5));       // 10000000 00000000
#define MAX_RUN_LENGTH 4095 // 2^12-1
#endif


class MoveRow{
    public:
#if MODE == 0 or MODE == 1 or MODE == 2
        MoveRow () {n = 0; id = 0; overflow_bits = 0;}
#endif
#if MODE == 3
        MoveRow () {n = 0; id = 0; offset = 0;}
#endif
        MoveRow(uint16_t n_, uint16_t offset_, uint64_t id_);
        void init(uint16_t n_, uint16_t offset_, uint64_t id_);
        friend std::ostream& operator<<(std::ostream& os, const MoveRow& mr);

        void set_n(uint16_t n_);
        void set_offset(uint16_t offset_);
        void set_id(uint64_t id_);
        void set_c(char c_, std::vector<uint64_t>& alphamap);

        uint16_t get_n() const;
        uint16_t get_offset() const;
        uint64_t get_id() const;
        char get_c() const;

        void set_overflow_n();
        void set_overflow_offset();
        void set_overflow_thresholds();
        bool is_overflow_n() const;
        bool is_overflow_offset() const;

#if MODE == 0 or MODE == 1 or MODE == 2
        uint8_t get_threshold_status(uint16_t i) const;
        void set_threshold_status(uint16_t i, uint8_t status);
        bool is_overflow_thresholds() const;
        uint16_t get_threshold() { return threshold; }
        void set_threshold(uint16_t t) { threshold = t; }
#endif

#if MODE == 1
        uint16_t get_next_up(uint32_t i) { return next_up[i]; }
        uint16_t get_next_down(uint32_t i) { return next_down[i]; }
        void set_next_up(uint32_t i, uint16_t t) { next_up[i] = t; }
        void set_next_down(uint32_t i, uint16_t t) { next_down[i] = t; }
#endif
        uint64_t row_size() {
#if MODE == 0
            return 12;
#endif
#if MODE == 1
            return 24;
#endif
#if MODE == 2
            return 12;
#endif
#if MODE == 3
            return 8;
#endif
        }
    private:
        uint32_t id; // bwt run after the LF-jump
        uint16_t n; // length of the run
        uint16_t offset; // offset of the bwt row head of the current run in the new run after the LF-jump

#if MODE == 0 or MODE == 1 or MODE == 2
        uint16_t threshold;
        uint8_t overflow_bits;
#endif
#if MODE == 0 or MODE == 1
        uint8_t thresholds_status; // Whether each threshold is at the boundary or it's a non-trivial value
#endif

#if MODE == 1
        // to store pointers for avoiding scanning
        uint16_t next_up[3];
        uint16_t next_down[3];
#endif
};

#if MODE == 0 or MODE == 1 or MODE == 2
inline uint8_t extract_value(uint8_t source, uint8_t mask, uint16_t offset_bits) {
    uint8_t res = (source & (~mask)) >> offset_bits;
    return res;
}
#endif
#if MODE == 3
inline uint16_t extract_value(uint16_t source, uint16_t mask, uint16_t offset_bits) {
    uint16_t res = (source & (~mask)) >> offset_bits;
    return res;
}
#endif

inline uint16_t MoveRow::get_n() const{
#if MODE == 0 or MODE == 1 or MODE == 2
    return n;
#endif
#if MODE == 3
    uint16_t res = static_cast<uint16_t>(extract_value(n, mask_n, 0));
    return res;
#endif
}

inline uint16_t MoveRow::get_offset() const{
#if MODE == 0 or MODE == 1 or MODE == 2
    return offset;
#endif
#if MODE == 3
    uint16_t res = static_cast<uint16_t>(extract_value(offset, mask_offset, 0));
    return res;
#endif
}

inline uint64_t MoveRow::get_id() const{
#if MODE == 0 or MODE == 1 or MODE == 2
    if (overflow_bits != 0) {
        uint64_t res = static_cast<uint64_t>(extract_value(overflow_bits, mask_id, 0));
#endif
#if MODE == 3
    if (offset >= 2^12) {
        uint64_t res = static_cast<uint64_t>(extract_value(offset, mask_id, 12));
#endif
        res = res << 32;
        uint64_t c = static_cast<uint64_t>(id);
        c = c | res;
        return c;
    } else {
        return static_cast<uint64_t>(id);
    }
}

inline char MoveRow::get_c() const{
#if MODE == 0 or MODE == 1 or MODE == 2
    return static_cast<char>(extract_value(thresholds_status, mask_c, 6));
#endif
#if MODE == 3
    return static_cast<char>(extract_value(n, mask_c, 12));
#endif
}

inline bool MoveRow::is_overflow_n() const{
#if MODE == 0 or MODE == 1 or MODE == 2
    uint8_t res = extract_value(overflow_bits, mask_overflow_n, 4);
#endif
#if MODE == 3
    std::cerr << "The length overflow should not occur in the compressed mode.\n";
    uint16_t res = 1;
    // uint16_t res = extract_value(n, mask_overflow_n, 14);
#endif
    return !static_cast<bool>(res);
}

inline bool MoveRow::is_overflow_offset() const{
#if MODE == 0 or MODE == 1 or MODE == 2
    uint8_t res = extract_value(overflow_bits, mask_overflow_offset, 5);
#endif
#if MODE == 3
    std::cerr << "The offset overflow should not occur in the compressed mode.\n";
    uint16_t res = 1;
    // uint16_t res = extract_value(n, mask_overflow_offset, 15);
#endif
    return !static_cast<bool>(res);
}

#if MODE == 0 or MODE == 1 or MODE == 2
inline uint8_t MoveRow::get_threshold_status(uint16_t i) const {
    const uint8_t mask_thresholds = static_cast<uint8_t>(~(((1U << 2) - 1) << i*2));
    uint8_t status = static_cast<uint8_t>((thresholds_status & (~mask_thresholds)) >> i*2);
    return status;
}

inline bool MoveRow::is_overflow_thresholds() const{
    uint8_t res = extract_value(overflow_bits, mask_overflow_thresholds, 6);
    return !static_cast<bool>(res);
}
#endif

#endif