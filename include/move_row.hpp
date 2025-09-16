#ifndef MOVE_ROW_HPP
#define MOVE_ROW_HPP

#include <iostream>
#include <vector>
#include <bitset>

#include "utils.hpp"

#if MODE == 0 or MODE == 1 or MODE == 4
const uint8_t mask_thresholds1 = static_cast<uint8_t>(~(((1U << 2) - 1) << 0));         // 00000011
const uint8_t mask_thresholds2 = static_cast<uint8_t>(~(((1U << 2) - 1) << 2));         // 00001100
const uint8_t mask_thresholds3 = static_cast<uint8_t>(~(((1U << 2) - 1) << 4));         // 00110000
const uint8_t mask_c = static_cast<uint8_t>(~(((1U << 2) - 1) << 6));                   // 11000000
const uint8_t mask_id =  static_cast<uint8_t>(~(((1U << 4) - 1) << 0));                 // 00001111
const uint8_t mask_overflow_n = static_cast<uint8_t>(~(((1U << 1) - 1) << 4));          // 00010000
const uint8_t mask_overflow_offset = static_cast<uint8_t>(~(((1U << 1) - 1) << 5));     // 00100000
const uint8_t mask_overflow_thresholds = static_cast<uint8_t>(~(((1U << 1) - 1) << 6)); // 01000000
#define MAX_RUN_LENGTH 65535 // 2^16 - 1
#endif

#if MODE == 3
#define SHIFT_C 13
#define SHIFT_ID 12
#define ID_SIG_BITS 4
#define LENGTH_BITS 12
#define C_BITS 3
const uint16_t mask_id =  static_cast<uint16_t>(~(((1U << ID_SIG_BITS) - 1) << SHIFT_ID)); // 11110000 00000000
const uint16_t mask_offset =  static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << 0));    // 00001111 11111111
const uint16_t mask_n =  static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << 0));         // 00001111 11111111
const uint16_t mask_c = static_cast<uint16_t>(~(((1U << C_BITS) - 1) << SHIFT_C));         // 11100000 00000000
#define MAX_RUN_LENGTH 4095 // 2^12-1
#endif

#if MODE == 6
#define SHIFT_C 13
#define SHIFT_ID 12
#define SHIFT_THRESHOLD_1 11
#define SHIFT_THRESHOLD_2 11
#define SHIFT_THRESHOLD_3 12
#define ID_SIG_BITS 4
#define LENGTH_BITS 11
#define C_BITS 3
const uint16_t mask_id = static_cast<uint16_t>(~(((1U << ID_SIG_BITS) - 1) << SHIFT_ID));         // 11110000 00000000
const uint16_t mask_offset = static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << 0));            // 00000111 11111111
const uint16_t mask_n = static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << 0));                 // 00000111 11111111
const uint16_t mask_c = static_cast<uint16_t>(~(((1U << C_BITS) - 1) << SHIFT_C));                // 11100000 00000000
const uint16_t mask_thresholds1 = static_cast<uint16_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_1)); // 00001000 00000000
const uint16_t mask_thresholds2 = static_cast<uint16_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_2)); // 00001000 00000000
const uint16_t mask_thresholds3 = static_cast<uint16_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_3)); // 00010000 00000000
#define MAX_RUN_LENGTH 2047 // 2^11-1
#endif

#if MODE == 2
#define SHIFT_ID1 10
#define SHIFT_ID2 14
#define SHIFT_ID1_RES 22
#define SHIFT_N 0
#define SHIFT_OFFSET 0
#define SHIFT_C 10
#define ID_SIG_BITS1 6
#define ID_SIG_BITS2 2
#define LENGTH_BITS 10
#define C_BITS 3
const uint16_t mask_id1 = static_cast<uint16_t>(~(((1U << ID_SIG_BITS1) - 1) << SHIFT_ID1));        // 11111100 00000000
const uint16_t mask_id2 = static_cast<uint16_t>(~(((1U << ID_SIG_BITS2) - 1) << SHIFT_ID2));        // 11000000 00000000
const uint16_t mask_offset =  static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << SHIFT_OFFSET));  // 00000011 11111111
const uint16_t mask_n =  static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << SHIFT_N));            // 00000011 11111111
const uint16_t mask_c =  static_cast<uint16_t>(~(((1U << C_BITS) - 1) << SHIFT_C));                 // 00011100 00000000
#define MAX_RUN_LENGTH  1023     // 2^10 - 1
#define BLOCK_SIZE      4194304  // 2^22 -- the actual value might be different
#define MAX_ALLOWED_BLOCKED_ID  16777215 // 2^24 - 1 -- the actual value might be different
#endif

#if MODE == 8
// Current implementation of MODE 8, does not allocate any bits in the offset for the id field.
// SHIFT_ID2, ID_SIG_BITS2, mask_id2 are only kept for consistency between MODE 2 and MODE 8,
// otherwise they can be removed for MODE 8.
#define SHIFT_ID1 10
#define SHIFT_ID2 16
#define SHIFT_ID1_RES 22
#define SHIFT_OFFSET 0
#define SHIFT_N 0
#define SHIFT_C 10
#define SHIFT_THRESHOLD_1 13
#define SHIFT_THRESHOLD_2 14
#define SHIFT_THRESHOLD_3 15
#define ID_SIG_BITS1 6
#define ID_SIG_BITS2 0
#define LENGTH_BITS 10
#define C_BITS 3
const uint16_t mask_id1 = static_cast<uint16_t>(~(((1U << ID_SIG_BITS1) - 1) << SHIFT_ID1));      // 11111100 00000000
const uint16_t mask_id2 = static_cast<uint16_t>(~(((1U << ID_SIG_BITS2) - 1) << SHIFT_ID2));      // 00000000 00000000
const uint16_t mask_offset = static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << SHIFT_OFFSET)); // 00000011 11111111
const uint16_t mask_n = static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << SHIFT_N));           // 00000011 11111111
const uint16_t mask_c = static_cast<uint16_t>(~(((1U << C_BITS) - 1) << SHIFT_C));                // 00011100 00000000
const uint16_t mask_thresholds1 = static_cast<uint16_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_1)); // 00100000 00000000
const uint16_t mask_thresholds2 = static_cast<uint16_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_2)); // 01000000 00000000
const uint16_t mask_thresholds3 = static_cast<uint16_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_3)); // 10000000 00000000
#define MAX_RUN_LENGTH          1023    // 2^10 - 1
#define BLOCK_SIZE              1048576 // 2^20 (the actual value might be different)
#define MAX_ALLOWED_BLOCKED_ID  4194303 // 2^22 - 1 (the actual value might be different)
#endif

#if MODE == 5
#define LENGTH_MASK_BITS 2
#define SHIFT_OFFSET 0
#define SHIFT_N 2
#define C_BITS 4
#define SHIFT_C 4
const uint8_t mask_offset =  static_cast<uint8_t>(~(((1U << LENGTH_MASK_BITS) - 1) << SHIFT_OFFSET));   // 00000011
const uint8_t mask_n =  static_cast<uint8_t>(~(((1U << LENGTH_MASK_BITS) - 1) << SHIFT_N));             // 00001100
const uint8_t mask_c = static_cast<uint8_t>(~(((1U << C_BITS) - 1) << SHIFT_C));                        // 11110000
#define MAX_RUN_LENGTH 1023 // 2^10-1
#endif

#if MODE == 7
#define LENGTH_MASK_BITS 1
#define SHIFT_OFFSET 0
#define SHIFT_N 1
#define SHIFT_THRESHOLD_1 5
#define SHIFT_THRESHOLD_2 6
#define SHIFT_THRESHOLD_3 7
#define C_BITS 3
#define SHIFT_C 2
const uint8_t mask_offset = static_cast<uint8_t>(~(((1U << LENGTH_MASK_BITS) - 1) << SHIFT_OFFSET)); // 00000001
const uint8_t mask_n = static_cast<uint8_t>(~(((1U << LENGTH_MASK_BITS) - 1) << SHIFT_N));           // 00000010
const uint8_t mask_c = static_cast<uint8_t>(~(((1U << C_BITS) - 1) << SHIFT_C));                     // 00011100
const uint8_t mask_thresholds1 = static_cast<uint8_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_1));      // 00100000
const uint8_t mask_thresholds2 = static_cast<uint8_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_2));      // 01000000
const uint8_t mask_thresholds3 = static_cast<uint8_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_3));      // 10000000
#define MAX_RUN_LENGTH 511 // 2^9-1
#endif

// #if MODE == 5 or MODE == 7
struct __attribute__((packed)) MoveTally {
    uint32_t right;
    uint8_t left;

    void set_value(uint64_t val) {
        if (val >= (static_cast<uint64_t>(1)<<40)) {
            std::cerr << "More than 40 bits are required for the id column.\n";
            std::cerr << val << "\n";
            exit(0);
        }
        left = 0;
        right = val;
        if (val >= (static_cast<uint64_t>(1) << 32)) {
            left = left | (val >> 32);
        }
    }

    uint64_t get() {
        uint64_t res = static_cast<uint64_t>(right);
        if (left == 0) {
            return res;
        } else {
            uint64_t left_64 = static_cast<uint64_t>(left);
            uint64_t left_shifted = left_64 << 32;
            res = res | left_shifted;
            return res;
        }
    }
};
// #endif

class __attribute__((packed)) MoveRow {
    public:
#if MODE == 0 or MODE == 1 or MODE == 4
        MoveRow () {n = 0; id = 0; overflow_bits = 0;}
#endif
#if MODE == 3 or MODE == 6 or MODE == 2 or MODE == 8
        MoveRow () {n = 0; id = 0; offset = 0;}
#endif
#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4 or MODE == 8 or MODE == 3 or MODE == 6
        MoveRow(uint16_t n_, uint16_t offset_, uint64_t id_);
        void init(uint16_t n_, uint16_t offset_, uint64_t id_);
#endif
#if MODE == 5 or MODE == 7
        MoveRow () {n = 0; offset = 0;}
        MoveRow(uint16_t n_, uint16_t offset_);
        void init(uint16_t n_, uint16_t offset_);
#endif
#if MODE == 2 or MODE == 8
        void print_all();
#endif
        friend std::ostream& operator<<(std::ostream& os, const MoveRow& mr);

        void set_n(uint16_t n_);
        void set_offset(uint16_t offset_);
        void set_c(char c_, std::vector<uint64_t>& alphamap);
        uint16_t get_n() const;
        uint16_t get_offset() const;
        char get_c() const;

#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4 or MODE == 8 or MODE == 3 or MODE == 6
        void set_id(uint64_t id_);
        uint64_t get_id() const;
#endif

#if MODE == 0 or MODE == 1 or MODE == 4
        uint8_t get_threshold_status(uint16_t i) const;
        void set_threshold_status(uint16_t i, uint8_t status);
        bool is_overflow_thresholds() const;
        uint16_t get_threshold() { return threshold; }
        void set_threshold(uint16_t t) { threshold = t; }

        void set_overflow_n();
        void set_overflow_offset();
        void set_overflow_thresholds();
        bool is_overflow_n() const;
        bool is_overflow_offset() const;
#endif

#if USE_NEXT_POINTERS
        uint16_t get_next_up(uint32_t i) { return next_up[i]; }
        uint16_t get_next_down(uint32_t i) { return next_down[i]; }
        void set_next_up(uint32_t i, uint16_t t) { next_up[i] = t; }
        void set_next_down(uint32_t i, uint16_t t) { next_down[i] = t; }
#endif

#if MODE == 8 or MODE == 7 or MODE == 6
        uint16_t get_threshold(uint16_t i) const;
        void set_threshold(uint16_t i, uint16_t value);
#endif
        uint64_t row_size() {
#if MODE == 0 or MODE == 4
            return 12;
#endif
#if MODE == 1
            return 24;
#endif
#if MODE == 3 or MODE == 6
            return 8;
#endif
#if MODE == 2 or MODE == 8
            return 6;
#endif
#if MODE == 5 or MODE == 7
            return 3;
#endif
        }
    // private:
#if MODE == 5 or MODE == 7
        uint8_t n; // length of the run
        uint8_t offset; // offset of the bwt row head of the current run in the new run after the LF-jump
        uint8_t c;
#endif

#if MODE == 2 or MODE == 8
        uint16_t id;        // The least significant bits of the bwt run after the LF-jump (distance from the block check point)
#endif
#if MODE == 0 or MODE == 1 or MODE == 4 or MODE == 3 or MODE == 6
        uint32_t id;        // The least significant bits of the bwt run after the LF-jump
#endif

#if COLOR_MODE == 1
        uint32_t color_id; // The color id of the run, only used in the colored-compressed mode
#endif

#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4 or MODE == 8 or MODE == 3 or MODE == 6
        uint16_t n;         // length of the run
        uint16_t offset;    // offset of the bwt row head of the current run in the new run after the LF-jump
#endif

#if MODE == 0 or MODE == 1 or MODE == 4
        uint16_t threshold;
        uint8_t overflow_bits;
        uint8_t thresholds_status; // Whether each threshold is at the boundary or it's a non-trivial value
#endif
#if MODE == 1
        // to store pointers for avoiding scanning
        uint16_t next_up[3];
        uint16_t next_down[3];
#endif
};

#if MODE == 0 or MODE == 1 or MODE == 4
inline uint8_t extract_value(uint8_t source, uint8_t mask, uint16_t shift) {
    uint8_t res = (source & (~mask)) >> shift;
    return res;
}

inline uint64_t MoveRow::get_id() const {
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

inline uint16_t MoveRow::get_n() const {
    return n;
}

inline uint16_t MoveRow::get_offset() const {
    return offset;
}

inline char MoveRow::get_c() const {
    return static_cast<char>(extract_value(thresholds_status, mask_c, 6));
}

inline bool MoveRow::is_overflow_n() const {
    uint8_t res = extract_value(overflow_bits, mask_overflow_n, 4);
    return !static_cast<bool>(res);
}

inline bool MoveRow::is_overflow_offset() const{
    uint8_t res = extract_value(overflow_bits, mask_overflow_offset, 5);
    return !static_cast<bool>(res);
}

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

#if MODE == 5 or MODE == 7
inline uint16_t MoveRow::get_n() const{
    uint16_t res = n;
    return res | (static_cast<uint16_t>((c & (~mask_n)) >> SHIFT_N) << 8);
}

inline uint16_t MoveRow::get_offset() const{
    uint16_t res = offset;
    return res | (static_cast<uint16_t>((c & (~mask_offset)) >> SHIFT_OFFSET) << 8);
}

inline char MoveRow::get_c() const{
    return static_cast<char>((c & (~mask_c)) >> SHIFT_C);
}
#endif

#if MODE == 2 or MODE == 8 or MODE == 3 or MODE == 6
inline uint16_t extract_value(uint16_t source, uint16_t mask, uint16_t shift) {
    uint16_t res = (source & (~mask)) >> shift;
    return res;
}
#endif

#if MODE == 3 or MODE == 6
inline uint64_t MoveRow::get_id() const {
    if (offset >= (1U << SHIFT_ID)) {
        uint64_t res = static_cast<uint64_t>(extract_value(offset, mask_id, SHIFT_ID));
        res = res << 32;
        uint64_t c = static_cast<uint64_t>(id);
        c = c | res;
        return c;
    } else {
        return static_cast<uint64_t>(id);
    }
}

inline uint16_t MoveRow::get_n() const {
    uint16_t res = static_cast<uint16_t>(extract_value(n, mask_n, 0));
    return res;
}

inline uint16_t MoveRow::get_offset() const {
    uint16_t res = static_cast<uint16_t>(extract_value(offset, mask_offset, 0));
    return res;
}

inline char MoveRow::get_c() const {
    return static_cast<char>(extract_value(n, mask_c, SHIFT_C));
}
#endif

#if MODE == 2 or MODE == 8
inline void MoveRow::print_all() {
    std::cerr << "id:\t" << std::bitset<16>(id) <<
                 "\nn:\t" << std::bitset<16>(n) <<
                 "\noffset:\t" << std::bitset<16>(offset) << "\n";
}

inline uint64_t MoveRow::get_id() const {
    if (n >= (1U << SHIFT_ID1) or offset >= (1U << SHIFT_ID2) ) {
        uint64_t res = 0;
        if (n >= (1U << SHIFT_ID1) ) {
            res = static_cast<uint64_t>(extract_value(n, mask_id1, SHIFT_ID1));
            res = res << 16;
        }
#if MODE == 2
        if (offset >= (1U << SHIFT_ID2) ) {
            uint64_t res2 = static_cast<uint64_t>(extract_value(offset, mask_id2, SHIFT_ID2));
            res2 = res2 << SHIFT_ID1_RES;
            res = res | res2;
        }
#endif
        uint64_t c = static_cast<uint64_t>(id);
        c = c | res;
        return c;
    } else {
        return static_cast<uint64_t>(id);
    }
}

inline uint16_t MoveRow::get_n() const{
    uint16_t res = static_cast<uint16_t>(extract_value(n, mask_n, SHIFT_N));
    return res;
}

inline uint16_t MoveRow::get_offset() const{
    uint16_t res = static_cast<uint16_t>(extract_value(offset, mask_offset, SHIFT_OFFSET));
    return res;
}

inline char MoveRow::get_c() const{
    return static_cast<char>(extract_value(offset, mask_c, SHIFT_C));
}
#endif

#if MODE == 6
inline uint16_t MoveRow::get_threshold(uint16_t i) const {
    switch (i) {
        case 0:
            return static_cast<uint16_t>((offset & (~mask_thresholds1)) >> SHIFT_THRESHOLD_1);
        case 1:
            return static_cast<uint16_t>((n & (~mask_thresholds2)) >> SHIFT_THRESHOLD_2);
        case 2:
            return static_cast<uint16_t>((n & (~mask_thresholds3)) >> SHIFT_THRESHOLD_3);
        default:
            std::cerr << "Only three thresholds exist per run: " << i << "\n";
            exit(0);
    }
}
#endif

#if MODE == 8
inline uint16_t MoveRow::get_threshold(uint16_t i) const {
    switch (i) {
        case 0:
            return static_cast<uint16_t>((offset & (~mask_thresholds1)) >> SHIFT_THRESHOLD_1);
        case 1:
            return static_cast<uint16_t>((offset & (~mask_thresholds2)) >> SHIFT_THRESHOLD_2);
        case 2:
            return static_cast<uint16_t>((offset & (~mask_thresholds3)) >> SHIFT_THRESHOLD_3);
        default:
            std::cerr << "Only three thresholds exist per run: " << i << "\n";
            exit(0);
    }
}
#endif

#if MODE == 7
inline uint16_t MoveRow::get_threshold(uint16_t i) const {
    switch (i) {
        case 0:
            return static_cast<uint16_t>((c & (~mask_thresholds1)) >> SHIFT_THRESHOLD_1);
        case 1:
            return static_cast<uint16_t>((c & (~mask_thresholds2)) >> SHIFT_THRESHOLD_2);
        case 2:
            return static_cast<uint16_t>((c & (~mask_thresholds3)) >> SHIFT_THRESHOLD_3);
        default:
            std::cerr << "Only three thresholds exist per run: " << i << "\n";
            exit(0);
    }
}
#endif

#endif // end of file