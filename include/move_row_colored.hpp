// This file defines a class MoveRowColored, which is just a duplicate of
// MoveRow (in move_row.hpp) with the addition of color_id
// It's only for the purpose of testing and not still used in the main version

#ifndef MOVE_ROW_COLORED_HPP
#define MOVE_ROW_COLORED_HPP

#include <iostream>
#include <vector>
#include <bitset>

#include "utils.hpp"
#include "move_row.hpp"

// removed all the MODE based definitions which already exist in move_row.hpp
// including the definitions of extract_value below

class __attribute__((packed)) MoveRowColored {
    public:
#if MOVI1_STYLE
        MoveRowColored () {n = 0; id = 0; overflow_bits = 0;}
#endif
#if REGULAR_MODES or BLOCKED_MODES
        MoveRowColored () {n = 0; id = 0; offset = 0;}
#endif
#if NO_SAMPPLED_ID
        MoveRowColored(uint16_t n_, uint16_t offset_, uint64_t id_);
        void init(uint16_t n_, uint16_t offset_, uint64_t id_);
#endif
#if TALLY_MODES
        MoveRowColored () {n = 0; offset = 0;}
        MoveRowColored(uint16_t n_, uint16_t offset_);
        void init(uint16_t n_, uint16_t offset_);
#endif
#if BLOCKED_MODES
        void print_all();
#endif
        friend std::ostream& operator<<(std::ostream& os, const MoveRowColored& mr);

        void set_n(uint16_t n_);
        void set_offset(uint16_t offset_);
        void set_c(char c_, std::vector<uint64_t>& alphamap);
        uint16_t get_n() const;
        uint16_t get_offset() const;
        char get_c() const;

#if NO_SAMPPLED_ID
        void set_id(uint64_t id_);
        uint64_t get_id() const;
#endif

#if MOVI1_STYLE
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

#if SPLIT_THRESHOLDS_TRUE
        uint16_t get_threshold(uint16_t i) const;
        void set_threshold(uint16_t i, uint16_t value);
#endif
        uint64_t row_size() {
#if LARGE_INDEX or SPLIT_INDEX
            return 12;
#endif
#if CONSTANT_INDEX
            return 24;
#endif
#if REGULAR_INDEX or REGULAR_THRESHOLDS_INDEX
            return 8;
#endif
#if BLOCKED_INDEX or BLOCKED_THRESHOLDS_INDEX
            return 6;
#endif
#if TALLY_INDEX or TALLY_THRESHOLDS_INDEX
            return 3;
#endif
        }

        uint32_t get_color_id() const { return color_id; }
        void set_color_id(uint32_t color_id_) { color_id = color_id_; }

    // private:
#if TALLY_MODES
        uint8_t n; // length of the run
        uint8_t offset; // offset of the bwt row head of the current run in the new run after the LF-jump
        uint8_t c;
#endif

#if BLOCKED_MODES
        uint16_t id;        // The least significant bits of the bwt run after the LF-jump (distance from the block check point)
#endif
#if NO_EXTRA_TABLE
        uint32_t id;        // The least significant bits of the bwt run after the LF-jump
#endif
        uint32_t color_id;  // The color id of the bwt run
#if NO_SAMPPLED_ID
        uint16_t n;         // length of the run
        uint16_t offset;    // offset of the bwt row head of the current run in the new run after the LF-jump
#endif
  
#if MOVI1_STYLE
        uint16_t threshold;
        uint8_t overflow_bits;
        uint8_t thresholds_status; // Whether each threshold is at the boundary or it's a non-trivial value
#endif
#if CONSTANT_INDEX
        // to store pointers for avoiding scanning
        uint16_t next_up[3];
        uint16_t next_down[3];
#endif
};

#if MOVI1_STYLE
// inline uint8_t extract_value(uint8_t source, uint8_t mask, uint16_t shift) {
//     uint8_t res = (source & (~mask)) >> shift;
//     return res;
// }

inline uint64_t MoveRowColored::get_id() const {
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

inline uint16_t MoveRowColored::get_n() const {
    return n;
}

inline uint16_t MoveRowColored::get_offset() const {
    return offset;
}

inline char MoveRowColored::get_c() const {
    return static_cast<char>(extract_value(thresholds_status, mask_c, 6));
}

inline bool MoveRowColored::is_overflow_n() const {
    uint8_t res = extract_value(overflow_bits, mask_overflow_n, 4);
    return !static_cast<bool>(res);
}

inline bool MoveRowColored::is_overflow_offset() const{
    uint8_t res = extract_value(overflow_bits, mask_overflow_offset, 5);
    return !static_cast<bool>(res);
}

inline uint8_t MoveRowColored::get_threshold_status(uint16_t i) const {
    const uint8_t mask_thresholds = static_cast<uint8_t>(~(((1U << 2) - 1) << i*2));
    uint8_t status = static_cast<uint8_t>((thresholds_status & (~mask_thresholds)) >> i*2);
    return status;
}

inline bool MoveRowColored::is_overflow_thresholds() const{
    uint8_t res = extract_value(overflow_bits, mask_overflow_thresholds, 6);
    return !static_cast<bool>(res);
}
#endif

#if TALLY_MODES
inline uint16_t MoveRowColored::get_n() const{
    uint16_t res = n;
    return res | (static_cast<uint16_t>((c & (~mask_n)) >> SHIFT_N) << 8);
}

inline uint16_t MoveRowColored::get_offset() const{
    uint16_t res = offset;
    return res | (static_cast<uint16_t>((c & (~mask_offset)) >> SHIFT_OFFSET) << 8);
}

inline char MoveRowColored::get_c() const{
    return static_cast<char>((c & (~mask_c)) >> SHIFT_C);
}
#endif

#if REGULAR_MODES or BLOCKED_MODES
// inline uint16_t extract_value(uint16_t source, uint16_t mask, uint16_t shift) {
//     uint16_t res = (source & (~mask)) >> shift;
//     return res;
// }
#endif

#if REGULAR_MODES
inline uint64_t MoveRowColored::get_id() const {
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

inline uint16_t MoveRowColored::get_n() const {
    uint16_t res = static_cast<uint16_t>(extract_value(n, mask_n, 0));
    return res;
}

inline uint16_t MoveRowColored::get_offset() const {
    uint16_t res = static_cast<uint16_t>(extract_value(offset, mask_offset, 0));
    return res;
}

inline char MoveRowColored::get_c() const {
    return static_cast<char>(extract_value(n, mask_c, SHIFT_C));
}
#endif

#if BLOCKED_MODES
inline void MoveRowColored::print_all() {
    INFO_MSG("id:\t" + std::bitset<16>(id) + " n:\t" + std::bitset<16>(n) + " offset:\t" + std::bitset<16>(offset));
}

inline uint64_t MoveRowColored::get_id() const {
    if (n >= (1U << SHIFT_ID1) or offset >= (1U << SHIFT_ID2) ) {
        uint64_t res = 0;
        if (n >= (1U << SHIFT_ID1) ) {
            res = static_cast<uint64_t>(extract_value(n, mask_id1, SHIFT_ID1));
            res = res << 16;
        }
        if (offset >= (1U << SHIFT_ID2) ) {
            uint64_t res2 = static_cast<uint64_t>(extract_value(offset, mask_id2, SHIFT_ID2));
            res2 = res2 << SHIFT_ID1_RES;
            res = res | res2;
        }
        uint64_t c = static_cast<uint64_t>(id);
        c = c | res;
        return c;
    } else {
        return static_cast<uint64_t>(id);
    }
}

inline uint16_t MoveRowColored::get_n() const{
    uint16_t res = static_cast<uint16_t>(extract_value(n, mask_n, SHIFT_N));
    return res;
}

inline uint16_t MoveRowColored::get_offset() const{
    uint16_t res = static_cast<uint16_t>(extract_value(offset, mask_offset, SHIFT_OFFSET));
    return res;
}

inline char MoveRowColored::get_c() const{
    return static_cast<char>(extract_value(offset, mask_c, SHIFT_C));
}
#endif

#if REGULAR_THRESHOLDS_INDEX
inline uint16_t MoveRowColored::get_threshold(uint16_t i) const {
    switch (i) {
        case 0:
            return static_cast<uint16_t>((offset & (~mask_thresholds1)) >> 11);
        case 1:
            return static_cast<uint16_t>((n & (~mask_thresholds2)) >> 11);
        case 2:
            return static_cast<uint16_t>((n & (~mask_thresholds3)) >> 12);
        default:
            throw std::runtime_error(ERROR_MSG("[MoveRowColored - get_threshold] Only three thresholds exist per run: " + std::to_string(i) + "\n"));
    }
}
#endif

#if BLOCKED_THRESHOLDS_INDEX
inline uint16_t MoveRowColored::get_threshold(uint16_t i) const {
    switch (i) {
        case 0:
            return static_cast<uint16_t>((offset & (~mask_thresholds1)) >> 12);
        case 1:
            return static_cast<uint16_t>((offset & (~mask_thresholds2)) >> 13);
        case 2:
            return static_cast<uint16_t>((offset & (~mask_thresholds3)) >> 14);
        default:
            throw std::runtime_error(ERROR_MSG("[MoveRowColored - get_threshold] Only three thresholds exist per run: " + std::to_string(i) + "\n"));
    }
}
#endif

#if TALLY_THRESHOLDS_INDEX
inline uint16_t MoveRowColored::get_threshold(uint16_t i) const {
    switch (i) {
        case 0:
            return static_cast<uint16_t>((c & (~mask_thresholds1)) >> 5);
        case 1:
            return static_cast<uint16_t>((c & (~mask_thresholds2)) >> 6);
        case 2:
            return static_cast<uint16_t>((c & (~mask_thresholds3)) >> 7);
        default:
            throw std::runtime_error(ERROR_MSG("[MoveRowColored - get_threshold] Only three thresholds exist per run: " + std::to_string(i) + "\n"));
    }
}
#endif

#endif // end of file