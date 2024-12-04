#ifndef __MOVE_ROW__
#define __MOVE_ROW__

#include <iostream>
#include <vector>
#include <bitset>

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
#define ID_SIG_BITS 4
#define LENGTH_BITS 11
#define C_BITS 3
const uint16_t mask_id =  static_cast<uint16_t>(~(((1U << ID_SIG_BITS) - 1) << SHIFT_ID)); // 11110000 00000000
const uint16_t mask_offset =  static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << 0));    // 00000111 11111111
const uint16_t mask_n =  static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << 0));         // 00000111 11111111
const uint16_t mask_c = static_cast<uint16_t>(~(((1U << C_BITS) - 1) << SHIFT_C));         // 11100000 00000000
const uint16_t mask_thresholds1 =  static_cast<uint16_t>(~(((1U << 1) - 1) << 11));        // 00001000 00000000
const uint16_t mask_thresholds2 =  static_cast<uint16_t>(~(((1U << 1) - 1) << 11));        // 00001000 00000000
const uint16_t mask_thresholds3 = static_cast<uint16_t>(~(((1U << 1) - 1) << 12));         // 00010000 00000000
#define MAX_RUN_LENGTH 2047 // 2^11-1
#endif

class MoveRow {
    public:
#if MODE == 0 or MODE == 1 or MODE == 4
        MoveRow () {n = 0; id = 0; overflow_bits = 0;}
#endif
#if MODE == 3 or MODE == 6
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

#if MODE == 1
        uint16_t get_next_up(uint32_t i) { return next_up[i]; }
        uint16_t get_next_down(uint32_t i) { return next_down[i]; }
        void set_next_up(uint32_t i, uint16_t t) { next_up[i] = t; }
        void set_next_down(uint32_t i, uint16_t t) { next_down[i] = t; }
#endif

#if MODE == 6
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
        }
    private:
        uint32_t id; // bwt run after the LF-jump
        uint16_t n; // length of the run
        uint16_t offset; // offset of the bwt row head of the current run in the new run after the LF-jump

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

#if MODE == 3 or MODE == 6
inline uint16_t extract_value(uint16_t source, uint16_t mask, uint16_t shift) {
    uint16_t res = (source & (~mask)) >> shift;
    return res;
}

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

#if MODE == 6
inline uint16_t MoveRow::get_threshold(uint16_t i) const {
    switch (i) {
        case 0:
            return static_cast<uint16_t>((offset & (~mask_thresholds1)) >> 11);
        case 1:
            return static_cast<uint16_t>((n & (~mask_thresholds2)) >> 11);
        case 2:
            return static_cast<uint16_t>((n & (~mask_thresholds3)) >> 12);
        default:
            std::cerr << "Only three thresholds exist per run: " << i << "\n";
            exit(0);
    }
}
#endif

#endif // end of file