#include <cstdint>
#include <limits>
#include <iostream>

#include "move_row.hpp"

#if TALLY_MODE
MoveRow::MoveRow(uint16_t n_, uint16_t offset_) {
    this->init(n_, offset_);
}
#else
MoveRow::MoveRow(uint16_t n_, uint16_t offset_, uint64_t id_) {
    this->init(n_, offset_, id_);
}
#endif

#if TALLY_MODE
void MoveRow::init(uint16_t n_, uint16_t offset_) {
#else
void MoveRow::init(uint16_t n_, uint16_t offset_, uint64_t id_) {
#if SPLIT_THRESHOLDS_FALSE
    // This should be set before calling other setters.
    overflow_bits = std::numeric_limits<uint8_t>::max();
#endif
    // TODO: In the blocked mode, the following (set_id) is not needed here (will be set later)
    this->set_id(id_);
    if (id_ != this->get_id()) {
        std::cerr << "The id setter or getter is not working properly.\n";
        std::cerr << id_ << " " << this->get_id() << "\n";
        exit(0);
    }
#endif
    this->set_n(n_);
    if (n_ != this->get_n()) {
        std::cerr << "The length setter or getter is not working properly.\n";
        std::cerr << n_ << " " << this->get_n() << "\n";
        exit(0);
    }
    this->set_offset(offset_);
    if (offset_ != this->get_offset()) {
        std::cerr << "The offset setter or getter is not working properly.\n";
        std::cerr << offset_ << " " << this->get_offset() << "\n";
        exit(0);
    }
}

std::ostream& operator<<(std::ostream& os, const MoveRow& mr)
{
#if TALLY_MODE
    os << "n:" << mr.get_n() <<  " offset: " << mr.get_offset();
#else
    os << "n:" << mr.get_n() <<  " offset: " << mr.get_offset() << " id:" << mr.get_id();
#endif
    return os;
}

#if SPLIT_THRESHOLDS_FALSE
void MoveRow::set_threshold_status(uint16_t i, uint8_t status) {
    const uint8_t mask_thresholds = static_cast<uint8_t>(~(((1U << 2) - 1) << i*2));
    thresholds_status = thresholds_status & mask_thresholds;
    thresholds_status = thresholds_status |  (status << i*2);
}

void MoveRow::set_overflow_thresholds() {
    overflow_bits = overflow_bits & mask_overflow_thresholds;
    overflow_bits = overflow_bits | (1U >> 6);
}

void MoveRow::set_n(uint16_t n_) {
    n = n_;
}

void MoveRow::set_overflow_n() {
    overflow_bits = overflow_bits & mask_overflow_n;
    overflow_bits = overflow_bits | (1U >> 4);
}

void MoveRow::set_offset(uint16_t offset_) {
    offset = offset_;
}

void MoveRow::set_overflow_offset() {
    overflow_bits = overflow_bits & mask_overflow_offset;
    overflow_bits = overflow_bits | (1U >> 5);
}


void MoveRow::set_id(uint64_t id_) {
    id = id_; // Store the least significant bits in the didicated id variable
    overflow_bits = overflow_bits & mask_id;
    overflow_bits = overflow_bits | (id_ >> 32);
}

void MoveRow::set_c(char c_, std::vector<uint64_t>& alphamap) {
    uint64_t c_64 = static_cast<uint64_t>(alphamap[c_]);
    thresholds_status = thresholds_status & mask_c;
    thresholds_status = thresholds_status | (c_64 << 6);
}
#endif

#if COMPACT_MODE
void MoveRow::set_id(uint64_t id_) {
    id = id_; // Store the least significant bits in the didicated id variable
    offset = offset & mask_id;
    offset = offset | ((id_ >> 32) << SHIFT_ID);
}

void MoveRow::set_n(uint16_t n_) {
    n = n & mask_n;
    if (n_ < (1U << LENGTH_BITS)) {
        n = n | n_;
    } else {
        std::cerr << "The length is greater than 2^" << LENGTH_BITS  << ": " << n_ << "\n";
        exit(0);
    }
}

void MoveRow::set_offset(uint16_t offset_) {
    offset = offset & mask_offset;
    if (offset_ < (1U << LENGTH_BITS)) {
        offset = offset | offset_;
    } else {
        std::cerr << "The offset is greater than 2^" << LENGTH_BITS  << ": " << offset_ << "\n";
        exit(0);
    }
}

void MoveRow::set_c(char c_, std::vector<uint64_t>& alphamap) {
    uint16_t c_16 = static_cast<uint16_t>(alphamap[c_]);
    n = n & mask_c;
    n = n | (c_16 << SHIFT_C);
}

#if USE_THRESHOLDS
void MoveRow::set_threshold(uint16_t i, uint16_t value) {
    if (value > 1) {
        std::cerr << "The theshold may be either 0 or 1: " << i << "\n";
        exit(0);
    }
    switch (i) {
        case 0:
            offset = offset & mask_thresholds1;
            offset = offset | (value << SHIFT_THRESHOLD_1);
            break;
        case 1:
            n = n & mask_thresholds2;
            n = n | (value << SHIFT_THRESHOLD_2);
            break;
        case 2:
            n = n & mask_thresholds3;
            n = n | (value << SHIFT_THRESHOLD_3);
            break;
        default:
            std::cerr << "Only three thresholds may be stored: " << i << "\n";
            exit(0);
    }
}
#endif
#endif

#if BLOCKED_MODE
void MoveRow::set_n(uint16_t n_) {
    if (n_ <= MAX_RUN_LENGTH) {
        n = n & mask_n;
        n = n | (n_ << SHIFT_N);
    } else {
        std::cerr << "The length is greater than 2^12: " << n_ << "\n";
        exit(0);
    }
}

void MoveRow::set_offset(uint16_t offset_) {
    if (offset_ <= MAX_RUN_LENGTH) {
        offset = offset & mask_offset;
        offset = offset | (offset_ << SHIFT_OFFSET);
    } else {
        std::cerr << "The offset is greater than 2^12: " << offset_ << "\n";
        exit(0);
    }
}

// example
// 00000000 00000000 00000000 11000011 11000011 11000011 00000000 00000000
// >> 16
// 00000000 00000000 00000000 00000000 00000000 11000000 11000011 11000011
// << 10
// 00000000 00000000 00000000 00000110 00000011 00001111 00001100 00000000
//                                                       11111100 00000000

// 00000000 00000000 00000000 11000011 11000011 11000011 00000000 00000000
// >> 22
// 00000000 00000000 00000000 00000000 00000000 00000011 00001111 00001111
// << 14
// 00000000 00000000 00000000 00000000 11000011 11000011 110000000 0000000
//                                                       11000000 00000000
void MoveRow::set_id(uint64_t id_) {
    id = id_; // Store the least significant bits in the didicated id variable

    // We don't have to do anything now, because it is checked before calling set_id
    // that the id_ is smaller than or equal to MAX_BLOCKED_ID
    // This function is also called during the first initialization of the row,
    // the id is greater than MAX_BLOCKED_ID then, but it is Ok, as it will be
    // set to the correct value later in the code.

    n = n & mask_id1;
    n = n | ((id_ >> 16) << SHIFT_ID1);
#if MODE == 2
    offset = offset & mask_id2;
    //Note: SHIFT_ID1_RES = 16 + ID_SIG_BITS1
    offset = offset | ((id_ >> (16 + ID_SIG_BITS1)) << SHIFT_ID2);
#endif
}

void MoveRow::set_c(char c_, std::vector<uint64_t>& alphamap) {
    uint16_t c_16 = static_cast<uint16_t>(alphamap[c_]);
    offset = offset & mask_c;
    offset = offset | (c_16 << SHIFT_C );
}

#if USE_THRESHOLDS
void MoveRow::set_threshold(uint16_t i, uint16_t value) {
    if (value > 1) {
        std::cerr << "The theshold may be either 0 or 1: " << i << "\n";
        exit(0);
    }
    switch (i) {
        case 0:
            offset = offset & mask_thresholds1;
            offset = offset | (value << SHIFT_THRESHOLD_1);
            break;
        case 1:
            offset = offset & mask_thresholds2;
            offset = offset | (value << SHIFT_THRESHOLD_2);
            break;
        case 2:
            offset = offset & mask_thresholds3;
            offset = offset | (value << SHIFT_THRESHOLD_3);
            break;
        default:
            std::cerr << "Only three thresholds may be stored: " << i << "\n";
            exit(0);
    }
}
#endif
#endif

#if TALLY_MODE
void MoveRow::set_n(uint16_t n_) {
    if (n_ <= MAX_RUN_LENGTH) {
        n = static_cast<uint8_t>(n_);
        c = c & mask_n;
        if (n_ >= 256) {
            uint8_t n_8 = (n_ >> 8) & static_cast<uint8_t>(MASK_BYTE);
            c = c | (n_8 << SHIFT_N);
#if BYTE
            uint8_t n_high = (n_ >> SHIFT_BYTE);
            l = l | (n_high << SHIFT_LN);
#endif
        }
    } else {
        std::cerr << "The length is greater than 2^12: " << n_ << "\n";
        exit(0);
    }
}

void MoveRow::set_offset(uint16_t offset_) {
    if (offset_ <= MAX_RUN_LENGTH) {
        offset = static_cast<uint8_t>(offset_);
        c = c & mask_offset;
        if (offset_ >= 256) {
            uint8_t offset_8 = (offset_ >> 8) & static_cast<uint8_t>(MASK_BYTE);
            c = c | (offset_8 << SHIFT_OFFSET);
#if BYTE
            uint8_t offset_high = (offset_ >> SHIFT_BYTE);
            l = l | (offset_high << SHIFT_LO);
#endif
        }
    }
    else {
        std::cerr << "The offset is greater than 2^12: " << offset_ << "\n";
        exit(0);
    }
}

void MoveRow::set_c(char c_, std::vector<uint64_t>& alphamap) {
    uint8_t c_8 = static_cast<uint8_t>(alphamap[c_]);
    c = c & mask_c;
    c = c | (c_8 << SHIFT_C);
}

#if USE_THRESHOLDS
void MoveRow::set_threshold(uint16_t i, uint16_t value) {
    if (value > 1) {
        std::cerr << "The theshold may be either 0 or 1: " << i << "\n";
        exit(0);
    }
    switch (i) {
        case 0:
            c = c & mask_thresholds1;
            c = c | (value << SHIFT_THRESHOLD_1);
            break;
        case 1:
            c = c & mask_thresholds2;
            c = c | (value << SHIFT_THRESHOLD_2);
            break;
        case 2:
            c = c & mask_thresholds3;
            c = c | (value << SHIFT_THRESHOLD_3);
            break;
        default:
            std::cerr << "Only three thresholds may be stored: " << i << "\n";
            exit(0);
    }
}
#endif
#endif