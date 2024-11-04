#include <cstdint>
#include <limits>
#include <iostream>

#include "move_row.hpp"

#if MODE == 0 or MODE == 1 or MDOE == 3 or MODE == 4 or MODE == 6
MoveRow::MoveRow(uint16_t n_, uint16_t offset_, uint64_t id_) {
    this->init(n_, offset_, id_);
}
#endif
#if MODE == 5 or MODE == 7
MoveRow::MoveRow(uint16_t n_, uint16_t offset_) {
    this->init(n_, offset_);
}
#endif

#if MODE == 0 or MODE == 1 or MODE == 3 or MODE == 4 or MODE == 6
void MoveRow::init(uint16_t n_, uint16_t offset_, uint64_t id_) {
    this->set_id(id_);
#endif
#if MODE == 5 or MODE == 7
void MoveRow::init(uint16_t n_, uint16_t offset_) {
#endif
#if MODE == 0 or MODE == 1 or MODE == 4
    overflow_bits = std::numeric_limits<uint8_t>::max();
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
#if MODE == 0 or MODE == 1 or MDOE == 3 or MODE == 4 or MODE == 6
    os << "n:" << mr.get_n() <<  " offset: " << mr.get_offset() << " id:" << mr.get_id();
#endif
#if MODE == 5 or MODE == 7
    os << "n:" << mr.get_n() <<  " offset: " << mr.get_offset();
#endif
    return os;
}

#if MODE == 0 or MODE == 1 or MODE == 4
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

#if MODE == 5 or MODE == 7
void MoveRow::set_n(uint16_t n_) {
    if (n_ <= MAX_RUN_LENGTH) {
        n = static_cast<uint8_t>(n_);
        c = c & mask_n;
        if (n_ >= 256) {
            uint8_t n_8 = n_ >> 8;
            c = c | (n_8 << SHIFT_N);
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
            uint8_t offset_8 = offset_ >> 8;
            c = c | (offset_8 << SHIFT_OFFSET);
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
#endif

#if MODE == 3 or MODE == 6
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
    offset = offset & mask_id2;
    //Note: SHIFT_ID1_RES = 16 + ID_SIG_BITS1
    offset = offset | ((id_ >> (16 + ID_SIG_BITS1)) << SHIFT_ID2);
}

void MoveRow::set_c(char c_, std::vector<uint64_t>& alphamap) {
    uint16_t c_16 = static_cast<uint16_t>(alphamap[c_]);
    offset = offset & mask_c;
    offset = offset | (c_16 << SHIFT_C );
}
#endif

#if MODE == 6
void MoveRow::set_threshold(uint16_t i, uint16_t value) {
    if (value > 1) {
        std::cerr << "The theshold may be either 0 or 1: " << i << "\n";
        exit(0);
    }
    switch (i) {
        case 0:
            offset = offset & mask_thresholds1;
            offset = offset | (value << 12);
            break;
        case 1:
            offset = offset & mask_thresholds2;
            offset = offset | (value << 13);
            break;
        case 2:
            offset = offset & mask_thresholds3;
            offset = offset | (value << 14);
            break;
        default:
            std::cerr << "Only three thresholds may be stored: " << i << "\n";
            exit(0);
    }
}
#endif

#if MODE == 7
void MoveRow::set_threshold(uint16_t i, uint16_t value) {
    if (value > 1) {
        std::cerr << "The theshold may be either 0 or 1: " << i << "\n";
        exit(0);
    }
    switch (i) {
        case 0:
            c = c & mask_thresholds1;
            c = c | (value << 5);
            break;
        case 1:
            c = c & mask_thresholds2;
            c = c | (value << 6);
            break;
        case 2:
            c = c & mask_thresholds3;
            c = c | (value << 7);
            break;
        default:
            std::cerr << "Only three thresholds may be stored: " << i << "\n";
            exit(0);
    }
}
#endif