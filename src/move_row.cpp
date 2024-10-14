#include <cstdint>
#include <limits>
#include <iostream>

#include "move_row.hpp"

#if MODE == 0 or MODE == 1 or MODE == 2 or MDOE == 3 or MODE == 4
MoveRow::MoveRow(uint16_t n_, uint16_t offset_, uint64_t id_) {
    this->init(n_, offset_, id_);
}
#endif
#if MODE == 5
MoveRow::MoveRow(uint16_t n_, uint16_t offset_) {
    this->init(n_, offset_);
}
#endif

#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 3 or MODE == 4
void MoveRow::init(uint16_t n_, uint16_t offset_, uint64_t id_) {
    this->set_id(id_);
#endif
#if MODE == 5
void MoveRow::init(uint16_t n_, uint16_t offset_) {
#endif
#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4
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
#if MODE == 0 or MODE == 1 or MODE == 2 or MDOE == 3 or MODE == 4
    os << "n:" << mr.get_n() <<  " offset: " << mr.get_offset() << " id:" << mr.get_id();
#endif
#if MODE == 5
    os << "n:" << mr.get_n() <<  " offset: " << mr.get_offset();
#endif
    return os;
}

#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4
void MoveRow::set_threshold_status(uint16_t i, uint8_t status) {
    const uint8_t mask_thresholds = static_cast<uint8_t>(~(((1U << 2) - 1) << i*2));
    thresholds_status = thresholds_status & mask_thresholds;
    thresholds_status = thresholds_status |  (status << i*2);
}

void MoveRow::set_overflow_thresholds() {
    overflow_bits = overflow_bits & mask_overflow_thresholds;
    overflow_bits = overflow_bits | (1 >> 6);
}
#endif

void MoveRow::set_n(uint16_t n_) {
#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4
    n = n_;
#endif
#if MODE == 3
    n = n & mask_n;
    if (n_ < MAX_RUN_LENGTH)
        n = n | n_;
    else {
        std::cerr << "The length is greater than 2^12: " << n_ << "\n";
        exit(0);
    }
#endif
#if MODE == 5
    n = static_cast<uint8_t>(n_);
    if (n_ <= MAX_RUN_LENGTH) {
        c = c & mask_n;
        if (n_ >= 256) {
            uint8_t n_8 = n_ >> 8;
            c = c | (n_8 << 2);
        }
    } else {
        std::cerr << "The length is greater than 2^12: " << n_ << "\n";
        exit(0);
    }
#endif
}

void MoveRow::set_overflow_n() {
#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4
    overflow_bits = overflow_bits & mask_overflow_n;
    overflow_bits = overflow_bits | (1 >> 4);
#endif
#if MODE == 3
    std::cerr << "The length overflow should not occur in the compressed mode.\n";
    exit(0);
#endif
}

void MoveRow::set_offset(uint16_t offset_) {
#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4
    offset = offset_;
#endif
#if MODE == 3
    offset = offset & mask_offset;
    if (offset_ < MAX_RUN_LENGTH)
        offset = offset | offset_;
    else {
        std::cerr << "The offset is greater than 2^12: " << offset_ << "\n";
        exit(0);
    }
#endif
#if MODE == 5
    offset = static_cast<uint8_t>(offset_);
    if (offset_ <= MAX_RUN_LENGTH) {
        uint8_t offset_8 = offset_ >> 8;
        if (offset_ >= 256) {
            c = c & mask_offset;
            c = c | offset_8;
        }
    }
    else {
        std::cerr << "The offset is greater than 2^12: " << offset_ << "\n";
        exit(0);
    }
#endif
}

void MoveRow::set_overflow_offset() {
#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4
    overflow_bits = overflow_bits & mask_overflow_offset;
    overflow_bits = overflow_bits | (1 >> 5);
#endif
#if MODE == 3
    std::cerr << "The offset overflow should not occur in the compressed mode.\n";
    exit(0);
#endif
}


#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 3 or MODE == 4
void MoveRow::set_id(uint64_t id_) {
    id = id_; // Store the least significant bits in the didicated id variable
#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4
    overflow_bits = overflow_bits & mask_id;
    overflow_bits = overflow_bits | (id_ >> 32);
#endif
#if MODE == 3
    offset = offset & mask_id;
    offset = offset | (id_ >> 32);
#endif
}
#endif

void MoveRow::set_c(char c_, std::vector<uint64_t>& alphamap) {
#if MODE == 0 or MODE == 1 or MODE == 2 or MODE == 4
    uint64_t c_64 = static_cast<uint64_t>(alphamap[c_]);
    thresholds_status = thresholds_status & mask_c;
    thresholds_status = thresholds_status | (c_64 << 6);
#endif
#if MODE == 3
    uint16_t c_16 = static_cast<uint16_t>(alphamap[c_]);
    n = n & mask_c;
    n = n | (c_16 << 12);
#endif
#if MODE == 5
    uint8_t c_8 = static_cast<uint8_t>(alphamap[c_]);
    c = c & mask_c;
    c = c | (c_8 << 4);
#endif
}