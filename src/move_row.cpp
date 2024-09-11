#include <cstdint>
#include <limits>
#include <iostream>

#include "move_row.hpp"

MoveRow::MoveRow(uint16_t n_, uint16_t offset_, uint64_t id_) {
    this->init(n_, offset_, id_);
}
void MoveRow::init(uint16_t n_, uint16_t offset_, uint64_t id_) {
#if MODE == 0 or MODE == 1 or MODE == 2
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
    this->set_id(id_);
}

std::ostream& operator<<(std::ostream& os, const MoveRow& mr)
{
    os << "n:" << mr.get_n() <<  " offset: " << mr.get_offset() << " id:" << mr.get_id();
    return os;
}

#if MODE == 0 or MODE == 1 or MODE == 2
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
#if MODE == 0 or MODE == 1 or MODE == 2
    n = n_;
#endif
#if MODE == 3
    if (n_ <= MAX_RUN_LENGTH) {
        n = n & mask_n;
        n = n | (n_ << shift_n);
    } else {
        std::cerr << "The length is greater than 2^12: " << n_ << "\n";
        exit(0);
    }
#endif
}

void MoveRow::set_overflow_n() {
#if MODE == 0 or MODE == 1 or MODE == 2
    overflow_bits = overflow_bits & mask_overflow_n;
    overflow_bits = overflow_bits | (1 >> 4);
#endif
#if MODE == 3
    std::cerr << "The length overflow should not occur in the compressed mode.\n";
    exit(0);
#endif
}

void MoveRow::set_offset(uint16_t offset_) {
#if MODE == 0 or MODE == 1 or MODE == 2
    offset = offset_;
#endif
#if MODE == 3
    if (offset_ <= MAX_RUN_LENGTH) {
        offset = offset & mask_offset;
        offset = offset | (offset_ << shift_offset);
    } else {
        std::cerr << "The offset is greater than 2^12: " << offset_ << "\n";
        exit(0);
    }
#endif
}

void MoveRow::set_overflow_offset() {
#if MODE == 0 or MODE == 1 or MODE == 2
    overflow_bits = overflow_bits & mask_overflow_offset;
    overflow_bits = overflow_bits | (1 >> 5);
#endif
#if MODE == 3
    std::cerr << "The offset overflow should not occur in the compressed mode.\n";
    exit(0);
#endif
}

void MoveRow::set_id(uint64_t id_) {
    id = id_; // Store the least significant bits in the didicated id variable
#if MODE == 0 or MODE == 1 or MODE == 2
    overflow_bits = overflow_bits & mask_id;
    overflow_bits = overflow_bits | (id_ >> 32);
#endif
#if MODE == 3
    n = n & mask_id1;
    n = n | ((id_ >> 16) << shift_id1);
    offset = offset & mask_id2;
    offset = offset | ((id_ >> 21) << shift_id2);
#endif
}

void MoveRow::set_c(char c_, std::vector<uint64_t>& alphamap) {
#if MODE == 0 or MODE == 1 or MODE == 2
    uint64_t c_64 = static_cast<uint64_t>(alphamap[c_]);
    thresholds_status = thresholds_status & mask_c;
    thresholds_status = thresholds_status | (c_64 << 6);
#endif
#if MODE == 3
    uint16_t c_16 = static_cast<uint16_t>(alphamap[c_]);
    offset = offset & mask_c;
    offset = offset | (c_16 << shift_c);
#endif
}