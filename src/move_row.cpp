#include <cstdint>
#include <limits>
#include <iostream>

#include "move_row.hpp"

MoveRow::MoveRow(uint16_t n_, uint16_t offset_, uint64_t id_) {
    this->init(n_, offset_, id_);
}
void MoveRow::init(uint16_t n_, uint16_t offset_, uint64_t id_) {
    overflow_bits = std::numeric_limits<uint16_t>::max();
    this->set_n(n_);
    this->set_offset(offset_);
    this->set_id(id_);
}

std::ostream& operator<<(std::ostream& os, const MoveRow& mr)
{
    os << "n:" << mr.get_n() <<  " offset: " << mr.get_offset() << " id:" << mr.get_id();
    return os;
}

void MoveRow::set_n(uint16_t n_) {
    n = n_;
}

void MoveRow::set_threshold_status(uint16_t i, uint8_t status) {
    const uint16_t mask_thresholds = static_cast<uint16_t>(~(((1U << 2) - 1) << i*2));
    thresholds_status = thresholds_status & mask_thresholds;
    thresholds_status = thresholds_status |  (status << i*2);
}

void MoveRow::set_overflow_n() {
    overflow_bits = overflow_bits & mask_overflow_n;
    overflow_bits = overflow_bits | (1 >> 10);
}

void MoveRow::set_offset(uint16_t offset_) {
    offset = offset_;
}

void MoveRow::set_overflow_offset() {
    overflow_bits = overflow_bits & mask_overflow_offset;
    overflow_bits = overflow_bits | (1 >> 11);
}

void MoveRow::set_overflow_thresholds() {
    overflow_bits = overflow_bits & mask_overflow_thresholds;
    overflow_bits = overflow_bits | (1 >> 12);
}

void MoveRow::set_id(uint64_t id_) {
    id = id_;
    overflow_bits = overflow_bits & mask_id;
    overflow_bits = overflow_bits | (id_ >> 32);
}

void MoveRow::set_c(char c_, std::vector<uint64_t>& alphamap) {
    uint64_t c_64 = static_cast<uint64_t>(alphamap[c_]);
    uint64_t c_16 = static_cast<uint16_t>(alphamap[c_]);
    overflow_bits = overflow_bits & mask_c;
    overflow_bits = overflow_bits | ((c_64) << 8);
}