#include <cstdint>
#include <limits>
#include <iostream>

#include "move_row.hpp"

MoveRow::MoveRow(uint64_t p_, uint16_t n_, uint64_t pp_, uint64_t offset_, uint64_t id_) {
    this->init(p_, n_, pp_, offset_, id_);
}

void MoveRow::init(uint64_t p_, uint16_t n_, uint64_t pp_, uint64_t offset_, uint64_t id_) {
    overflow_bits = std::numeric_limits<uint32_t>::max();
    this->set_p(p_);
    this->set_n(n_);
    this->set_pp(pp_);
    this->set_offset(offset_);
    this->set_id(id_);
}

void MoveRow::init(uint64_t p_, uint16_t n_, uint64_t pp_, uint64_t offset_, uint64_t id_, char c_) {
    overflow_bits = std::numeric_limits<uint32_t>::max();
    this->set_p(p_);
    this->set_n(n_);
    this->set_pp(pp_);
    this->set_offset(offset_);
    this->set_id(id_);
    this->set_c(c_);
}

std::ostream& operator<<(std::ostream& os, const MoveRow& mr)
{
    os << "p:" << mr.get_p() << " n:" << mr.get_n() << " pp:" << mr.get_pp() << " offset: " << mr.get_offset() << " id:" << mr.get_id();
    return os;
}
void MoveRow::set_p(uint64_t p_) {
    p = p_;
    overflow_bits = overflow_bits & mask_p;
    overflow_bits = overflow_bits | (p_ >> 32);
}

void MoveRow::set_n(uint16_t n_) {
    n = n_;
}

void MoveRow::set_offset(uint64_t offset_) {
    offset = offset_;
}

void MoveRow::set_pp(uint64_t pp_) {
    pp = pp_;
    overflow_bits = overflow_bits & mask_pp;
    overflow_bits = overflow_bits | ((pp_ >> 32) << 8);
}

void MoveRow::set_id(uint64_t id_) {
    id = id_;
    overflow_bits = overflow_bits & mask_id;
    overflow_bits = overflow_bits | ((id_ >> 32) << 16);
}

void MoveRow::set_c(char c_) {
    uint64_t c_64 = static_cast<uint64_t>(c_);
    overflow_bits = overflow_bits & mask_c;
    overflow_bits = overflow_bits | ((c_64) << 24);
}