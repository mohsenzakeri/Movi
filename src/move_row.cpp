#include <cstdint>
#include <limits>

#include "move_row.hpp"

MoveRow::MoveRow(uint64_t p_, uint16_t n_, uint64_t pp_, uint64_t id_) {
    this->init(p_, n_, pp_, id_);
}

void MoveRow::init(uint64_t p_, uint16_t n_, uint64_t pp_, uint64_t id_) {
    overflow_bits = std::numeric_limits<uint32_t>::max();
    this->set_p(p_);
    this->set_n(n_);
    this->set_pp(pp_);
    this->set_id(id_);
}

void MoveRow::set_p(uint64_t p_) {
    p = p_;
    overflow_bits = overflow_bits & mask_p;
    overflow_bits = overflow_bits | (p_ >> 32);
}

void MoveRow::set_n(uint16_t n_) {
    n = n_;
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