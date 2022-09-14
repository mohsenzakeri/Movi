#include <cstdint>
#include <limits>

#include "move_row.hpp"

uint32_t mask_p = ~((1U << 8) - 1);
uint32_t mask_pp = ~(((1U << 8) - 1) << 8);
uint32_t mask_id = ~(((1U << 8) - 1) << 16);

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

uint64_t MoveRow::get_p() {
    uint32_t a = overflow_bits & (~mask_p);
    uint64_t b = a;
    b = (b << 32);
    uint64_t c = p;
    c = c | b;
    return c;
}

uint64_t MoveRow::get_pp() {
    uint32_t a = (overflow_bits & (~mask_pp)) >> 8;
    uint64_t b = a;
    b = (b << 32);
    uint64_t c = pp;
    c = c | b;
    return c;
}

uint64_t MoveRow::get_id() {
    uint32_t a = (overflow_bits & (~mask_id)) >> 16;
    uint64_t b = a;
    b = (b << 32);
    uint64_t c = id;
    c = c | b;
    return c;
}