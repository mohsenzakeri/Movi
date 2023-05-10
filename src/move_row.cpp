#include <cstdint>
#include <limits>
#include <iostream>

#include "move_row.hpp"

MoveRow::MoveRow(uint16_t n_, uint16_t offset_, uint64_t id_) {
    this->init(n_, offset_, id_);
}
void MoveRow::init(uint16_t n_, uint16_t offset_, uint64_t id_) {
    overflow_bits = std::numeric_limits<uint16_t>::max();
    // this->set_p(p_);
    // this->set_pp(pp_);
    this->set_n(n_);
    this->set_offset(offset_);
    this->set_id(id_);
}

/*void MoveRow::init(uint16_t n_, uint16_t offset_, uint64_t id_, char c_) {
    overflow_bits = std::numeric_limits<uint16_t>::max();
    // this->set_p(p_);
    // this->set_pp(pp_);
    this->set_n(n_);
    this->set_offset(offset_);
    this->set_id(id_);
    this->set_c(c_);
}*/

std::ostream& operator<<(std::ostream& os, const MoveRow& mr)
{
    // os << "p:" << mr.get_p() << " n:" << mr.get_n() << " pp:" << mr.get_pp() << " offset: " << mr.get_offset() << " id:" << mr.get_id();
    os << "n:" << mr.get_n() <<  " offset: " << mr.get_offset() << " id:" << mr.get_id();
    return os;
}

/*void MoveRow::set_p(uint64_t p_) {
    p = p_;
    overflow_bits = overflow_bits & mask_p;
    overflow_bits = overflow_bits | (p_ >> 32);
}*/

void MoveRow::set_n(uint16_t n_) {
    n = n_;
}

void MoveRow::set_overflow_n() {
    // std::cerr<< std::bitset<16>(overflow_bits)<<"\n";
    overflow_bits = overflow_bits & mask_overflow_n;
    overflow_bits = overflow_bits | (1 >> 10);
    // std::cerr<< std::bitset<16>(overflow_bits)<<"\n";
}

void MoveRow::set_offset(uint16_t offset_) {
    offset = offset_;
}

void MoveRow::set_overflow_offset() {
    // std::cerr<< std::bitset<16>(overflow_bits)<<"\n";
    overflow_bits = overflow_bits & mask_overflow_offset;
    overflow_bits = overflow_bits | (1 >> 11);
    // std::cerr<< std::bitset<16>(overflow_bits)<<"\n";
}

void MoveRow::set_overflow_thresholds() {
    // std::cerr<< std::bitset<16>(overflow_bits)<<"\n";
    overflow_bits = overflow_bits & mask_overflow_thresholds;
    overflow_bits = overflow_bits | (1 >> 12);
    // std::cerr<< std::bitset<16>(overflow_bits)<<"\n";
}


/* void MoveRow::set_pp(uint64_t pp_) {
    pp = pp_;
    overflow_bits = overflow_bits & mask_pp;
    overflow_bits = overflow_bits | ((pp_ >> 32) << 8);
}*/

void MoveRow::set_id(uint64_t id_) {
    id = id_;
    overflow_bits = overflow_bits & mask_id;
    overflow_bits = overflow_bits | (id_ >> 32);
}

void MoveRow::set_c(char c_, std::vector<uint64_t>& alphamap) {
    // std::cerr<< c_ << "\n";
    // std::cerr<< "overflow_bits: " << std::bitset<16>(overflow_bits) << "\n";
    uint64_t c_64 = static_cast<uint64_t>(alphamap[c_]);
    uint64_t c_16 = static_cast<uint16_t>(alphamap[c_]);
    // std::cerr<< "c_16: " << std::bitset<16>(c_16) << "\n";
    // std::cerr<< "c_64: " << std::bitset<64>(c_64) << "\n";
    // std::cerr<< "mask_c: " << std::bitset<16>(mask_c) << "\n";
    overflow_bits = overflow_bits & mask_c;
    // std::cerr<< "overflow_bits: " << std::bitset<16>(overflow_bits) << "\n";
    overflow_bits = overflow_bits | ((c_64) << 8);
    // std::cerr<< "overflow_bits: " << std::bitset<16>(overflow_bits) << "\n";
}