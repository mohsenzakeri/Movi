#ifndef __UTILS__
#define __UTILS__

#include <sys/stat.h> 

#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>

#include "move_query.hpp"

char complement(char c);

std::string reverse_complement(std::string& fw);

std::string reverse_complement_from_pos(MoveQuery& mq_fw, int32_t pos_on_r, uint64_t match_len);

std::string number_to_kmer(size_t j, size_t m, std::vector<unsigned char>& alphabet, std::vector<uint64_t>& alphamap);

uint64_t kmer_to_number(size_t k, std::string& r, int32_t pos, std::vector<uint64_t>& alphamap, bool rc = false) ;

#endif