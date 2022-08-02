#include "move_structure.hpp"

MoveStructure::MoveStructure(char* input_file) {
    build();
}

uint32_t MoveStructure::LF(uint32_t row_number) {
    return 0;
}

void MoveStructure::build() {
    length = bwt_string.length();
    std::cerr<<"length: " << length << "\n";
}