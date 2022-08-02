#include <iostream>

#include "move_structure.hpp"

int main(int argc, char* argv[]) {
    std::string command = argv[1];
    if (command == "build") {
        std::cerr<<"The move structure is being built.\n";
        MoveStructure mv_(argv[2]);
        std::cerr<<"The move structure is successfully built!\n";
    }
}