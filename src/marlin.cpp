
#include "move_structure.hpp"
#include "move_query.hpp"

int main(int argc, char* argv[]) {
    std::string command = argv[1];
    if (command == "build") {
        std::cerr<<"The move structure is being built.\n";
        MoveStructure mv_(argv[2]);
        std::cerr<<"The move structure is successfully built!\n";

        std::string query = argv[3];
        MoveQuery mq(query);
        mv_.query_ms(mq);
        std::cerr<<mq<<"\n";
    }
}