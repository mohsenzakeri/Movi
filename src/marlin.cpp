
#include "move_structure.hpp"
#include "move_query.hpp"

int main(int argc, char* argv[]) {
    std::string command = argv[1];
    if (command == "build") {
        std::cerr<<"The move structure is being built.\n";
        bool verbose = (argc > 4 and std::string(argv[4]) == "verbose");
        MoveStructure mv_(argv[2], verbose);
        std::cerr<<"The move structure is successfully built!\n";
        std::cerr<<"The original string is:\n" << mv_.reconstruct() << "\n";
        // mv_.seralize(argv[3]);
        // exit(0);

        std::string query = argv[3];
        MoveQuery mq(query);
        bool random_jump = true;
        mv_.query_ms(mq, random_jump);
        std::cerr<<mq<<"\n";
    }
}