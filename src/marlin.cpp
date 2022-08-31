
#include "move_structure.hpp"
#include "move_query.hpp"

int main(int argc, char* argv[]) {
    std::string command = argv[1];
    if (command == "build") {
        std::cerr<<"The move structure is being built.\n";
        bool mode = argv[2] == "2bits" ? true : false;
        bool verbose = (argc > 5 and std::string(argv[5]) == "verbose");
        MoveStructure mv_(argv[3], mode, verbose);
        std::cerr<<"The move structure is successfully built!\n";
        // mv_.reconstruct();
        // std::cerr<<"The original string is reconstructed.\n";
        // std::cerr<<"The original string is:\n" << mv_.reconstruct() << "\n";
        mv_.seralize(argv[4]);
        std::cerr<<"The move structure is successfully stored at ./" << argv[4] << "/\n";
    } else if (command == "query") {
        bool verbose = (argc > 4 and std::string(argv[4]) == "verbose");
        MoveStructure mv_(verbose);
        mv_.deseralize(argv[2]);
        std::cerr<< "The move structure is read from the file successfully.\n";
        // std::cerr<<"The original string is: " << mv_.reconstruct() << "\n";

        std::string query = argv[3];
        MoveQuery mq(query);
        bool random_jump = true;
        mv_.query_ms(mq, random_jump);
        std::cerr<<mq<<"\n";
    }
}