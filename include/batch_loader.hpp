// The following code is borrowed from SPUMONI by Omar Y. Ahmed with minor modifications
// https://github.com/oma219/spumoni/blob/main/include/batch_loader.hpp

/*
  * File: batch_loader.hpp 
  * Description: Header file for batch_loader.cpp
  *             
  * Start Date: April 21, 2022
  *
  * Note: The BatchLoader class is inspired by the implementation
  *       in the Kraken2 repository at this link 
  *       (https://github.com/DerrickWood/kraken2/blob/master/src/seqreader.h). It uses
  *       different logic to determine the batch size to load balance between
  *       multiple reader threads.
  */

#ifndef BATCH_LOADER_H
#define BATCH_LOADER_H

#include <sstream> 
#include <fstream>


#define FATAL_ERROR(...) do {std::fprintf(stderr, "\n\033[31mError: \033[0m"); std::fprintf(stderr, __VA_ARGS__);\
                              std::fprintf(stderr, "\n\n"); std::exit(1);} while(0)
#define ASSERT(condition, msg) do {if (!condition){std::fprintf(stderr, "Assertion Failed: %s\n", msg); \
                                                   std::exit(1);}} while(0)

enum query_input_type {FA, FQ, NOT_CLEAR};

struct Read {
    std::string header_line; // name of read
    std::string id; // only header until first whitespace
    std::string seq; // sequence
    std::string qual; // qualities if needed

    query_input_type read_format;
};

class BatchLoader {
public:
    BatchLoader();
    ~BatchLoader() {};

    bool loadBatch(std::ifstream& input, size_t num_bases, size_t min_number_of_reads);
    bool grabNextRead(Read& curr_read);
    
private:
    query_input_type input_format;
    std::stringstream batch_stream; // will hold entire batch
    std::string batch_buffer; // will be used for reading
};

#endif /* End of BATCH_LOADER_H */