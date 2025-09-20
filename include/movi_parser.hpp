#ifndef MOVI_PARSER_HPP
#define MOVI_PARSER_HPP

#include <iostream>

#include <cxxopts.hpp>

#include "utils.hpp"
#include "movi_options.hpp"

bool parse_command(int argc, char** argv, MoviOptions& movi_options, bool supress_messages = false);

#endif