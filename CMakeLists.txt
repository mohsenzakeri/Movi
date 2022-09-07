cmake_minimum_required(VERSION 3.0)

project(marlin)

include_directories("include")
link_directories("/home/mzakeri1/lib/")

include_directories("/home/mzakeri1/include/")

find_package(ZLIB)

set(move_lib_srcs src/move_structure.cpp)
add_library(move STATIC ${move_lib_srcs})

add_executable(marlin src/marlin.cpp)
target_link_libraries(marlin move sdsl ZLIB::ZLIB)