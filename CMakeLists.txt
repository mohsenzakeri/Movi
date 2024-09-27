include(FetchContent)

cmake_minimum_required(VERSION 3.0)

project(movi)

#set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED On)
#set(CMAKE_CXX_EXTENSIONS Off)

include_directories("include")

find_package(ZLIB)

# set(move_lib_srcs src/move_structure.cpp)
# add_library(move_structure STATIC ${move_lib_srcs})

# set(move_row_lib_srcs src/move_row.cpp)
# add_library(move_row STATIC ${move_row_lib_srcs})

#----------------------------------------
# Add sdsl
#----------------------------------------
FetchContent_Declare(
  sdsl
  # GIT_REPOSITORY https://github.com/elarielcl/sdsl-lite.git
  GIT_REPOSITORY https://github.com/simongog/sdsl-lite
)

FetchContent_GetProperties(sdsl)
if(NOT sdsl_POPULATED)
  FetchContent_Populate(sdsl)

  set(GENERATE_DOC OFF CACHE BOOL "Do not generate doxygen for sdsl-lite")

  add_subdirectory(${sdsl_SOURCE_DIR} ${sdsl_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# SET(COMPILE_FLAGS "-O2 -DMODE=0")
# debug:
# SET(COMPILE_FLAGS "-g -DMODE=0")

## SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${COMPILE_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14")
# compile with address sanitizer
# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)

# UBSAN_OPTIONS=silence_unsigned_overflow=1
# add_compile_options(-fsanitize-recover=unsigned-integer-overflow)
# add_link_options(-fsanitize-recover=unsigned-integer-overflow)

add_executable(movi-default src/movi.cpp src/move_structure.cpp src/move_row.cpp src/read_processor.cpp)
target_link_libraries(movi-default sdsl ZLIB::ZLIB)
target_compile_options(movi-default PUBLIC -O2 -DMODE=0)

add_executable(movi-constant src/movi.cpp src/move_structure.cpp src/move_row.cpp src/read_processor.cpp)
target_link_libraries(movi-constant sdsl ZLIB::ZLIB)
target_compile_options(movi-constant PUBLIC -O2 -DMODE=1)

add_executable(movi-compact src/movi.cpp src/move_structure.cpp src/move_row.cpp src/read_processor.cpp)
target_link_libraries(movi-compact sdsl ZLIB::ZLIB)
target_compile_options(movi-compact PUBLIC -O2 -DMODE=3)

add_executable(movi-split src/movi.cpp src/move_structure.cpp src/move_row.cpp src/read_processor.cpp)
target_link_libraries(movi-split sdsl ZLIB::ZLIB)
target_compile_options(movi-split PUBLIC -O2 -DMODE=4)

add_executable(prepare_ref src/prepare_ref.cpp)
target_link_libraries(prepare_ref ZLIB::ZLIB)
