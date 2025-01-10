include(FetchContent)

cmake_minimum_required(VERSION 3.0)

project(movi)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
#set(CMAKE_CXX_EXTENSIONS Off)

find_package(ZLIB)

include_directories("include")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

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

include(ExternalProject)

set(EXTERNAL_REPOS_DIR ${CMAKE_BINARY_DIR}/external_repos)

if (NOT EXISTS ${EXTERNAL_REPOS_DIR}/r-permute-build/test/src/run_constructor)
  ExternalProject_Add(
      pfp-thresholds
      GIT_REPOSITORY https://github.com/mohsenzakeri/pfp-thresholds
      GIT_TAG probing
      PREFIX ${EXTERNAL_REPOS_DIR}
      SOURCE_DIR ${EXTERNAL_REPOS_DIR}/pfp-thresholds
      BINARY_DIR ${EXTERNAL_REPOS_DIR}/pfp-thresholds-build
      STAMP_DIR ${EXTERNAL_REPOS_DIR}/pfp-thresholds-stamp
      TMP_DIR ${EXTERNAL_REPOS_DIR}/pfp-thresholds-tmp
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_REPOS_DIR}/pfp-thresholds/install
      BUILD_ALWAYS OFF
      INSTALL_COMMAND make install
  )
endif()

if (NOT EXISTS ${EXTERNAL_REPOS_DIR}/pfp-thresholds-build/pfp_thresholds)
  ExternalProject_Add(
      r-permute
      GIT_REPOSITORY https://github.com/drnatebrown/r-permute
      GIT_TAG main
      PREFIX ${EXTERNAL_REPOS_DIR}
      SOURCE_DIR ${EXTERNAL_REPOS_DIR}/r-permute
      BINARY_DIR ${EXTERNAL_REPOS_DIR}/r-permute-build
      STAMP_DIR ${EXTERNAL_REPOS_DIR}/r-permute-stamp
      TMP_DIR ${EXTERNAL_REPOS_DIR}/r-permute-tmp
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_REPOS_DIR}/r-permute/install
      BUILD_ALWAYS OFF
      INSTALL_COMMAND make install
  )
endif()

# SET(COMPILE_FLAGS "-O2 -DMODE=0")
# debug:
# SET(COMPILE_FLAGS "-g -DMODE=0")

## SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${COMPILE_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17")
# compile with address sanitizer
# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)

# UBSAN_OPTIONS=silence_unsigned_overflow=1
# add_compile_options(-fsanitize-recover=unsigned-integer-overflow)
# add_link_options(-fsanitize-recover=unsigned-integer-overflow)

add_subdirectory(src)
