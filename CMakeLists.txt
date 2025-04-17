cmake_minimum_required(VERSION 3.0)
project(movi)

# C++ standard settings
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17")

# Output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

# Dependencies
include_directories("include")

include(FetchContent)

find_package(OpenMP REQUIRED)

#----------------------------------------
# Add zlib
#----------------------------------------
# find_package(ZLIB)
FetchContent_Declare(
  zlib
  GIT_REPOSITORY https://github.com/madler/zlib.git
)

FetchContent_GetProperties(zlib)
if(NOT zlib_POPULATED)
  FetchContent_Populate(zlib)
  add_subdirectory(${zlib_SOURCE_DIR} ${zlib_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

#----------------------------------------
# Add sdsl
#----------------------------------------
FetchContent_Declare(
  sdsl
  GIT_REPOSITORY https://github.com/simongog/sdsl-lite
)

FetchContent_GetProperties(sdsl)
if(NOT sdsl_POPULATED)
  FetchContent_Populate(sdsl)
  set(GENERATE_DOC OFF CACHE BOOL "Do not generate doxygen for sdsl-lite")
  add_subdirectory(${sdsl_SOURCE_DIR} ${sdsl_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# External repositories
include(ExternalProject)
set(EXTERNAL_REPOS_DIR ${CMAKE_BINARY_DIR}/external_repos)

#----------------------------------------
# Add pfp-thresholds
#----------------------------------------
if (NOT EXISTS ${EXTERNAL_REPOS_DIR}/pfp-thresholds-build/pfp_thresholds)
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

#----------------------------------------
# Add r-permute
#----------------------------------------
if (NOT EXISTS ${EXTERNAL_REPOS_DIR}/r-permute-build/test/src/run_constructor)
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

add_subdirectory(src)