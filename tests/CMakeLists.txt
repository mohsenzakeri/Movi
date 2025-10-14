# Add Catch2 (header-only)
include(FetchContent)
FetchContent_Declare(
  catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.5.2 # or latest
)
FetchContent_MakeAvailable(catch2)

# Add test executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests_bin)

add_executable(basics-tests test_basics.cpp test_utils.cpp)
target_link_libraries(basics-tests PRIVATE movi_lib_test Catch2::Catch2WithMain)
target_include_directories(basics-tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(basics-tests PRIVATE -Wno-deprecated-declarations)
target_compile_definitions(basics-tests PRIVATE BINARY_DIR="${CMAKE_BINARY_DIR}" MODE=0 COLOR_MODE=0)

add_executable(build-tests test_build.cpp test_utils.cpp)
target_link_libraries(build-tests PRIVATE movi_lib_test Catch2::Catch2WithMain)
target_include_directories(build-tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(build-tests PRIVATE -Wno-deprecated-declarations)
target_compile_definitions(build-tests PRIVATE BINARY_DIR="${CMAKE_BINARY_DIR}" MODE=0 COLOR_MODE=0)

add_executable(pml-tests test_pml.cpp test_utils.cpp)
target_link_libraries(pml-tests PRIVATE movi_lib_test Catch2::Catch2WithMain)
target_include_directories(pml-tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(pml-tests PRIVATE -Wno-deprecated-declarations)
target_compile_definitions(pml-tests PRIVATE BINARY_DIR="${CMAKE_BINARY_DIR}" MODE=0 COLOR_MODE=0)

add_executable(classification-tests test_classification.cpp test_utils.cpp)
target_link_libraries(classification-tests PRIVATE movi_lib_test Catch2::Catch2WithMain)
target_include_directories(classification-tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(classification-tests PRIVATE -Wno-deprecated-declarations)
target_compile_definitions(classification-tests PRIVATE BINARY_DIR="${CMAKE_BINARY_DIR}" MODE=0 COLOR_MODE=0)

# CTest registration removed - using direct Catch2 execution instead
# add_test(NAME basics-tests COMMAND basics-tests)
# add_test(NAME build-tests COMMAND build-tests)
# add_test(NAME pml-tests COMMAND pml-tests)

# Old Bidirectional tests
add_executable(test_bidirectional test_bidirectional.cpp)
target_link_libraries(test_bidirectional movi_lib_test sdsl zlibstatic OpenMP::OpenMP_CXX hclust)
target_include_directories(test_bidirectional PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(test_bidirectional PUBLIC -DMODE=3 -DCOLOR_MODE=0)

add_executable(test_bidirectional_blocked test_bidirectional.cpp)
target_link_libraries(test_bidirectional_blocked movi_lib_test sdsl zlibstatic OpenMP::OpenMP_CXX hclust)
target_include_directories(test_bidirectional_blocked PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(test_bidirectional_blocked PUBLIC -DMODE=2 -DCOLOR_MODE=0)