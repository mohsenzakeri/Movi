# Add Catch2 (header-only)
include(FetchContent)
FetchContent_Declare(
  catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.5.2 # or latest
)
FetchContent_MakeAvailable(catch2)

# Add test executables
add_executable(basics-tests test_basics.cpp test_utils.cpp)
target_link_libraries(basics-tests PRIVATE movi_lib_test Catch2::Catch2WithMain)
target_include_directories(basics-tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(basics-tests PRIVATE -Wno-deprecated-declarations)
target_compile_definitions(basics-tests PRIVATE BINARY_DIR="${CMAKE_BINARY_DIR}")

add_executable(build-tests test_build.cpp test_utils.cpp)
target_link_libraries(build-tests PRIVATE movi_lib_test Catch2::Catch2WithMain)
target_include_directories(build-tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(build-tests PRIVATE -Wno-deprecated-declarations)
target_compile_definitions(build-tests PRIVATE BINARY_DIR="${CMAKE_BINARY_DIR}")

add_executable(pml-tests test_pml.cpp test_utils.cpp)
target_link_libraries(pml-tests PRIVATE movi_lib_test Catch2::Catch2WithMain)
target_include_directories(pml-tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(pml-tests PRIVATE -Wno-deprecated-declarations)
target_compile_definitions(pml-tests PRIVATE BINARY_DIR="${CMAKE_BINARY_DIR}")

# CTest registration removed - using direct Catch2 execution instead
# add_test(NAME basics-tests COMMAND basics-tests)
# add_test(NAME build-tests COMMAND build-tests)
# add_test(NAME pml-tests COMMAND pml-tests)