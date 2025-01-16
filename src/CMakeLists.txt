# release
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
# debug:
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# compile with address sanitizer
# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)

add_executable(movi-large movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-large sdsl zlibstatic)
target_compile_options(movi-large PUBLIC -DMODE=0)

add_executable(movi-constant movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-constant sdsl zlibstatic)
target_compile_options(movi-constant PUBLIC -DMODE=1)

add_executable(movi-split movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-split sdsl zlibstatic)
target_compile_options(movi-split PUBLIC -DMODE=4)

add_executable(movi-regular movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-regular sdsl zlibstatic)
target_compile_options(movi-regular PUBLIC -DMODE=3)

add_executable(movi-regular-thresholds movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-regular-thresholds sdsl zlibstatic)
target_compile_options(movi-regular-thresholds PUBLIC -DMODE=6)

add_executable(movi-blocked movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-blocked sdsl zlibstatic)
target_compile_options(movi-blocked PUBLIC -DMODE=2)

add_executable(movi-blocked-thresholds movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-blocked-thresholds sdsl zlibstatic)
target_compile_options(movi-blocked-thresholds PUBLIC -DMODE=8)

add_executable(movi-tally movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-tally sdsl zlibstatic)
target_compile_options(movi-tally PUBLIC -DMODE=5)

add_executable(movi-tally-thresholds movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-tally-thresholds sdsl zlibstatic)
target_compile_options(movi-tally-thresholds PUBLIC -DMODE=7)

add_executable(movi movi_launcher.cpp)

add_executable(prepare_ref prepare_ref.cpp)
target_link_libraries(prepare_ref zlibstatic)

if (BUILD_TEST STREQUAL "1")
  add_executable(test_bidirectional test_bidirectional.cpp move_structure.cpp move_row.cpp utils.cpp sequitur.cpp)
  target_link_libraries(test_bidirectional sdsl zlibstatic)
  target_compile_options(test_bidirectional PUBLIC -DMODE=3)

  add_executable(test_bidirectional test_bidirectional.cpp move_structure.cpp move_row.cpp utils.cpp sequitur.cpp)
  target_link_libraries(test_bidirectional sdsl zlibstatic)
  target_compile_options(test_bidirectional_blocked PUBLIC -DMODE=2)
endif ()