add_executable(movi-large movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-large sdsl ZLIB::ZLIB)
target_compile_options(movi-large PUBLIC -O2 -DMODE=0)

add_executable(movi-constant movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-constant sdsl ZLIB::ZLIB)
target_compile_options(movi-constant PUBLIC -O2 -DMODE=1)

add_executable(movi-split movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-split sdsl ZLIB::ZLIB)
target_compile_options(movi-split PUBLIC -O2 -DMODE=4)

add_executable(movi-regular movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-regular sdsl ZLIB::ZLIB)
target_compile_options(movi-regular PUBLIC -O2 -DMODE=3)

add_executable(movi-regular-thresholds movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-regular-thresholds sdsl ZLIB::ZLIB)
target_compile_options(movi-regular-thresholds PUBLIC -O2 -DMODE=6)

add_executable(movi-blocked movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-blocked sdsl ZLIB::ZLIB)
target_compile_options(movi-blocked PUBLIC -O2 -DMODE=2)

add_executable(movi-blocked-thresholds movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-blocked-thresholds sdsl ZLIB::ZLIB)
target_compile_options(movi-blocked-thresholds PUBLIC -O2 -DMODE=8)

add_executable(movi-tally movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-tally sdsl ZLIB::ZLIB)
target_compile_options(movi-tally PUBLIC -O2 -DMODE=5)

add_executable(movi-tally-thresholds movi.cpp move_structure.cpp move_row.cpp read_processor.cpp utils.cpp sequitur.cpp)
target_link_libraries(movi-tally-thresholds sdsl ZLIB::ZLIB)
target_compile_options(movi-tally-thresholds PUBLIC -O2 -DMODE=7)

add_executable(movi movi_launcher.cpp)

add_executable(prepare_ref prepare_ref.cpp)
target_link_libraries(prepare_ref ZLIB::ZLIB)

if (BUILD_TEST STREQUAL "1")
  add_executable(test_bidirectional test_bidirectional.cpp move_structure.cpp move_row.cpp utils.cpp sequitur.cpp)
  target_link_libraries(test_bidirectional sdsl ZLIB::ZLIB)
  target_compile_options(test_bidirectional PUBLIC -O2 -DMODE=3)

  add_executable(test_bidirectional test_bidirectional.cpp move_structure.cpp move_row.cpp utils.cpp sequitur.cpp)
  target_link_libraries(test_bidirectional sdsl ZLIB::ZLIB)
  target_compile_options(test_bidirectional_blocked PUBLIC -O2 -DMODE=2)
endif ()
