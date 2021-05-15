cmake_minimum_required(VERSION 3.17)
project(motif_finder)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(motif_finder
        create_stat_model.cpp
        create_stat_model.h
        main.cpp)

target_link_libraries(motif_finder
        lib::stat
        lib::iterative_finder
        )

