cmake_minimum_required(VERSION 3.17)
project(motif_finder_cpu CXX)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(motif_finder_cpu
        main.cpp)

target_link_libraries(motif_finder_cpu
        lib::stat
        lib::cpu
        lib::iterative_finder
        )

