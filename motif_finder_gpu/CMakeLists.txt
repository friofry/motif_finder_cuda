cmake_minimum_required(VERSION 3.18)
project(motif_finder_gpu CUDA CXX)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(motif_finder_gpu
        main.cpp)

target_link_libraries(motif_finder_gpu
        lib::stat
        lib::gpu
        lib::iterative_finder
        )

