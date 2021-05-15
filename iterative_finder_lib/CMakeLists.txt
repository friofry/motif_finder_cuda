cmake_minimum_required(VERSION 3.17)
project(iterative_finder_lib CXX)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_library(iterative_finder_lib
        argo_cuda_params.h
        argo_cuda_params_reader.cpp
        argo_cuda_params_reader.h
        find_motifs_iterative.cpp
        find_motifs_iterative.h
        )

target_link_libraries(iterative_finder_lib
        lib::common
        lib::stat
        lib::cpu
        )

add_library(lib::iterative_finder ALIAS iterative_finder_lib)

target_include_directories(iterative_finder_lib PUBLIC ${PROJECT_SOURCE_DIR}/.)