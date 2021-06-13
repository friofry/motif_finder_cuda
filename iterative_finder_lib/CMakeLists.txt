cmake_minimum_required(VERSION 3.18)
project(iterative_finder_lib CXX)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_library(iterative_finder_lib
        argo_cuda_params.h
        argo_cuda_params_reader.cpp
        argo_cuda_params_reader.h
        iterative_finder.cpp
        iterative_finder.h
        create_stat_model.h
        create_stat_model.cpp
        find_most_important_motif.cpp
        find_most_important_motif.h
        motif_data.h
        filter_shifted_motifs.cpp
        filter_shifted_motifs.h
        motif_data_utils.cpp
        motif_data_utils.h)

target_link_libraries(iterative_finder_lib
        lib::common
        lib::conversion
        lib::stat
        lib::cpu
        )

add_library(lib::iterative_finder ALIAS iterative_finder_lib)

target_include_directories(iterative_finder_lib PUBLIC ${PROJECT_SOURCE_DIR}/.)