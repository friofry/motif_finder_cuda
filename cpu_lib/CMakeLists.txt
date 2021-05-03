cmake_minimum_required(VERSION 3.17)
project(motif_finder_cpu CXX)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_library(motif_finder_cpu
        motif_finder_cpu_external.h
        motif_finder_cpu_external.cpp
        i_range_motif_finder.h
        motif_finder_cpu_internal.h
        motif_finder_cpu_internal.cpp
        cpu_range_motif_finder.cpp
        cpu_range_motif_finder.h
        sequence_hashes_params.h)

add_library(lib::motif_finder_cpu ALIAS motif_finder_cpu)

target_link_libraries(motif_finder_cpu
                      lib::common
                      lib::conversion)

target_include_directories(motif_finder_cpu PUBLIC ${PROJECT_SOURCE_DIR}/.)