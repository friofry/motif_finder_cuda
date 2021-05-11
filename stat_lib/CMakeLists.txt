cmake_minimum_required(VERSION 3.17)
project(stat_lib CXX)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_library(stat_lib
        probability.cpp
        probability.h
        bernulli_stat_model.cpp
        bernulli_stat_model.h
        external_generator_cpu.cpp
        external_generator_cpu.h
        ideal_stat_model.cpp
        ideal_stat_model.h
        markov_stat_model.cpp
        markov_stat_model.h
        stat_model.cpp
        stat_model.h)

add_library(lib::stat ALIAS stat_lib)

target_link_libraries(stat_lib
        lib::common
        lib::conversion)

target_include_directories(stat_lib PUBLIC ${PROJECT_SOURCE_DIR}/.)