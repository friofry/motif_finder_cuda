cmake_minimum_required(VERSION 3.17)

project(benchmark_common_lib CXX)
set(CMAKE_CXX_STANDARD 14)

set(SOURCES
        benchmark_common.h
        benchmark_common.cpp
        parse_helper.h
        parse_helper.cpp
        sequence_generator.h
        sequence_generator.cpp
        )

add_library(benchmark_common_lib ${SOURCES})
add_library(lib::benchmark_common ALIAS benchmark_common_lib)

target_include_directories(benchmark_common_lib PUBLIC ${PROJECT_SOURCE_DIR}/.)
target_include_directories(benchmark_common_lib PRIVATE .)

target_link_libraries(benchmark_common_lib
        lib::common
        lib::conversion
        lib::stat
        3rdparty::json
        3rdparty::cxxopts
        )
