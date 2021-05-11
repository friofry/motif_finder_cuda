cmake_minimum_required(VERSION 3.17)
project(benchmark_gpu CUDA CXX)

set(CMAKE_CXX_STANDARD 14)

add_executable(benchmark_gpu
               cmd_parser.h
               cmd_parser.cpp
               main.cpp)

target_link_libraries(benchmark_gpu
        lib::common
        lib::conversion
        lib::stat
        lib::benchmark_common
        lib::gpu
        3rdparty::cxxopts
        3rdparty::json
        )
