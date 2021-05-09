cmake_minimum_required(VERSION 3.17)
project(benchmark_cpu CXX)

set(CMAKE_CXX_STANDARD 14)

add_executable(benchmark_cpu
        main.cpp)

target_link_libraries(benchmark_cpu
        lib::cpu
        lib::occurrence_algorithm_lib
        3rdparty::cxxopts
        )
