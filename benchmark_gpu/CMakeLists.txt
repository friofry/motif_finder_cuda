cmake_minimum_required(VERSION 3.17)
project(benchmark_gpu CXX)

set(CMAKE_CXX_STANDARD 14)

add_executable(benchmark_gpu
               main.cpp)

target_link_libraries(benchmark_gpu
        lib::common
        lib::gpu
        3rdparty::cxxopts
        )
