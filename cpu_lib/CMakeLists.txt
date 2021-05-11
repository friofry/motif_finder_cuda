cmake_minimum_required(VERSION 3.17)
project(cpu_lib CXX)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_library(cpu_lib
        cpu_common_algorithm.cpp
        cpu_common_algorithm.h
        external_cpu_algorithm.cpp
        external_cpu_algorithm.h
        internal_cpu_algorithm.cpp
        internal_cpu_algorithm.h)

add_library(lib::cpu ALIAS cpu_lib)

target_link_libraries(cpu_lib
                      lib::common
                      lib::conversion)

target_include_directories(cpu_lib PUBLIC ${PROJECT_SOURCE_DIR}/.)