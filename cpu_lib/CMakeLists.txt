cmake_minimum_required(VERSION 3.17)
project(cpu_lib CXX)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_library(cpu_lib
        CpuCommonAlgorithm.cpp
        CpuCommonAlgorithm.h
        ExternalCpuAlgorithm.cpp
        ExternalCpuAlgorithm.h
        InternalCpuAlgorithm.cpp
        InternalCpuAlgorithm.h)

add_library(lib::cpu ALIAS cpu_lib)

target_link_libraries(cpu_lib
                      lib::common
                      lib::conversion
                      lib::occurrence_algorithm_lib)

target_include_directories(cpu_lib PUBLIC ${PROJECT_SOURCE_DIR}/.)