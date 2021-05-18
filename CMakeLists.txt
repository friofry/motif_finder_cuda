cmake_minimum_required(VERSION 3.17)
project(motif_finder_cpu)
set(CMAKE_CXX_STANDARD 14)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG=1)
endif()

if(MSVC)
    add_compile_options("/W4")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
    add_compile_options("-Wall" "-Wextra" "-Werror")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
message("Generated with config types: ${CMAKE_BUILD_TYPE}")

add_subdirectory(3rdparty/cxxopts)
add_subdirectory(3rdparty/json)

add_subdirectory(common_lib)
add_subdirectory(conversion_lib)
add_subdirectory(benchmark_common_lib)
add_subdirectory(stat_lib)
add_subdirectory(cpu_lib)
add_subdirectory(iterative_finder_lib)
add_subdirectory(benchmark_cpu)

find_package( CUDA 10.0 QUIET )
if(NOT CUDA_FOUND)
    message("CUDA 10.0 not found")
else()
    add_subdirectory(gpu_lib)
    add_subdirectory(benchmark_gpu)
    add_subdirectory(motif_finder_gpu)
endif()

add_subdirectory(motif_finder_cpu)
