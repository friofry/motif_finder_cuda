cmake_minimum_required(VERSION 3.18)
project(gpu_lib LANGUAGES CUDA)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 14)
include_directories( .  ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

set(ARCH "30;35;37;50;52;60;61;70")

set(SOURCES
        finder_kernel.cu
        finder_kernel.cuh
        internal_gpu_algorithm.cpp
        internal_gpu_algorithm.h
        external_gpu_algorithm.cpp
        external_gpu_algorithm.h
        gpu_cuda_params.h
        gpu_info.cpp
        gpu_info.h
        gpu_memory.cu
        gpu_memory.h)

add_library(${PROJECT_NAME} ${SOURCES})
add_library(lib::gpu ALIAS ${PROJECT_NAME})

set_target_properties(
        ${PROJECT_NAME}
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES "${ARCH}")


target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${COMMON_NVCC_FLAGS}>)
target_link_libraries(${PROJECT_NAME}
        lib::common
        lib::conversion)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/.)
target_include_directories(${PROJECT_NAME} PRIVATE .)


