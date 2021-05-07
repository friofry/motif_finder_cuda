cmake_minimum_required(VERSION 3.17)
project(motif_finder_gpu LANGUAGES CUDA)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 14)
include_directories(.)
set(COMMON_NVCC_FLAGS -Xptxas -O3)

set(ARCH "30 35 37 50 52 60 61 70 80 86")

set(SOURCES
        finder_kernel.cu
        finder_kernel.cuh
        gpu_cuda_params.h
        gpu_info.cpp
        gpu_info.h
        gpu_memory.cu
        gpu_memory.h
        gpu_range_motif_finder_external.cpp
        gpu_range_motif_finder_external.h
        gpu_range_motif_finder_internal.cpp
        gpu_range_motif_finder_internal.h
        motif_finder_gpu_external.h
        motif_finder_gpu_external.cpp
        motif_finder_gpu_internal.cpp
        motif_finder_gpu_internal.h
        sequence_hashes_params.h)

add_library(${PROJECT_NAME} ${SOURCES})
add_library(lib::gpu ALIAS ${PROJECT_NAME})

set_target_properties(
        ${PROJECT_NAME}
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES ${ARCH})


target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${COMMON_NVCC_FLAGS}>)
target_link_libraries(${PROJECT_NAME} lib::common)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/.)
target_include_directories(${PROJECT_NAME} PRIVATE .)


