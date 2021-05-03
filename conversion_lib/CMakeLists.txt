cmake_minimum_required(VERSION 3.17)
project(conversion_lib CXX)
set(CMAKE_CXX_STANDARD 14)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)

set(COMMON_CXX_FLAGS -c -m64 -D_USE_MATH_DEFINES)

set(SOURCES
        letter_conversions.h
        letter_conversions.cpp
        sequences_to_hashes.h
        sequences_to_hashes.cpp
        hash_conversions.h
        hash_conversions.cpp
        )

add_library(conversion_lib ${SOURCES})
add_library(lib::conversion ALIAS conversion_lib)

target_compile_options(conversion_lib PUBLIC ${COMMON_CXX_FLAGS})
target_link_libraries(conversion_lib Threads::Threads)
target_link_libraries(conversion_lib lib::common)

target_include_directories(conversion_lib PUBLIC ${PROJECT_SOURCE_DIR}/.)
target_include_directories(conversion_lib PRIVATE .)
