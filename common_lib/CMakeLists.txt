cmake_minimum_required(VERSION 3.17)
project(common_lib CXX)
set(CMAKE_CXX_STANDARD 14)


set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)

set(COMMON_CXX_FLAGS -c -m64 -D_USE_MATH_DEFINES)

set(SOURCES
        sequence_hashes.h
        gmps.cpp
        gmps.h
        safe_counter.h
        timer.h
        config.h
        fst_reader.cpp
        fst_reader.h
        run_parallel.cpp
        run_parallel.h)


add_library(common_lib ${SOURCES})
add_library(lib::common ALIAS common_lib)

target_compile_options(common_lib PUBLIC ${COMMON_CXX_FLAGS})
target_link_libraries(common_lib Threads::Threads)

target_include_directories(common_lib PUBLIC ${PROJECT_SOURCE_DIR}/.)
target_include_directories(common_lib PRIVATE .)

