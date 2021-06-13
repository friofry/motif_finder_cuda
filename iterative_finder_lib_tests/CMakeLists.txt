cmake_minimum_required(VERSION 3.14)
project(iterative_finder_lib_tests)

# GoogleTest requires at least C++11
set(CMAKE_CXX_STANDARD 11)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
        iterative_finder_lib_tests
        filter_shifted_motifs_test.cpp
)
target_link_libraries(
        iterative_finder_lib_tests
        gtest_main
        lib::conversion
        lib::iterative_finder
)
include(GoogleTest)
gtest_discover_tests(iterative_finder_lib_tests)