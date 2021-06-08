cmake_minimum_required(VERSION 3.18)
project(cxxopts)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_library(cxxopts INTERFACE)
add_library(3rdparty::cxxopts ALIAS cxxopts)

target_include_directories( cxxopts INTERFACE ${PROJECT_SOURCE_DIR}/.)
