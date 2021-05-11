cmake_minimum_required(VERSION 3.17)
project(json)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_library(json INTERFACE)
add_library(3rdparty::json ALIAS json)

target_include_directories( json INTERFACE ${PROJECT_SOURCE_DIR}/.)
