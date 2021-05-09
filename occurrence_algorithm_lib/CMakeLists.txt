cmake_minimum_required(VERSION 3.17)
project(occurrence_algorithm_lib CXX)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_library(occurrence_algorithm_lib
        IAlgorithm.h
        IAlgorithmFactory.h
        SequenceHashesParams.h
        MotifOccurrenceParallelEvaluator.h
        MotifOccurrenceParallelEvaluator.cpp)

add_library(lib::occurrence_algorithm_lib ALIAS occurrence_algorithm_lib)

target_link_libraries(occurrence_algorithm_lib
                      lib::common
                      lib::conversion)

target_include_directories(occurrence_algorithm_lib PUBLIC ${PROJECT_SOURCE_DIR}/.)