cmake_minimum_required(VERSION 3.20)
project(SLAE_project)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(SLAE_WITH_TESTS "Build tests" OFF)

add_subdirectory(src)
add_subdirectory(time_measure)

if(SLAE_WITH_TESTS)
    add_subdirectory(tests)
endif ()