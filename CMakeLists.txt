cmake_minimum_required(VERSION 3.20)
project(SLAE_project)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(time_measure)