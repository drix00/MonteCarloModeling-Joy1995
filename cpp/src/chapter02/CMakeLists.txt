cmake_minimum_required(VERSION 3.6)
project(chapter02)

set(CMAKE_CXX_STANDARD 14)

message("${CMAKE_CXX_FLAGS}")

set(SOURCE_FILES random_walk_2D.h random_walk_2D.cpp)
add_executable(random_walk_2D ${SOURCE_FILES})
