cmake_minimum_required(VERSION 3.10)
project(EP1)

add_library(algelin include/algelin.cpp)

add_executable(ep1 ep1.cpp)

target_link_libraries(ep1 PRIVATE algelin)
