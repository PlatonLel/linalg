cmake_minimum_required(VERSION 3.10)

project(Linalg)

add_library(Linalg Source/Linalg.cpp)

target_include_directories(Linalg PUBLIC include/)
