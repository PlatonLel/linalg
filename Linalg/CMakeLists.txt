cmake_minimum_required(VERSION 3.10)

project(Matrix)

add_library(Matrix Source/Linalg.cpp)

target_include_directories(Matrix PUBLIC include/)
