cmake_minimum_required(VERSION 3.10)

project(Linalg)

add_subdirectory(complex)

add_library(Linalg INTERFACE)

target_include_directories(Linalg INTERFACE Linalg/Include)

target_link_libraries(Linalg INTERFACE Complex)
