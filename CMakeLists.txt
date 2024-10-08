cmake_minimum_required(VERSION 3.10)

project(testMatrix)

add_subdirectory(Linalg)

add_executable(proga Source/main.cpp)

target_link_libraries(proga PUBLIC Linalg)
