cmake_minimum_required(VERSION 3.10)
project(testM)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(Linalg)

# Добавляем исполняемый файл для запуска тестов
add_executable(tests Source/main.cpp)

target_link_libraries(tests PUBLIC
        Linalg
)
