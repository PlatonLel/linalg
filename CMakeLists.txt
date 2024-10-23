cmake_minimum_required(VERSION 3.10)
project(testM)


set(CMAKE_CXX_STANDARD 17)

# Добавляем библиотеку linalg и gtest
add_subdirectory(Linalg)
add_subdirectory(gtest/googletest)

# Добавляем исполняемый файл для запуска тестов
add_executable(run_tests Source/main.cpp Source/test.cpp)

# Линкуем с Google Test и linalg
target_link_libraries(run_tests PRIVATE
        gtest
        gtest_main
        Linalg
)
