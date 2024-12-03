cmake_minimum_required(VERSION 3.10)
project(testM)
set(CMAKE_CXX_STANDARD 17)

# Добавляем библиотеку linalg и gtest
add_subdirectory(Linalg)
#add_subdirectory(gtest/googletest)

# Добавляем исполняемый файл для запуска тестов
add_executable(tests Source/main.cpp)

# Линкуем с Google Test и linalg
target_link_libraries(tests PUBLIC
        Linalg
)
#``````````````````
#cmake_minimum_required(VERSION 3.10)
#
#set(CMAKE_CXX_STANDARD 17)
#project(testVector)
#
#add_subdirectory(vector)
#
#add_subdirectory(googletest)
#
#add_executable(proga source/test.cpp)
#
#target_link_libraries(proga PUBLIC vector GTest::gtest_main)