cmake_minimum_required(VERSION 3.10)

project(testMatrix)

# Определяем стандарт C++
set(CMAKE_CXX_STANDARD 17)

# Добавляем библиотеку Linalg
add_subdirectory(Linalg)

# Добавляем основной исполняемый файл
add_executable(proga Source/main.cpp)
target_link_libraries(proga PUBLIC Linalg)

# Включаем тестирование
enable_testing()

# Загружаем Google Test
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/heads/main.zip
        DOWNLOAD_EXTRACT_TIMESTAMP true
)

# Собираем Google Test
FetchContent_MakeAvailable(googletest)

# Добавляем тестовый исполняемый файл
add_executable(testM Source/test.cpp)

# Линкуем тестовый файл с библиотеками Google Test и Linalg
target_link_libraries(testM gtest_main Linalg)

include_directories(D:/mine/Laboratory 1/cmake-build-debug/_deps/googletest-src/googletest/include)

# Регистрируем тест
add_test(NAME example_test COMMAND testM)
