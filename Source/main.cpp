#include <gtest/gtest.h>
#include <Linalg.h>

int main(int argc, char **argv) {
//
//    linalg::Matrix m = {{1,2},{1,2}};
//    linalg::Matrix m2 = {{1,2},{1,2},{1,2}};
//    m2 = std::move(m);
//    std::cout << m.empty();
//    std::cout << m2;

    // инициализация Google Test
    testing::InitGoogleTest(&argc, argv);

    // запуск всех тестов
    return RUN_ALL_TESTS();

}