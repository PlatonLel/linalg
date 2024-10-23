#include <gtest/gtest.h>
#include <Linalg.h>

int main(int argc, char **argv) {

    linalg::Matrix m = {{1,2},{1,2}};
    linalg::Matrix m2 = {{1,2},{1,2},{1,2}};
    linalg::Matrix m3 = {};
    m2 = std::move(m);
    std::cout << m3.empty() << m3.rows()<<m3.columns();
    std::cout << m2;

//    // инициализация Google Test
//    testing::InitGoogleTest(&argc, argv);
//
//    // запуск всех тестов
//    return RUN_ALL_TESTS();

}