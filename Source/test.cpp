#include "../Include/test.h"

TEST_F(Matrix_Test, TestMatrixAddition) {
    // Тестируем сложение матриц
    Linalg::Matrix result_1 = *mat_sum1 + *mat_sum2;
    Linalg::Matrix result_2 = *mat_sum2 + *mat_sum1;

    // Проверяем корректность сложения
    EXPECT_DOUBLE_EQ(result_1(0,0), 10.0);
    EXPECT_DOUBLE_EQ(result_1(0,1), 10.0);
    EXPECT_DOUBLE_EQ(result_1(0,2), 10.0);
    EXPECT_DOUBLE_EQ(result_1(1,0), 10.0);
    EXPECT_DOUBLE_EQ(result_1(1,1), 10.0);
    EXPECT_DOUBLE_EQ(result_1(1,2), 10.0);
    EXPECT_DOUBLE_EQ(result_1(2,0), 10.0);
    EXPECT_DOUBLE_EQ(result_1(2,1), 10.0);
    EXPECT_DOUBLE_EQ(result_1(2,2), 10.0);

    EXPECT_DOUBLE_EQ(result_2(0,0), 10.0);
    EXPECT_DOUBLE_EQ(result_2(0,1), 10.0);
    EXPECT_DOUBLE_EQ(result_2(0,2), 10.0);
    EXPECT_DOUBLE_EQ(result_2(1,0), 10.0);
    EXPECT_DOUBLE_EQ(result_2(1,1), 10.0);
    EXPECT_DOUBLE_EQ(result_2(1,2), 10.0);
    EXPECT_DOUBLE_EQ(result_2(2,0), 10.0);
    EXPECT_DOUBLE_EQ(result_2(2,1), 10.0);
    EXPECT_DOUBLE_EQ(result_2(2,2), 10.0);
}

TEST_F(Matrix_Test, TestMatrixPower) {
    // Тестируем сложение матриц
    Linalg::Matrix result = Linalg::power(*mat_pow_e,10);

    // Проверяем корректность сложения
    EXPECT_DOUBLE_EQ(result(0,0), 1.0);
    EXPECT_DOUBLE_EQ(result(0,1), 0.0);
    EXPECT_DOUBLE_EQ(result(0,2), 0.0);
    EXPECT_DOUBLE_EQ(result(1,0), 0.0);
    EXPECT_DOUBLE_EQ(result(1,1), 1.0);
    EXPECT_DOUBLE_EQ(result(1,2), 0.0);
    EXPECT_DOUBLE_EQ(result(2,0), 0.0);
    EXPECT_DOUBLE_EQ(result(2,1), 0.0);
    EXPECT_DOUBLE_EQ(result(2,2), 1.0);
}

