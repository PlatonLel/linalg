#pragma once
#include <gtest/gtest.h>
#include "Linalg.h"

class Matrix_Test : public ::testing::Test {
protected:
    Linalg::Matrix* mat_sum1;  // Указатели на объекты Matrix
    Linalg::Matrix* mat_sum2;
    Linalg::Matrix* mat_pow_e;
//    Linalg::Matrix* mat_sum2;
//    Linalg::Matrix* mat_sum2;
//    Linalg::Matrix* mat_sum2;
//    Linalg::Matrix* mat_sum2;
//    Linalg::Matrix* mat_sum2;

    // SetUp вызывается перед каждым тестом
    void SetUp() override {
        mat_sum1 = new Linalg::Matrix({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}});  // Инициализация матриц
        mat_sum2 = new Linalg::Matrix({{9.0, 8.0, 7.0}, {6.0, 5.0, 4.0}, {3.0, 2.0, 1.0}});
        mat_pow_e = new Linalg::Matrix({{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
    }

    // TearDown вызывается после каждого теста
    void TearDown() override {
        delete mat_sum1;  // Освобождаем ресурсы
        delete mat_sum2;
        delete mat_pow_e;
    }
};