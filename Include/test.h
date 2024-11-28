#pragma once
#include <gtest/gtest.h>
#include "Linalg.hpp"


// класс для тестирования, от которого будут наследовать остальные
class Matrix_Test : public ::testing::Test {
};

// класс для тестирования возведения в степень
class Matrix_test_pow : public ::Matrix_Test {
protected:
    linalg::Matrix* mat_pow_e;
    linalg::Matrix* mat_pow_1;
    linalg::Matrix* mat_pow_2;

    void SetUp() override {
        mat_pow_e = new linalg::Matrix({{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
        mat_pow_1 = new linalg::Matrix({{{1.0, 2.0, 3.0}, {3.0, 1.0, 4.0}, {6.0, 7.0, 19.0}}});
        mat_pow_2 = new linalg::Matrix({{0.9, 2.0, 4.0}, {4.0, 5.0, 6.0}, {1.0, 2.0, 3.0} });
    }

    void TearDown() override {
        delete mat_pow_e;
        delete mat_pow_2;
        delete mat_pow_1;
    }
};

//класс для тестирования суммирования и вычитания
class Matrix_test_addition_and_subtraction : public ::Matrix_Test {
protected:
    linalg::Matrix* mat_sum1;
    linalg::Matrix* mat_sum2;


    void SetUp() override {
        mat_sum1 = new linalg::Matrix({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}});
        mat_sum2 = new linalg::Matrix({{9.0, 8.0, 7.0}, {6.0, 5.0, 4.0}, {3.0, 2.0, 1.0}});
    }

    void TearDown() override {
        delete mat_sum1;
        delete mat_sum2;
    }};

//класс для тестирования умножения
class Matrix_test_multiplication : public ::Matrix_Test {
protected:
    linalg::Matrix* mat_mult1;
    linalg::Matrix* mat_mult2;

    void SetUp() override {
        mat_mult1 = new linalg::Matrix({{1.0, 2.0}, {3.0, 4.0}});
        mat_mult2 = new linalg::Matrix({{2.0, 0.0}, {1.0, 2.0}});
    }

    void TearDown() override {
        delete mat_mult1;
        delete mat_mult2;
    }
};

//транспонирование
class Matrix_test_transpose : public ::Matrix_Test {
protected:
    linalg::Matrix* mat_transpose;
    linalg::Matrix* mat_transpose_2;

    void SetUp() override {
        mat_transpose = new linalg::Matrix({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
        mat_transpose_2 = new linalg::Matrix({{1.0, 2.0, 3.0}, {3.0, 1.0, 4.0}, {6.0, 7.0, 19.0}});
    }

    void TearDown() override {
        delete mat_transpose;
        delete mat_transpose_2;
    }
};

//обратная матрица тест
class Matrix_test_inversion : public ::Matrix_Test {
protected:
    linalg::Matrix* mat_invertible;
    linalg::Matrix* mat_singular;

    void SetUp() override {
        mat_invertible = new linalg::Matrix({{4.0, 7.0}, {2.0, 6.0}});
        mat_singular = new linalg::Matrix({{1.0, 2.0}, {2.0, 4.0}});
    }

    void TearDown() override {
        delete mat_invertible;
        delete mat_singular;
    }
};

//умножение на число
class Matrix_test_scalar_multiplication : public ::Matrix_Test {
protected:
    linalg::Matrix* mat_scalar;

    void SetUp() override {
        mat_scalar = new linalg::Matrix({{1.0, 2.0}, {3.0, 4.0}});
    }

    void TearDown() override {
        delete mat_scalar;
    }
};

//определитель
class Matrix_test_determinant : public ::Matrix_Test {
protected:
    linalg::Matrix* mat_determinant;

    void SetUp() override {
        mat_determinant = new linalg::Matrix({{1.0, 2.0}, {3.0, 4.0}});
    }

    void TearDown() override {
        delete mat_determinant;
    }
};

//матрицы с большими числами
class Matrix_test_large_numbers : public ::Matrix_Test {
protected:
    linalg::Matrix* mat_large;
    linalg::Matrix* mat_small;

    void SetUp() override {
        mat_large = new linalg::Matrix({{1e18, 2e18}, {3e18, 4e18}});
        mat_small = new linalg::Matrix({{1e-18, 2e-18}, {3e-18, 4e-18}});
    }

    void TearDown() override {
        delete mat_large;
        delete mat_small;
    }
};

//след
class Matrix_test_trace : public ::Matrix_Test {
protected:
    linalg::Matrix* mat_square;
    linalg::Matrix* mat_non_square;

    void SetUp() override {
        mat_square = new linalg::Matrix({{1.0, 0.0}, {0.0, 1.0}});
        mat_non_square = new linalg::Matrix({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});  // Неквадратная
    }

    void TearDown() override {
        delete mat_square;
        delete mat_non_square;
    }
};


// сложение с нулевой матрицей
class Matrix_test_zero_matrix_addition : public ::Matrix_Test {
protected:
    linalg::Matrix* mat1;
    linalg::Matrix* zero_matrix;

    void SetUp() override {
        mat1 = new linalg::Matrix({{1.0, 2.0}, {3.0, 4.0}});
        zero_matrix = new linalg::Matrix({{0.0, 0.0}, {0.0, 0.0}});
    }

    void TearDown() override {
        delete mat1;
        delete zero_matrix;
    }
};

//транспонирование не квадратной
class Matrix_test_non_square_transpose : public ::Matrix_Test {
protected:
    linalg::Matrix* mat_non_square;

    void SetUp() override {
        mat_non_square = new linalg::Matrix({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
    }

    void TearDown() override {
        delete mat_non_square;
    }
};

//для доступа вне индексов
class Matrix_test_out_of_bounds : public ::Matrix_Test {
protected:
    linalg::Matrix* mat_bounds;

    void SetUp() override {
        mat_bounds = new linalg::Matrix({{1.0, 2.0}, {3.0, 4.0}});
    }

    void TearDown() override {
        delete mat_bounds;
    }
};

//умножение не квадратных матриц
class Matrix_test_multiplication_identity_non_square : public ::Matrix_Test {
protected:
    linalg::Matrix* mat_non_square;
    linalg::Matrix* identity_non_square;

    void SetUp() override {
        mat_non_square = new linalg::Matrix({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
        identity_non_square = new linalg::Matrix({{1.0, 0.0}, {0.0, 1.0}, {0.0, 0.0}});
    }

    void TearDown() override {
        delete mat_non_square;
        delete identity_non_square;
    }
};


//умножение пустых матриц
class Matrix_test_empty_multiplication : public ::Matrix_Test {
protected:
    linalg::Matrix* empty_matrix1;
    linalg::Matrix* empty_matrix2;

    void SetUp() override {
        empty_matrix1 = new linalg::Matrix(0, 0);
        empty_matrix2 = new linalg::Matrix(0, 0);
    }

    void TearDown() override {
        delete empty_matrix1;
        delete empty_matrix2;
    }
};

//полупустая матрица
class Matrix_test_sparse_matrices : public ::Matrix_Test {
protected:
    linalg::Matrix* sparse_matrix;

    void SetUp() override {
        sparse_matrix = new linalg::Matrix({{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}});
    }

    void TearDown() override {
        delete sparse_matrix;
    }
};

//полностью единичная матрица
class Matrix_test_all_ones : public ::Matrix_Test {
protected:
    linalg::Matrix* all_ones_matrix;

    void SetUp() override {
        all_ones_matrix = new linalg::Matrix({{1.0, 1.0}, {1.0, 1.0}});
    }

    void TearDown() override {
        delete all_ones_matrix;
    }
};

//одна строка
class Matrix_test_single_row : public ::Matrix_Test {
protected:
    linalg::Matrix* single_row_matrix;

    void SetUp() override {
        single_row_matrix = new linalg::Matrix({{1.0, 2.0, 3.0}});
    }

    void TearDown() override {
        delete single_row_matrix;
    }
};

//очень маленькая разница
class Matrix_test_floating_point_precision : public ::Matrix_Test {
protected:
    linalg::Matrix* mat1;
    linalg::Matrix* mat2;

    void SetUp() override {
        mat1 = new linalg::Matrix({{1.0000001, 2.0000001}, {3.0000001, 4.0000001}});
        mat2 = new linalg::Matrix({{1.0000002, 2.0000002}, {3.0000002, 4.0000002}});
    }

    void TearDown() override {
        delete mat1;
        delete mat2;
    }
};

//решейп
class Matrix_test_reshape : public ::Matrix_Test {
protected:
    linalg::Matrix* mat_modify;

    void SetUp() override {
        mat_modify = new linalg::Matrix({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
    }

    void TearDown() override {
        delete mat_modify;
    }
};

//умножение нулевой матрицы
class Matrix_test_zero_matrix_multiplication : public ::Matrix_Test {
protected:
    linalg::Matrix* zero_matrix;
    linalg::Matrix* mat;

    void SetUp() override {
        zero_matrix = new linalg::Matrix({{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});
        mat = new linalg::Matrix({{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}});
    }

    void TearDown() override {
        delete zero_matrix;
        delete mat;
    }
};

class Matrix_test_norm : public ::Matrix_Test {
protected:
    linalg::Matrix* mat;

    void SetUp() override {
         mat = new linalg::Matrix({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}});
    }

    void TearDown() override {
        delete mat;
    };

};







