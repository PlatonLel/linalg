#include "../Include/test.h"

// сложение тест
TEST_F(Matrix_test_addition_and_subtraction, TestMatrixAddition) {

    EXPECT_EQ(*mat_sum1 + *mat_sum2, linalg::Matrix({{10.0, 10.0, 10.0}, {10.0, 10.0, 10.0}, {10.0, 10.0, 10.0}}));
    EXPECT_EQ(*mat_sum2 + *mat_sum1, linalg::Matrix({{10.0, 10.0, 10.0}, {10.0, 10.0, 10.0}, {10.0, 10.0, 10.0}}));
}

// вычитание тест
TEST_F(Matrix_test_addition_and_subtraction, TestMatrixSubtraction) {

    EXPECT_EQ(*mat_sum1 - *mat_sum2, linalg::Matrix({{-8.0, -6.0, -4.0}, {-2.0, 0.0, 2.0}, {4.0, 6.0, 8.0}}));
    EXPECT_EQ(*mat_sum2 - *mat_sum1, linalg::Matrix({{8.0, 6.0, 4.0}, {2.0, 0.0, -2.0}, {-4.0, -6.0, -8.0}}));
    EXPECT_EQ(*mat_sum2 - *mat_sum2, linalg::Matrix({{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}));
}

// возведение в степень тест
TEST_F(Matrix_test_pow, TestMatrixPower) {

    EXPECT_EQ(linalg::power(*mat_pow_e, 10), *mat_pow_e);
    EXPECT_EQ(linalg::power(*mat_pow_e, -10), *mat_pow_e);
    EXPECT_EQ(linalg::power(power(*mat_pow_1, 2), 3), linalg::power(power(*mat_pow_1, 3), 2));
    EXPECT_EQ(linalg::power(*mat_pow_1, -1) * linalg::power(*mat_pow_2, -1),
              linalg::power(*mat_pow_2 * *mat_pow_1, -1));
}

// умножение тест
TEST_F(Matrix_test_multiplication, MatrixMultiplication) {

    EXPECT_EQ(*mat_mult1 * *mat_mult2, linalg::Matrix({{4.0, 4.0}, {10.0, 8.0}}));
}

// транспонирование тест
TEST_F(Matrix_test_transpose, MatrixTranspose) {

    EXPECT_EQ(linalg::transpose(*mat_transpose), linalg::Matrix({{1.0, 4.0}, {2.0, 5.0}, {3.0, 6.0}}));
    EXPECT_EQ(linalg::power(transpose(*mat_transpose_2), -1), linalg::transpose(power(*mat_transpose_2, -1)));
}

// обращение матрицы тест
TEST_F(Matrix_test_inversion, InvertibleMatrix) {

    EXPECT_EQ(linalg::invert(*mat_invertible), linalg::Matrix({{0.6, -0.7}, {-0.2, 0.4}}));
}

TEST_F(Matrix_test_inversion, SingularMatrix) {
    // проверка обращения вырожденной матрицы
    EXPECT_THROW(linalg::invert(*mat_singular), linalg::Singular_matrix);
}

// умножение на скаляр тест
TEST_F(Matrix_test_scalar_multiplication, ScalarMultiplication) {

    EXPECT_EQ(*mat_scalar * 2.0, linalg::Matrix({{2.0, 4.0}, {6.0, 8.0}}));
    EXPECT_EQ(*mat_scalar * 0.0, linalg::Matrix({{0.0, 0.0}, {0.0, 0.0}}));
}

// определитель тест
TEST_F(Matrix_test_determinant, MatrixDeterminant) {
    EXPECT_DOUBLE_EQ(mat_determinant->det(), -2.0);
}

// умножение с большими числами тест
TEST_F(Matrix_test_large_numbers, LargeNumberMatrix) {
    EXPECT_EQ(*mat_large * 2.0, linalg::Matrix({{2e18, 4e18}, {6e18, 8e18}}));
}

// умножение с маленькими числами тест
TEST_F(Matrix_test_large_numbers, SmallNumberMatrix) {

    EXPECT_EQ(*mat_small * 2.0, linalg::Matrix({{2e-18, 4e-18}, {6e-18, 8e-18}}));
}

// след матрицы тест
TEST_F(Matrix_test_trace, SquareMatrixTrace) {

    EXPECT_DOUBLE_EQ(mat_square->trace(), 2.0);
}

TEST_F(Matrix_test_trace, NonSquareMatrixTrace) {
    // проверка следа для неквадратной матрицы
    EXPECT_THROW(mat_non_square->trace(), linalg::Wrong_matrix_size);
}

// сложение с нулевой матрицей тест
TEST_F(Matrix_test_zero_matrix_addition, ZeroMatrixAddition) {

    EXPECT_EQ(*mat1 + *zero_matrix, *mat1);
}

// выход за пределы тест
TEST_F(Matrix_test_out_of_bounds, OutOfBoundsAccess) {

    EXPECT_THROW((*mat_bounds)(3, 0), linalg::Wrong_matrix_size);
}

// умножение неквадратной матрицы на единичную тест
TEST_F(Matrix_test_multiplication_identity_non_square, NonSquareMatrixMultiplication) {

    EXPECT_EQ(*mat_non_square * *identity_non_square, linalg::Matrix({{1.0, 2.0}, {4.0, 5.0}}));
}

// умножение пустых матриц тест
TEST_F(Matrix_test_empty_multiplication, EmptyMatrixMultiplication) {

    EXPECT_NO_THROW(*empty_matrix1 * *empty_matrix2);
}

// операции с разреженной матрицей тест
TEST_F(Matrix_test_sparse_matrices, SparseMatrixOperations) {

    EXPECT_EQ(*sparse_matrix * 2.0, linalg::Matrix({{0.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 0.0}}));
}

// операции с матрицей из единиц тест
TEST_F(Matrix_test_all_ones, AllOnesMatrixOperations) {

    EXPECT_EQ(*all_ones_matrix * 2.0, linalg::Matrix({{2.0, 2.0}, {2.0, 2.0}}));
}

// транспонирование строки тест
TEST_F(Matrix_test_single_row, SingleRowMatrix) {

    linalg::Matrix expected_result = std::initializer_list<std::initializer_list<double>>{{1.0}, {2.0}, {3.0}};
    EXPECT_EQ(transpose(*single_row_matrix), expected_result);
}

TEST_F(Matrix_test_single_row, SingleColumnMatrix) {
    // проверка транспонирования строки
    EXPECT_EQ(transpose(*single_row_matrix), linalg::Matrix({1.0, 2.0, 3.0}));
}

// точность вычислений тест
TEST_F(Matrix_test_floating_point_precision, FloatingPointPrecision) {

    linalg::Matrix result = *mat1 - *mat2;
    linalg::Matrix expected_result({{-1e-7, -1e-7}, {-1e-7, -1e-7}});
    EXPECT_NEAR(result(0, 0), expected_result(0, 0), 1e-7);
}

// изменение размеров матрицы тест
TEST_F(Matrix_test_reshape, RowColumnModification) {

    mat_modify->reshape(3, 2);
    EXPECT_EQ((*mat_modify), linalg::Matrix({{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}}));
}

// умножение на нулевую матрицу тест
TEST_F(Matrix_test_zero_matrix_multiplication, ZeroMatrixMultiplication) {

    EXPECT_EQ(*mat * *zero_matrix, linalg::Matrix({{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}));
}

//неверная инициализация
TEST_F(Matrix_Test, WrongInitialization) {

    EXPECT_THROW(linalg::Matrix({{1, 2, 3}, {1, 2}}), linalg::Wrong_matrix_size);
}

// вычисление нормы
TEST_F(Matrix_test_norm, TestNorm) {

    EXPECT_NEAR(mat->norm(), 16.88194, 1e-5);
}




