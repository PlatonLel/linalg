#include "../Include/test.h"

// сложение тест
TEST_F(Matrix_test_addition_and_subtraction, TestMatrixAddition) {

    EXPECT_EQ(*mat_sum1 + *mat_sum2, Linalg::Matrix({{10.0,10.0,10.0},{10.0,10.0,10.0},{10.0,10.0,10.0}}));
    EXPECT_EQ(*mat_sum2 + *mat_sum1, Linalg::Matrix({{10.0,10.0,10.0},{10.0,10.0,10.0},{10.0,10.0,10.0}}));
}

//вычитание тест
TEST_F(Matrix_test_addition_and_subtraction, TestMatrixSubtraction) {

    EXPECT_EQ(*mat_sum1 - *mat_sum2, Linalg::Matrix({{-8.0,-6.0,-4.0},{-2.0,0.0,2.0},{4.0,6.0,8.0}}));
    EXPECT_EQ(*mat_sum2 - *mat_sum1, Linalg::Matrix({{8.0,6.0,4.0},{2.0,0.0,-2.0},{-4.0,-6.0,-8.0}}));
    EXPECT_EQ(*mat_sum2 - *mat_sum2, Linalg::Matrix({{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}}));
}

//возведение в степень тест
TEST_F(Matrix_test_pow, TestMatrixPower) {

    Linalg::Matrix result_pow_n = {{239.0/300.0,161.0/900.0,-31.0/180.0},{113.0/300.0,587.0/900.0,-37.0/180.0},{-5.0/12.0,-11.0/36.0,5.0/36.0}};

    EXPECT_EQ(Linalg::power(*mat_pow_e,10), *mat_pow_e);
    EXPECT_EQ(Linalg::power(*mat_pow_e,-10), *mat_pow_e);
    EXPECT_EQ(Linalg::power(*mat_pow,2), Linalg::Matrix({{25.0, 25.0, 68.0}, {30.0, 35.0, 89.0}, {141.0, 152.0, 407.0}}));
    EXPECT_EQ(Linalg::power(*mat_pow,-2), result_pow_n);
}

TEST_F(Matrix_test_multiplication, MatrixMultiplication) {

    EXPECT_EQ(*mat_mult1 * *mat_mult2, Linalg::Matrix({{4.0, 4.0}, {10.0, 8.0}}));
}

TEST_F(Matrix_test_transpose, MatrixTranspose) {

    EXPECT_EQ(Linalg::transpose(*mat_transpose), Linalg::Matrix({{1.0, 4.0}, {2.0, 5.0}, {3.0, 6.0}}));
}

TEST_F(Matrix_test_inversion, InvertibleMatrix) {

    EXPECT_EQ(Linalg::invert(*mat_invertible), Linalg::Matrix({{0.6, -0.7}, {-0.2, 0.4}}));
}

TEST_F(Matrix_test_inversion, SingularMatrix) {

    EXPECT_THROW(Linalg::invert(*mat_singular), Singular_matrix);
}

TEST_F(Matrix_test_scalar_multiplication, ScalarMultiplication) {

    EXPECT_EQ(*mat_scalar * 2.0, Linalg::Matrix({{2.0, 4.0}, {6.0, 8.0}}));
    EXPECT_EQ(*mat_scalar * 0.0, Linalg::Matrix({{0.0, 0.0}, {0.0, 0.0}}));
}

TEST_F(Matrix_test_determinant, MatrixDeterminant) {

    EXPECT_DOUBLE_EQ(mat_determinant->det(), -2.0);
}

TEST_F(Matrix_test_large_numbers, LargeNumberMatrix) {

    EXPECT_EQ(*mat_large * 2.0, Linalg::Matrix({{2e18, 4e18}, {6e18, 8e18}}));
}

TEST_F(Matrix_test_large_numbers, SmallNumberMatrix) {

    EXPECT_EQ(*mat_small * 2.0, Linalg::Matrix({{2e-18, 4e-18}, {6e-18, 8e-18}}));
}

TEST_F(Matrix_test_trace, SquareMatrixTrace) {
    EXPECT_DOUBLE_EQ(mat_square->trace(), 2.0);
}

TEST_F(Matrix_test_trace, NonSquareMatrixTrace) {
    EXPECT_THROW(mat_non_square->trace(), Wrong_matrix_size);
}


TEST_F(Matrix_test_zero_matrix_addition, ZeroMatrixAddition) {

    EXPECT_EQ(*mat1 + *zero_matrix, *mat1);
}

TEST_F(Matrix_test_out_of_bounds, OutOfBoundsAccess) {

    EXPECT_THROW((*mat_bounds)(3, 0), Wrong_matrix_size);
}

TEST_F(Matrix_test_multiplication_identity_non_square, NonSquareMatrixMultiplication) {

    EXPECT_EQ(*mat_non_square * *identity_non_square, Linalg::Matrix({{1.0, 2.0}, {4.0, 5.0}}));
}

TEST_F(Matrix_test_empty_multiplication, EmptyMatrixMultiplication) {

    EXPECT_NO_THROW(*empty_matrix1 * *empty_matrix2);
}

TEST_F(Matrix_test_sparse_matrices, SparseMatrixOperations) {

    EXPECT_EQ(*sparse_matrix * 2.0, Linalg::Matrix({{0.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 0.0}}));
}

TEST_F(Matrix_test_all_ones, AllOnesMatrixOperations) {

    EXPECT_EQ(*all_ones_matrix * 2.0, Linalg::Matrix({{2.0, 2.0}, {2.0, 2.0}}));
}

TEST_F(Matrix_test_single_row_or_column, SingleRowMatrix) {
    Linalg::Matrix expected_result = std::initializer_list<std::initializer_list<double>>{{1.0}, {2.0}, {3.0}};
    EXPECT_EQ(transpose(*single_row_matrix), expected_result);
}

TEST_F(Matrix_test_single_row_or_column, SingleColumnMatrix) {

    EXPECT_EQ(transpose(*single_column_matrix), Linalg::Matrix({{1.0, 2.0, 3.0}}));
}

TEST_F(Matrix_test_floating_point_precision, FloatingPointPrecision) {

    Linalg::Matrix result = *mat1 - *mat2;
    Linalg::Matrix expected_result({{-1e-7, -1e-7}, {-1e-7, -1e-7}});
    EXPECT_NEAR(result(0, 0), expected_result(0, 0), 1e-7);
}

TEST_F(Matrix_test_reshape, RowColumnModification) {

    mat_modify->reshape(3, 2);
    EXPECT_EQ((*mat_modify), Linalg::Matrix({{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}}));
}


TEST_F(Matrix_test_zero_matrix_multiplication, ZeroMatrixMultiplication) {

    EXPECT_EQ(*mat * *zero_matrix , Linalg::Matrix({{0.0,0.0,0.0},{0.0,0.0,0.0}, {0.0,0.0,0.0}}));
}

