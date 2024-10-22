#include "../Include/test.h"

// сложение тест
TEST_F(Matrix_test_addition_and_subtraction, TestMatrixAddition) {
    // проверка сложения двух матриц
    EXPECT_EQ(*mat_sum1 + *mat_sum2, Linalg::Matrix({{10.0,10.0,10.0},{10.0,10.0,10.0},{10.0,10.0,10.0}}));  // сложение mat_sum1 и mat_sum2
    EXPECT_EQ(*mat_sum2 + *mat_sum1, Linalg::Matrix({{10.0,10.0,10.0},{10.0,10.0,10.0},{10.0,10.0,10.0}}));  // обратное сложение
}

// вычитание тест
TEST_F(Matrix_test_addition_and_subtraction, TestMatrixSubtraction) {
    // проверка вычитания матриц
    EXPECT_EQ(*mat_sum1 - *mat_sum2, Linalg::Matrix({{-8.0,-6.0,-4.0},{-2.0,0.0,2.0},{4.0,6.0,8.0}}));  // вычитание mat_sum2 из mat_sum1
    EXPECT_EQ(*mat_sum2 - *mat_sum1, Linalg::Matrix({{8.0,6.0,4.0},{2.0,0.0,-2.0},{-4.0,-6.0,-8.0}}));  // обратное вычитание
    EXPECT_EQ(*mat_sum2 - *mat_sum2, Linalg::Matrix({{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}}));      // вычитание одинаковых матриц
}

// возведение в степень тест
TEST_F(Matrix_test_pow, TestMatrixPower) {
    // проверка возведения в степень

    EXPECT_EQ(Linalg::power(*mat_pow_e,10), *mat_pow_e);
    EXPECT_EQ(Linalg::power(*mat_pow_e,-10), *mat_pow_e);
    EXPECT_EQ(Linalg::power(power(*mat_pow_1,2),3), Linalg::power(power(*mat_pow_1,3),2));
    EXPECT_EQ(Linalg::power(*mat_pow_1,-1)*Linalg::power(*mat_pow_2,-1),
              Linalg::power(*mat_pow_2 * *mat_pow_1,-1));
}

// умножение тест
TEST_F(Matrix_test_multiplication, MatrixMultiplication) {
    // проверка умножения матриц
    EXPECT_EQ(*mat_mult1 * *mat_mult2, Linalg::Matrix({{4.0, 4.0}, {10.0, 8.0}}));  // умножение двух матриц
}

// транспонирование тест
TEST_F(Matrix_test_transpose, MatrixTranspose) {
    // проверка транспонирования
    EXPECT_EQ(Linalg::transpose(*mat_transpose), Linalg::Matrix({{1.0, 4.0}, {2.0, 5.0}, {3.0, 6.0}}));
    EXPECT_EQ(Linalg::power(transpose(*mat_transpose),-1), Linalg::transpose(power(*mat_transpose,-1)));
}

// обращение матрицы тест
TEST_F(Matrix_test_inversion, InvertibleMatrix) {
    // проверка обращения обратимой матрицы
    EXPECT_EQ(Linalg::invert(*mat_invertible), Linalg::Matrix({{0.6, -0.7}, {-0.2, 0.4}}));  // обращение матрицы
}

TEST_F(Matrix_test_inversion, SingularMatrix) {
    // проверка обращения вырожденной матрицы
    EXPECT_THROW(Linalg::invert(*mat_singular), Linalg::Singular_matrix);  // исключение для сингулярной матрицы
}

// умножение на скаляр тест
TEST_F(Matrix_test_scalar_multiplication, ScalarMultiplication) {
    // проверка умножения на скаляр
    EXPECT_EQ(*mat_scalar * 2.0, Linalg::Matrix({{2.0, 4.0}, {6.0, 8.0}}));  // умножение на 2
    EXPECT_EQ(*mat_scalar * 0.0, Linalg::Matrix({{0.0, 0.0}, {0.0, 0.0}}));  // умножение на 0
}

// определитель тест
TEST_F(Matrix_test_determinant, MatrixDeterminant) {
    // проверка вычисления определителя
    EXPECT_DOUBLE_EQ(mat_determinant->det(), -2.0);  // вычисление определителя
}

// умножение с большими числами тест
TEST_F(Matrix_test_large_numbers, LargeNumberMatrix) {
    // проверка умножения с большими числами
    EXPECT_EQ(*mat_large * 2.0, Linalg::Matrix({{2e18, 4e18}, {6e18, 8e18}}));  // умножение матрицы с большими числами
}

// умножение с маленькими числами тест
TEST_F(Matrix_test_large_numbers, SmallNumberMatrix) {
    // проверка умножения с маленькими числами
    EXPECT_EQ(*mat_small * 2.0, Linalg::Matrix({{2e-18, 4e-18}, {6e-18, 8e-18}}));  // умножение матрицы с маленькими числами
}

// след матрицы тест
TEST_F(Matrix_test_trace, SquareMatrixTrace) {
    // проверка вычисления следа квадратной матрицы
    EXPECT_DOUBLE_EQ(mat_square->trace(), 2.0);  // вычисление следа квадратной матрицы
}

TEST_F(Matrix_test_trace, NonSquareMatrixTrace) {
    // проверка следа для неквадратной матрицы
    EXPECT_THROW(mat_non_square->trace(), Linalg::Wrong_matrix_size);  // исключение для неквадратной матрицы
}

// сложение с нулевой матрицей тест
TEST_F(Matrix_test_zero_matrix_addition, ZeroMatrixAddition) {
    // проверка сложения с нулевой матрицей
    EXPECT_EQ(*mat1 + *zero_matrix, *mat1);  // сложение матрицы с нулевой
}

// выход за пределы тест
TEST_F(Matrix_test_out_of_bounds, OutOfBoundsAccess) {
    // проверка выхода за пределы матрицы
    EXPECT_THROW((*mat_bounds)(3, 0), Linalg::Wrong_matrix_size);  // исключение для доступа за пределы
}

// умножение неквадратной матрицы на единичную тест
TEST_F(Matrix_test_multiplication_identity_non_square, NonSquareMatrixMultiplication) {
    // проверка умножения неквадратной матрицы на единичную
    EXPECT_EQ(*mat_non_square * *identity_non_square, Linalg::Matrix({{1.0, 2.0}, {4.0, 5.0}}));  // умножение неквадратной матрицы
}

// умножение пустых матриц тест
TEST_F(Matrix_test_empty_multiplication, EmptyMatrixMultiplication) {
    // проверка умножения пустых матриц
    EXPECT_NO_THROW(*empty_matrix1 * *empty_matrix2);  // умножение пустых матриц не выбрасывает исключение
}

// операции с разреженной матрицей тест
TEST_F(Matrix_test_sparse_matrices, SparseMatrixOperations) {
    // проверка операций с разреженной матрицей
    EXPECT_EQ(*sparse_matrix * 2.0, Linalg::Matrix({{0.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 0.0}}));  // умножение разреженной матрицы
}

// операции с матрицей из единиц тест
TEST_F(Matrix_test_all_ones, AllOnesMatrixOperations) {
    // проверка операций с матрицей из единиц
    EXPECT_EQ(*all_ones_matrix * 2.0, Linalg::Matrix({{2.0, 2.0}, {2.0, 2.0}}));  // умножение матрицы из единиц
}

// транспонирование строки тест
TEST_F(Matrix_test_single_row, SingleRowMatrix) {
    Linalg::Matrix expected_result = std::initializer_list<std::initializer_list<double>>{{1.0}, {2.0}, {3.0}};
    EXPECT_EQ(transpose(*single_row_matrix), expected_result);
}

TEST_F(Matrix_test_single_row, SingleColumnMatrix) {
    // проверка транспонирования строки
    EXPECT_EQ(transpose(*single_row_matrix), Linalg::Matrix({1.0, 2.0, 3.0}));
}

// точность вычислений тест
TEST_F(Matrix_test_floating_point_precision, FloatingPointPrecision) {
    // проверка точности вычислений
    Linalg::Matrix result = *mat1 - *mat2;
    Linalg::Matrix expected_result({{-1e-7, -1e-7}, {-1e-7, -1e-7}});
    EXPECT_NEAR(result(0, 0), expected_result(0, 0), 1e-7);
}

// изменение размеров матрицы тест
TEST_F(Matrix_test_reshape, RowColumnModification) {
    // проверка изменения количества строк/столбцов
    mat_modify->reshape(3, 2);
    EXPECT_EQ((*mat_modify), Linalg::Matrix({{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}}));  // изменение размеров матрицы
}

// умножение на нулевую матрицу тест
TEST_F(Matrix_test_zero_matrix_multiplication, ZeroMatrixMultiplication) {
    // проверка умножения на нулевую матрицу
    EXPECT_EQ(*mat * *zero_matrix, Linalg::Matrix({{0.0,0.0,0.0},{0.0,0.0,0.0}, {0.0,0.0,0.0}}));  // результат - нулевая матрица
}

//неверная инициализация
TEST_F(Matrix_Test, WrongInitialization) {
    EXPECT_THROW(Linalg::Matrix({{1,2,3},{1,2}}),Linalg::Wrong_matrix_size);
}

// вычисление нормы
TEST_F(Matrix_test_norm, TestNorm) {
    EXPECT_NEAR(mat->norm(), 16.88194, 1e-5);
}




