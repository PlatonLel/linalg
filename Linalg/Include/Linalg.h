#pragma once
#include <cmath> //std::sqrt, std::abs
#include <utility> //std::swap, std::move
#include <initializer_list> //std::initializer_list
#include <limits> //std::numeric_limits
#include <sstream>   // std::ostringstream
#include <iomanip>   // std::setw
#include <algorithm> // std::max
#include <iostream> //std::cout


namespace linalg {
    template <typename T>
    class Matrix {
    public:
        Matrix() noexcept : m_rows{0}, m_columns{0}, m_ptr{nullptr} {}
        Matrix(const size_t& rows);
        Matrix(const size_t& rows, const size_t& columns);
        template<typename Y>
        Matrix(const Matrix<Y>& m) { copy_constructor(m); }
        Matrix(const Matrix& v) { copy_constructor(v); }
        Matrix(Matrix&& m) noexcept { swap(m); }

        Matrix(std::initializer_list<std::initializer_list<T>> m);

        Matrix(std::initializer_list<T> lst);

        ~Matrix() { delete[] m_ptr; }

        const T* begin() const noexcept { return m_ptr; }
        const T* end() const noexcept { return m_ptr + m_size; }
        T* begin() noexcept { return m_ptr; }
        T* end() noexcept { return m_ptr + m_size; }
        size_t size() const noexcept { return m_size; }
        size_t capacity() const noexcept { return m_capacity; }
        const size_t& rows() const noexcept { return m_rows; }

        const size_t& columns() const noexcept { return m_columns; }

        double *get_ptr() const noexcept { return m_ptr; }

        void swap_rows(const size_t& row1,const size_t& row2) noexcept;

        void swap(Matrix& m) noexcept;

        bool empty() const noexcept { return ((m_ptr == nullptr)||(m_columns==0 && m_rows==0)); }

        void reshape(const size_t& rows,const size_t& new_columns);

        double norm() const;

        double trace() const;

        double det() const;
//метод гаусс возвращает количество раз, когда строки менялись местами, он нужен для det(), но можно вызвать для любой матрицы
        size_t gauss() noexcept;
        template <typename Y>
        Matrix &operator=(Matrix<Y>&& m) noexcept;

        template <typename Y>
        Matrix &operator=(const Matrix<Y>& m);

        Matrix operator+(const Matrix& m) const;

        Matrix &operator+=(const Matrix& m) ;
        template <typename Y>
        Matrix &operator+=(const Matrix<Y>& m) ;

        Matrix operator-(const Matrix& m) const;

        Matrix &operator-=(const Matrix& m);

        Matrix &operator*=(const Matrix& m);

        Matrix &operator*=(const double& v) noexcept;

        bool operator==(const Matrix& m) const;

        bool operator!=(const Matrix& m) const;

        double &operator()(const size_t& m_row, const size_t& m_column);

        double operator()(const size_t& m_row, const size_t& m_column) const;

        const T& operator[](size_t i) const noexcept { return m_ptr[i]; }
        T& operator[](size_t i) noexcept { return m_ptr[i]; }
    private:
        template<typename Y>
        void copy_constructor(const Matrix<Y>& m);
    private:
        size_t m_rows = 0;
        size_t m_columns = 0;
        size_t m_size = 0;
        size_t m_capacity = 0;
        T* m_ptr = nullptr;
    };
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix);

    template <typename T>
    void swap(Matrix<T>& m1, Matrix<T>& m2) { m1.swap(m2); }
//
//    Matrix operator*(const double& v, const Matrix& m);
//
//    Matrix operator*( const Matrix& m, const double& v);
//
//    Matrix operator*(const Matrix& m1, const Matrix& m2);
//
//    std::ostream& operator<<(std::ostream& os, const Matrix& m);
//
//    Matrix power(const Matrix& m,const int& p);
//
//    Matrix concatenate(const Matrix& m_left, const Matrix& m_right);
//
//    Matrix transpose(const Matrix& m);
//
//    Matrix invert(const Matrix& m);
//// ниже находятся 4 функции, которые вычисляют элементы матриц L и U соответственно в разложении матрицы, которая передается в invert
////а также проихводят прямую подстановку единичного вектора и обратную подстановку вектора, который является результатом прямой
////ф
//    double get_sum_L(size_t& m_row, size_t& m_column, Matrix& m_L, Matrix& m_U);
//
//    double get_sum_U(size_t& m_row, size_t& m_column, Matrix& m_L, Matrix& m_U);
//
//    Matrix forward_substitution(const Matrix& m_L, const Matrix& b);
//
//    Matrix backward_substitution(const Matrix& m_U, const Matrix& y);
////класс исключений, от которого будут наследовать остальные
    class Matrix_exception {
    public:
        std::string what() {return description;}
    protected:
        std::string description;
    };
//в зависимости от параметра задается значение description
    class Wrong_matrix_size: public ::linalg::Matrix_exception {
    public:
        Wrong_matrix_size(size_t p);
    };

    class Empty_matrix: public ::linalg::Matrix_exception {
    public:
        Empty_matrix(size_t p);
    };

    class Singular_matrix: public ::linalg::Matrix_exception {
    public:
        Singular_matrix(size_t p);
    };
}

#include<Linalg.hpp>