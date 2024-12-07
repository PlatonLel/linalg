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
        Matrix() noexcept {}
        Matrix(const size_t& rows);
        Matrix(const size_t& rows, const size_t& columns);
        template<typename Y>
        Matrix(const Matrix<Y>& m) { copy_constructor(m); }
        Matrix(const Matrix& m) { copy_constructor(m); }

        Matrix(Matrix&& m) noexcept { swap(m); }
        template<typename Y>
        Matrix(std::initializer_list<std::initializer_list<Y>> m);

        template<typename Y>
        Matrix(std::initializer_list<Y> m);

        ~Matrix() noexcept;

        const T* begin() const noexcept { return m_ptr; }
        const T* end() const noexcept { return m_ptr + m_size; }
        T* begin() noexcept { return m_ptr; }
        T* end() noexcept { return m_ptr + m_size; }
        const size_t size() const noexcept { return m_size; }
        const size_t capacity() const noexcept { return m_capacity; }
        const size_t& rows() const noexcept { return m_rows; }

        const size_t& columns() const noexcept { return m_columns; }

        T* get_ptr() const noexcept { return m_ptr; }

        void swap_rows(const size_t& row1,const size_t& row2) noexcept;

        void swap(Matrix& m) noexcept;

        bool empty() const noexcept { return ((m_ptr == nullptr)||(m_columns==0 && m_rows==0)); }

        void reshape(const size_t& rows,const size_t& new_columns);

        T norm() const;

        T trace() const;

        T det() const;
//метод гаусс возвращает количество раз, когда строки менялись местами, он нужен для det(), но можно вызвать для любой матрицы
        size_t gauss() noexcept;
        Matrix& operator=(Matrix&& m) noexcept;

        template <typename Y>
        Matrix& operator=(const Matrix<Y>& m);

        Matrix &operator=(const Matrix& m);
        template <typename Y>
        Matrix operator+(const Matrix<Y>& m) const;

        template <typename Y>
        Matrix& operator+=(const Matrix<Y>& m) ;

        template <typename Y>
        Matrix operator-(const Matrix<Y>& m) const;

        template <typename Y>
        Matrix& operator-=(const Matrix<Y>& m);

        template <typename Y>
        Matrix &operator*=(const Matrix<Y>& m);

        template <typename Y>
        Matrix &operator*=(const Y& v) noexcept;

        template <typename Y>
        bool operator==(const Matrix<Y>& m) const;

        bool operator==(const Matrix<double>& m) const;

        bool operator!=(const Matrix<double>& m) const;

        template <typename Y>
        bool operator!=(const Matrix<Y>& m) const;

        T &operator()(const size_t& m_row, const size_t& m_column);

        T operator()(const size_t& m_row, const size_t& m_column) const;

        const T& operator[](size_t i) const noexcept { return m_ptr[i]; }
        T& operator[](size_t i) noexcept { return m_ptr[i]; }

        void shrink_to_fit();

        void clear() noexcept;

        void reserve(size_t n);
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
    template <typename T>
    Matrix<T> operator*(const T& v, const Matrix<T>& m);
    template <typename T>
    Matrix<T> operator*( const Matrix<T>& m, const T& v);
    template <typename T>
    Matrix<T> operator*(const Matrix<T>& m1, const Matrix<T>& m2);

    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Matrix<T>& m);

    template <typename T>
    Matrix<T> power(const Matrix<T>& m,const int& p);

    template <typename T>
    Matrix<T> concatenate(const Matrix<T>& m_left, const Matrix<T>& m_right);

    template <typename T>
    Matrix<T> transpose(const Matrix<T>& m);

    template <typename T>
    Matrix<T> invert(const Matrix<T>& m);
//// ниже находятся 4 функции, которые вычисляют элементы матриц L и U соответственно в разложении матрицы, которая передается в invert
////а также проихводят прямую подстановку единичного вектора и обратную подстановку вектора, который является результатом прямой
////ф
    template <typename T>
    T get_sum_L(size_t& m_row, size_t& m_column, Matrix<T>& m_L, Matrix<T>& m_U);

    template <typename T>
    T get_sum_U(size_t& m_row, size_t& m_column, Matrix<T>& m_L, Matrix<T>& m_U);

    template <typename T>
    Matrix<T> forward_substitution(const Matrix<T>& m_L, const Matrix<T>& b);

    template <typename T>
    Matrix<T> backward_substitution(const Matrix<T>& m_U, const Matrix<T>& y);
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

#include"Linalg.hpp"