#pragma once
#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>


namespace Linalg {
    class Matrix {
    public:
        Matrix() noexcept : m_rows{0}, m_columns{0}, m_ptr{nullptr} {}

        Matrix(const size_t& rows) : m_rows{rows}, m_columns{1} { m_ptr = new double[rows]; }

        Matrix(const size_t& rows, const size_t& columns) : m_rows{rows}, m_columns{columns} {m_ptr = new double[columns * rows];}

        Matrix(Matrix&& m) noexcept;

        Matrix(const Matrix& m);

        Matrix(std::initializer_list<std::initializer_list<double>> m);

        Matrix(std::initializer_list<double> m);

        ~Matrix() { delete[] m_ptr; }

        const size_t& get_rows() const noexcept { return m_rows; }

        const size_t& get_columns() const noexcept { return m_columns; }

        double *get_ptr() const noexcept { return m_ptr; }

        void swap_rows(const size_t& row1,const size_t& row2);

        bool empty() const noexcept { return ((m_ptr == nullptr)||(m_columns==0 && m_rows==0)); }

        void reshape(const size_t& new_m_columns,const size_t& new_m_rows);

        double norm() const;

        double trace() const;

        double det() const;

        void print() const;

        size_t gauss();

        Matrix &operator=(const Matrix& m);

        Matrix &operator=(Matrix&& m) noexcept;

        Matrix operator+(const Matrix& m) const;

        Matrix &operator+=(const Matrix& m);

        Matrix operator-(const Matrix& m) const;

        Matrix &operator-=(const Matrix& m);

        Matrix operator*(const Matrix& m) const;

        Matrix &operator*=(const Matrix& m);

        Matrix operator+(Matrix&& m) const;

        Matrix &operator+=(Matrix&& m);

        Matrix operator-(Matrix&& m) const;

        Matrix &operator-=(Matrix&& m);

        Matrix operator*(Matrix&& m) const;

        Matrix &operator*=(Matrix&& m);

        Matrix operator*(const double& v);

        Matrix &operator*=(const double& v);

        bool operator==(const Matrix& m) const;

        bool operator==(Matrix&& m) const;

        bool operator!=(const Matrix& m) const;

        bool operator!=(Matrix&& m) const;

        double &operator()(const size_t& m_row, const size_t& m_column);

        double operator()(const size_t& m_row, const size_t& m_column) const;

//        double rank() const;
    private:
        size_t m_rows;
        size_t m_columns;
        double *m_ptr;
    };

    std::ostream& operator<<(std::ostream& os, const Matrix& matrix);

    Matrix operator*(const double& v, const Matrix& m);

    Matrix operator*(const double& v, Matrix&& m);

    std::ostream& operator<<(std::ostream& os, const Matrix& m);

    std::ostream& operator<<(std::ostream& os, Matrix&& m);

    Matrix power(const Matrix& m, int& power);

    Matrix power(Matrix&& m, int& power);

    Matrix power(const Matrix& m, int&& power);

    Matrix power(Matrix&& m, int&& power);

    Matrix concatenate(const Matrix& m_left, const Matrix& m_right);

    Matrix concatenate(Matrix&& m_left, const Matrix& m_right);

    Matrix concatenate(Matrix&& m_left, Matrix&& m_right);

    Matrix concatenate(const Matrix& m_left, Matrix&& m_right);

    Matrix transpose(const Matrix& m);

    Matrix transpose(Matrix&& m);

    Matrix invert(const Matrix& m);

    Matrix invert(Matrix&& m);

    double get_sum_L(size_t& m_row, size_t& m_column, Matrix& m_L, Matrix& m_U);

    double get_sum_U(size_t& m_row, size_t& m_column, Matrix& m_L, Matrix& m_U);

    Matrix forward_substitution(const Matrix& m_L, const Matrix& b);

    Matrix backward_substitution(const Matrix& m_U, const Matrix& y);
}

class Wrong_matrix_size {
public:
    Wrong_matrix_size() {};
};

class Empty_matrix {
public:
    Empty_matrix() {};
};

class Singular_matrix {
public:
    Singular_matrix() {};
};