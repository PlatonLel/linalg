#pragma once
#include <cstdio>
#include <cmath>
#include <iostream>

namespace Linalg {
    class Matrix {
    public:
        Matrix(): m_rows{0}, m_columns{0}, m_ptr{nullptr} {}
        Matrix(size_t& rows): m_rows{rows}, m_columns{0} {m_ptr = new double[rows];}
        Matrix(size_t& rows, size_t& columns): m_rows{rows}, m_columns{columns}{m_ptr = new double[columns*rows];}
        Matrix(size_t&& rows): m_rows{rows}, m_columns{0} {m_ptr = new double[rows];}
        Matrix(size_t&& rows, size_t&& columns): m_rows{rows}, m_columns{columns}{m_ptr = new double[columns*rows];}
        Matrix(const Matrix& m);
        Matrix(size_t& rows, size_t& columns, std::initializer_list<double> m);
        Matrix(std::initializer_list<std::initializer_list<double>> m);
        ~Matrix() {delete[] m_ptr;}
        size_t get_rows() const {return m_rows;}
        size_t get_columns() const {return m_columns;}
        double* get_ptr() const {return m_ptr;}
        bool empty() const {return (m_ptr == nullptr);}
        void reshape(size_t new_m_columns, size_t new_m_rows);
        const double& operator[](size_t i) const { return m_ptr[i]; }
        double& operator[](size_t i) { return m_ptr[i]; }
        Matrix& operator = (const Matrix& v);
        double norm();
        Matrix& race();
        double det();
        Matrix& rank();
    private:
        size_t m_rows;
        size_t m_columns;
        double* m_ptr;
    };

}
