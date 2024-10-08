#pragma once
#include <cstdio>
#include <iostream>

namespace Linalg {
    class Matrix {
    public:
        Matrix(): m_rows{0}, m_columns{0}, m_ptr{nullptr} {}
        Matrix(size_t& rows): m_rows{rows}, m_columns{0} {m_ptr = new double[rows];}
        Matrix(size_t& rows, size_t& columns): m_rows{rows}, m_columns{columns}{m_ptr = new double[columns*rows];}
        Matrix(size_t&& rows): m_rows{rows}, m_columns{0} {m_ptr = new double[rows];}
        Matrix(size_t&& rows, size_t&& columns): m_rows{rows}, m_columns{columns}{m_ptr = new double[columns*rows];}
        ~Matrix() {delete[] m_ptr;}
        size_t get_rows() const {return m_rows;}
        size_t get_columns() const {return m_columns;}
        double* get_ptr() const {return m_ptr;}
        bool empty() const {return (m_columns == 0 & m_rows == 0);}
        void reshape(size_t new_m_columns, size_t new_m_rows);
    private:
        size_t m_rows;
        size_t m_columns;
        double* m_ptr;
    };

}
