#include "Linalg.h"
void Linalg::Matrix::reshape(size_t new_m_rows, size_t new_m_columns) {
    if (new_m_columns * new_m_rows == m_columns * m_rows) {
        m_columns = new_m_columns;
        m_rows = new_m_rows;
    }
}
Linalg::Matrix::Matrix(const Matrix& m) {
    *this = m;
}
Linalg::Matrix& Linalg::Matrix::operator = (const Matrix& m) {
    if (m_columns*m_rows != m.m_columns * m.m_rows) {
        delete[] m_ptr;
        m_ptr = new double[m.m_rows*m.m_columns];
        m_rows = m.m_rows;
        m_columns = m.m_columns;
    }
    for (size_t i = 0; i < (m.m_columns*m.m_rows); ++i) {
        m_ptr[i] = m.m_ptr[i];
    }
    return *this;
}

double Linalg::Matrix::norm() const {
    if (!empty()) {
        double m_norm=0;
        for (size_t i=0; i<(m_rows*m_columns); ++i) {
            m_norm += (m_ptr[i])*(m_ptr[i]);
        }
        return std::sqrt(m_norm);
    }
    return 0.52;
    //доделать
}

double Linalg::Matrix::det() const {
    if (m_rows != m_columns){
        return 52;
    }
    double m_det=0;
    for (size_t i=0; i<=(m_rows); ++i) {m_det += std::pow(-1,i)*m_ptr[m_rows*i + i];}
    return m_det;
    //доделать
}


Linalg::Matrix::Matrix(std::initializer_list<std::initializer_list<double>> m) {
    m_rows = m.size();
    m_columns = m.begin()->size();
    m_ptr = new double[m_rows*m_columns];
    size_t i=0;
    for (const auto& row : m) {
        for (double value : row) {
            m_ptr[i] = value;
            ++i;
        }
    }
    //доделать
}
double  Linalg::Matrix::trace() const {
    if (m_rows != m_columns){
        return 53;
    }
    double m_trace=0;
    for (size_t i=0; i<=(m_rows); ++i) {m_trace += m_ptr[m_rows*i + i];}
    return m_trace;
};
Linalg::Matrix Linalg::Matrix::operator+(const Matrix& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    Matrix m_sum(m_rows, m_columns);
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_sum.m_ptr[i] = m.m_ptr[i] + m_ptr[i];
    }
    return m_sum;
    //доделать
};
//double rank() const {
//
//};