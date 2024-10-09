#include "Linalg.h"
void Linalg::Matrix::reshape(size_t new_m_rows, size_t new_m_columns) {
    if (new_m_columns * new_m_rows == m_columns * m_rows) {
        m_columns = new_m_columns;
        m_rows = new_m_rows;
    }
}
Linalg::Matrix::Matrix(const Matrix& m) {
    m_ptr = new double[m.m_rows * m.m_columns];
    for (size_t i=0; i<(m.m_rows*m.m_columns); ++i) {
        m_ptr[i] = m.m_ptr[i];
    }
    m_rows = m.m_rows;
    m_columns = m.m_columns;
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
};

double Linalg::Matrix::norm() {
    if (!empty()) {
        double m_norm=0;
        for (size_t i=0; i<(m_rows*m_columns); ++i) {
            m_norm+=(m_ptr[i])*(m_ptr[i]);
        }
        return m_norm;
    }
};
//Matrix& race();
double Linalg::Matrix::det() {
    if (m_rows != m_columns){
        return 52;
    }
    double m_det=0;
    for (size_t i=0; i<=(m_rows); ++i) {m_det += m_ptr[m_rows*i + i];}
    return m_det;
};
Linalg::Matrix::Matrix(size_t& rows, size_t& columns, std::initializer_list<double> m): m_rows{rows}, m_columns{columns}{
        if (m.size()!= m_rows*m_columns) {
            m_ptr= nullptr;
            m_rows = 0;
            m_columns = 0;
        }
        m_ptr = new double[m_rows*m_columns];
        size_t i = 0;
        for (double value : m) {
            m_ptr[i++] = value;
        }
};
//Matrix& rank();