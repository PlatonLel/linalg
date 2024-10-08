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