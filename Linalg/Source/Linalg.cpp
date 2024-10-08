#include "Linalg.h"
void Linalg::Matrix::reshape(size_t new_m_columns, size_t new_m_rows) {
    if (new_m_columns * new_m_rows == m_columns * m_rows) {
        m_columns = new_m_columns;
        m_rows = new_m_rows;
    }
}