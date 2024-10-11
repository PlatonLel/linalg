#include "Linalg.h"

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

Linalg::Matrix::Matrix(std::initializer_list<double> m) {
    m_columns = 1;
    m_rows = m.size();
    m_ptr = new double[m_rows];
    size_t i=0;
    for (double value : m) {
        m_ptr[i] = value;
        ++i;
    }
}

Linalg::Matrix::Matrix(Linalg::Matrix&& m) noexcept {
    m_ptr = m.m_ptr;
    m.m_ptr = nullptr;
    m_columns = m.m_columns;
    m_rows = m.m_rows;
}

Linalg::Matrix::Matrix(const Matrix& m) {
    m_ptr = new double[m.m_columns*m.m_rows];
    for (size_t i = 0; i<(m.m_columns*m.m_rows); ++i) {
        m_ptr[i] = m.m_ptr[i];
    }
    m_columns = m.m_columns;
    m_rows = m.m_rows;
}

Linalg::Matrix Linalg::Matrix::operator*(const Matrix& m) const {
    if (m_columns != m.m_rows) {
        return {};
    }
    Matrix m_return(m_rows, m.m_columns);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m.m_columns; ++j) {
            m_return.m_ptr[i * m.m_columns + j] = 0;
            for (size_t k = 0; k < m_columns; ++k) {
                m_return.m_ptr[i * m.m_columns + j] += m_ptr[i * m_columns + k] * m.m_ptr[k * m.m_columns + j];
            }
        }
    }

    return m_return;
}

Linalg::Matrix Linalg::Matrix::operator*(Matrix&& m) const {
    if (m_columns != m.m_rows) {
        return {};
    }
    Matrix m_return(m_rows, m.m_columns);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m.m_columns; ++j) {
            m_return.m_ptr[i * m.m_columns + j] = 0;
            for (size_t k = 0; k < m_columns; ++k) {
                m_return.m_ptr[i * m.m_columns + j] += m_ptr[i * m_columns + k] * m.m_ptr[k * m.m_columns + j];
            }
        }
    }

    return m_return;
}

Linalg::Matrix& Linalg::Matrix::operator*=(const Matrix& m) {
    if (m_columns != m.m_rows) {
        return *this;
    }
    Matrix m_return(m_rows, m.m_columns);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m.m_columns; ++j) {
            m_return.m_ptr[i * m.m_columns + j] = 0;
            for (size_t k = 0; k < m_columns; ++k) {
                m_return.m_ptr[i * m.m_columns + j] += m_ptr[i * m_columns + k] * m.m_ptr[k * m.m_columns + j];
            }
        }
    }
    *this = std::move(m_return);
    return *this;
}

Linalg::Matrix& Linalg::Matrix::operator*=(Matrix&& m) {
    if (m_columns != m.m_rows) {
        return *this;
    }
    Matrix m_return(m_rows, m.m_columns);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m.m_columns; ++j) {
            m_return.m_ptr[i * m.m_columns + j] = 0;
            for (size_t k = 0; k < m_columns; ++k) {
                m_return.m_ptr[i * m.m_columns + j] += m_ptr[i * m_columns + k] * m.m_ptr[k * m.m_columns + j];
            }
        }
    }
    *this = std::move(m_return);
    return *this;
}

Linalg::Matrix Linalg::Matrix::operator*(double v) {
    Matrix m_return(*this);
    for (size_t i=0; i<m_columns*m_rows; ++i) {
        m_return.m_ptr[i] *= v;
    }
    return m_return;
};

Linalg::Matrix Linalg::operator*(double v, Linalg::Matrix& m) {
    Linalg::Matrix m_return(m);
    for (size_t i=0; i<m.get_columns()*m.get_rows(); ++i) {
        m_return.get_ptr()[i] *= v;
    }
    return m_return;
};

Linalg::Matrix Linalg::operator*(double v, Linalg::Matrix&& m) {
    Linalg::Matrix m_return = std::move(m);
    for (size_t i=0; i<m.get_columns()*m.get_rows(); ++i) {
        m_return.get_ptr()[i] *= v;
    }
    return m_return;
};


Linalg::Matrix& Linalg::Matrix::operator*=(double v) {
    for (size_t i=0; i<m_columns*m_rows; ++i) {
        m_ptr[i] *= v;
    }
    return *this;
};


Linalg::Matrix& Linalg::Matrix::operator=(const Matrix& m) {
    if (this == &m) {
        return *this;
    }
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

Linalg::Matrix& Linalg::Matrix::operator=(Matrix&& m) {
    if (this == &m) {
        return *this;
    }
    if (m_columns * m_rows != m.m_columns * m.m_rows) {
        delete[] m_ptr;
        m_ptr = new double[m.m_rows * m.m_columns];
        m_rows = m.m_rows;
        m_columns = m.m_columns;
    }
    for (size_t i = 0; i < (m.m_columns * m.m_rows); ++i) {
        m_ptr[i] = m.m_ptr[i];
    }
    return *this;
}

Linalg::Matrix Linalg::Matrix::operator+(const Matrix& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    Matrix m_sum(m_rows, m_columns);
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_sum.m_ptr[i] = m.m_ptr[i] + m_ptr[i];
    }
    return m_sum;
    //доделать
};

Linalg::Matrix Linalg::Matrix::operator+(Matrix&& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    Matrix m_sum(m_rows, m_columns);
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_sum.m_ptr[i] = m.m_ptr[i] + m_ptr[i];
    }
    return m_sum;
    //доделать
};

Linalg::Matrix& Linalg::Matrix::operator+=(const Matrix& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] += m.m_ptr[i];
    }
    return *this;
    //доделать
};

Linalg::Matrix& Linalg::Matrix::operator+=(Matrix&& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] += m.m_ptr[i];
    }
    return *this;
    //доделать
};

Linalg::Matrix Linalg::Matrix::operator-(const Matrix& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    Matrix m_sum(m_rows, m_columns);
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_sum.m_ptr[i] = m.m_ptr[i] - m_ptr[i];
    }
    return m_sum;
    //доделать
};

Linalg::Matrix Linalg::Matrix::operator-(Matrix&& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    Matrix m_sum(m_rows, m_columns);
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_sum.m_ptr[i] = m.m_ptr[i] - m_ptr[i];
    }
    return m_sum;
    //доделать
};

Linalg::Matrix& Linalg::Matrix::operator-=(const Matrix& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] -= m.m_ptr[i];
    }
    return *this;
    //доделать
};
Linalg::Matrix& Linalg::Matrix::operator-=(Matrix&& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] -= m.m_ptr[i];
    }
    return *this;
    //доделать
};

bool Linalg::Matrix::operator==(const Matrix& m) const {
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) { return false;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (m_ptr[i]!=m.m_ptr[i]) { return false;}
    }
    return true;
};

bool Linalg::Matrix::operator==(Matrix&& m) const {
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) { return false;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (m_ptr[i]!=m.m_ptr[i]) { return false;}
    }
    return true;
};

bool Linalg::Matrix::operator!=(const Matrix& m) const {
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) { return true;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (m_ptr[i]!=m.m_ptr[i]) { return true;}
    }
    return false;
};

bool Linalg::Matrix::operator!=(Matrix&& m) const {
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) { return true;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (m_ptr[i]!=m.m_ptr[i]) { return true;}
    }
    return false;
};

//std::ostream& operator<<(std::ostream& os, const Linalg::Matrix& m) {
//    for (size_t i = 0; i < m.get_rows(); ++i) {
//        for (size_t j = 0; j < m.get_columns(); ++j) {
//            os << m.get_ptr()[i * m.get_columns() + j] << " ";  // Выводим элемент
//        }
//        os << '\n';  // Переход на новую строку после каждого ряда
//    }
//    return os;
//}

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

void Linalg::Matrix::print() const {
    for (size_t i=0; i<(m_columns*m_rows); ++i) {
        std::cout << m_ptr[i] << ", ";
    }
};

void Linalg::Matrix::reshape(size_t new_m_rows, size_t new_m_columns) {
    if (new_m_columns * new_m_rows == m_columns * m_rows) {
        m_columns = new_m_columns;
        m_rows = new_m_rows;
    }
}

size_t Linalg::Matrix::find_max_column_element(size_t column) const {
    if (column >= m_columns) {
        return 0;
    }

    double max_value = std::abs(m_ptr[column]);
    size_t max_row = 0;

    for (size_t row = 1; row < m_rows; ++row) {
        double element = std::abs(m_ptr[row * m_columns + column]);
        if (element > max_value) {
            max_value = element;
            max_row = row;
        }
    }

    return max_row;
}


void Linalg::Matrix::swap_rows(size_t row1, size_t row2) {
    if (row1 >= m_rows || row2 >= m_rows) {
        return;
    }
    for (size_t column = 0; column < m_columns; ++column) {
        std::swap(m_ptr[row1*m_columns + column], m_ptr[row2*m_columns + column]);
    }
}

double Linalg::Matrix::det() const {
    if (m_columns != m_rows){
        return 0.521;
    }
    Matrix m = *this;
    size_t swap_counter = 0;
    double det_value = 1;
    for (size_t i = 0; i < m_rows-1; ++i) {
        size_t max_column_element = this->find_max_column_element(i);
        if (i != max_column_element) {
            m.swap_rows(i, max_column_element);
            ++swap_counter;
        }
        if (std::fabs(m.m_ptr[i * m_columns + i]) < 1e-9) {
            return 0;
        }

        for (size_t j = i + 1; j < m_rows; ++j) {
            double multiplication = -m.m_ptr[j*m_columns + i]/m.m_ptr[i*m_columns + i];

            for (size_t k = i; k < m_rows; ++k) {
                m.m_ptr[j*m_columns + k] += m.m_ptr[i*m_columns + k]*multiplication;
            }
        }
    }
    for (size_t l=0; l<m_rows; ++l) {
        det_value *= m.m_ptr[m_columns*l + l];
    }
    det_value *= (swap_counter % 2 == 0 ? 1 : -1);
    return det_value;
}


double  Linalg::Matrix::trace() const {
    if (m_rows != m_columns){
        return 53;
    }
    double m_trace=0;
    for (size_t i=0; i<=(m_rows); ++i) {m_trace += m_ptr[m_rows*i + i];}
    return m_trace;
};
//double rank() const {
//
//};