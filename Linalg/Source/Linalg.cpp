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
}

Linalg::Matrix Linalg::operator*(double v, Linalg::Matrix& m) {
    Linalg::Matrix m_return(m);
    for (size_t i=0; i<m.get_columns()*m.get_rows(); ++i) {
        m_return.get_ptr()[i] *= v;
    }
    return m_return;
}

Linalg::Matrix Linalg::operator*(double v, Linalg::Matrix&& m) {
    Linalg::Matrix m_return = std::move(m);
    for (size_t i=0; i<m.get_columns()*m.get_rows(); ++i) {
        m_return.get_ptr()[i] *= v;
    }
    return m_return;
}

Linalg::Matrix& Linalg::Matrix::operator*=(double v) {
    for (size_t i=0; i<m_columns*m_rows; ++i) {
        m_ptr[i] *= v;
    }
    return *this;
}

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

Linalg::Matrix& Linalg::Matrix::operator=(Matrix&& m) noexcept {
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
}

Linalg::Matrix Linalg::Matrix::operator+(Matrix&& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    Matrix m_sum(m_rows, m_columns);
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_sum.m_ptr[i] = m.m_ptr[i] + m_ptr[i];
    }
    return m_sum;
    //доделать
}

Linalg::Matrix& Linalg::Matrix::operator+=(const Matrix& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] += m.m_ptr[i];
    }
    return *this;
    //доделать
}

Linalg::Matrix& Linalg::Matrix::operator+=(Matrix&& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] += m.m_ptr[i];
    }
    return *this;
    //доделать
}

Linalg::Matrix Linalg::Matrix::operator-(const Matrix& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    Matrix m_sum(m_rows, m_columns);
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_sum.m_ptr[i] = m.m_ptr[i] - m_ptr[i];
    }
    return m_sum;
    //доделать
}

Linalg::Matrix Linalg::Matrix::operator-(Matrix&& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    Matrix m_sum(m_rows, m_columns);
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_sum.m_ptr[i] = m.m_ptr[i] - m_ptr[i];
    }
    return m_sum;
    //доделать
}

Linalg::Matrix& Linalg::Matrix::operator-=(const Matrix& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] -= m.m_ptr[i];
    }
    return *this;
    //доделать
}

Linalg::Matrix& Linalg::Matrix::operator-=(Matrix&& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {return *this;}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] -= m.m_ptr[i];
    }
    return *this;
    //доделать
}

bool Linalg::Matrix::operator==(const Matrix& m) const {
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) { return false;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (m_ptr[i]!=m.m_ptr[i]) { return false;}
    }
    return true;
}

bool Linalg::Matrix::operator==(Matrix&& m) const {
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) { return false;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (m_ptr[i]!=m.m_ptr[i]) { return false;}
    }
    return true;
}

bool Linalg::Matrix::operator!=(const Matrix& m) const {
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) { return true;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (m_ptr[i]!=m.m_ptr[i]) { return true;}
    }
    return false;
}

bool Linalg::Matrix::operator!=(Matrix&& m) const {
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) { return true;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (m_ptr[i]!=m.m_ptr[i]) { return true;}
    }
    return false;
}

double& Linalg::Matrix::operator()(const size_t& m_row, const size_t& m_column) {
    if (m_row >= m_rows || m_column >= m_columns) {
        std::cerr << "Invalid row or column!";
        exit(EXIT_FAILURE);
    }
    return m_ptr[m_row * m_columns + m_column];
}

double Linalg::Matrix::operator()(const size_t& m_row, const size_t& m_column) const {
    if (m_row >= m_rows || m_column >= m_columns) {
        std::cerr << "Invalid row or column!";
        exit(EXIT_FAILURE);
    }
    return m_ptr[m_row * m_columns + m_column];
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

void Linalg::Matrix::print() const {
    for (size_t i=0; i<(m_columns*m_rows); ++i) {
        std::cout << m_ptr[i] << ", ";
    }
}

void Linalg::Matrix::reshape(size_t& new_m_rows, size_t& new_m_columns) {
    if (new_m_columns * new_m_rows == m_columns * m_rows) {
        m_columns = new_m_columns;
        m_rows = new_m_rows;
    }
}

void Linalg::Matrix::swap_rows(size_t& row1, size_t& row2) {
    if (row1 >= m_rows || row2 >= m_rows) {return;}
    if (row1 == row2) {return;}
    for (size_t column = 0; column < m_columns; ++column) {
        std::swap(m_ptr[row1*m_columns + column], m_ptr[row2*m_columns + column]);
    }
}

double Linalg::Matrix::det() const {
    if (m_columns != m_rows){
        return 0.52;
    }
    Matrix m = *this;
    size_t swap_counter = m.gauss();
    double det_value = 1;

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
}

size_t Linalg::Matrix::gauss() {
    size_t swap_counter = 0;
    size_t lead_element = 0;

    for (size_t i = 0; i < m_rows; ++i) {
        if (lead_element >= m_columns) {
            break;
        }

        size_t j = i;

        while (m_ptr[j * m_columns + lead_element] == 0) {
            ++j;
            if (j == m_rows) {
                j = i;
                ++lead_element;
                if (lead_element == m_columns) {
                    return swap_counter;
                }
            }
        }

        this->swap_rows(j, i);
        if (j != i) {
            ++swap_counter;
        }

        for (size_t l = i + 1; l < m_rows; ++l) {
            double coefficient = m_ptr[l * m_columns + lead_element] / m_ptr[i * m_columns + lead_element];
            for (size_t s = 0; s < m_columns; ++s) {
                m_ptr[l * m_columns + s] -= coefficient * m_ptr[i * m_columns + s];
            }
        }
        ++lead_element;
    }
    return swap_counter;
}

Linalg::Matrix Linalg::power(const Matrix& m, size_t power) {
    if (m.get_columns() != m.get_rows()) {
        return m;
    }

    Matrix m_return(m.get_rows(), m.get_columns());
    for (size_t i = 0; i < m.get_rows(); ++i) {
        m_return.get_ptr()[i*m_return.get_columns() + i] = 1;
    }

    Matrix matrix_for_power = m;

    while (power > 0) {
        if (power % 2 == 1) {
            m_return *= matrix_for_power;
        }
        matrix_for_power *= matrix_for_power;
        power /= 2;
    }

    return m_return;
}

Linalg::Matrix Linalg::power(Matrix&& m, size_t power) {
    if (m.get_columns() != m.get_rows()) {
        return m;
    }

    Matrix m_return(m.get_rows(), m.get_columns());
    for (size_t i = 0; i < m.get_rows(); ++i) {
        m_return.get_ptr()[i * m_return.get_columns() + i] = 1;
    }

    Matrix matrix_for_power = m;

    while (power > 0) {
        if (power % 2 == 1) {
            m_return *= matrix_for_power;
        }
        matrix_for_power *= matrix_for_power;
        power /= 2;
    }

    return m_return;
}

Linalg::Matrix Linalg::concatenate(const Linalg::Matrix &m_left, const Linalg::Matrix &m_right) {
    if (m_left.get_rows() != m_right.get_rows()) {
        return {};
    }

    Matrix m_return(m_left.get_rows(), m_left.get_columns() + m_right.get_columns());

    for (size_t i = 0; i < m_left.get_rows(); ++i) {
        for (size_t j = 0; j < m_left.get_columns(); ++j) {
            m_return(i, j) = m_left(i, j);
        }
        for (size_t j = 0; j < m_right.get_columns(); ++j) {
            m_return(i, m_left.get_columns() + j) = m_right(i, j);
        }
    }

    return m_return;
}

Linalg::Matrix Linalg::concatenate(Linalg::Matrix &&m_left, const Linalg::Matrix &m_right) {
    if (m_left.get_rows() != m_right.get_rows()) {
        return {};
    }

    Matrix m_return(m_left.get_rows(), m_left.get_columns() + m_right.get_columns());

    for (size_t i = 0; i < m_left.get_rows(); ++i) {
        for (size_t j = 0; j < m_left.get_columns(); ++j) {
            m_return(i, j) = m_left(i, j);
        }
        for (size_t j = 0; j < m_right.get_columns(); ++j) {
            m_return(i, m_left.get_columns() + j) = m_right(i, j);
        }
    }

    return m_return;
}

Linalg::Matrix Linalg::concatenate(Linalg::Matrix &&m_left, Linalg::Matrix &&m_right) {
    if (m_left.get_rows() != m_right.get_rows()) {
        return {};
    }

    Matrix m_return(m_left.get_rows(), m_left.get_columns() + m_right.get_columns());

    for (size_t i = 0; i < m_left.get_rows(); ++i) {
        for (size_t j = 0; j < m_left.get_columns(); ++j) {
            m_return(i, j) = m_left(i, j);
        }
        for (size_t j = 0; j < m_right.get_columns(); ++j) {
            m_return(i, m_left.get_columns() + j) = m_right(i, j);
        }
    }

    return m_return;
}

Linalg::Matrix Linalg::concatenate(const Linalg::Matrix &m_left, Linalg::Matrix &&m_right) {
    if (m_left.get_rows() != m_right.get_rows()) {
        return {};
    }

    Matrix m_return(m_left.get_rows(), m_left.get_columns() + m_right.get_columns());

    for (size_t i = 0; i < m_left.get_rows(); ++i) {
        for (size_t j = 0; j < m_left.get_columns(); ++j) {
            m_return(i, j) = m_left(i, j);
        }
        for (size_t j = 0; j < m_right.get_columns(); ++j) {
            m_return(i, m_left.get_columns() + j) = m_right(i, j);
        }
    }

    return m_return;
}

Linalg::Matrix Linalg::transpose(const Matrix& m) {
    Matrix m_return(m.get_columns(), m.get_rows());
    for (size_t i = 0; i < m.get_rows(); ++i) {
        for (size_t j =0; j < m.get_columns(); ++j) {
            m_return(i,j) = m(j,i);
        }
    }
    return m_return;
}

Linalg::Matrix Linalg::transpose(Matrix&& m) {
    Matrix m_return(m.get_columns(), m.get_rows());
    for (size_t i = 0; i < m.get_rows(); ++i) {
        for (size_t j =0; j < m.get_columns(); ++j) {
            m_return(i,j) = m(j,i);
        }
    }
    return m_return;
}

