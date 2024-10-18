#include "Linalg.h"

double eps = std::numeric_limits<double>::epsilon();

Linalg::Matrix::Matrix(std::initializer_list<std::initializer_list<double>> m)

{
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
        throw Wrong_matrix_size{};
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
        throw Wrong_matrix_size{};
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
        throw Wrong_matrix_size{};
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
        throw Wrong_matrix_size{};
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

Linalg::Matrix Linalg::Matrix::operator*(const double& v) {
    Matrix m_return(*this);
    for (size_t i=0; i<m_columns*m_rows; ++i) {
        m_return.m_ptr[i] *= v;
    }
    return m_return;
}

Linalg::Matrix Linalg::operator*(const double& v, const Matrix& m) {
    Linalg::Matrix m_return(m);
    for (size_t i=0; i<m.get_columns()*m.get_rows(); ++i) {
        m_return.get_ptr()[i] *= v;
    }
    return m_return;
}

Linalg::Matrix Linalg::operator*(const double& v, Matrix&& m) {
    Linalg::Matrix m_return = std::move(m);
    for (size_t i=0; i<m.get_columns()*m.get_rows(); ++i) {
        m_return.get_ptr()[i] *= v;
    }
    return m_return;
}

Linalg::Matrix& Linalg::Matrix::operator*=(const double& v) {
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
    }
    for (size_t i = 0; i < (m.m_columns*m.m_rows); ++i) {
        m_ptr[i] = m.m_ptr[i];
    }
    m_rows = m.m_rows;
    m_columns = m.m_columns;
    return *this;
}

Linalg::Matrix& Linalg::Matrix::operator=(Matrix&& m) noexcept {
    if (this == &m) {
        return *this;
    }
    delete[] m_ptr;
    m_ptr = m.m_ptr;
    m_rows = m.m_rows;
    m_columns = m.m_columns;
    m.m_ptr = nullptr;
    return *this;
}

Linalg::Matrix Linalg::Matrix::operator+(const Matrix& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {throw Wrong_matrix_size{};}
    Matrix m_sum(m_rows, m_columns);
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_sum.m_ptr[i] = m.m_ptr[i] + m_ptr[i];
    }
    return m_sum;
}

Linalg::Matrix Linalg::Matrix::operator+(Matrix&& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {throw Wrong_matrix_size{};}
    Matrix m_sum(m_rows, m_columns);
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_sum.m_ptr[i] = m.m_ptr[i] + m_ptr[i];
    }
    return m_sum;
}

Linalg::Matrix& Linalg::Matrix::operator+=(const Matrix& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {throw Wrong_matrix_size{};}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] += m.m_ptr[i];
    }
    return *this;
}

Linalg::Matrix& Linalg::Matrix::operator+=(Matrix&& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {throw Wrong_matrix_size{};}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] += m.m_ptr[i];
    }
    return *this;
}

Linalg::Matrix Linalg::Matrix::operator-(const Matrix& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {throw Wrong_matrix_size{};}
    Matrix m_sum(m_rows, m_columns);
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_sum.m_ptr[i] = m.m_ptr[i] - m_ptr[i];
    }
    return m_sum;
}

Linalg::Matrix Linalg::Matrix::operator-(Matrix&& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {throw Wrong_matrix_size{};}
    Matrix m_sum(m_rows, m_columns);
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_sum.m_ptr[i] = m.m_ptr[i] - m_ptr[i];
    }
    return m_sum;
    //доделать
}

Linalg::Matrix& Linalg::Matrix::operator-=(const Matrix& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {throw Wrong_matrix_size{};}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] -= m.m_ptr[i];
    }
    return *this;
    //доделать
}

Linalg::Matrix& Linalg::Matrix::operator-=(Matrix&& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {throw Wrong_matrix_size{};}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] -= m.m_ptr[i];
    }
    return *this;
    //доделать
}

bool Linalg::Matrix::operator==(const Matrix& m) const {
    if (this == &m) {return true;}
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) {return false;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (std::abs(m_ptr[i]-m.m_ptr[i]) >= eps) {return false;}
    }
    return true;
}

bool Linalg::Matrix::operator==(Matrix&& m) const {
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) { return false;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (std::abs(m_ptr[i]-m.m_ptr[i]) >= eps) { return false;}
    }
    return true;
}

bool Linalg::Matrix::operator!=(const Matrix& m) const {
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) { return true;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (std::abs(m_ptr[i]-m.m_ptr[i]) >= eps) { return true;}
    }
    return false;
}

bool Linalg::Matrix::operator!=(Matrix&& m) const {
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) { return true;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (std::abs(m_ptr[i]-m.m_ptr[i]) >= eps) { return true;}
    }
    return false;
}

double& Linalg::Matrix::operator()(const size_t& m_row, const size_t& m_column) {
    if (m_row >= m_rows || m_column >= m_columns) {
        throw Wrong_matrix_size{};
    }
    return m_ptr[m_row * m_columns + m_column];
}

double Linalg::Matrix::operator()(const size_t& m_row, const size_t& m_column) const {
    if (m_row >= m_rows || m_column >= m_columns) {
        throw Wrong_matrix_size{};
    }
    return m_ptr[m_row * m_columns + m_column];
}

double Linalg::Matrix::norm() const {
    if (empty()) {
        throw Empty_matrix{};
    }
    double m_norm=0;
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_norm += (m_ptr[i])*(m_ptr[i]);
    }
    return std::sqrt(m_norm);
}

//void Linalg::Matrix::print() const {
//    for (size_t i=0; i<(m_columns*m_rows); ++i) {
//        std::cout << m_ptr[i] << ", ";
//    }
//}

std::ostream& Linalg::operator<<(std::ostream& os, const Matrix& m) {
    if(m.empty()) {throw Empty_matrix{};}

    for (size_t i=0; i<(m.get_rows()*m.get_columns()); ++i) {
        os << m.get_ptr()[i] << " ";
    }

    return os;
}

void Linalg::Matrix::reshape(const size_t& new_m_rows, const size_t& new_m_columns) {
    if (new_m_columns * new_m_rows != m_columns * m_rows) {
        throw Wrong_matrix_size{};
    }
    m_columns = new_m_columns;
    m_rows = new_m_rows;
}

void Linalg::Matrix::swap_rows(const size_t& row1, const size_t& row2) {
    if (row1 >= m_rows || row2 >= m_rows) {throw Wrong_matrix_size{};}
    if (row1 == row2) {return;}
    for (size_t column = 0; column < m_columns; ++column) {
        std::swap(m_ptr[row1*m_columns + column], m_ptr[row2*m_columns + column]);
    }
}

double Linalg::Matrix::det() const {
    if (empty()) {
        throw Empty_matrix{};
    }
    if (m_columns != m_rows){
        throw Wrong_matrix_size{};
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
    if (empty()) {
        throw Empty_matrix{};
    }

    if (m_rows != m_columns){
        throw Wrong_matrix_size{};
    }
    double m_trace=0;
    for (size_t i=0; i<=(m_rows); ++i) {m_trace += m_ptr[m_rows*i + i];}
    return m_trace;
}

size_t Linalg::Matrix::gauss() {
    if (empty()) {
        throw Empty_matrix{};
    }

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

Linalg::Matrix Linalg::power(const Matrix& m, int& power) {
    if (m.empty()) {
        throw Empty_matrix{};
    }
    if (m.get_columns() != m.get_rows()) {
        throw Wrong_matrix_size{};
    }

    Matrix m_return(m.get_rows(), m.get_columns());
    for (size_t i = 0; i < m.get_rows(); ++i) {
        m_return(i,i) = 1;
    }
    Matrix matrix_for_power;
    if (power>0) {matrix_for_power = m;}

    if (power==0) {return m_return;}

    if (power<0) {
        matrix_for_power = invert(m);
        power = std::abs(power);
    }

    while (power > 0) {
        if (power % 2 == 1) {
            m_return *= matrix_for_power;
        }
        matrix_for_power *= matrix_for_power;
        power /= 2;
    }

    for (size_t i=0; i<m.get_rows()*m.get_columns(); ++i) {
        if(std::abs(m_return.get_ptr()[i])<eps) {m_return.get_ptr()[i]=0;}
    }


    return m_return;
}

Linalg::Matrix Linalg::power(Matrix&& m, int& power) {
    if (m.empty()) {
        throw Empty_matrix{};
    }
    if (m.get_columns() != m.get_rows()) {
        throw Wrong_matrix_size{};
    }

    Matrix m_return(m.get_rows(), m.get_columns());
    for (size_t i = 0; i < m.get_rows(); ++i) {
        m_return(i,i) = 1;
    }
    Matrix matrix_for_power;
    if (power>0) {matrix_for_power = m;}

    if (power==0) {return m_return;}

    if (power<0) {
        matrix_for_power = invert(m);
        power = std::abs(power);
    }

    while (power > 0) {
        if (power % 2 == 1) {
            m_return *= matrix_for_power;
        }
        matrix_for_power *= matrix_for_power;
        power /= 2;
    }

    for (size_t i=0; i<m.get_rows()*m.get_columns(); ++i) {
        if(std::abs(m_return.get_ptr()[i])<eps) {m_return.get_ptr()[i]=0;}
    }


    return m_return;
}

Linalg::Matrix Linalg::power(const Matrix& m, int&& power) {
    if (m.empty()) {
        throw Empty_matrix{};
    }
    if (m.get_columns() != m.get_rows()) {
        throw Wrong_matrix_size{};
    }

    Matrix m_return(m.get_rows(), m.get_columns());
    for (size_t i = 0; i < m.get_rows(); ++i) {
        m_return(i,i) = 1;
    }
    Matrix matrix_for_power;
    if (power>0) {matrix_for_power = m;}

    if (power==0) {return m_return;}

    if (power<0) {
        matrix_for_power = invert(m);
        power = std::abs(power);
    }

    while (power > 0) {
        if (power % 2 == 1) {
            m_return *= matrix_for_power;
        }
        matrix_for_power *= matrix_for_power;
        power /= 2;
    }

    for (size_t i=0; i<m.get_rows()*m.get_columns(); ++i) {
        if(std::abs(m_return.get_ptr()[i])<eps) {m_return.get_ptr()[i]=0;}
    }

    return m_return;
}

Linalg::Matrix Linalg::power(Matrix&& m, int&& power) {
    if (m.empty()) {
        throw Empty_matrix{};
    }
    if (m.get_columns() != m.get_rows()) {
        throw Wrong_matrix_size{};
    }

    Matrix m_return(m.get_rows(), m.get_columns());
    for (size_t i = 0; i < m.get_rows(); ++i) {
        m_return(i,i) = 1;
    }
    Matrix matrix_for_power;
    if (power>0) {matrix_for_power = m;}

    if (power==0) {return m_return;}

    if (power<0) {
        matrix_for_power = invert(m);
        power = std::abs(power);
    }

    while (power > 0) {
        if (power % 2 == 1) {
            m_return *= matrix_for_power;
        }
        matrix_for_power *= matrix_for_power;
        power /= 2;
    }

    for (size_t i=0; i<m.get_rows()*m.get_columns(); ++i) {
        if(std::abs(m_return.get_ptr()[i])<eps) {m_return.get_ptr()[i]=0;}
    }

    return m_return;
}

Linalg::Matrix Linalg::concatenate(const Linalg::Matrix& m_left, const Linalg::Matrix& m_right) {
    if (m_left.empty()||m_right.empty()) {
        throw Empty_matrix{};
    }

    if (m_left.get_rows() != m_right.get_rows()) {
        throw Wrong_matrix_size{};
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
    if (m_left.empty()||m_right.empty()) {
        throw Empty_matrix{};
    }
    if (m_left.get_rows() != m_right.get_rows()) {
        throw Wrong_matrix_size{};;
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
    if (m_left.empty()||m_right.empty()) {
        throw Empty_matrix{};
    }
    if (m_left.get_rows() != m_right.get_rows()) {
        throw Wrong_matrix_size{};;
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
    if (m_left.empty()||m_right.empty()) {
        throw Empty_matrix{};
    }
    if (m_left.get_rows() != m_right.get_rows()) {
        throw Wrong_matrix_size{};;
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
    if (m.empty()) {
        throw Empty_matrix{};
    }
    Matrix m_return(m.get_columns(), m.get_rows());
    for (size_t i = 0; i < m.get_rows(); ++i) {
        for (size_t j =0; j < m.get_columns(); ++j) {
            m_return(i,j) = m(j,i);
        }
    }
    return m_return;
}

Linalg::Matrix Linalg::transpose(Matrix&& m) {
    if (m.empty()) {
        throw Empty_matrix{};
    }
    Matrix m_return(m.get_columns(), m.get_rows());
    for (size_t i = 0; i < m.get_rows(); ++i) {
        for (size_t j =0; j < m.get_columns(); ++j) {
            m_return(i,j) = m(j,i);
        }
    }
    return m_return;
}

Linalg::Matrix Linalg::invert(const Matrix& m) {
    if (m.get_rows() != m.get_columns()) { throw Wrong_matrix_size{}; }

    size_t n = m.get_rows();
    Matrix m_L(n, n);
    Matrix m_U(n, n);

    for (size_t i=0; i<n; ++i) {
        m_L.get_ptr()[i] = 0;
        m_U.get_ptr()[i] = 0;
    }

    // Инициализация m_L как единичной матрицы
    for (size_t i = 0; i < n; ++i) {
        m_L(i,i) = 1;
    }

    // Инициализация первой строки и первого столбца
    for (size_t j = 0; j < n; j++) {
        m_U(0,j) = m(0,j);
        if (m_U(0,0) == 0) {throw Singular_matrix{};}
        m_L(j,0) = m(j,0) / m_U(0,0);
    }

    // Заполнение матриц m_L и m_U
    for (size_t i = 1; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            m_U(i,j) = m(i,j) - get_sum_U(i,j,m_L,m_U);
            if (m_U(i,i) == 0) {throw Singular_matrix{};}
            m_L(j,i) = (m(j,i) - get_sum_L(i, j, m_L, m_U)) / m_U(i,i);
        }
    }

    Matrix m_return(n, n);
    Matrix e(n);
    for (size_t i=0; i<n; ++i) {e(i,0) = 0;}
    // Обратный ход
    for (size_t j = 0; j < n; ++j) {
        e(j,0) = 1;

        Matrix y = forward_substitution(m_L, e);
        Matrix x = backward_substitution(m_U, y);

        for (size_t i = 0; i < n; ++i) {
            m_return(i, j) = x(i,0);
        }
        e(j,0) = 0;
    }

    return m_return;
}


double Linalg::get_sum_U(size_t& m_row, size_t& m_column, Matrix& m_L, Matrix& m_U) {
    double result = 0;

    for (size_t i=0; i<m_column; ++i) {
        if (std::abs(m_L(m_row,i)) <= eps) {m_L(m_row,i) = 0;}
        if (std::abs(m_U(i,m_column)) <= eps) {m_U(i,m_column) = 0;}
        result += m_L(m_row,i)*m_U(i,m_column);
    }
    return result;
}

double Linalg::get_sum_L(size_t& m_row, size_t& m_column, Matrix& m_L, Matrix& m_U) {
    double result = 0;

    for (size_t i=0; i<m_column; ++i) {
        if (std::abs(m_L(m_column,i)) <= eps) {m_L(m_column,i) = 0;}
        if (std::abs(m_U(i,m_row)) <= eps) {m_U(i,m_row) = 0;}
        result += m_L(m_column,i)*m_U(i,m_row);
    }
    return result;
}
Linalg::Matrix Linalg::backward_substitution(const Matrix& m_U, const Matrix& y) {
    size_t n = m_U.get_rows();
    Matrix x(n);

    for (int i = n - 1; i >= 0; --i) {
        x(i,0) = y(i,0);
        for (size_t j = i + 1; j < n; ++j) {
            x(i,0) -= m_U(i,j) * x(j,0);
        }
        x(i,0) /= m_U(i,i);
    }

    return x;
}
Linalg::Matrix Linalg::forward_substitution(const Matrix& m_L, const Matrix& b) {
    size_t n = m_L.get_rows();
    Matrix y(n);

    for (size_t i = 0; i < n; ++i) {
        y(i,0) = b(i,0);
        for (size_t j = 0; j < i; ++j) {
            y(i,0) -= m_L(i,j) * y(j,0);
        }
        y(i,0) /= m_L(i,i);
    }

    return y;
}


