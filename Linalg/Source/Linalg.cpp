#include "Linalg.h"

double eps = std::numeric_limits<double>::epsilon();

Linalg::Matrix::Matrix(std::initializer_list<std::initializer_list<double>> m)
{
    //проверка на правильный ввод
    for (const auto& row : m) {
        if (row.size()!=m.begin()->size()) {
            throw Wrong_matrix_size(0);
        }
    }
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
//перемещающий конструктор
Linalg::Matrix::Matrix(Matrix&& m) noexcept {
    m_ptr = m.m_ptr;
    m.m_ptr = nullptr;
    m_columns = m.m_columns;
    m_rows = m.m_rows;
}
//копирующий
Linalg::Matrix::Matrix(const Matrix& m) {
    m_ptr = new double[m.m_columns*m.m_rows];
    for (size_t i = 0; i<(m.m_columns*m.m_rows); ++i) {
        m_ptr[i] = m.m_ptr[i];
    }
    m_columns = m.m_columns;
    m_rows = m.m_rows;
}

Linalg::Matrix Linalg::operator*(const Matrix& m1, const Matrix& m2) {
    if (m1.get_columns() != m2.get_rows()) {
        throw Wrong_matrix_size(1);
    }

    Matrix m_return(m1);
    m_return*=m2;

    return m_return;
}

Linalg::Matrix& Linalg::Matrix::operator*=(const Matrix& m) {
    if (m_columns != m.m_rows) {
        throw Wrong_matrix_size(2);
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

Linalg::Matrix Linalg::operator*(const Matrix& m, const double& v) {
    Matrix m_return(m);
    m_return*=v;
    return m_return;
}

Linalg::Matrix Linalg::operator*(const double& v, const Matrix& m) {
    Matrix m_return(m);
    m_return*=v;
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

Linalg::Matrix Linalg::Matrix::operator+(const Matrix& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {throw Wrong_matrix_size(3);}
    Matrix m_sum(*this);
    m_sum+=m;
    return m_sum;
}

Linalg::Matrix& Linalg::Matrix::operator+=(const Matrix& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {throw Wrong_matrix_size(4);}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] += m.m_ptr[i];
    }
    return *this;
}

Linalg::Matrix Linalg::Matrix::operator-(const Matrix& m) const {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {throw Wrong_matrix_size(5);}
    Matrix m_sum(*this);
    m_sum-=m;
    return m_sum;
}

Linalg::Matrix& Linalg::Matrix::operator-=(const Matrix& m) {
    if (m_rows != m.m_rows|| m_columns != m.m_columns) {throw Wrong_matrix_size(6);}
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] -= m.m_ptr[i];
    }
    return *this;
}

bool Linalg::Matrix::operator==(const Matrix& m) const {
    if (this == &m) {return true;}
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) {return false;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (std::abs(m_ptr[i]-m.m_ptr[i]) >= eps*100) {return false;}
    }
    return true;
}

bool Linalg::Matrix::operator!=(const Matrix& m) const {
    if (m_columns!=m.m_columns||m_rows!=m.m_rows) { return true;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (std::abs(m_ptr[i]-m.m_ptr[i]) >= eps*100) { return true;}
    }
    return false;
}

double& Linalg::Matrix::operator()(const size_t& m_row, const size_t& m_column) {
    if (m_row >= m_rows || m_column >= m_columns) {
        throw Wrong_matrix_size(7);
    }
    return m_ptr[m_row * m_columns + m_column];
}

double Linalg::Matrix::operator()(const size_t& m_row, const size_t& m_column) const {
    if (m_row >= m_rows || m_column >= m_columns) {
        throw Wrong_matrix_size(7);
    }
    return m_ptr[m_row * m_columns + m_column];
}

std::ostream& Linalg::operator<<(std::ostream& os, const Matrix& m) {
    if (m.empty()) {
        throw Empty_matrix(6);
    }
    size_t max_width_first_column = 0;
    size_t max_width = 0;
//проходимся по элементам и вычисляем максимальную длину в 1 столбце и в остальных
    for (size_t i = 0; i < m.get_rows(); ++i) {
        std::ostringstream temp_1;
        temp_1 << m(i, 0);
        max_width_first_column = std::max(max_width_first_column, temp_1.str().length());
        for (size_t j = 0; j < m.get_columns(); ++j) {
            std::ostringstream temp_2;
            temp_2 << m(i, j);
            max_width = std::max(max_width, temp_2.str().length());
        }
    }

    for (size_t i = 0; i < m.get_rows(); ++i) {
        os << "|";
        for (size_t j = 0; j < m.get_columns(); ++j) {
            if (j==0) {
                os << std::setw(max_width_first_column) << m(i, j) << " ";
            }
            else {
                os << std::setw(max_width) << m(i, j) << " ";
            }
        }
        os << "|\n";
    }
    return os;
}

double Linalg::Matrix::norm() const {
    if (empty()) {
        throw Empty_matrix(0);
    }
    double norm=0;
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        norm += (m_ptr[i])*(m_ptr[i]);
    }
    return std::sqrt(norm);
}

void Linalg::Matrix::reshape(const size_t& new_m_rows, const size_t& new_m_columns) {
    if (new_m_columns * new_m_rows != m_columns * m_rows) {
        throw Wrong_matrix_size(8);
    }
    m_columns = new_m_columns;
    m_rows = new_m_rows;
}

void Linalg::Matrix::swap_rows(const size_t& row1, const size_t& row2) noexcept {
    if (row1 == row2) {return;}
    //по одному элементу меняем местами
    for (size_t column = 0; column < m_columns; ++column) {
        std::swap(m_ptr[row1*m_columns + column], m_ptr[row2*m_columns + column]);
    }
}

double Linalg::Matrix::det() const {
    if (empty()) {
        throw Empty_matrix(1);
    }
    if (m_columns != m_rows){
        throw Wrong_matrix_size(9);
    }
    Matrix m = *this;
    //вызываем метод гаусса, который поменяет матрицу и вернет кол-во свапнутых строк
    size_t swap_counter = m.gauss();
    double det_value = 1;

    for (size_t l=0; l<m_rows; ++l) {
        det_value *= m.m_ptr[m_columns*l + l];
    }
    det_value *= (swap_counter % 2 == 0 ? 1 : -1);//если количество свапов четно, то 1, иначе -1
    return det_value;
}

double  Linalg::Matrix::trace() const {
    if (empty()) {
        throw Empty_matrix(2);
    }

    if (m_rows != m_columns){
        throw Wrong_matrix_size(10);
    }
    double m_trace=0;
    for (size_t i=0; i<m_rows; ++i) {
        m_trace += m_ptr[m_rows*i + i];
    }
    return m_trace;
}
//только прямой ход
size_t Linalg::Matrix::gauss() noexcept{
    size_t swap_counter = 0;
    size_t lead_element = 0;

    for (size_t i = 0; i < m_rows; ++i) {
        //если вышли за пределы кол ва столбцов прерываем
        if (lead_element >= m_columns) {
            break;
        }

        size_t j = i;

        while (m_ptr[j * m_columns + lead_element] == 0) {  // если элемент 0
            ++j;  // следующая строка
            if (j == m_rows) {  // если дошли до конца столбца
                j = i;          // возвращаемся к исходной строке
                ++lead_element; // переходим к следующему столбцу
                if (lead_element == m_columns) {  // если столбцы закончились
                    return swap_counter;
                }
            }
        }

        if (j != i) {
            this->swap_rows(j, i);
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

Linalg::Matrix Linalg::power(const Matrix& m, const int& p) {
    int power = p;
    if (m.empty()) {
        throw Empty_matrix(3);
    }
    if (m.get_columns() != m.get_rows()) {
        throw Wrong_matrix_size(11);
    }

    Matrix m_return(m.get_rows(), m.get_columns());

    for (size_t i = 0; i < m.get_rows(); ++i) {
        for (size_t j = 0; j < m.get_columns(); ++j) {
            if (i == j) {
                m_return(i,i) = 1;
            } else {
                m_return(i,j) = 0;
            }
        }
    }

    if (power == 0) {
        return m_return;
    }

    Matrix matrix_for_power;
    if (power > 0) {
        matrix_for_power = m;
    }

    if (power < 0) {
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

    for (size_t i = 0; i < m.get_rows() * m.get_columns(); ++i) {
        if (std::abs(m_return.get_ptr()[i]) < eps) {
            m_return.get_ptr()[i] = 0;
        }
    }

    return m_return;
}

Linalg::Matrix Linalg::concatenate(const Matrix& m_left, const Matrix& m_right) {
    if (m_left.empty()||m_right.empty()) {
        throw Empty_matrix(4);
    }

    if (m_left.get_rows() != m_right.get_rows()) {
        throw Wrong_matrix_size(12);
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
        throw Empty_matrix(5);
    }
    Matrix m_return(m.get_columns(), m.get_rows());
    for (size_t i = 0; i < m.get_rows(); ++i) {
        for (size_t j =0; j < m.get_columns(); ++j) {
            m_return(j,i) = m(i,j);
        }
    }
    return m_return;
}


Linalg::Matrix Linalg::invert(const Matrix& m) {
    if (m.get_rows() != m.get_columns()) {
        throw Wrong_matrix_size(13);
    }

    size_t n = m.get_rows();
    Matrix m_L(n, n);//тк реализация с помощью разложения на L и U матрицы
    Matrix m_U(n, n);

    for (size_t i=0; i<n*n; ++i) {
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
        if (std::abs(m_U(0,0)) < eps) {
            throw Singular_matrix(0);
        }
        m_L(j,0) = m(j,0) / m_U(0,0);
    }

    // Заполнение матриц m_L и m_U
    for (size_t i = 1; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            m_U(i,j) = m(i,j) - get_sum_U(i,j,m_L,m_U); //заполняем матрицу по формуле
            if (std::abs(m_U(i,i)) < eps) {
                throw Singular_matrix(0);
            }
            //заполняем матрицу по формуле
            m_L(j,i) = (m(j,i) - get_sum_L(i, j, m_L, m_U)) / m_U(i,i);
        }
    }

    Matrix m_return(n, n);
    Matrix e(n);
    for (size_t i=0; i<n; ++i) {
        e(i,0) = 0;
    }
    for (size_t j = 0; j < n; ++j) {
        e(j,0) = 1;
        //подставляем сначала вектор с единицей в j строке
        Matrix y = forward_substitution(m_L, e);
        //затем подставляем полученный после прямой подстановки вектор и получаем в итоге столбец матрицы
        Matrix x = backward_substitution(m_U, y);
        //тут заполняем
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
        if (std::abs(m_L(m_row,i)) <= eps) {
            m_L(m_row,i) = 0;
        }

        if (std::abs(m_U(i,m_column)) <= eps) {
            m_U(i,m_column) = 0;
        }

        result += m_L(m_row,i)*m_U(i,m_column);
    }
    return result;
}

double Linalg::get_sum_L(size_t& m_row, size_t& m_column, Matrix& m_L, Matrix& m_U) {
    double result = 0;

    for (size_t i=0; i<m_column; ++i) {
        if (std::abs(m_L(m_column,i)) <= eps) {
            m_L(m_column,i) = 0;
        }

        if (std::abs(m_U(i,m_row)) <= eps) {
            m_U(i,m_row) = 0;
        }

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

Linalg::Wrong_matrix_size::Wrong_matrix_size(size_t p) {
    switch (p) {
        case 0: {
            description = "Wrong initializer_list in constructor of Matrix";
            break;
        }
        case 1: {
            description = "Wrong matrix sizes for multiplication";
            break;
        }
        case 2: {
            description = "Wrong matrix sizes for multiplication and assignment";
            break;
        }
        case 3: {
            description = "Wrong matrix sizes for addition";
            break;
        }
        case 4: {
            description = "Wrong matrix sizes for addition and assignment";
            break;
        }
        case 5: {
            description = "Wrong matrix sizes for subtraction";
            break;
        }
        case 6: {
            description = "Wrong matrix sizes for subtraction and assignment";
            break;
        }
        case 7: {
            description = "Wrong matrix size index access";
            break;
        }
        case 8: {
            description = "Wrong matrix sizes for reshape";
            break;
        }
        case 9: {
            description = "Not square matrix for determinant";
            break;
        }
        case 10: {
            description = "Not square matrix for trace";
            break;
        }
        case 11: {
            description = "Not square matrix for power";
            break;
        }
        case 12: {
            description = "Wrong number of rows in right matrix for concatenate";
            break;
        }
        case 13: {
            description = "Not square matrix for invert";
            break;
        }
    }
}

Linalg::Empty_matrix::Empty_matrix(size_t p) {
    switch (p) {
        case 0: {
            description = "Empty matrix in norm";
            break;
        }
        case 1: {
            description = "Empty matrix in det";
            break;
        }
        case 2: {
            description = "Empty matrix in trace";
            break;
        }
        case 3: {
            description = "Empty matrix in power";
            break;
        }
        case 4: {
            description = "Empty matrix in concatenate";
            break;
        }
        case 5: {
            description = "Empty matrix in transpose";
            break;
        }
        case 6: {
            description = "Empty matrix in cout";
            break;
        }
    }
}

Linalg::Singular_matrix::Singular_matrix(size_t p) {
    switch (p) {
        case 0: {
            description = "Singular matrix in invert";
            break;
        }
    }
}

