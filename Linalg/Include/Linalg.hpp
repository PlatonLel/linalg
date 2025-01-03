#include <algorithm>
#include <iomanip>

double eps = std::numeric_limits<double>::epsilon();

template <typename T>
template <typename Y>
void linalg::Matrix<T>::copy_constructor(const Matrix<Y>& m)  {
    m_ptr = reinterpret_cast<T*>(operator new(m.size() * sizeof(T)));
    T* ptr = m_ptr;
    for (const Y* ptr2 = m.begin(); ptr2 != m.end(); ++ptr2, ++ptr) {
        new(ptr) T( *ptr2 );
    }
    m_size = m.size();
    m_columns = m.columns();
    m_rows = m.rows();
    m_capacity = m.capacity();
}

template <typename T>
template <typename Y>
linalg::Matrix<T>::Matrix(std::initializer_list<std::initializer_list<Y>> m) {
    if (m.size() == 0) {
        return;
    }

    size_t cols = m.begin()->size();
    for (const auto& row : m) {
        if (row.size() != cols) {
            throw Wrong_matrix_size(0);
        }
    }

    m_rows = m.size();
    m_columns = cols;
    m_size = m_rows * m_columns;
    m_capacity = m_size;
    m_ptr = reinterpret_cast<T*>(operator new(m_size * sizeof(T)));

    T* ptr = m_ptr;
    for (const auto& row : m) {
        for (const auto& el : row) {
            new (ptr) T(static_cast<T>(el));
            ++ptr;
        }
    }
}


template <typename T>
template <typename Y>
linalg::Matrix<T>::Matrix(std::initializer_list<Y> m) {

    if (m.size() == 0) {
        return;
    }
    m_ptr = reinterpret_cast<T*>(operator new(m.size() * sizeof(T)));
    T* ptr = m_ptr;
    for (const T &el: m) {
        new(ptr) T(static_cast<T>(el));
        ++ptr;
    }
    m_columns = m.size();
    m_rows = 1;
    m_size = m_columns*m_rows;
    m_capacity = m_size;
}

template <typename T>
linalg::Matrix<T>::Matrix(const size_t& rows) {
    if (rows == 0) {
        return;
    }
    m_ptr = reinterpret_cast<T*>(operator new(rows * sizeof(T)));
    for (T* ptr = m_ptr; ptr < m_ptr + rows;++ptr) {
        new(ptr) T{};
    }
    m_rows = rows;
    m_columns = 1;
    m_size = rows;
    m_capacity = rows;
}

template <typename T>
linalg::Matrix<T>::Matrix(const size_t& rows, const size_t& columns) {
    if (rows == 0 || columns == 0) {
        return;
    }
    m_ptr = reinterpret_cast<T*>(operator new(rows * columns * sizeof(T)));
    for (T* ptr = m_ptr; ptr < m_ptr + rows;++ptr) {
        new(ptr) T{};
    }
    m_rows = rows;
    m_columns = columns;
    m_size = rows*columns;
    m_capacity = m_size;
}

template <typename T>
template <typename Y>
linalg::Matrix<T>& linalg::Matrix<T>::operator+=(const Matrix<Y>& m){
    if (m_columns!=m.columns() || m_rows!=m.rows()) {
        throw Wrong_matrix_size(3);
    }
    for (size_t i=0;i<m_size;++i){
        m_ptr[i]+=static_cast<T>(m[i]);
    };
    return *this;
}

template <typename T>
linalg::Matrix<T> linalg::operator*(const Matrix<T>& m1, const Matrix<T>& m2) {
    if (m1.columns() != m2.rows()) {
        throw Wrong_matrix_size(1);
    }
    Matrix m_return(m1);
    m_return*=m2;

    return m_return;
}

template <typename T>
template <typename Y>
linalg::Matrix<T>& linalg::Matrix<T>::operator*=(const Matrix<Y>& m) {
    if (m_columns != m.rows()) {
        throw Wrong_matrix_size(2);
    }
    Matrix m_return(m_rows, m.columns());
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m.columns(); ++j) {
            m_return[i * m.columns() + j] = 0;
            for (size_t k = 0; k < m_columns; ++k) {
                m_return[i * m.columns() + j] += static_cast<T>(m_ptr[i * m_columns + k]) * static_cast<T>(m[k * m.columns() + j]);
            }
        }
    }
    for (size_t i=0; i<m_return.columns()*m_return.rows(); ++i) {
        if(std::abs(m_return[i])<eps*1000) {
            m_return[i]=0;
        }
    }
    *this = std::move(m_return);
    return *this;
}
template <typename T>
linalg::Matrix<T> linalg::operator*(const Matrix<T>& m, const double& v) {
    Matrix m_return(m);
    m_return*=v;
    return m_return;
}
template <typename T>
linalg::Matrix<T> linalg::operator*(const double& v, const Matrix<T>& m) {
    Matrix m_return(m);
    m_return*=v;
    return m_return;
}

template <typename T>
template <typename Y>
linalg::Matrix<T>& linalg::Matrix<T>::operator*=(const Y& v) noexcept {
    for (size_t i=0; i<m_columns*m_rows; ++i) {
        m_ptr[i] *= v;
    }
    for (size_t i=0; i<m_columns*m_rows; ++i) {
        if(std::abs(m_ptr[i])<eps*1000) {
            m_ptr[i]=0;
        }
    }
    return *this;
}

template <typename T>
void linalg::Matrix<T>::swap(Matrix<T>& m) noexcept {
    std::swap(m_size, m.m_size);
    std::swap(m_ptr, m.m_ptr);
    std::swap(m_capacity, m.m_capacity);
    std::swap(m_columns, m.m_columns);
    std::swap(m_rows, m.m_rows);
}

template <typename T>
linalg::Matrix<T>& linalg::Matrix<T>::operator=(Matrix<T>&& m) noexcept {
    swap(m);
    return *this;
}

template <typename T>
template <typename Y>
linalg::Matrix<T>& linalg::Matrix<T>::operator=(const Matrix<Y>& m) {
    return *this = Matrix<T>(m);
//    if (m_capacity < m.size()) {
//        return *this = Matrix{m};
//    }
//
//    size_t i = 0;
//    for (; i < std::min(m_size, m.size()); ++i) {
//        m_ptr[i] = static_cast<T>(m[i]);
//    }
//
//    if (m_size < m.size()) {
//        for (; i < m.size(); ++i) {
//            new (m_ptr + i) T(static_cast<T>(m[i]));
//        }
//    } else {
//        for (; i < m_size; ++i) {
//            m_ptr[i].~T();
//        }
//    }
//
//    // Обновляем метаданные
//    m_rows = m.rows();
//    m_columns = m.columns();
//    m_size = m.size();
//
//    return *this;
}

template <typename T>
linalg::Matrix<T>& linalg::Matrix<T>::operator=(const Matrix<T>& m) {
    if (*this == m) return *this;
    if (m_capacity < m.size()) {
        return *this = Matrix{m};
    }

    size_t i = 0;
    for (; i < std::min(m_size, m.size()); ++i) {
        m_ptr[i] = m[i];
    }

    if (m_size < m.size()) {
        for (; i < m.size(); ++i) {
            new (m_ptr + i) T(m[i]);
        }
    } else {
        for (; i < m_size; ++i) {
            m_ptr[i].~T();
        }
    }
    m_rows = m.rows();
    m_columns = m.columns();
    m_size = m.size();
    return *this;
}


template <typename T>
template <typename Y>
linalg::Matrix<T> linalg::Matrix<T>::operator+(const Matrix<Y>& m) const {
    if (m_rows != m.rows()|| m_columns != m.columns()) {
        throw Wrong_matrix_size(3);
    }
    Matrix m_return(*this);
    m_return+=m;
    return m_return;
}


template <typename T>
template <typename Y>
linalg::Matrix<T> linalg::Matrix<T>::operator-(const Matrix<Y>& m) const {
    if (m_rows != m.rows()|| m_columns != m.columns()) {
        throw Wrong_matrix_size(5);
    }
    Matrix m_return(*this);
    m_return-=m;
    return m_return;
}

template <typename T>
template <typename Y>
linalg::Matrix<T>& linalg::Matrix<T>::operator-=(const Matrix<Y>& m) {
    if (m_rows != m.rows()|| m_columns != m.columns()) {
        throw Wrong_matrix_size(6);
    }
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        m_ptr[i] -= static_cast<T>(m[i]);
    }
    return *this;
}

template <typename T>
template <typename Y>
bool linalg::Matrix<T>::operator==(const Matrix<Y>& m) const {
    if constexpr (std::is_same_v<T, Y>) {
        if (this == &m) {
            return true;
        }
    }
    if (m_columns != m.columns() || m_rows != m.rows()) {
        return false;
    }
    for (size_t i = 0; i < m_rows * m_columns; ++i) {
        if (m_ptr[i] != static_cast<T>(m[i])) {
            return false;
        }
    }
    return true;
}

template <typename T>
template <typename Y>
bool linalg::Matrix<T>::operator!=(const Matrix<Y>& m) const {
    return !(*this == m);
}


template <typename T>
bool linalg::Matrix<T>::operator==(const Matrix<double>& m) const {
    if (m_columns!=m.columns()||m_rows!=m.rows()) {return false;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (std::abs(m_ptr[i] - m[i]) <= eps*100) {return true;}
    }
    return true;
}

template <typename T>
bool linalg::Matrix<T>::operator!=(const Matrix<double>& m) const {
    if (m_columns!=m.columns()||m_rows!=m.rows()) {return true;}
    for (size_t i=0; i<m_rows*m_columns; ++i) {
        if (std::abs(m_ptr[i] - m[i]) <= eps*100) {return false;}
    }
    return true;
}

template <typename T>
T& linalg::Matrix<T>::operator()(const size_t& m_row, const size_t& m_column) {
    if (m_row >= m_rows || m_column >= m_columns) {
        throw Wrong_matrix_size(7);
    }
    return m_ptr[m_row * m_columns + m_column];
}

template <typename T>
T linalg::Matrix<T>::operator()(const size_t& m_row, const size_t& m_column) const {
    if (m_row >= m_rows || m_column >= m_columns) {
        throw Wrong_matrix_size(7);
    }
    return m_ptr[m_row * m_columns + m_column];
}
//
template <typename T>
std::ostream& linalg::operator<<(std::ostream& os, const Matrix<T>& m) {
    if (m.empty()) {
        throw Empty_matrix(6);
    }

    size_t max_width_first_column = 0;
    size_t max_width = 0;
    std::ostringstream temp;

    for (size_t i = 0; i < m.rows(); ++i) {
        for (size_t j = 0; j < m.columns(); ++j) {
            temp.str("");
            temp.clear();
            temp << m(i, j);

            if (j == 0) {
                max_width_first_column = std::max(max_width_first_column, temp.str().length());
            }
            max_width = std::max(max_width, temp.str().length());
        }
    }

    std::ostringstream sout;
    for (size_t i = 0; i < m.rows(); ++i) {
        sout << "|";
        for (size_t j = 0; j < m.columns(); ++j) {
            temp.str("");
            temp.clear();
            temp << m(i, j);

            if (j == 0) {
                sout << std::setw(max_width_first_column) << temp.str() << " ";
            } else {
                sout << std::setw(max_width) << temp.str() << (j == m.columns() - 1 ? "" : " ");
            }
        }
        sout << "|\n";
    }

    os << sout.str();
    return os;
}

template <typename T>
T linalg::Matrix<T>::norm() const {
    if (empty()) {
        throw Empty_matrix(0);
    }
    T norm=0;
    for (size_t i=0; i<(m_rows*m_columns); ++i) {
        norm += (m_ptr[i])*(m_ptr[i]);
    }
    return std::sqrt(norm);
}

template <typename T>
void linalg::Matrix<T>::reshape(const size_t& new_rows, const size_t& new_columns) {
    if(this->empty()) {
        return;
    }
    Matrix<T> m_temp = *this;
    *this = Matrix<T>(new_rows,new_columns);
    if (m_capacity<new_columns*new_rows) {
        *this = m_temp;
        m_columns = new_columns;
        m_rows = new_rows;
        m_size = m_rows*m_columns;
        return;
    }
    for (size_t i = 0; i < this->m_size; ++i) {
        m_ptr[i] = m_temp[i];
    }
}

template <typename T>
void linalg::Matrix<T>::swap_rows(const size_t& row_1, const size_t& row_2) noexcept {
    if (row_1 == row_2) {return;}
    //по одному элементу меняем местами
    for (size_t column = 0; column < m_columns; ++column) {
        std::swap(m_ptr[row_1*m_columns + column], m_ptr[row_2*m_columns + column]);
    }
}

template <typename T>
T linalg::Matrix<T>::det() const {
    if (empty()) {
        throw Empty_matrix(1);
    }
    if (m_columns != m_rows){
        throw Wrong_matrix_size(9);
    }
    Matrix m = *this;
    //вызываем метод гаусса, который поменяет матрицу и вернет кол-во свапнутых строк
    size_t swap_counter = m.gauss();
    T det_value = 1;

    for (size_t l=0; l<m_rows; ++l) {
        det_value *= m.m_ptr[m_columns*l + l];
    }
    det_value *= (swap_counter % 2 == 0 ? 1 : -1);//если количество свапов четно, то 1, иначе -1
    return det_value;
}

template <typename T>
T  linalg::Matrix<T>::trace() const {
    if (empty()) {
        throw Empty_matrix(2);
    }

    if (m_rows != m_columns){
        throw Wrong_matrix_size(10);
    }
    T m_trace=0;
    for (size_t i=0; i<m_rows; ++i) {
        m_trace += m_ptr[m_rows*i + i];
    }
    return m_trace;
}

//только прямой ход
template <typename T>
size_t linalg::Matrix<T>::gauss() noexcept{
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
            T coefficient = m_ptr[l * m_columns + lead_element] / m_ptr[i * m_columns + lead_element];
            for (size_t s = 0; s < m_columns; ++s) {
                m_ptr[l * m_columns + s] -= coefficient * m_ptr[i * m_columns + s];
            }
        }
        ++lead_element;
    }
    return swap_counter;
}

template <typename T>
linalg::Matrix<T> linalg::power(const Matrix<T>& m, const int& p) {
    int power = p;
    if (m.empty()) {
        throw Empty_matrix(3);
    }
    if (m.columns() != m.rows()) {
        throw Wrong_matrix_size(11);
    }

    Matrix<T> m_return(m.rows(), m.columns());

    for (size_t i = 0; i < m.rows(); ++i) {
        for (size_t j = 0; j < m.columns(); ++j) {
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

    Matrix matrix_for_power(m);

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

    for (size_t i = 0; i < m.rows() * m.columns(); ++i) {
        if (std::abs(m_return.get_ptr()[i]) < eps) {
            m_return.get_ptr()[i] = 0;
        }
    }

    return m_return;
}

template <typename T>
linalg::Matrix<T> linalg::concatenate(const Matrix<T>& m_left, const Matrix<T>& m_right) {
    if (m_left.empty()||m_right.empty()) {
        throw Empty_matrix(4);
    }

    if (m_left.rows() != m_right.rows()) {
        throw Wrong_matrix_size(12);
    }

    Matrix<T> m_return(m_left.rows(), m_left.columns() + m_right.columns());

    for (size_t i = 0; i < m_left.rows(); ++i) {
        for (size_t j = 0; j < m_left.columns(); ++j) {
            m_return(i, j) = m_left(i, j);
        }
        for (size_t j = 0; j < m_right.columns(); ++j) {
            m_return(i, m_left.columns() + j) = m_right(i, j);
        }
    }

    return m_return;
}

template <typename T>
linalg::Matrix<T> linalg::transpose(const Matrix<T>& m) {
    if (m.empty()) {
        throw Empty_matrix(5);
    }
    Matrix<T> m_return(m.columns(), m.rows());
    for (size_t i = 0; i < m.rows(); ++i) {
        for (size_t j =0; j < m.columns(); ++j) {
            m_return(j,i) = m(i,j);
        }
    }
    return m_return;
}

template <typename T>
linalg::Matrix<T> linalg::invert(const Matrix<T>& m) {
    if (m.rows() != m.columns()) {
        throw Wrong_matrix_size(13);
    }

    size_t n = m.rows();
    Matrix<T> m_L(n, n);//тк реализация с помощью разложения на L и U матрицы
    Matrix<T> m_U(n, n);

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

    Matrix<T> m_return(n, n);
    Matrix<T> e(n);
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
    for (size_t i=0; i< m_return.columns() * m_return.rows(); ++i) {
        if(std::abs(m_return.get_ptr()[i])<eps*10000) {
            m_return.get_ptr()[i]=0;
        }
    }

    return m_return;
}


template <typename T>
T linalg::get_sum_U(size_t& m_row, size_t& m_column, Matrix<T>& m_L, Matrix<T>& m_U) {
    T result = 0;

    for (size_t i=0; i<m_column; ++i) {
        if (std::abs(m_L(m_row,i)) <= eps*1000) {
            m_L(m_row,i) = 0;
        }

        if (std::abs(m_U(i,m_column)) <= eps*1000) {
            m_U(i,m_column) = 0;
        }

        result += m_L(m_row,i)*m_U(i,m_column);
    }
    return result;
}
//
template <typename T>
T linalg::get_sum_L(size_t& m_row, size_t& m_column, Matrix<T>& m_L, Matrix<T>& m_U) {
    T result = 0;

    for (size_t i=0; i<m_column; ++i) {
        if (std::abs(m_L(m_column,i)) <= eps*1000) {
            m_L(m_column,i) = 0;
        }

        if (std::abs(m_U(i,m_row)) <= eps*1000) {
            m_U(i,m_row) = 0;
        }

        result += m_L(m_column,i)*m_U(i,m_row);
    }
    return result;
}
//
template <typename T>
linalg::Matrix<T> linalg::backward_substitution(const Matrix<T>& m_U, const Matrix<T>& y) {
    size_t n = m_U.rows();
    Matrix<T> x(n);

    for (int i = n - 1; i >= 0; --i) {
        x(i,0) = y(i,0);
        for (size_t j = i + 1; j < n; ++j) {
            x(i,0) -= m_U(i,j) * x(j,0);
        }
        x(i,0) /= m_U(i,i);
    }

    return x;
}
//
template <typename T>
linalg::Matrix<T> linalg::forward_substitution(const Matrix<T>& m_L, const Matrix<T>& b) {
    size_t n = m_L.rows();
    Matrix<T> y(n);

    for (size_t i = 0; i < n; ++i) {
        y(i,0) = b(i,0);
        for (size_t j = 0; j < i; ++j) {
            y(i,0) -= m_L(i,j) * y(j,0);
        }
        y(i,0) /= m_L(i,i);
    }

    return y;
}
//
linalg::Wrong_matrix_size::Wrong_matrix_size(size_t p) {
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

linalg::Empty_matrix::Empty_matrix(size_t p) {
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

linalg::Singular_matrix::Singular_matrix(size_t p) {
    switch (p) {
        case 0: {
            description = "Singular matrix in invert";
            break;
        }
    }
}

template <typename T>
linalg::Matrix<T>::~Matrix() noexcept {
    for (T* ptr = m_ptr; ptr < m_ptr + m_size; ++ptr) {
        ptr->~T();
    }
    operator delete(reinterpret_cast<void*>(m_ptr));
}

template <typename T>
void linalg::Matrix<T>::shrink_to_fit(){
    if (m_size == m_capacity) {
        return;
    }
    *this = Matrix{*this};
}

template <typename T>
void linalg::Matrix<T>::clear() noexcept {
    for (T* ptr = m_ptr; ptr < m_ptr + m_size; ++ptr) {
        ptr->~T();
    }
    m_size = 0;
    m_capacity = 0;
}

template <typename T>
void linalg::Matrix<T>::reserve(size_t n) {
    this->clear();
    m_ptr = reinterpret_cast<T*>(operator new(n * sizeof(T)));
    m_capacity = n;
}

Complex linalg::parse_complex(const std::string& complex_str) {
    size_t i_pos = complex_str.find('i');
    if (i_pos == std::string::npos) {
        throw std::invalid_argument("Invalid complex number format: " + complex_str);
    }

    size_t plus_pos = complex_str.find('+');
    size_t minus_pos = complex_str.find('-', 1);

    double real = 0.0;
    double img = 0.0;

    try {
        if (plus_pos != std::string::npos) {
            real = std::stod(complex_str.substr(0, plus_pos));
            img = std::stod(complex_str.substr(plus_pos + 1, i_pos - plus_pos - 1));
        } else if (minus_pos != std::string::npos) {
            real = std::stod(complex_str.substr(0, minus_pos));
            img = std::stod(complex_str.substr(minus_pos, i_pos - minus_pos));
        } else {
            throw std::invalid_argument("Invalid complex number format: " + complex_str);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error parsing complex number: " << complex_str << "\n";
        throw;
    }

    return {real, img};
}

linalg::Matrix<Complex> linalg::load_matrix(const char* file_name) {
    std::ifstream file(file_name);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file");
    }

    std::string line;
    size_t rows = 0;
    size_t cols = 0;

    while (std::getline(file, line)) {
        line.erase(std::remove(line.begin(), line.end(), '|'), line.end());
        std::istringstream row_stream(line);
        size_t col_count = 0;
        std::string complex_str;

        while (row_stream >> complex_str) {
            col_count++;
        }

        if (cols == 0) {
            cols = col_count;
        } else if (cols != col_count) {
            throw std::runtime_error("Row size mismatch");
        }

        rows++;
    }

    file.clear();
    file.seekg(0);

    linalg::Matrix<Complex> matrix(rows, cols);

    size_t row = 0;
    while (std::getline(file, line)) {
        line.erase(std::remove(line.begin(), line.end(), '|'), line.end());
        std::istringstream row_stream(line);
        std::string complex_str;
        size_t col = 0;

        while (row_stream >> complex_str) {
            Complex complex_num = parse_complex(complex_str);
            matrix(row, col) = complex_num;
            col++;
        }
        row++;
    }

    return matrix;
}
