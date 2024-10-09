#include "../Linalg/Include/Linalg.h"

int main() {
    Linalg::Matrix m;
    const Linalg::Matrix m2(1);
    Linalg::Matrix m5 = {{3,4}, {5,6}};
    m5[3] = 5;
    std::cout << m5.norm() << "\n" << m5.get_columns() << "\n" << m5.get_ptr() << "\n" << m5.det();
}