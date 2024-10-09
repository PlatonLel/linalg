#include "../Linalg/Include/Linalg.h"

int main() {
    Linalg::Matrix m;
    const Linalg::Matrix m2(1);
    Linalg::Matrix m3(3,4);
    Linalg::Matrix m4(m3);
    Linalg::Matrix m5 = {{3,4}, {5,6}};
    std::cout << m5.norm() << "\n" << m5.get_columns() << "\n" << m3.get_ptr() << "\n" << m5.det();
}