#include "../Linalg/Include/Linalg.h"

int main() {
    Linalg::Matrix m;
    const Linalg::Matrix m2(1);
    Linalg::Matrix m5 = {{3,4}, {5,6}};
    Linalg::Matrix m6 = {{2,3}, {7,6}};
    m6 += m5;
    std::cout << m5.norm() << "\n" << m5.trace() << "\n" << m5.get_ptr() << "\n" << m6.det();
}