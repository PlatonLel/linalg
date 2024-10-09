#include "../Linalg/Include/Linalg.h"

int main() {
    Linalg::Matrix m;
    const Linalg::Matrix m2(1);
    const Linalg::Matrix m5 = {{3,4}, {5,6}};
    std::cout << m5.norm() << "\n" << m5.trace() << "\n" << m5.get_ptr() << "\n" << m5.det();
}