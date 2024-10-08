#include "../Linalg/Include/Linalg.h"

int main() {
    Linalg::Matrix m;
    const Linalg::Matrix m2(1);
    Linalg::Matrix m3(3,4);
    Linalg::Matrix m3(3,4);
    std::cout << m3.get_rows() << "\n" << m3.get_columns() << "\n" << m3.get_ptr() << "\n" << m2.empty();
}