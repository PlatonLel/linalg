#include "../Linalg/Include/Linalg.h"

int main() {
    Linalg::Matrix m;
    const Linalg::Matrix m2(1);
    Linalg::Matrix m5 = {{1,2}, {1,2}};
    Linalg::Matrix m6 = {{12,634,51,64},{76,82,71,75}};

    Linalg::Matrix m7 = Linalg::concatenate({{1,2},{3,4}},m5);
//    m7.det();
//    Linalg::Matrix m7(m5*m6);
    m7.print();
//    std::cout << m5.det();
}