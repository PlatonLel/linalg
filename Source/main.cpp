#include "../Linalg/Include/Linalg.h"

int main() {
    Linalg::Matrix m;
    const Linalg::Matrix m2(1);
    const Linalg::Matrix m5 = {{1,3}, {3,5}};
    Linalg::Matrix m6 = {{12,634,51,64},{76,82,71,75}};
//    m5 *= 2.6 ;
//    Linalg::Matrix m7(m5*m6);
//    m7.det();
//    Linalg::Matrix m7(m5*m6);
//    m7.print();
    std::cout << m5.det();
}