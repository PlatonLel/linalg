#include "../Linalg/Include/Linalg.h"

int main() {
    Linalg::Matrix m;
    const Linalg::Matrix m2(1);
    Linalg::Matrix m5 = {{1,2}, {3,4}};
    Linalg::Matrix m6 = {{5,6},{7,8}};
//    m5 *= 2.6 ;
//    m5.print();
//    Linalg::Matrix m7(m5*m6);
//    m7.print();
    std::cout << m5(2,2);
}