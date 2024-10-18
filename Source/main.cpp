#include "../Linalg/Include/Linalg.h"

int main() {
//    Linalg::Matrix m;
//    const Linalg::Matrix m2(1);
    Linalg::Matrix m5 = {{1,0,0}, {0,1,0},{0,0,1}};
//    Linalg::Matrix m6 = Linalg::Matrix({1,2}) * Linalg::Matrix({{3,4}, {5,6}});
    Linalg::Matrix m7 = power(m5,0);
    std::cout << m7;
}