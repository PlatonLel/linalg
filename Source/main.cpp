#include <Linalg.h>

int main() {
    linalg::Matrix<size_t> m = {1,3,1,3};
    linalg::Matrix<int> m2(8,8);
    linalg::Matrix<double> m3 = {2,3,1,3};
    linalg::Matrix<double> m4 = {2,3,1,3,5,6,4,5,7};
    m4.reshape(2,2);
    std::cout << m4;
    std::cout << m4.capacity();
}