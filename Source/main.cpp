#include <Linalg.h>

int main() {
    linalg::Matrix<size_t> m = {1,3,1,3};
    linalg::Matrix m2 = {{1,3},{1,3}};
    linalg::Matrix<double> m3 = {2,3,1,3};
    linalg::Matrix<double> m4 = {2,3,1,3};
    m4.reserve(6);
    std::cout << m4.size() << " " << m4.capacity();
    m4.clear();
    std::cout << m4.size() << " " << m4.capacity();
    m4.shrink_to_fit();
    std::cout << m4.size() << " " << m4.capacity();
}