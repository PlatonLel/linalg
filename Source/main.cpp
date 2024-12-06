#include <Linalg.h>

int main() {
    linalg::Matrix<size_t> m = {1,3,1,3};
    linalg::Matrix m2 = {{1,3},{1,3}};
    linalg::Matrix<double> m3 = {2,3,1,3};
    std::cout << (m2 != m2);
}