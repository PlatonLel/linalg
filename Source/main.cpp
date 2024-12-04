//#include <gtest/gtest.h>
#include <Linalg.h>

int main() {
    linalg::Matrix<char> m(3);
    linalg::Matrix m2 = {1,2,2};
    m2 = linalg::Matrix<int> {1,2,3};
    std::cout << m2;
}