//#include <gtest/gtest.h>
#include <Linalg.h>

int main() {
    linalg::Matrix m = {1.3,2.6,1.3,2.9};
    linalg::Matrix<int> m2 = {{1,3},{1,3}};
    m2 = m;
    std::cout << m2;
}