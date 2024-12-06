//#include <gtest/gtest.h>
#include <Linalg.h>

int main() {
    linalg::Matrix m = {1.3,2.6,1.3,2.9};
    linalg::Matrix m2 = {{1,3},{1,3}};
    linalg::Matrix<int> m3 = {{1,3},{1,3}};
    std::cout << (m3*=m2);
}