//#include <gtest/gtest.h>
#include <Linalg.h>

int main() {
    linalg::Matrix<char> m(3);
    linalg::Matrix m2 = {1,2,2};
    std::cout << (m+=m2);
}