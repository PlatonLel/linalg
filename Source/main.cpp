//#include <gtest/gtest.h>
#include <Linalg.h>

int main() {
    linalg::Matrix m = {{1,2},{1,2}};
    linalg::Matrix m2 = {{1,3},{1,3}};
    m2 *= 4;
//    m2*=m;
    std::cout << m2;
}