#include "../Linalg/Include/Linalg.h"

int main() {

    Linalg::Matrix m5 = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
    Linalg::Matrix m8 = Linalg::Matrix({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}) + Linalg::Matrix({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
    std::cout << m8;
}