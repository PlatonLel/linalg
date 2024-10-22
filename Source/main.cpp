#include "../Linalg/Include/Linalg.h"

int main() {

    Linalg::Matrix m5 = {{3.0, 2.0}, {5.0, 7.0}, {11.0, 13.0}};
    Linalg::Matrix m8(0);
    m5.gauss();
    std::cout << m5;
}