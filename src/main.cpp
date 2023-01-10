#include <iostream>
#include "matrix.hpp"

int main(void) {
    Matrix<double> matrix = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 1.0, 9.0}
    };
    auto det = matrix.determinant();
    std::cout << "Determinant: " << det << "\n";
    auto matrix1 = matrix * matrix;
    auto [rows, cols] = matrix1.size();
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            std::cout << matrix1[i][j] << " ";
        }
        std::cout << "\n";
    }
}
