#ifndef TECHNOPARK_CPP_TASK_2_MATRIX_H
#define TECHNOPARK_CPP_TASK_2_MATRIX_H

#include <cstddef>

template<class T>
class Matrix {
private:
    T** data = nullptr;
    size_t rows = 0;
    size_t cols = 0;
    void swap(Matrix& other) noexcept;
public:
    Matrix() = delete;
    Matrix(size_t rows, size_t cols);
    ~Matrix();
    Matrix& operator=(const Matrix& other);
    Matrix(Matrix& other);
    Matrix& operator=(Matrix&& other);
};


#endif //TECHNOPARK_CPP_TASK_2_MATRIX_H
