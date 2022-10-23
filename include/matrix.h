#ifndef TECHNOPARK_CPP_TASK_2_MATRIX_H
#define TECHNOPARK_CPP_TASK_2_MATRIX_H

#include <cstddef>

template<class T>
class Matrix {
private:
    T** data = nullptr;
    size_t rows = 0;
    size_t cols = 0;
    friend void swap(Matrix& a, Matrix& b) noexcept;
public:
    Matrix() = delete;
    Matrix(size_t rows, size_t cols);
    Matrix(T** data, size_t rows, size_t cols);
    ~Matrix();
    Matrix(const Matrix& other);
    Matrix(Matrix&& other) noexcept;
    Matrix& operator=(const Matrix& other);
    Matrix& operator=(Matrix&& other) noexcept;

    Matrix transposed();
    Matrix inversed();
    T determinant();

    Matrix operator-(const Matrix& other);
    Matrix operator+(const Matrix& other);
    Matrix operator*(const Matrix& other);

    friend Matrix operator*(const Matrix& a, const T& b);
    Matrix operator*(const T& other);
    Matrix& operator*=(const T& other);
};



#endif //TECHNOPARK_CPP_TASK_2_MATRIX_H
