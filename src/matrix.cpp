#include <cassert>
#include <algorithm>

#include "matrix.h"

template<class T>
Matrix<T>::Matrix(size_t _rows, size_t _cols): rows(_rows), cols(_cols) {
    assert(rows > 0);
    assert(cols > 0);
    data = new T*[rows];
    for (size_t i = 0 ; i < rows; ++i)
        data[i] = new T[cols]();

}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& other) {
    Matrix(other.rows, other.cols);
    for (size_t i = 0; i < rows; ++i)
        std::copy(other.data[i], other.data[i] + cols, data[i]);
}

template<class T>
Matrix<T>::~Matrix() {
    for (size_t i = 0; i < rows; ++i)
        delete[] data[i], data[i] = nullptr;
    delete data;
    data = nullptr;
    rows = 0;
    cols = 0;
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) {
    if (this == &other)
        return *this;
    Matrix<T> tmp(other);
    swap(*this, tmp);
    return *this;
}

template<class T>
void swap(Matrix<T>& a, Matrix<T> &b) noexcept {
    std::swap(a.data, b.data);
    std::swap(a.rows, b.rows);
    std::swap(a.cols, b.cols);
}

template<class T>
Matrix<T>& Matrix<T>::operator=(Matrix&& other) noexcept {
    if (this == &other)
        return *this;
    ~Matrix();
    swap(*this, other);
}

template<class T>
Matrix<T>::Matrix(Matrix<T> &&other) noexcept {
    if (this == &other)
        return;
    ~Matrix();
    swap(*this, other);
}

template<class T>
Matrix<T>::Matrix(T **_data, size_t rows, size_t cols) {
    Matrix(rows, cols);
    for (size_t i = 0; i < rows; ++i)
        std::copy(_data[i], _data[i] + cols, data[i]);
}

template<class T>
Matrix<T> Matrix<T>::transposed() {
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            result.data[j][i] = data[i][j];
    return result;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix &other) {
    assert(rows == other.rows);
    assert(cols == other.cols);
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            result[i][j] = data[i][j] - other.data[i][j];
    return result;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix &other) {
    assert(rows == other.rows);
    assert(cols == other.cols);
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            result[i][j] = data[i][j] + other.data[i][j];
    return result;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix &other) {
    assert(rows == other.rows);
    assert(cols == other.cols);
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            result[i][j] = data[i][j] * other.data[i][j];
    return result;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const T &other) {
    Matrix result(*this);
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            result[i][j] = data[i][j] * other;
    return result;
}

template<class T>
Matrix<T>& Matrix<T>::operator*=(const T& other) {
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; i < cols; ++j)
            data[i][j] = data[i][j] * other;
    return *this;
}

template<class T>
Matrix<T> operator*(const Matrix<T> &a, const T &b) {
    Matrix result(a.rows, a.cols);
    for (size_t i = 0; i < a.rows; ++i)
        for (size_t j = 0; j < a.cols; ++j)
            result.data[i][j] = a.data[i][j] * b;
    return result;
}
