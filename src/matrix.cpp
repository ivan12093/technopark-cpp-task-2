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
Matrix<T>::Matrix(Matrix<T>& other) {
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
    swap(tmp);
    return *this;
}

template<class T>
void Matrix<T>::swap(Matrix<T> &other) noexcept {
    std::swap(data, other.data);
    std::swap(rows, other.rows);
    std::swap(cols, other.cols);
}

