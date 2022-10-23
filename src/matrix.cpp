#include <cassert>
#include <algorithm>

#include "matrix.hpp"

template<class T>
Matrix<T>::Matrix(size_t _rows, size_t _cols): rows(_rows), cols(_cols) {
    assert(rows != 0);
    assert(cols != 0);
    data = Vector<Vector<T>>(rows);
    for (size_t i = 0 ; i < rows; ++i)
        data[i] = Vector<T>(cols);
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& other) {
    Matrix(other.rows, other.cols);
    for (size_t i = 0; i < rows; ++i)
        std::copy(other.data[i], other.data[i] + cols, data[i]);
}

template<class T>
Matrix<T>::~Matrix() {
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
void Matrix<T>::swap(Matrix<T> &other) noexcept {
    std::swap(data, other.data);
    std::swap(rows, other.rows);
    std::swap(cols, other.cols);
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
Matrix<T>::Matrix(T **_data, size_t _rows, size_t _cols) {
    assert(_rows != 0);
    assert(_rows != 0);
    data = Vector<Vector<T>>(rows);
    rows = _rows;
    cols = _cols;
    for (size_t i = 0; i < rows; ++i)
        data[i] = Vector<T>(_data[i], cols);
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

template<class T>
T Matrix<T>::determinant() {

}

template<class T>
Vector<T> Matrix<T>::getDiag() const {
    size_t elems_in_diag = std::min(rows, cols);
    Vector<T> result(elems_in_diag);
    for (size_t i = 0; i < elems_in_diag; ++i)
        result[i] = data[i][i];
}

template<class T>
Matrix<T>::Matrix(Vector<Vector<T>> vec_of_vec): Matrix(vec_of_vec.size(), vec_of_vec[0].size()) {
    data = vec_of_vec;
}

template<class T>
Matrix<T>::Matrix(Vector<T> *vectors, size_t size) {
    assert(size != 0);
    rows = size;
    cols = vectors[0].size();
    data = Vector<Vector<T>>(size);
    for (size_t i = 0; i < rows; ++i)
        data[i] = vectors[i];
}

template<class T>
Matrix<T>::Matrix(std::initializer_list<Vector<T>> initializerList) {
    static_assert(initializerList.size() != 0);
    rows = initializerList.size();
    cols = initializerList.begin()->size();
    data = Vector<Vector<T>>(rows);
    for (auto it_src = initializerList.begin(),
              it_dst = data.begin(); it_src < initializerList.end(); ++it_src, ++it_dst)
        *it_dst = *it_src;
}

template<class T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> initializerList) {
    static_assert(initializerList.size() != 0);
    rows = initializerList.size();
    cols = initializerList.begin()->size();
    static_assert(cols != 0);
    data = Vector<Vector<T>>(rows);
    for (auto it_src = initializerList.begin(),
              it_dst = data.begin(); it_src < initializerList.end(); ++it_src, ++it_dst) {
        static_assert(it_src->size() == cols);
        *it_dst = *it_src;
    }
}

template<class T>
Vector<T> Matrix<T>::getColumn(size_t idx) const {
    Vector<T> result(rows, Column);
    for (size_t i = 0; i < rows; ++i)
        result[i] = data[i][idx];
    return result;
}

template<class T>
Vector<T> Matrix<T>::operator[](size_t idx) const {
    return data[idx];
}

template<class T>
Vector<T>& Matrix<T>::operator[](size_t idx) {
    return data[idx];
}

template<class T>
std::pair<size_t, size_t> Matrix<T>::size() const {
    return std::make_pair(rows, cols);
}
