#ifndef TECHNOPARK_CPP_TASK_2_MATRIX_HPP
#define TECHNOPARK_CPP_TASK_2_MATRIX_HPP

#include <utility>
using std::size_t;

#include "vector.hpp"

template<class T>
class Matrix {
private:
    Vector<Vector<T>> data;
    size_t rows = 0;
    size_t cols = 0;
public:
    Matrix() = default;
    Matrix(size_t _rows, size_t _cols);
    Matrix(T** _data, size_t _rows, size_t _cols);
    explicit Matrix(Vector<Vector<T>> vec_of_vec);
    Matrix(std::initializer_list<Vector<T>> initializerList);
    Matrix(std::initializer_list<std::initializer_list<T>> initializerList);

    ~Matrix();
    Matrix(const Matrix& other);
    Matrix(Matrix&& other) noexcept;
    Matrix& operator=(const Matrix& other);
    Matrix& operator=(Matrix&& other) noexcept;

    Matrix transposed();
    //Matrix inversed();
    T determinant() const;

    Matrix operator-(const Matrix& other);
    Matrix operator+(const Matrix& other);
    Matrix operator*(const Matrix& other);

    template<class U>
    friend Matrix<U> operator*(const Matrix<U>& a, const U& b);
    Matrix operator*(const T& other);
    Matrix& operator*=(const T& other);

    Vector<T> getDiag() const;
    Vector<T> getColumn(size_t idx) const;
    Vector<T> operator[](size_t idx) const;
    Vector<T>& operator[](size_t idx);

    std::pair<size_t, size_t> size() const;

    void swap(Matrix& other) noexcept;
};

#include <cassert>
#include <algorithm>

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
    swap(*this, other);
}

template<class T>
Matrix<T>::Matrix(Matrix<T> &&other) noexcept {
    if (this == &other)
        return;
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
T Matrix<T>::determinant() const {
    assert(rows == cols);
    T n1;
    T n2;
    T det(1);
    T total(1);
    size_t index;
    size_t _size = size().first;

    Vector<T> temp(_size + 1);
    for (size_t i = 0; i < _size; ++i) {
        index = i;
        while (index < _size && data[index][i] == T(0))
            ++index;
        if (index == _size)
            continue;
        if (index != i) {
            for (size_t j = 0; j < _size; ++j)
                std::swap(data[index][j], data[i][j]);
            det = det * ((index - i) % 2 == 0 ? T(1) : T(-1));
        }
        for (size_t j = 0; j < _size; ++j)
            temp[j] = data[i][j];
        for (size_t j = i + 1; j < _size; ++j) {
            n1 = temp[i];
            n2 = data[j][i];
            for (size_t k = 0; k < _size; ++k)
                data[j][k] = (n1 * data[j][k]) - (n2 * temp[k]);
            total *= n1;
        }
    }
    for (size_t i = 0; i < _size; ++i)
        det *= data[i][i];
    return det / total;
}

template<class T>
Vector<T> Matrix<T>::getDiag() const {
    size_t elems_in_diag = std::min(rows, cols);
    Vector<T> result(elems_in_diag);
    for (size_t i = 0; i < elems_in_diag; ++i)
        result[i] = data[i][i];
}

template<class T>
Matrix<T>::Matrix(Vector<Vector<T>> vec_of_vec) {
    if (vec_of_vec[0].get_format() == Row) {
        data = vec_of_vec;
        rows = vec_of_vec.size();
        cols = vec_of_vec[0].size();
        return;
    }
    rows = vec_of_vec[0].size();
    cols = vec_of_vec.size();
    data = Vector<Vector<T>>(rows);
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; i < rows; ++j)
            data[i][j] = vec_of_vec[j][i];
}

template<class T>
Matrix<T>::Matrix(std::initializer_list<Vector<T>> initializerList) {
    assert(initializerList.size() != 0);
    rows = initializerList.size();
    cols = initializerList.begin()->size();
    data = Vector<Vector<T>>(rows);
    for (auto it_src = initializerList.begin(),
                 it_dst = data.begin(); it_src < initializerList.end(); ++it_src, ++it_dst)
        *it_dst = *it_src;
}

template<class T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> initializerList) {
    assert(initializerList.size() != 0);
    rows = initializerList.size();
    cols = initializerList.begin()->size();
    assert(cols != 0);
    data = Vector<Vector<T>>(rows);
    auto it_dst = data.begin();
    for (auto it_src = initializerList.begin(); it_src < initializerList.end(); ++it_src, ++it_dst) {
        assert(it_src->size() == cols);
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

#endif //TECHNOPARK_CPP_TASK_2_MATRIX_HPP
