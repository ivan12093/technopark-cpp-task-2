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
    T _determinant(const Matrix& mtx) const;
    Matrix downsized(size_t row, size_t col) const;
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

    Matrix& operator-=(const Matrix& other);
    Matrix& operator+=(const Matrix& other);
    Matrix& operator*=(const Matrix& other);

    Matrix operator-(const T& other);
    Matrix operator+(const T& other);
    Matrix operator*(const T& other);

    template<class U>
    friend Matrix<U> operator*(const U& a, const Matrix<U>& b);
    template<class U>
    friend Matrix<U> operator-(const U& a, const Matrix<U>& b);
    template<class U>
    friend Matrix<U> operator+(const U& a, const Matrix<U>& b);

    Matrix& operator*=(const T& other);
    Matrix& operator+=(const T& other);
    Matrix& operator-=(const T& other);

    Matrix operator*(const Vector<T>& other);
    template<class U>
    friend Matrix<U> operator*(const Vector<U>& a, const Matrix<U>& b);

    Matrix dot(const Matrix& other);

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
    rows = other.rows;
    cols = other.cols;
    data = other.data;
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
    swap(tmp);
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
    swap(other);
    return *this;
}

template<class T>
Matrix<T>::Matrix(Matrix<T> &&other) noexcept {
    if (this == &other)
        return;
    swap(other);
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
    Matrix result(cols, rows);
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
        result[i] = data[i] - other.data[i];
    return result;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix &other) {
    assert(rows == other.rows);
    assert(cols == other.cols);
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i)
        result[i] = data[i] + other.data[i];
    return result;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix &other) {
    assert(rows == other.rows);
    assert(cols == other.cols);
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i)
        result[i] = data[i] * other.data[i];
    return result;
}

template<class T>
Matrix<T>& Matrix<T>::operator*=(const T& other) {
    for (size_t i = 0; i < rows; ++i)
        data[i] = data[i] * other;
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
T Matrix<T>::_determinant(const Matrix& mtx) const {
    assert(mtx.rows == mtx.cols);
    if (mtx.rows == 1)
        return mtx[0][0];
    if (mtx.rows == 2)
        return mtx[0][0] * mtx[1][1] - mtx[0][1] * mtx[1][0];
    T result(0);
    for (size_t i = 0; i < mtx.rows; ++i) {
        T mul = i % 2 == 0 ? 1 : -1;
        result += mul * mtx.data[i][0] * _determinant(mtx.downsized(i, 0));
    }
    return result;
}

template<class T>
T Matrix<T>::determinant() const {
    return _determinant(*this);
}

template<class T>
Vector<T> Matrix<T>::getDiag() const {
    size_t elems_in_diag = std::min(rows, cols);
    Vector<T> result(elems_in_diag);
    for (size_t i = 0; i < elems_in_diag; ++i)
        result[i] = data[i][i];
    return result;
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
    auto it_dst = data.begin();
    for (auto it_src = initializerList.begin(); it_src < initializerList.end(); ++it_src, ++it_dst)
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

template<class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix &other) {
    for (size_t i = 0; i < rows; ++i)
        data[i] -= other;
    return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix &other) {
    for (size_t i = 0; i < rows; ++i)
        data[i] += other;
    return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator*=(const Matrix &other) {
    for (size_t i = 0; i < rows; ++i)
        data[i] *= other;
    return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const T &other) {
    Matrix result(*this);
    for (size_t i = 0; i < rows; ++i)
        result.data[i] -= other;
    return result;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const T &other) {
    Matrix result(*this);
    for (size_t i = 0; i < rows; ++i)
        result.data[i] += other;
    return result;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const T &other) {
    Matrix result(*this);
    for (size_t i = 0; i < rows; ++i)
        result.data[i] *= other;
    return result;
}

template<class U>
Matrix<U> operator*(const U &a, const Matrix<U> &b) {
    Matrix result(b);
    for (size_t i = 0; i < b.rows; ++i)
        result.data[i] *= a;
    return result;
}

template<class U>
Matrix<U> operator-(const U &a, const Matrix<U> &b) {
    Matrix result(b);
    for (size_t i = 0; i < b.rows; ++i)
        result.data[i] *= -1, result.data[i] += a;
    return result;
}

template<class U>
Matrix<U> operator+(const U &a, const Matrix<U> &b) {
    Matrix result(b);
    for (size_t i = 0; i < b.rows; ++i)
        result.data[i] += a;
    return result;
}

template<class T>
Matrix<T>& Matrix<T>::operator+=(const T &other) {
    for (size_t i = 0; i < rows; ++i)
        data[i] += other;
    return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator-=(const T &other) {
    for (size_t i = 0; i < rows; ++i)
        data[i] -= other;
    return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Vector<T> &other) {
    if (other.get_format() == Row) {
        assert(cols == 1);
        Matrix result(rows, other.size());
        for (size_t i = 0; i < result.rows; ++i) {
            for (size_t j = 0; j < result.cols; ++j) {
                result[i][j] += data[i][0] * other[j];
            }
        }
        return result;
    }
    assert(cols == other.size());
    Matrix result(rows, 1);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < other.size(); ++j)
            result.data[i][0] += data[i][j] * other[j];
    }
    return result;
}

template<class U>
Matrix<U> operator*(const Vector<U> &a, const Matrix<U> &b) {
    if (a.get_format() == Row) {
        assert(a.size() == b.rows);
        Matrix<U> result(1, b.cols);
        for (size_t i = 0; i < result.cols; ++i)
            for (size_t j = 0; j < b.rows; ++j)
                result[0][i] += a[j] * b[j][i];
        return result;
    }
    assert(b.rows == 1);
    Matrix<U> result(a.size(), b.cols);
    for (size_t i = 0; i < result.rows; ++i)
        for (size_t j = 0; j < result.cols; ++j)
            result[i][j] += a[i] * b[0][j];
    return result;
}

template<class T>
Matrix<T> Matrix<T>::dot(const Matrix &other) {
    assert(cols == other.rows);
    Matrix result(rows, other.cols);
    for (size_t i = 0; i < result.rows; ++i)
        for (size_t j = 0; j < result.cols; ++j) {
            for (size_t k = 0; k < cols; ++k)
                result[i][j] += data[i][k] * other[k][j];
        }
    return result;
}

template<class T>
Matrix<T> Matrix<T>::downsized(size_t row, size_t col) const {
    Matrix res(*this);
    assert(res.rows != 0);
    for (size_t i = row; i < rows - 1; ++i)
        res.data[i] = res.data[i + 1];
    --res.rows;
    assert(res.rows != 0);
    for (size_t j = col; j < cols - 1; ++j)
        for (size_t i = 0; i < res.rows; ++i)
            res.data[i][j] = res.data[i][j + 1];
    --res.cols;
    return res;
}

#endif //TECHNOPARK_CPP_TASK_2_MATRIX_HPP
