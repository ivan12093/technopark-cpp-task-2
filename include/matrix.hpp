#ifndef TECHNOPARK_CPP_TASK_2_MATRIX_HPP
#define TECHNOPARK_CPP_TASK_2_MATRIX_HPP

#include <utility>
#include <cassert>
#include <algorithm>
#include <cfloat>

#include "vector.hpp"

// T must support the multiplication, addition, subtraction operator
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

    Matrix transposed() const;
    Matrix inversed() const;
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

template<class T>
Matrix<T>::Matrix(size_t _rows, size_t _cols): rows(_rows), cols(_cols) {
    if (_rows == 0 || _cols == 0)
        throw std::runtime_error("Matrix can't be zero size");
    data = Vector<Vector<T>>(rows);
    for (size_t i = 0 ; i < rows; ++i)
        data[i] = Vector<T>(cols);
}

template<class T>
Matrix<T>::Matrix(T **_data, size_t _rows, size_t _cols) {
    if (_rows == 0 || _cols == 0)
        throw std::runtime_error("Matrix can't be zero size");
    data = Vector<Vector<T>>(rows);
    rows = _rows;
    cols = _cols;
    for (size_t i = 0; i < rows; ++i)
        data[i] = Vector<T>(_data[i], cols);
}

template<class T>
Matrix<T>::Matrix(std::initializer_list<Vector<T>> initializerList) {
    if (initializerList.size() == 0)
        throw std::runtime_error("Matrix can't be zero size");
    rows = initializerList.size();
    cols = initializerList.begin()->size();
    data = Vector<Vector<T>>(rows);
    auto it_dst = data.begin();
    for (auto it_src = initializerList.begin(); it_src < initializerList.end(); ++it_src, ++it_dst)
        *it_dst = *it_src;
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
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> initializerList) {
    if (initializerList.size() == 0)
        throw std::runtime_error("Matrix can't be zero size");
    rows = initializerList.size();
    cols = initializerList.begin()->size();
    if (cols == 0)
        throw std::runtime_error("Matrix can't be zero size");
    data = Vector<Vector<T>>(rows);
    auto it_dst = data.begin();
    for (auto it_src = initializerList.begin(); it_src < initializerList.end(); ++it_src, ++it_dst) {
        if (it_src->size() != cols)
            throw std::runtime_error("The rows in the matrix cannot be of different sizes");
        *it_dst = *it_src;
    }
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& other) {
    rows = other.rows;
    cols = other.cols;
    data = other.data;
}

template<class T>
Matrix<T>::Matrix(Matrix<T> &&other) noexcept {
    if (this == &other)
        return;
    swap(other);
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
Matrix<T>& Matrix<T>::operator=(Matrix&& other) noexcept {
    if (this == &other)
        return *this;
    swap(other);
    return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix &other) {
    if (rows != other.rows || cols != other.cols)
        throw std::runtime_error("Operands must be of the same dimension");
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i)
        result[i] = data[i] - other.data[i];
    return result;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix &other) {
    if (rows != other.rows || cols != other.cols)
        throw std::runtime_error("Operands must be of the same dimension");
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i)
        result[i] = data[i] + other.data[i];
    return result;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix &other) {
    if (rows != other.rows || cols != other.cols)
        throw std::runtime_error("Operands must be of the same dimension");
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
        result[i] -= other;
    return result;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const T &other) {
    Matrix result(*this);
    for (size_t i = 0; i < rows; ++i)
        result[i] += other;
    return result;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const T &other) {
    Matrix result(*this);
    for (size_t i = 0; i < rows; ++i)
        result[i] *= other;
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
        if (cols != 1)
            throw std::runtime_error("When multiplying by a vector-column there "
                                     "should be one column in the matrix");
        Matrix result(rows, other.size());
        for (size_t i = 0; i < result.rows; ++i) {
            for (size_t j = 0; j < result.cols; ++j) {
                result[i][j] += data[i][0] * other[j];
            }
        }
        return result;
    }
    if (cols != other.size())
        throw std::runtime_error("The number of columns in the matrix doesn't "
                                 "match the size of the vector");
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
        if (a.size() != b.rows)
            throw std::runtime_error("The number of rows in the matrix doesn't "
                                     "match the size of the vector");
        Matrix<U> result(1, b.cols);
        for (size_t i = 0; i < result.cols; ++i)
            for (size_t j = 0; j < b.rows; ++j)
                result[0][i] += a[j] * b[j][i];
        return result;
    }
    if (b.rows != 1)
        throw std::runtime_error("When multiplying by a vector-row there "
                                 "should be one row in the matrix");
    Matrix<U> result(a.size(), b.cols);
    for (size_t i = 0; i < result.rows; ++i)
        for (size_t j = 0; j < result.cols; ++j)
            result[i][j] += a[i] * b[0][j];
    return result;
}

template<class T>
Matrix<T> Matrix<T>::dot(const Matrix &other) {
    if (cols != other.rows)
        throw std::runtime_error("The number of columns of the left operand "
                                 "differs from the number of rows of the "
                                 "right operand");
    Matrix result(rows, other.cols);
    for (size_t i = 0; i < result.rows; ++i)
        for (size_t j = 0; j < result.cols; ++j) {
            for (size_t k = 0; k < cols; ++k)
                result[i][j] += data[i][k] * other[k][j];
        }
    return result;
}

template<class T>
T Matrix<T>::determinant() const {
    if (rows != cols)
        throw std::runtime_error("To get the determinant of matrix, the matrix "
                                 "must be square");
    if (rows == 1)
        return data[0][0];
    if (rows == 2)
        return data[0][0] * data[1][1] - data[0][1] * data[1][0];
    if (rows == 3)
        return data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) -
        data[0][1] * (data[1][0] * data[2][2] - data[2][0] * data[1][2]) +
        data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]);
    Matrix<T> tmp(*this);
    T det(1);
    for (size_t i = 0; i < rows; ++i) {
        size_t k = i;
        for (size_t j = i + 1; j < rows; ++j)
            if (abs(tmp[j][i]) > abs(tmp[k][i]))
                k = j;
        if (abs(tmp[k][i]) < DBL_EPSILON) {
            det = 0;
            break;
        }
        tmp[i].swap(tmp[k]);
        if (i != k)
            det = -det;
        det *= tmp[i][i];
        for (size_t j = i + 1; j < rows; ++j)
            tmp.data[i][j] /= tmp[i][i];
        for (size_t j = 0; j < rows; ++j)
            if (j != i && abs(tmp[j][i]) > DBL_EPSILON)
                for (size_t q = i + 1; q < rows; ++q)
                    tmp[j][q] -= tmp[i][q] * tmp[j][i];
    }
    return det;
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
Vector<T> Matrix<T>::getColumn(size_t idx) const {
    Vector<T> result(rows, Column);
    for (size_t i = 0; i < rows; ++i)
        result[i] = data[i][idx];
    return result;
}

template<class T>
Matrix<T> Matrix<T>::inversed() const {
    if (rows != cols)
        throw std::runtime_error("To get the inverse matrix, the matrix must "
                                 "be square");
    Matrix result(rows, rows);
    Matrix tmp(*this);
    for (size_t i = 0; i < rows; ++i)
        result[i][i] = 1;

    Matrix<T> concat(rows, 2 * rows);
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < rows; ++j) {
            concat[i][j] = tmp[i][j];
            concat[i][j + rows] = result[i][j];
        }

    for (size_t k = 0; k < rows; ++k) {
        for (size_t i = 0; i < 2 * rows; ++i)
            concat[k][i] /= tmp[k][k];
        for (size_t i = k + 1; i < rows; ++i) {
            T K = concat[i][k] / concat[k][k];
            for (size_t j = 0; j < 2 * rows; ++j)
                concat[i][j] -= concat[k][j] * K;
        }
        for (size_t i = 0; i < rows; ++i)
            for (size_t j = 0; j < rows; ++j)
                tmp[i][j] = concat[i][j];
    }
    int rows_int = static_cast<int>(rows);
    for (int k = rows_int - 1; k >= 0; --k) {
        for (int i = 2 * rows_int - 1; i >= 0; --i)
            concat[k][i] /= tmp[k][k];
        for (int i = k - 1; i >= 0; --i) {
            T K = concat[i][k] / concat[k][k];
            std::cerr << k << " " << concat[k][k] << "\n";
            for (int j = 2 * rows_int - 1; j >= 0; --j)
                concat[i][j] -= concat[k][j] * K;
        }
    }
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < rows; ++j)
            result[i][j] = concat[i][j + rows];
    return result;
}

template<class T>
Matrix<T> Matrix<T>::transposed() const {
    Matrix result(cols, rows);
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            result.data[j][i] = data[i][j];
    return result;
}

template<class T>
void Matrix<T>::swap(Matrix<T> &other) noexcept {
    std::swap(data, other.data);
    std::swap(rows, other.rows);
    std::swap(cols, other.cols);
}

#endif //TECHNOPARK_CPP_TASK_2_MATRIX_HPP
