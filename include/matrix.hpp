#ifndef TECHNOPARK_CPP_TASK_2_MATRIX_HPP
#define TECHNOPARK_CPP_TASK_2_MATRIX_HPP

#include <utility>
using std::size_t;

#include "vector.hpp"

template<class T>
class Matrix {
private:
    Vector<Vector<T>> data = nullptr;
    size_t rows = 0;
    size_t cols = 0;
public:
    Matrix() = delete;
    Matrix(size_t _rows, size_t _cols);
    Matrix(T** _data, size_t _rows, size_t _cols);
    explicit Matrix(Vector<Vector<T>> vec_of_vec);
    Matrix(Vector<T>* vectors, size_t size);
    Matrix(std::initializer_list<Vector<T>> initializerList);
    Matrix(std::initializer_list<std::initializer_list<T>> initializerList);

    ~Matrix();
    Matrix(const Matrix& other);
    Matrix(Matrix&& other) noexcept;
    Matrix& operator=(const Matrix& other);
    Matrix& operator=(Matrix&& other) noexcept;

    Matrix transposed();
    //Matrix inversed();
    T determinant();

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

#endif //TECHNOPARK_CPP_TASK_2_MATRIX_HPP
