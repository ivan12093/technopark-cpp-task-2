#include <algorithm>
#include <cassert>

#include "vector.hpp"

template<class T, Format format>
Vector<T, format>::Vector(T *data, size_t size): Vector(size) {
    std::copy(data, data + size, _data);
}

template<class T, Format format>
Vector<T, format>::Vector(size_t size): _size(size) {
    _data = new T[size]();
}

template<class T, Format format>
Vector<T, format>::~Vector() {
    delete[] _data;
    _data = nullptr;
    _size = 0;
}

template<class T, Format format>
Vector<T, format>::Vector(std::initializer_list<T> initializerList): Vector(initializerList.size()) {
    std::copy(initializerList.begin(), initializerList.end(), _data);
}

template<class T, Format format>
Vector<T, format>::Vector(const Vector &other): Vector(other._size) {
    std::copy(other.begin(), other.end(), _data);
}

template<class T, Format format>
T *Vector<T, format>::begin() {
    return _data;
}

template<class T, Format format>
T *Vector<T, format>::end() {
    return _data + _size;
}

template<class T, Format format>
void Vector<T, format>::swap(Vector &other) noexcept {
    std::swap(_data, other.data);
    std::swap(_size, other._size);
}

template<class T, Format format>
Vector<T, format>::Vector(Vector &&other) noexcept {
    if (this != &other)
        swap(other);
}

template<class T, Format format>
Vector<T, format>& Vector<T, format>::operator=(const Vector &other) {
    if (this == &other)
        return *this;
    Vector<T, format> tmp(other);
    swap(tmp);
    return *this;
}

template<class T, Format format>
Vector<T, format>& Vector<T, format>::operator=(Vector &&other) noexcept {
    if (this != &other)
        swap(other);
    return *this;
}

template<class T, Format format>
Vector<T, format> Vector<T, format>::operator-(const Vector &other) const {
    assert(_size == other.size);
    Vector<T, format> result(size());
    for (size_t i = 0; i < size(); ++i)
        result[i] = _data[i] - other[i];
    return result;
}

template<class T, Format format>
Vector<T, format> Vector<T, format>::operator+(const Vector &other) const {
    assert(size() == other.size());
    Vector<T, format> result(size());
    for (size_t i = 0; i < size(); ++i)
        result[i] = _data[i] + other[i];
    return result;
}

template<class T, Format format>
Vector<T, format> Vector<T, format>::operator*(const Vector &other) const {
    assert(size() == other.size());
    Vector<T, format> result(size());
    for (size_t i = 0; i < _size; ++i)
        result[i] = _data[i] * other[i];
    return result;
}

template<class T, Format format>
Vector<T, format> operator*(const T& a, const Vector<T, format>& b) {
    Vector<T, format> result(b.size());
    for (size_t i = 0; i < b.size(); ++i)
        result[i] = b[i] * a;
    return result;
}

template<class T, Format format>
Vector<T, format> Vector<T, format>::operator*(const T &elem) const {
    Vector<T, format> result(size());
    for (size_t i = 0; i < size(); ++i)
        result[i] = _data[i] * elem;
    return result;
}

template<class T, Format format>
Vector<T, format>& Vector<T, format>::operator*=(const T &elem) {
    for (size_t i = 0; i < size(); ++i)
        _data[i] = _data[i] * elem;
    return *this;
}

template<class T, Format format>
Vector<T, format> Vector<T, format>::operator-(const T &elem) const {
    Vector<T, format> result(size());
    for (size_t i = 0; i < size(); ++i)
        result[i] = _data[i] - elem;
    return result;
}

template<class T, Format format>
Vector<T, format> Vector<T, format>::operator+(const T &elem) const {
    Vector<T, format> result(size());
    for (size_t i = 0; i < size(); ++i)
        result[i] = _data[i] + elem;
    return result;
}

template<class T, Format format>
Vector<T, format>& Vector<T, format>::operator-=(const Vector &other) {
    for (size_t i = 0; i < size(); ++i)
        _data[i] = _data[i] * other[i];
    return *this;
}

template<class T, Format format>
Vector<T, format>& Vector<T, format>::operator+=(const Vector &other) {
    for (size_t i = 0; i < size(); ++i)
        _data[i] = _data[i] + other[i];
    return *this;
}

template<class T, Format format>
Vector<T, format>& Vector<T, format>::operator*=(const Vector &other) {
    for (size_t i = 0; i < size(); ++i)
        _data[i] = _data[i] * other[i];
    return *this;
}

template<class T, Format format>
Vector<T, format>& Vector<T, format>::operator-=(const T &elem) {
    for (size_t i = 0; i < size(); ++i)
        _data[i] = _data[i] - elem;
    return *this;
}

template<class T, Format format>
Vector<T, format>& Vector<T, format>::operator+=(const T &elem) {
    for (size_t i = 0; i < size(); ++i)
        _data[i] = _data[i] + elem;
    return *this;
}

template<class T, Format format>
Vector<T, format> operator+(const T& a, const Vector<T, format>& b) {
    Vector<T, format> result(b.size());
    for (size_t i = 0; i < b.size(); ++i)
        result[i] = b[i] + a;
    return result;
}

template<class T, Format format>
Vector<T, format> operator-(const T& a, const Vector<T, format>& b) {
    Vector<T, format> result(b.size());
    for (size_t i = 0; i < b.size(); ++i)
        result[i] = b[i] - a;
    return result;
}

template<class T, Format format>
T Vector<T, format>::operator[](size_t idx) const {
    return _data[idx];
}

template<class T, Format format>
T& Vector<T, format>::operator[](size_t idx) {
    return _data[idx];
}

template<class T, Format format>
size_t Vector<T, format>::size() const {
    return _size;
}

template<class T, Format format>
void Vector<T, format>::set_format(Format new_format) {
}


