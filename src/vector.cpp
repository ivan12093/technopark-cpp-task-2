#include <algorithm>
#include <cassert>

#include "vector.hpp"

template<class T>
Vector<T>::Vector(T* data, size_t size, Format format): Vector(size, format) {
    assert(size != 0);
    std::copy(data, data + size, _data);
}

template<class T>
Vector<T>::Vector(size_t size, Format format): _size(size), _format(format) {
    assert(size != 0);
    _data = new T[size]();
}

template<class T>
Vector<T>::~Vector() {
    delete[] _data;
    _data = nullptr;
    _size = 0;
}

template<class T>
Vector<T>::Vector(std::initializer_list<T> initializerList, Format format): Vector(initializerList.size(), format) {
    static_assert(initializerList.size() != 0);
    std::copy(initializerList.begin(), initializerList.end(), _data);
}

template<class T>
Vector<T>::Vector(const Vector &other): Vector(other._size, other._format) {
    std::copy(other.begin(), other.end(), _data);
}

template<class T>
T* Vector<T>::begin() {
    return _data;
}

template<class T>
T* Vector<T>::end() {
    return _data + _size;
}

template<class T>
void Vector<T>::swap(Vector &other) noexcept {
    std::swap(_data, other.data);
    std::swap(_size, other._size);
    std::swap(_format, other._format);
}

template<class T>
Vector<T>::Vector(Vector &&other) noexcept {
    if (this != &other)
        swap(other);
}

template<class T>
Vector<T>& Vector<T>::operator=(const Vector &other) {
    if (this == &other)
        return *this;
    Vector<T> tmp(other);
    swap(tmp);
    return *this;
}

template<class T>
Vector<T>& Vector<T>::operator=(Vector &&other) noexcept {
    if (this != &other)
        swap(other);
    return *this;
}

template<class T>
Vector<T> Vector<T>::operator-(const Vector &other) const {
    assert(size() == other.size());
    assert(get_format() == other.get_format());
    Vector<T> result(size(), _format);
    for (size_t i = 0; i < size(); ++i)
        result[i] = _data[i] - other[i];
    return result;
}

template<class T>
Vector<T> Vector<T>::operator+(const Vector &other) const {
    assert(size() == other.size());
    assert(get_format() == other.get_format());
    Vector<T> result(size(), _format);
    for (size_t i = 0; i < size(); ++i)
        result[i] = _data[i] + other[i];
    return result;
}

template<class T>
Vector<T> Vector<T>::operator*(const Vector &other) const {
    assert(size() == other.size());
    assert(get_format() == other.get_format());
    Vector<T> result(size(), _format);
    for (size_t i = 0; i < _size; ++i)
        result[i] = _data[i] * other[i];
    return result;
}

template<class T, class Format>
Vector<T> operator*(const T& a, const Vector<T>& b) {
    Vector<T> result(b.size(), b.get_format());
    for (size_t i = 0; i < b.size(); ++i)
        result[i] = b[i] * a;
    return result;
}

template<class T>
Vector<T> Vector<T>::operator*(const T &elem) const {
    Vector<T> result(size(), get_format());
    for (size_t i = 0; i < size(); ++i)
        result[i] = _data[i] * elem;
    return result;
}

template<class T>
Vector<T>& Vector<T>::operator*=(const T &elem) {
    for (size_t i = 0; i < size(); ++i)
        _data[i] = _data[i] * elem;
    return *this;
}

template<class T>
Vector<T> Vector<T>::operator-(const T &elem) const {
    Vector<T> result(size(), get_format());
    for (size_t i = 0; i < size(); ++i)
        result[i] = _data[i] - elem;
    return result;
}

template<class T>
Vector<T> Vector<T>::operator+(const T &elem) const {
    Vector<T> result(size(), get_format());
    for (size_t i = 0; i < size(); ++i)
        result[i] = _data[i] + elem;
    return result;
}

template<class T>
Vector<T>& Vector<T>::operator-=(const Vector &other) {
    assert(size() == other.size());
    assert(get_format() == other.get_format());
    for (size_t i = 0; i < size(); ++i)
        _data[i] = _data[i] * other[i];
    return *this;
}

template<class T>
Vector<T>& Vector<T>::operator+=(const Vector &other) {
    assert(size() == other.size());
    assert(get_format() == other.get_format());
    for (size_t i = 0; i < size(); ++i)
        _data[i] = _data[i] + other[i];
    return *this;
}

template<class T>
Vector<T>& Vector<T>::operator*=(const Vector &other) {
    assert(size() == other.size());
    assert(get_format() == other.get_format());
    for (size_t i = 0; i < size(); ++i)
        _data[i] = _data[i] * other[i];
    return *this;
}

template<class T>
Vector<T>& Vector<T>::operator-=(const T &elem) {
    for (size_t i = 0; i < size(); ++i)
        _data[i] = _data[i] - elem;
    return *this;
}

template<class T>
Vector<T>& Vector<T>::operator+=(const T &elem) {
    for (size_t i = 0; i < size(); ++i)
        _data[i] = _data[i] + elem;
    return *this;
}

template<class T>
Vector<T> operator+(const T& a, const Vector<T>& b) {
    Vector<T> result(b.size(), b.get_format());
    for (size_t i = 0; i < b.size(); ++i)
        result[i] = b[i] + a;
    return result;
}

template<class T>
Vector<T> operator-(const T& a, const Vector<T>& b) {
    Vector<T> result(b.size(), b.get_format());
    for (size_t i = 0; i < b.size(); ++i)
        result[i] = b[i] - a;
    return result;
}

template<class T>
T Vector<T>::operator[](size_t idx) const {
    return _data[idx];
}

template<class T>
T& Vector<T>::operator[](size_t idx) {
    return _data[idx];
}

template<class T>
size_t Vector<T>::size() const {
    return _size;
}

template<class T>
void Vector<T>::set_format(Format format) {
    _format = format;
}

template<class T>
Format Vector<T>::get_format() {
    return _format;
}
