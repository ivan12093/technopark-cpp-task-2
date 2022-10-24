#ifndef TECHNOPARK_CPP_TASK_2_VECTOR_HPP
#define TECHNOPARK_CPP_TASK_2_VECTOR_HPP

#include <utility>
using std::size_t;

using Format = enum Format_{Row, Column};

template<class T>
class Vector {
private:
    T* _data = nullptr;
    size_t _size = 0;
    Format _format = Row;
public:
    Vector() = default;
    explicit Vector(size_t size, Format format = Row);
    Vector(T* data, size_t size, Format format = Row);
    Vector(std::initializer_list<T> initializerList, Format = Row);
    Vector(const Vector& other);
    ~Vector();
    Vector(Vector&& other) noexcept;
    Vector& operator=(const Vector& other);
    Vector& operator=(Vector&& other) noexcept;

    Vector operator-(const Vector& other) const;
    Vector operator+(const Vector& other) const;
    Vector operator*(const Vector& other) const;

    Vector& operator-=(const Vector& other);
    Vector& operator+=(const Vector& other);
    Vector& operator*=(const Vector& other);

    Vector operator*(const T& elem) const;
    Vector operator-(const T& elem) const;
    Vector operator+(const T& elem) const;

    Vector& operator-=(const T& elem);
    Vector& operator+=(const T& elem);
    Vector& operator*=(const T& elem);

    template<class U>
    friend Vector<U> operator*(const U& a, const Vector<U>& b);
    template<class U>
    friend Vector<U> operator+(const U& a, const Vector<U>& b);
    template<class U>
    friend Vector<U> operator-(const U& a, const Vector<U>& b);

    T operator[](size_t idx) const;
    T& operator[](size_t idx);

    size_t size() const;
    void set_format(Format format);
    Format get_format() const;
    void swap(Vector& other) noexcept;

    T* begin() const;
    T* end() const;
};

#include <algorithm>
#include <cassert>

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
    assert(initializerList.size() != 0);
    std::copy(initializerList.begin(), initializerList.end(), _data);
}

template<class T>
Vector<T>::Vector(const Vector &other): Vector(other._size, other._format) {
    std::copy(other.begin(), other.end(), _data);
}

template<class T>
T* Vector<T>::begin() const {
    return _data;
}

template<class T>
T* Vector<T>::end() const {
    return _data + _size;
}

template<class T>
void Vector<T>::swap(Vector &other) noexcept {
std::swap(_data, other._data);
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
        _data[i] = _data[i] - other[i];
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
        result[i] = a - b[i];
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
Format Vector<T>::get_format() const {
    return _format;
}


#endif //TECHNOPARK_CPP_TASK_2_VECTOR_HPP
