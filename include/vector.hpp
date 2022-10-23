#ifndef TECHNOPARK_CPP_TASK_2_VECTOR_HPP
#define TECHNOPARK_CPP_TASK_2_VECTOR_HPP

using std::size_t;
#include <utility>

using Format = enum Format_{Row, Column};

template<class T>
class Vector {
private:
    T* _data = nullptr;
    size_t _size = 0;
    Format _format = Row;
public:
    Vector() = delete;
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
    Format get_format();
    void swap(Vector& other) noexcept;

    T* begin();
    T* end();
};

#endif //TECHNOPARK_CPP_TASK_2_VECTOR_HPP
