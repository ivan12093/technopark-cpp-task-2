#ifndef TECHNOPARK_CPP_TASK_2_VECTOR_HPP
#define TECHNOPARK_CPP_TASK_2_VECTOR_HPP

using std::size_t;
#include <utility>

using Format = enum {kRow, kColumn};

template<class T, Format format = kRow>
class Vector {
private:
    T* _data = nullptr;
    size_t _size = 0;
public:
    Vector() = delete;
    explicit Vector(size_t size);
    Vector(T* data, size_t size);
    Vector(std::initializer_list<T> initializerList);
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

    friend Vector operator*(const T& a, const Vector& b);
    friend Vector operator+(const T& a, const Vector& b);
    friend Vector operator-(const T& a, const Vector& b);

    T operator[](size_t idx) const;
    T& operator[](size_t idx);

    void set_format(Format new_format);

    size_t size() const;
    void swap(Vector& other) noexcept;

    using iterator = T*;
    using const_iterator = const T*;
    iterator begin();
    iterator end();
};

#endif //TECHNOPARK_CPP_TASK_2_VECTOR_HPP
