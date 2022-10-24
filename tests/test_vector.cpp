#include <gtest/gtest.h>

#include "vector.hpp"

TEST(Build_Vector, From_size) {
    Vector<int> vec(4);
    for (size_t i = 0; i < 4; ++i)
        EXPECT_EQ(vec[i], 0);
}

TEST(Build_Vector, From_raw_data) {
    int a[] = {1, 2, 3, 4};
    Vector<int> vec(a, 4);
    for (size_t i = 0; i < vec.size(); ++i)
        EXPECT_EQ(vec[i], a[i]);
}

TEST(Build_Vector, From_initializer_list) {
    int a[] = {4, 3, 2, 1, 0};
    Vector<int> vec = {4, 3, 2, 1, 0};
    for (size_t i = 0; i < vec.size(); ++i)
        EXPECT_EQ(vec[i], a[i]);
}

TEST(Copy_Vector, From_copy_constructor) {
    Vector<int> vec1 = {1, 2, 3, 4, 5};
    Vector<int> vec2(vec1);
    for (size_t i = 0; i < vec1.size(); ++i)
        EXPECT_EQ(vec1[i], vec2[i]);
}

TEST(Copy_Vector, From_assign_operator) {
    Vector<int> vec1 = {1, 2, 3, 4, 5};
    Vector<int> vec2 = vec1;
    for (size_t i = 0; i < vec1.size(); ++i)
        EXPECT_EQ(vec1[i], vec2[i]);
}

TEST(Move_Vector, From_constructor) {
    int res[] = {1, 2, 3, 4, 5};
    Vector<int> vec1 = {1, 2, 3, 4, 5};
    Vector<int> vec2(std::move(vec1));
    EXPECT_EQ(vec2.size(), 5);
    for (size_t i = 0; i < vec2.size(); ++i)
        EXPECT_EQ(vec2[i], res[i]);
}

TEST(Move_Vector, From_assign_operator) {
    int res[] = {1, 2, 3, 4, 5};
    Vector<int> vec1 = {1, 2, 3, 4, 5};
    Vector<int> vec2 = std::move(vec1);
    EXPECT_EQ(vec2.size(), 5);
    for (size_t i = 0; i < vec2.size(); ++i)
        EXPECT_EQ(vec2[i], res[i]);
}

TEST(Vector_math, With_other_vectors) {
    Vector<int> vec1 = {1, 2, 3, 4, 5};
    Vector<int> vec2 = {3, 4, 5, 6, 7};
    Vector<int> vec3 = vec2 - vec1; // vec3 = {2, 2, 2, 2, 2}
    Vector<int> res = {2, 2, 2, 2, 2};
    for (size_t i = 0; i < vec3.size(); ++i)
        EXPECT_EQ(vec3[i], res[i]);
    vec2 -= vec1; // vec2 = {2, 2, 2, 2, 2}
    for (size_t i = 0; i < vec2.size(); ++i)
        EXPECT_EQ(vec2[i], res[i]);
    vec2 += vec1; // vec2 {3, 4, 5, 6, 7}
    res = {3, 4, 5, 6, 7};
    for (size_t i = 0; i < vec2.size(); ++i)
        EXPECT_EQ(vec2[i], res[i]);
    Vector<int> vec4 = {1, 0, 0, 0, 0};
    Vector<int> vec5 = vec4 * vec2; // {3, 0, 0, 0, 0}
    res = {3, 0, 0, 0, 0};
    for (size_t i = 0; i < vec5.size(); ++i)
        EXPECT_EQ(vec5[i], res[i]);
    vec5 *= vec1; // {3, 0, 0, 0, 0}
    for (size_t i = 0; i < vec5.size(); ++i)
        EXPECT_EQ(vec5[i], res[i]);
    vec5 = vec1 + vec2; // {4, 6, 8, 10, 12}
    res = {4, 6, 8, 10, 12};
    for (size_t i = 0; i < vec5.size(); ++i)
        EXPECT_EQ(vec5[i], res[i]);
}

TEST(Vector_math, With_scalars) {
    Vector<int> vec = {1, 2, 3, 4, 5};
    vec -= 1;
    Vector<int> res = {0, 1, 2, 3, 4};
    for (size_t i = 0; i < vec.size(); ++i)
        EXPECT_EQ(vec[i], res[i]);
    vec += 1;
    res = {1, 2, 3, 4, 5};
    for (size_t i = 0; i < vec.size(); ++i)
        EXPECT_EQ(vec[i], res[i]);
    vec *= 2;
    res = {2, 4, 6, 8, 10};
    for (size_t i = 0; i < vec.size(); ++i)
        EXPECT_EQ(vec[i], res[i]);
    Vector<int> vec1 = vec - 2;
    res = {0, 2, 4, 6, 8};
    for (size_t i = 0; i < vec1.size(); ++i)
        EXPECT_EQ(vec1[i], res[i]);
    vec1 = vec + 2;
    res = {4, 6, 8, 10, 12};
    for (size_t i = 0; i < vec.size(); ++i)
        EXPECT_EQ(vec1[i], res[i]);
    vec1 = vec * 2;
    res = {4, 8, 12, 16, 20};
    for (size_t i = 0; i < vec.size(); ++i)
        EXPECT_EQ(vec1[i], res[i]);
}

TEST(Vector_math, With_scalars_rightside) {
    Vector<int> vec = {1, 2, 3, 4, 5};
    Vector<int> vec1 = 2 - vec;
    Vector<int> res = {1, 0, -1, -2, -3};
    for (size_t i = 0; i < vec1.size(); ++i)
        EXPECT_EQ(vec1[i], res[i]);
    vec1 = 2 + vec;
    res = {3, 4, 5, 6, 7};
    for (size_t i = 0; i < vec.size(); ++i)
        EXPECT_EQ(vec1[i], res[i]);
    vec1 = vec * 2;
    res = {2, 4, 6, 8, 10};
    for (size_t i = 0; i < vec.size(); ++i)
        EXPECT_EQ(vec1[i], res[i]);
}
