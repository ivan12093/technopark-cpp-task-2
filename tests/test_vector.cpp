#include <gtest/gtest.h>

#include "vector.hpp"

TEST(Build, From_raw_data) {
    int a[] = {1, 2, 3, 4};
    Vector<int> vec(a, 4);
    for (size_t i = 0; i < vec.size(); ++i)
        EXPECT_EQ(vec[i], a[i]);
}

TEST(Build, From_initializer_list) {
    int a[] = {4, 3, 2, 1, 0};
    Vector<int> vec = {4, 3, 2, 1, 0};
    for (size_t i = 0; i < vec.size(); ++i)
        EXPECT_EQ(vec[i], a[i]);
}

