#include <gtest/gtest.h>

#include "matrix.hpp"

#define EXPECT_EQ_MATRIX(a, b) ASSERT_EQ(a.size().first, b.size().first); \
ASSERT_EQ(a.size().second, b.size().second);                          \
for (size_t i = 0; i < a.size().first; ++i)                               \
for (size_t j = 0; j < a.size().second; ++j)                              \
EXPECT_EQ(a[i][j], b[i][j]);

#define EXPECT_EQ_MATRIX_DBL(a, b) ASSERT_EQ(a.size().first, b.size().first); \
ASSERT_EQ(a.size().second, b.size().second);                          \
for (size_t i = 0; i < a.size().first; ++i)                               \
for (size_t j = 0; j < a.size().second; ++j)                              \
EXPECT_DOUBLE_EQ(a[i][j], b[i][j]);

#define EXPECT_EQ_VECTOR(a, b) ASSERT_EQ(a.size(), b.size()); \
EXPECT_EQ(a.get_format(), b.get_format());                                \
for (size_t i = 0; i < a.size(); ++i)                         \
EXPECT_EQ(a[i], b[i]);

TEST(Build_Matrix, From_size) {
    Matrix<int> mtx(10, 9);
    EXPECT_EQ(mtx.size().first, 10);
    EXPECT_EQ(mtx.size().second, 9);
    for (size_t i = 0; i < mtx.size().first; ++i)
        for (size_t j = 0; j < mtx.size().second; ++j)
            EXPECT_EQ(mtx[i][j], 0);
}

TEST(Build_Matrix, From_initializer_list) {
    Matrix<int> matrix = {
            {1, 2, 3},
            {1, 2, 3},
            {1, 2, 3}
    };
    int mtx[3][3] = {
            {1, 2, 3},
            {1, 2, 3},
            {1, 2, 3}
    };

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_EQ(matrix[i][j], mtx[i][j]);
}

TEST(Build_Matrix, From_vec_of_vec) {
    Vector<Vector<int>> vec = {
            {1, 2, 3},
            {1, 2, 3},
            {1, 2, 3}
    };
    Matrix<int> matrix(vec);
    ASSERT_EQ(matrix.size().first, 3);
    ASSERT_EQ(matrix.size().second, 3);
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_EQ(matrix[i][j], vec[i][j]);
}

TEST(Build_Matrix, From_initializer_of_vec) {
    int mtx[3][3] = {
            {1, 2, 3},
            {1, 2, 3},
            {1, 2, 3}
    };
    Matrix<int> matrix({
                Vector({1, 2, 3}),
                Vector({1, 2, 3}),
                Vector({1, 2, 3})
    });
    ASSERT_EQ(matrix.size().first, 3);
    ASSERT_EQ(matrix.size().second, 3);
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_EQ(matrix[i][j], mtx[i][j]);
}

TEST(Build_Matrix, From_other_matrix) {
    Matrix<int> matrix({Vector({1, 2, 3}),
                        Vector({1, 2, 3}),
                        Vector({1, 2, 3})
                       });
    Matrix<int> matrix1(matrix);
    EXPECT_EQ_MATRIX(matrix, matrix1);
}

TEST(Build_Matrix, From_other_matrix_assign) {
    Matrix<int> matrix({Vector({1, 2, 3}),
                        Vector({1, 2, 3}),
                        Vector({1, 2, 3})
                       });
    Matrix<int> matrix1 = matrix;
    EXPECT_EQ_MATRIX(matrix, matrix1);
}

TEST(Build_Matrix, From_move_matrix) {
    Matrix<int> matrix({Vector({1, 2, 3}),
                        Vector({1, 2, 3}),
                        Vector({1, 2, 3})
                       });
    Matrix<int> res(matrix);
    Matrix<int> matrix1(std::move(matrix));
    EXPECT_EQ_MATRIX(res, matrix1);
}

TEST(Build_Matrix, From_move_assign_matrix) {
    Matrix<int> matrix({Vector({1, 2, 3}),
                        Vector({1, 2, 3}),
                        Vector({1, 2, 3})
                       });
    Matrix<int> res(matrix);
    Matrix<int> matrix1 = std::move(matrix);
    EXPECT_EQ_MATRIX(res, matrix1);
}

TEST(Matrix_functionality, Transpose_matrix) {
    Matrix<int> res = {
            {1, 2, 3},
            {1, 2, 3}
    };
    Matrix<int> a = {
            {1, 1},
            {2, 2},
            {3, 3}
    };
    Matrix<int> b = a.transposed();
    EXPECT_EQ_MATRIX(res, b);
}

TEST(Matrix_functionality, Math_operations_with_scalars) {
    Matrix<int> a = {
            {1, 2, 3},
            {1, 2, 3}
    };
    a -= 1;
    Matrix<int> res = {
            {0, 1, 2},
            {0, 1, 2}
    };
    EXPECT_EQ_MATRIX(a, res);
    a += 2;
    res = {
            {2, 3, 4},
            {2, 3, 4}
    };
    EXPECT_EQ_MATRIX(a, res);
    a *= 2;
    res = {
            {4, 6, 8},
            {4, 6, 8}
    };

    Matrix<int> b = a - 4;
    res = {
            {0, 2, 4},
            {0, 2, 4}
    };
    EXPECT_EQ_MATRIX(b, res);

    b = a + 2;
    res = {
            {6, 8, 10},
            {6, 8, 10}
    };
    EXPECT_EQ_MATRIX(b, res);

    b = a * 2;
    res = {
            {8, 12, 16},
            {8, 12, 16}
    };
    EXPECT_EQ_MATRIX(b, res);
}

TEST(Matrix_functionality, Math_operations_with_vectors) {
    Vector<int> vec = {1, 2, 3};
    Matrix<int> mtx = {
            {1, 2},
            {3, 4},
            {5, 6}
    };
    Matrix<int> mul_mtx = vec * mtx;
    Matrix<int> res = {
            {22, 28}
    };
    EXPECT_EQ_MATRIX(res, mul_mtx);

    vec.set_format(Column);
    mtx = {
            {1, 2, 3}
    };

    mul_mtx = mtx * vec;

    res = { {14} };
    EXPECT_EQ_MATRIX(res, mul_mtx);
}

TEST(Matrix_functionality, Math_operations_with_matrix) {
    Matrix<int> mtx1 = {
            {1, 2},
            {3, 4},
            {5, 6}
    };
    Matrix<int> mtx2 = {
            {1, 2, 3},
            {4, 5, 6}
    };

    Matrix<int> mtx3 = mtx1.dot(mtx2);
    Matrix<int> res = {
            {9, 12, 15},
            {19, 26, 33},
            {29, 40, 51}
    };
    EXPECT_EQ_MATRIX(mtx3, res);

    Matrix<int> mtx4 = mtx1.transposed();
    Matrix<int> mtx5 = mtx2.transposed();
    EXPECT_EQ(mtx5.size().first, 3);
    EXPECT_EQ(mtx5.size().second, 2);
    Matrix<int> mtx6 = mtx4.dot(mtx5);
    res = {
            {22, 49},
            {28, 64}
    };
    EXPECT_EQ_MATRIX(mtx6, res);
}

TEST(Matrix_functionality, Get_Projections) {
    Matrix<int> mtx = {
            {1, 2, 3},
            {4, 5, 6},
            {7, 8, 9}
    };
    Vector<int> vec = mtx.getDiag();
    Vector<int> res = {1, 5, 9};
    EXPECT_EQ_VECTOR(vec, res);
    vec = mtx.getColumn(1);
    res = {2, 5, 8};
    res.set_format(Column);
    EXPECT_EQ_VECTOR(vec, res);
    mtx[1][1] = 4;
    EXPECT_EQ_VECTOR(vec, res);
    vec[1] = 0;
    EXPECT_EQ(mtx[1][1], 4);
    vec = mtx[1];
    res = {4, 4, 6};
    EXPECT_EQ_VECTOR(vec, res);
    mtx[1][1] = 5;
    EXPECT_EQ_VECTOR(vec, res);
    vec[1] = 0;
    EXPECT_EQ(mtx[1][1], 5);
}


TEST(Matrix_Determinant, Default) {
    Matrix<double> matrix = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 1.0, 9.0}
    };
    double det = matrix.determinant();
    EXPECT_DOUBLE_EQ(det, -42);
}

TEST(Matrix_Inversed, Default) {
    Matrix<double> matrix = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 1.0, 9.0}
    };
    Matrix<double> res = {
            {-13. / 14, 5. / 14, 1. / 14},
            {-1. / 7, 2. / 7, -1. / 7},
            {31. / 42, -13. / 42, 1. / 14}
    };
    auto inversed = matrix.inversed();
    EXPECT_EQ_MATRIX_DBL(inversed, res);
}
