#include <gtest/gtest.h>
#include <primitives/matrix_from_vector.h>
#include <primitives/vector_from_vector.h>

TEST(matrix_from_vector, multiply_matrix_on_vector) {
    Vector<double> b{{1, 2, 3}, 3};
    Matrix<double> A{{1, 2, 3, 4, 5, 6, 7, 8, 10}, 3, 3};
    Vector<double> x_res = A * b;
    Vector<double> x_exp{{14, 32, 53}, 3};
    ASSERT_EQ(x_res == x_exp, 1);
}

TEST(vector_from_vector, diada_vectors) {
    const Vector<double> b{{1, 2, 3}, 3};

    Matrix<double> res = diada_for_vector(b);
    Matrix<double> exp{{1, 2, 3, 2, 4, 6, 3, 6, 9}, 3, 3};
    ASSERT_EQ(res == exp, 1);
}