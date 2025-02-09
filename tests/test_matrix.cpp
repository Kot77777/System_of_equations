#include <gtest/gtest.h>
#include <primitives/matrix.h>

TEST(matrix, determinant) {
    Matrix<3, 3, double> matrix{{1, 2, 3, 4, 5, 6, 7, 8, 10}};

    ASSERT_NEAR(determinant(matrix), -3, 10e-15);
}

TEST(matrix, dot_matrix_on_vector) {
    Vector<3, double> b{{1, 2, 3,}};
    Matrix<3, 3, double> A{{1, 2, 3, 4, 5, 6, 7, 8, 10}};
    Vector<3, double> x_res = A * b;
    Vector<3, double> x_exp{{14, 32, 53}};

    ASSERT_EQ(x_res == x_exp, 1);
}

TEST(vector, dot_vector_on_number) {
    Vector<3, double> b{{1, 2, 3,}};
    const double c = 3;
    Vector<3, double> x_res_1 = c * b;
    Vector<3, double> x_res_2 = b * c;
    Vector<3, double> x_exp{{3, 6, 9}};

    ASSERT_EQ(x_res_1 == x_exp, 1);
    ASSERT_EQ(x_res_2 == x_exp, 1);
}

TEST(vector, dot_vector_on_vector) {
    Vector<3, double> b{{1, 2, 3,}};
    Vector_transp<3, double> c = transposed(b);
    const double res = c * b;

    ASSERT_NEAR(res, 14, 1e-15);
}

TEST(vector, sum_vector_and_vector) {
    Vector<3, double> b{{1, 2, 3,}};
    Vector<3, double> d{{2, 2, 3,}};
    Vector<3, double> res = b + d;
    Vector<3, double> exp{{3, 4, 6}};

    ASSERT_EQ(res == exp, 1);
}