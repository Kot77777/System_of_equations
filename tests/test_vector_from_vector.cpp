#include <gtest/gtest.h>
#include <primitives/vector_from_vector.h>

TEST(vector_from_vector, operator_equal) {
    const Vector<double> b{{1, 2, 3}, 3};
    const Vector<double> c{{1, 2, 3}, 3};

    const Vector<double> b1{{1, 2, 3}, 3};
    const Vector<double> c1{{1, 9, 3}, 3};

    ASSERT_EQ(b == c, 1);
    ASSERT_EQ(b1 == c1, 0);
}

TEST(vector_from_vector, sum_vectors) {
    const Vector<double> b{{1, 2, 3}, 3};
    const Vector<double> c{{1, 2, 3}, 3};

    const Vector<double> x_res = c + b;
    const Vector<double> x_exp{{2, 4, 6}, 3};
    ASSERT_EQ(x_res == x_exp, 1);
}

TEST(vector_from_vector, differ_vectors) {
    const Vector<double> b{{1, 2, 3}, 3};
    const Vector<double> c{{1, 2, 3}, 3};

    const Vector<double> x_res = c - b;
    const Vector<double> x_exp{{0, 0, 0}, 3};
    ASSERT_EQ(x_res == x_exp, 1);
}

TEST(vector_from_vector, multiply_vector_on_number) {
    const Vector<double> b{{1, 2, 3}, 3};
    double c = 3.;

    const Vector<double> x_res = c * b;
    const Vector<double> x_exp{{3, 6, 9}, 3};
    ASSERT_EQ(x_res == x_exp, 1);
}

TEST(vector_from_vector, dot_vectors) {
    const Vector<double> b{{1, 2, 3}, 3};
    const Vector<double> c{{4, 2, 8}, 3};

    const double res = dot(b, c);
    const double exp = 32;
    ASSERT_NEAR(res, exp, 10e-15);
}