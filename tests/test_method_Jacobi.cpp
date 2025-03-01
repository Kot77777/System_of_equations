#include <gtest/gtest.h>
#include "solution_SLAE/method_Jacobi.h"

TEST(method_Jacobi, method_Jacobi) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 8.0;
    DOK[{1, 1}] = 10.0;
    DOK[{2, 2}] = 6.0;
    DOK[{0, 1}] = 4.0;
    DOK[{1, 2}] = 4.0;
    DOK[{2, 1}] = 2.0;

    const CSR_matrix<double> matrix_1{DOK, 3, 3};
    Vector<double> b{{1., 2., 4.}, 3};
    const double eps = 10e-15;

    const Vector<double> res = method_Jacobi(matrix_1, b, eps);
    const Vector<double> exp{{17./104, -1./13, 9./13}, 3};

    for(std::size_t i = 0; i < 3; ++i) {
        ASSERT_NEAR(res(i), exp(i), eps);
    }

}
