#include <gtest/gtest.h>
#include <solution_SLAE/method_Gauss_Seidel.h>

TEST(method_Gauss_Seidel, method_Gauss_Seidel_1) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 8.0;
    DOK[{1, 1}] = 10.0;
    DOK[{2, 2}] = 6.0;
    DOK[{0, 1}] = 2.0;
    DOK[{1, 0}] = 2.0;

    const CSR_matrix<double> matrix_1{DOK, 3, 3};
    Vector<double> b{{1., 2., 4.}, 3};
    const double eps = 10e-15;

    const Vector<double> res = method_Gauss_Seidel(matrix_1, b, 1000, eps);
    const Vector<double> exp{{3./38, 7./38, 2./3}, 3};

    for(std::size_t i = 0; i < 3; ++i) {
        ASSERT_NEAR(res(i), exp(i), eps);
    }
}

TEST(method_Gauss_Seidel, method_Gauss_Seidel_2) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 28.0;
    DOK[{1, 1}] = 10.0;
    DOK[{2, 2}] = 6.0;
    DOK[{0, 1}] = 12.0;
    DOK[{1, 0}] = 2.0;
    DOK[{3, 3}] = 2.0;

    const CSR_matrix<double> matrix_1{DOK, 4, 4};
    Vector<double> b{{1., 2., 4., 8.}, 4};
    const double eps = 10e-15;

    const Vector<double> res = method_Gauss_Seidel(matrix_1, b, 1000, eps);
    const Vector<double> exp{{-7./128, 27./128, 2./3, 4.}, 4};

    for(std::size_t i = 0; i < 3; ++i) {
        ASSERT_NEAR(res(i), exp(i), eps);
    }
}
