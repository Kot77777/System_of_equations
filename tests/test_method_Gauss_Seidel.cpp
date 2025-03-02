#include <gtest/gtest.h>
#include <solution_SLAE/method_Gauss_Seidel.h>

TEST(method_Gauss_Seidel, method_Gauss_Seidel_1) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 10.0;
    DOK[{1, 1}] = 9.0;
    DOK[{1, 2}] = 1.0;
    DOK[{2, 1}] = 2.0;
    DOK[{2, 2}] = 8.0;
    DOK[{2, 9}] = 3.0;
    DOK[{3, 0}] = 1.0;
    DOK[{3, 3}] = 7.0;
    DOK[{4, 2}] = 3.0;
    DOK[{4, 4}] = 12.0;
    DOK[{4, 7}] = 4.0;
    DOK[{5, 1}] = 5.0;
    DOK[{5, 5}] = 11.0;
    DOK[{6, 2}] = 6.0;
    DOK[{6, 6}] = 15.0;
    DOK[{7, 3}] = 1.0;
    DOK[{7, 4}] = 2.0;
    DOK[{7, 7}] = 14.0;
    DOK[{8, 0}] = 4.0;
    DOK[{8, 2}] = 5.0;
    DOK[{8, 3}] = 3.0;
    DOK[{8, 8}] = 13.0;
    DOK[{9, 1}] = 6.0;
    DOK[{9, 4}] = 1.0;
    DOK[{9, 6}] = 2.0;
    DOK[{9, 8}] = 3.0;
    DOK[{9, 9}] = 9.0;

    const CSR_matrix<double> matrix_1{DOK, 10, 10};
    Vector<double> b{{1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}, 10};
    const Vector<double> x_0{10};
    const double eps = 10e-15;

    const Vector<double> res = method_Gauss_Seidel(matrix_1, b, x_0, 1000, eps);
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
    const Vector<double> b{{1., 2., 4., 8.}, 4};
    const Vector<double> x_0{10};
    const double eps = 10e-15;

    const Vector<double> res = method_Gauss_Seidel(matrix_1, b, x_0, 1000, eps);
    const Vector<double> exp{{-7. / 128, 27. / 128, 2. / 3, 4.}, 4};

    for (std::size_t i = 0; i < 3; ++i) {
        ASSERT_NEAR(res(i), exp(i), eps);
    }
}
