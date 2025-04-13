#include <gtest/gtest.h>
#include <solution_SLAE/method_symmetriz_Gauss_Seidel.h>

TEST(method_Gauss_Seidel, method_Gauss_Seidel_1) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 6.0;
    DOK[{1, 1}] = 5.0;
    DOK[{2, 2}] = 4.0;
    DOK[{3, 3}] = 3.0;
    DOK[{4, 4}] = 4.0;
    DOK[{0, 1}] = 2.0;
    DOK[{1, 0}] = 2.0;
    DOK[{0, 2}] = 1.0;
    DOK[{2, 0}] = 1.0;
    DOK[{1, 2}] = 1.0;
    DOK[{2, 1}] = 1.0;
    DOK[{1, 3}] = 1.0;
    DOK[{3, 1}] = 1.0;
    DOK[{2, 4}] = 1.0;
    DOK[{4, 2}] = 1.0;
    DOK[{3, 4}] = 2.0;
    DOK[{4, 3}] = 2.0;

    const CSR_matrix<double> A{DOK, 5, 5};
    const Vector<double> b{{1., 2., 4., 8., 10.}, 5};
    const double eps = 10e-15;
    const Vector<double> x_0{5};

    const Vector<double> res = method_symmetriz_Gauss_Seidel(A, b, x_0, 1000, eps);
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
    const Vector<double> x_0{4};
    const double eps = 10e-15;

    const Vector<double> res = method_symmetriz_Gauss_Seidel(matrix_1, b, x_0, 1000, eps);
    const Vector<double> exp{{-7. / 128, 27. / 128, 2. / 3, 4.}, 4};

    for (std::size_t i = 0; i < 3; ++i) {
        ASSERT_NEAR(res(i), exp(i), eps);
    }
}

TEST(method_Gauss_Seidel, method_Gauss_Seidel_3) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 28.0;
    DOK[{1, 1}] = 10.0;
    DOK[{2, 2}] = 6.0;
    DOK[{0, 1}] = 2.0;
    DOK[{1, 0}] = 2.0;

    const CSR_matrix<double> A{DOK, 3, 3};
    const Vector<double> b{{1., 2., 4.}, 3};
    const Vector<double> x_0{3};
    const double eps = 1e-15;

    const Vector<double> res = method_symmetriz_Gauss_Seidel(A, b, x_0, 1000, eps);
}
