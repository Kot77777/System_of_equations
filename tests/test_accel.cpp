#include <gtest/gtest.h>

#include "algorithms/lamda_max_for_accel.h"
#include "algorithms/accel_for_symm_method.h"
#include "solution_SLAE/method_symmetriz_Gauss_Seidel.h"
#include "solution_SLAE/method_Jacobi.h"

TEST(accel, lamda_max) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 28.0;
    DOK[{1, 1}] = 10.0;
    DOK[{2, 2}] = 6.0;
    DOK[{0, 1}] = 2.0;
    DOK[{1, 0}] = 2.0;

    const CSR_matrix<double> A{DOK, 3, 3};
    const Vector<double> b{{1., 2., 4.}, 3};
    const Vector<double> x_0{4};
    const double eps = 1e-15;

    const auto symm_Gauss_Seidel = [&A, eps](const Vector<double> &b, const Vector<double> &x_0) {
        return method_symmetriz_Gauss_Seidel(A, b, x_0, 1, eps);
    };

    const double rho_symm_GS = lamda_max_for_accel(symm_Gauss_Seidel, b, 1000);

    ASSERT_NEAR(rho_symm_GS, 1./70, 1e-14);
}

TEST(accel, accel_for_symm_GS_1) {
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

    const auto symm_Gauss_Seidel = [&A, eps](const Vector<double> &b, const Vector<double> &x_0) {
        return method_symmetriz_Gauss_Seidel(A, b, x_0, 1, eps);
    };

    const double rho_symm_GS = lamda_max_for_accel(symm_Gauss_Seidel, b, 1000);

    const Vector<double> res = accel(A, x_0, b, rho_symm_GS, symm_Gauss_Seidel, 1000, eps);
}

TEST(accel, accel_for_symm_GS_2) {
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
    const double eps = 1e-15;
    const Vector<double> x_0{5};

    const auto symm_Gauss_Seidel = [&A, eps](const Vector<double> &b, const Vector<double> &x_0) {
        return method_symmetriz_Gauss_Seidel(A, b, x_0, 1, eps);
    };

    const double rho_symm_GS = lamda_max_for_accel(symm_Gauss_Seidel, b, 1000);

    const Vector<double> res = accel(A, x_0, b, rho_symm_GS, symm_Gauss_Seidel, 1000, eps);
}

TEST(accel, accel_for_symm_GS_3) {
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
    const double eps = 1e-15;
    const Vector<double> x_0{5};

    const auto Jacobi = [&A, eps](const Vector<double> &b, const Vector<double> &x_0) {
        return method_Jacobi(A, b, x_0, 1, eps);
    };

    const double rho_Jacobi = lamda_max_for_accel(Jacobi, b, 1000);

    const Vector<double> res = accel(A, x_0, b, rho_Jacobi, Jacobi, 1000, eps);
}