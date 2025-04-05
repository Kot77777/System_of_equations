#include <gtest/gtest.h>
#include "solution_SLAE/method_SOR.h"
#include "algorithms/lamda_max_for_SOR.h"

TEST(method_SOR, method_Gauss_SOR) {
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
    const double w = 1;
    const Vector<double> x_0{5};

    const Vector<double> res = method_SOR(A, b, x_0, 1000, eps, w);
}

TEST(method_SOR, lamda_max) {
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
    const Vector<double> r_0{5};
    const double mu = lamda_max_for_SOR(A, b, 1000);
    const double mu_p = mu / (1.0 + sqrt(1.0 - mu*mu));
    const double w = 1 + mu_p * mu_p;

    const Vector<double> res = method_SOR(A, b, r_0, 1000, eps, w);

}


