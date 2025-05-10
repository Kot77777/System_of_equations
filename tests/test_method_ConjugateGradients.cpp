#include <gtest/gtest.h>

#include "solution_SLAE/method_ConjugateGradients.h"

TEST(method_ConjugateGradients, method_ConjugateGradients) {
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

    const Vector<double> res = method_ConjugateGradients(A, b, x_0, 1000, eps);
}

TEST(method_ConjugateGradients, method_ConjugateGradients2) {
    std::map<std::array<std::size_t, 2>, double> DOK{

                {{0, 0}, 10.0},
                {{0, 1}, 3.0},
                {{0, 2}, 6.0},
                {{1, 0}, 3.0},
                {{1, 1}, 5.0},
                {{1, 2}, 1.0},
                {{2, 0}, 6.0},
                {{2, 1}, 1.0},
                {{2, 2}, 8.0}
    };


    const CSR_matrix<double> A{DOK, 3, 3};
    const Vector<double> b{{2., 2., 2.}, 3};
    const double eps = 1e-15;
    const Vector<double> x_0{3};

    const Vector<double> res = method_ConjugateGradients(A, b, x_0, 1000, eps);
}

