#include <gtest/gtest.h>
#include "solution_SLAE/method_GMRES.h"

#include "primitives/Arnoldi.h"
#include "primitives/Col_Major_Matrix.h"
#include "primitives/Upper_Triangular_Matrix.h"

TEST(GMRES, Upper_Triangular_Matrix) {
    Upper_Triangular_Matrix<double> R{{1., 2., 3., 4., 5., 6., 7., 8., 9., 10.}, 4};
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            std::cout << R(i, j) << " ";
        }
        std::cout << '\n';
    }
}

TEST(GMRES, Col_Major_Matrix) {
    Col_Major_Matrix<double> V{{1., 2., 3., 4., 5., 6., 7., 8.}, 4, 2};

    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 5.0;
    DOK[{1, 1}] = 4.0;
    DOK[{2, 2}] = 2.0;
    DOK[{0, 1}] = -1.0;
    DOK[{1, 0}] = -1.0;
    DOK[{1, 2}] = -1.0;
    DOK[{2, 1}] = -1.0;

    const CSR_matrix<double> A{DOK, 3, 3};
    const Vector<double> b{{0, 1., 0}, 3};
    const Vector<double> x_0{3};
    const double eps = 1e-12;


    Vector<double> res = method_GMRES(A, b, x_0, 1000, eps);
}

TEST(GMRES, method) {
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


    Vector<double> res = method_GMRES(A, b, x_0, 1000, eps);
}
