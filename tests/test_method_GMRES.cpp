#include <gtest/gtest.h>
#include "solution_SLAE/method_GMRES.h"

#include "primitives/Arnoldi.h"
#include "primitives/Col_Major_Matrix.h"
#include "primitives/Upper_Triangular_Matrix.h"

TEST(GMRES, Upper_Triangular_Matrix) {
    Upper_Triangular_Matrix<double> R{{1., 2., 3., 4., 5., 6., 7., 8., 9., 10.}, 4};
    std::cout << R(2,3);
}

TEST(GMRES, Col_Major_Matrix) {
    Col_Major_Matrix<double> V{{1., 2., 3., 4., 5., 6., 7., 8., 9., 10.}, 4};

    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 4.0;
    DOK[{1, 1}] = 4.0;
    DOK[{2, 2}] = 4.0;
    DOK[{0, 1}] = -1.0;
    DOK[{1, 0}] = -1.0;
    DOK[{0, 3}] = -1.0;
    DOK[{3, 0}] = -1.0;
    DOK[{1, 2}] = -1.0;
    DOK[{2, 1}] = -1.0;
    DOK[{2, 3}] = -1.0;
    DOK[{3, 2}] = -1.0;
    DOK[{3, 3}] = 4.0;

    const CSR_matrix<double> A{DOK, 4, 4};
    const Vector<double> b{{4., 2., 4., 3.}, 4};
    const Vector<double> x_0{4};
    const double eps = 1e-12;


    Vector<double> res = method_GMRES(A, b, x_0, 1000, eps);
}