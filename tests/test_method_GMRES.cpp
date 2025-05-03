#include <gtest/gtest.h>
#include "primitives/Arnoldi.h"
#include "primitives/Col_Major_Matrix.h"
#include "primitives/Upper_Triangular_Matrix.h"

TEST(GMRES, Upper_Triangular_Matrix) {
    Upper_Triangular_Matrix<double> R{{1., 2., 3., 4., 5., 6., 7., 8., 9., 10.}, 4};
    std::cout << R(2,3);
}

TEST(GMRES, Col_Major_Matrix) {
    Col_Major_Matrix<double> V{{1., 2., 3., 4., 5., 6., 7., 8., 9., 10.}, 4};
    std::span<double> v = V.col(1);

}