#include <gtest/gtest.h>
#include "primitives/CSR-matrix.h"
#include "primitives/matrix.h"

TEST(CSR_matrix, constructor) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 1.0;
    DOK[{2, 1}] = 1.0;
    DOK[{1, 2}] = 4.0;
    DOK[{0, 3}] = 3.0;
    DOK[{0, 1}] = 2.0;
    DOK[{2, 3}] = 11.0;

    const CSR_matrix<double> matrix_1{DOK};
    const Matrix<3, 4, double> A{{1, 2, 0, 3, 0, 0, 4, 0, 0, 1, 0, 11}};

    for(std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            ASSERT_NEAR(matrix_1(i, j), A(i, j), 10e-15);
        }
    }
}

TEST(CSR_matrix, dot_CSR_matrix_on_vector) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 1.0;
    DOK[{2, 1}] = 1.0;
    DOK[{1, 2}] = 4.0;
    DOK[{0, 3}] = 3.0;
    DOK[{0, 1}] = 2.0;
    DOK[{2, 3}] = 11.0;

    const CSR_matrix<double> matrix_1{DOK};
    const Matrix<3, 4, double> A{{1, 2, 0, 3, 0, 0, 4, 0, 0, 1, 0, 11}};
    const Vector<4, double> vec{{2, 3, 4, 5}};

    const Vector<3, double> res = A * vec;
    const Vector<3, double> exp{{23, 16, 58}};

    ASSERT_EQ(res == exp, 1);
}
