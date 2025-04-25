#include <gtest/gtest.h>
#include "primitives/CSR_matrix.h"
#include "primitives/matrix_array.h"
#include "primitives/matrix_from_vector.h"

TEST(CSR_matrix, constructor) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 1.0;
    DOK[{2, 1}] = 1.0;
    DOK[{1, 2}] = 4.0;
    DOK[{0, 3}] = 3.0;
    DOK[{0, 1}] = 2.0;
    DOK[{2, 3}] = 11.0;

    const CSR_matrix<double> matrix_1{DOK, 3, 4};
    const Matrix_arr<3, 4, double> A{{1, 2, 0, 3, 0, 0, 4, 0, 0, 1, 0, 11}};

    for(std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            ASSERT_NEAR(matrix_1(i, j), A(i, j), 10e-15);
        }
    }
}

TEST(CSR_matrix, dot_CSR_matrix_on_vector_1) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 1.0;
    DOK[{2, 1}] = 1.0;
    DOK[{1, 2}] = 4.0;
    DOK[{0, 3}] = 3.0;
    DOK[{0, 1}] = 2.0;
    DOK[{2, 3}] = 11.0;

    const CSR_matrix<double> matrix_1{DOK, 3, 4};
    const Matrix<double> A{{1, 2, 0, 3, 0, 0, 4, 0, 0, 1, 0, 11}, 3, 4};
    const Vector<double> vec{{2, 3, 4, 5}, 4};

    const Vector<double> res1 = matrix_1 * vec;
    const Vector<double> res2 = A * vec;

    const Vector<double> exp{{23, 16, 58}, 3};

    ASSERT_EQ(res1 == exp, 1);
    ASSERT_EQ(res2 == exp, 1);
}

TEST(CSR_matrix, dot_CSR_matrix_on_vector_2) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 1.0;
    DOK[{2, 1}] = 1.0;
    DOK[{1, 2}] = 4.0;
    DOK[{0, 3}] = 3.0;
    DOK[{0, 1}] = 2.0;
    DOK[{2, 3}] = 11.0;

    const CSR_matrix<double> matrix_1{DOK, 4, 4};
    const Matrix<double> A{{1, 2, 0, 3, 0, 0, 4, 0, 0, 1, 0, 11, 0, 0, 0, 0}, 4, 4};
    const Vector<double> vec{{2, 3, 4, 5}, 4};

    const Vector<double> res1 = matrix_1 * vec;
    const Vector<double> res2 = A * vec;

    const Vector<double> exp{{23, 16, 58, 0}, 4};

    ASSERT_EQ(res1 == exp, 1);
    ASSERT_EQ(res2 == exp, 1);
}

TEST(CSR_matrix, dot_transpose_CSR_matrix_on_vector_3) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 1.0;
    DOK[{2, 1}] = 1.0;
    DOK[{1, 2}] = 4.0;
    DOK[{0, 3}] = 3.0;
    DOK[{0, 1}] = 2.0;
    DOK[{2, 3}] = 11.0;
    DOK[{3, 2}] = 2.0;

    const CSR_matrix<double> matrix_1{DOK, 4, 4};
    const Vector<double> vec{{2, 3, 4, 5}, 4};

    const Vector<double> res1 = transpose_multiply(matrix_1, vec);

    const Vector<double> exp{{2, 8, 22, 50}, 4};

    ASSERT_EQ(res1 == exp, 1);
}