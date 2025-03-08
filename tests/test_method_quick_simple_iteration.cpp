#include <gtest/gtest.h>
#include "solution_SLAE/method_quick_simple_iteration.h"
#include "algorithms/lamda_max.h"
#include "algorithms/get_polynom_roots.h"
#include "algorithms/permutation.h"

TEST(method_quick_simple_iteration, lamda_max) {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 1.0;
    DOK[{0, 1}] = 2.0;
    DOK[{1, 0}] = 4.0;
    DOK[{1, 1}] = 5.0;
    DOK[{1, 2}] = 6.0;
    DOK[{2, 2}] = 9.0;


    const CSR_matrix<double> matrix_1{DOK, 3, 3};
    const Vector<double> r_0{{1., 1., 1.}, 3};
    const double eps = 10e-15;
    const double lamd_max = lamda_max(matrix_1, r_0, eps);

    ASSERT_NEAR(lamd_max, 9, 10e-14);
}

TEST(method_quick_simple_iteration, method_quick_simple_iteration_1) {
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
    const std::size_t n = 64;
    const double lamd_min = 2.;
    const double lamd_max = lamda_max(matrix_1, b, eps);
    const Vector<double> t{get_polynom_roots<double>(n, (lamd_max - lamd_min) / 2, (lamd_max + lamd_min) / 2)};
    const Vector<std::size_t> perm = permutation(Vector<size_t>{{0, 1}, 2}, n);

    const Vector<double> res = method_quick_simple_iteration(matrix_1, b, x_0, t, perm, n, 1000, eps);
    const Vector<double> exp{{-7. / 128, 27. / 128, 2. / 3, 4.}, 4};

    for (std::size_t i = 0; i < 3; ++i) {
        ASSERT_NEAR(res(i), exp(i), eps);
    }
}

TEST(method_quick_simple_iteration, method_quick_simple_iteration_2) {
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
    const Vector<double> b{{1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}, 10};
    const Vector<double> x_0{10};
    const double eps = 10e-15;
    const std::size_t n = 128;
    const double lamd_min = 6.438447187191168;
    const double lamd_max = lamda_max(matrix_1, b, eps);
    const Vector<double> t{get_polynom_roots<double>(n, (lamd_max - lamd_min) / 2, (lamd_max + lamd_min) / 2)};
    const Vector<std::size_t> perm = permutation(Vector<size_t>{{0, 1}, 2}, n);

    const Vector<double> res = method_quick_simple_iteration(matrix_1, b, x_0, t, perm, n, 1000, eps);

}
