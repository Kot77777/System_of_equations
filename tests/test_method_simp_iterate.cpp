#include <gtest/gtest.h>
#include <solution_SLAE/method_simple_iteration.h>
#include "algorithms/lamda_max.h"
#include <fstream>

TEST(method_simple_iteration, method_simple_iteration_1) {
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
    const double lamd_min = 2.;
    const double lamd_max = lamda_max(matrix_1, b, 1000);

    const Vector<double> res = method_simple_iteration(matrix_1, b, x_0, 2. / (lamd_max + lamd_min),1000, eps);
    const Vector<double> exp{{-7. / 128, 27. / 128, 2. / 3, 4.}, 4};

    for (std::size_t i = 0; i < 3; ++i) {
        ASSERT_NEAR(res(i), exp(i), eps);
    }
}

TEST(method_simple_iteration, method_simple_iteration_2) {
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
    const double eps = 1e-13;
    const Vector<double> x_0{3};
    const double lamd_min = 2.251;
    const double lamd_max = lamda_max(A, b, 1000);

    const Vector<double> res = method_simple_iteration(A, b, x_0, 2. / (lamd_max + lamd_min),1000, eps);
}
// TEST(method_simple_iteration, method_simple_iteration_2) {
//     std::map<std::array<std::size_t, 2>, double> DOK;
//
//     DOK[{0, 0}] = 89.0;
//     DOK[{1, 1}] = 50.0;
//     DOK[{2, 2}] = 95.0;
//     DOK[{3, 3}] = 93.0;
//     DOK[{0, 2}] = 2.0;
//     DOK[{0, 3}] = 2.0;
//     DOK[{3, 0}] = 9.0;
//     DOK[{3, 2}] = 9.0;
//     DOK[{0, 1}] = 5.0;
//
//     const CSR_matrix<double> matrix_1{DOK, 4, 4};
//     const Vector<double> b{{1., 1., 1., 1.}, 4};
//     const Vector<double> x_0{4};
//     const double eps = 10e-12;
//     const double lamd_max = lamda_max(matrix_1, b, 1000);
//     std::ofstream data_measure("data_measure.csv");
//     data_measure << "n_iterate" << "," << "tau" << '\n';
//     const std::size_t N = 100;
//     std::size_t i = 1;
//     while(i < N) {
//         std::size_t res = method_simple_iteration(matrix_1, b, x_0, 2. / (i*10), 10000, eps);
//         data_measure << res << (2. / (i*10)) << '\n';
//         i+=1;
//     }
//}
