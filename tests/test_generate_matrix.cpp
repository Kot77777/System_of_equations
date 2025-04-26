#include <gtest/gtest.h>
#include "algorithms/generate_matrix_of_elliptic_equation.h"

TEST(generate_matrix, generate_matrix) {
    const std::size_t N = 2;
    CSR_matrix<double> A = generate_matrix<double>(N);

    for (std::size_t i = 0; i < N * N; ++i) {
        for(std::size_t j = 0; j < N * N; ++j) {
            std::cout << A(i, j) << ' ';
        }
        std::cout << '\n';
    }
}