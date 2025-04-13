#ifndef GENERATE_MATRIX_OF_ELLIPTIC_EQUATION_H
#define GENERATE_MATRIX_OF_ELLIPTIC_EQUATION_H
#include <array>

template<typename T>
CSR_matrix<T> generate_matrix(const std::size_t N) {
    Vector<T> b{N * N};
    std::map<std::array<std::size_t, 2>, T> DOK;

    const T h = 1. / (N + 1);
    const T h_kv = h * h;

    for (std::size_t i = 0; i < N * N; ++i) {
        DOK[{i, i}] = -4. / h_kv;
    }

    for (std::size_t i = 0; i < N * N - (N - 1); i += N) {
        for (std::size_t j = 0; j < N - 1; ++j) {
            DOK[{i + j, i + j + 1}] = 1. / h_kv;
            DOK[{i + j + 1, i + j}] = 1. / h_kv;
        }
    }

    for (std::size_t i = 0; i < N * N - N; ++i) {
        DOK[{i, i + N}] = 1. / h_kv;
        DOK[{i + N, i}] = 1. / h_kv;
    }

    return CSR_matrix<T>{DOK, N*N, N*N};
}

#endif //GENERATE_MATRIX_OF_ELLIPTIC_EQUATION_H
