#ifndef SOLUTION_SLAE_METHOD_SYMMETRIZ_GAUSS_SEIDEL_H
#define SOLUTION_SLAE_METHOD_SYMMETRIZ_GAUSS_SEIDEL_H
#include "primitives/CSR_matrix.h"
#include "primitives/vector_from_vector.h"
#include "algorithms/stop_cond.h"

template<typename T>
Vector<T> method_symmetriz_Gauss_Seidel(const CSR_matrix<T> &A, const Vector<T> &b, const Vector<T> &x_0, const std::size_t N_iter, const T eps) {
    Vector<T> x_i{x_0};
    std::size_t count{}, n_iter{};
    T sum_1{}, sum_2{};
    while (!cond_stop(A, b, x_i, eps) && n_iter < N_iter) {
        for (std::size_t i = 0; i != A.N(); ++i) {
            sum_1 = 0;
            sum_2 = 0;
            const std::size_t k_start = A.rows(i);
            const std::size_t k_end = A.rows(i + 1);
            for (std::size_t k = k_start; A.cols(k) < i; ++k) {
                sum_1 += A.values(k) * x_i(A.cols(k));
                count += 1;
            }
            for (std::size_t k = k_start + count + 1; k < k_end; ++k) {
                sum_2 += A.values(k) * x_i(A.cols(k));
            }
            x_i(i) = (b(i) - sum_1 - sum_2) / A(i, i);
            count = 0;
        }

        for (std::size_t i = A.N() - 1; i != -1; --i) {
            sum_1 = 0;
            sum_2 = 0;
            const std::size_t k_start = A.rows(i);
            const std::size_t k_end = A.rows(i + 1);
            for (std::size_t k = k_start; A.cols(k) < i; ++k) {
                sum_1 += A.values(k) * x_i(A.cols(k));
                count += 1;
            }
            for (std::size_t k = k_start + count + 1; k < k_end; ++k) {
                sum_2 += A.values(k) * x_i(A.cols(k));
            }
            x_i(i) = (b(i) - sum_1 - sum_2) / A(i, i);
            count = 0;
        }
        n_iter += 1;
    }
    return x_i;
}

#endif //SOLUTION_SLAE_METHOD_SYMMETRIZ_GAUSS_SEIDEL_H
