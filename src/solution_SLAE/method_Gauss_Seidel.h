#ifndef SOLUTION_SLAE_METHOD_GAUSS_SEIDEL_H
#define SOLUTION_SLAE_METHOD_GAUSS_SEIDEL_H
#include "primitives/CSR_matrix.h"
#include "primitives/vector_from_vector.h"
#include "algorithms/stop_cond.h"

template<typename T>
Vector<T> method_Gauss_Seidel(const CSR_matrix<T> &A, const Vector<T> &b, const std::size_t N_iter, const T eps) {
    Vector<T> x_i{b.N()};
    Vector<T> x_i_next{b.N()};
    std::size_t count{}, n_iter{};
    T d_ii{};
    while (!cond_stop(A, b, x_i, eps) && n_iter < N_iter) {
        for (std::size_t i = 0; i != A.N(); ++i) {
            const std::size_t k_start = A.rows(i);
            const std::size_t k_end = A.rows(i + 1);
            for (std::size_t k = k_start; A.cols(k) < i; ++k) {
                x_i_next(i) += A.values(k) * x_i_next(A.cols(k));
                count += 1;
            }
            d_ii = A.values(k_start + count);
            for (std::size_t k = k_start + count + 1; k < k_end; ++k) {
                x_i_next(i) += A.values(k) * x_i(A.cols(k));
            }
            x_i_next(i) = (b(i) - x_i_next(i)) / d_ii;
            count = 0;
        }
        x_i = x_i_next;
        x_i_next.clean();
        n_iter += 1;
    }
    return x_i;
}

#endif //SOLUTION_SLAE_METHOD_GAUSS_SEIDEL_H
