#ifndef SOLUTION_SLAE_METHOD_SIMPLE_ITERATION_H
#define SOLUTION_SLAE_METHOD_SIMPLE_ITERATION_H
#include "primitives/vector_from_vector.h"
#include "primitives/CSR_matrix.h"
#include "algorithms/stop_cond.h"

template<typename T>
Vector<T> method_simple_iteration(const CSR_matrix<T> &A, const Vector<T> &b, const Vector<T> &x_0, const T t,
                                  const std::size_t N_iter,
                                  const T eps) {
    Vector<T> x_i{x_0};
    std::size_t n_iter{};
    while (!cond_stop(A, b, x_i, eps) && n_iter < N_iter) {
        x_i = x_i - t * (A * x_i - b);
        n_iter += 1;
    }
    return x_i;
}

#endif //SOLUTION_SLAE_METHOD_SIMPLE_ITERATION_H
