#ifndef SOLUTION_SLAE_METHOD_STEEPEST_GRADIENT_DESCENT_H
#define SOLUTION_SLAE_METHOD_STEEPEST_GRADIENT_DESCENT_H
#include "primitives/CSR_matrix.h"
#include "algorithms/stop_cond.h"

template<typename T>
Vector<T> method_steepest_gradient_descent(const CSR_matrix<T> &A, const Vector<T> &b, const Vector<T> &x_0,
                                           const std::size_t N_iter, const T eps) {
    Vector<T> x_i{x_0};
    Vector<T> r_i = A * x_i - b;
    T alpha_i = dot(r_i, r_i) / dot(r_i, A * r_i);
    std::size_t n_iter{};
    while (!cond_stop(A, b, x_i, eps) && n_iter < N_iter) {
        x_i = x_i - alpha_i * r_i;
        n_iter += 1;
        r_i = A * x_i - b;
        alpha_i = dot(r_i, r_i) / dot(r_i, A * r_i);
    }
    return x_i;
}

#endif //SOLUTION_SLAE_METHOD_STEEPEST_GRADIENT_DESCENT_H
