#ifndef SOLUTION_SLAE_METHOD_CONJUGATE_GRADIENT_H
#define SOLUTION_SLAE_METHOD_CONJUGATE_GRADIENT_H
#include "primitives/vector_from_vector.h"
#include "primitives/CSR_matrix.h"

template<typename T>
Vector<T> method_ConjugateGradients(const CSR_matrix<T> &A, const Vector<T> &b, const Vector<T> &x_0,
                                    const std::size_t N_iter, const T eps) {
    Vector<T> r_i = A * x_0 - b;
    Vector<T> r_i_next{r_i};
    Vector<T> d_i{r_i};
    Vector<T> x_i{x_0};
    T alpha_i, betta_i{};
    std::size_t n_iter{};

    while (r_i.norm() >= eps && n_iter < N_iter) {
        alpha_i = dot(r_i, r_i) / dot(d_i, A * d_i);
        x_i = x_i - alpha_i * d_i;
        r_i_next = A * x_i - b;

        betta_i = dot(r_i_next, r_i_next) / dot(r_i, r_i);
        d_i = r_i_next + betta_i * d_i;
        r_i = r_i_next;

        n_iter += 1;
    }

    return x_i;
}

#endif //SOLUTION_SLAE_METHOD_CONJUGATE_GRADIENT_H
