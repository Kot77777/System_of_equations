#ifndef SOLUTION_SLAE_METHOD_CGS_H
#define SOLUTION_SLAE_METHOD_CGS_H

#include "primitives/vector_from_vector.h"
#include "primitives/CSR_matrix.h"

template<typename T>
Vector<T> method_CGS(const CSR_matrix<T> &A, const Vector<T> &b, const Vector<T> &x_0,
                      const std::size_t N_iter, const T eps) {
    Vector<T> x_i{x_0};
    Vector<T> r_0 = A * x_0 - b;
    Vector<T> r_i{r_0}, r_i_next{r_0};
    Vector<T> d_i{r_i}, u_i{r_i}, q_i{x_0.N()};
    T alpha_i{}, betta_i{};
    std::size_t n_iter{};

    while (r_i.norm() >= eps && n_iter < N_iter) {
        alpha_i = dot(r_0, r_i) / dot(r_0, A * d_i);
        q_i = u_i - alpha_i * (A * d_i);
        x_i = x_i - alpha_i * (u_i + q_i);
        r_i_next = A * x_i - b;

        betta_i = dot(r_0, r_i_next) / dot(r_0, r_i);
        u_i = r_i_next + betta_i * q_i;
        d_i = u_i + betta_i * (q_i + betta_i * d_i);

        r_i = r_i_next;

        n_iter += 1;
    }

    return x_i;
}

#endif //SOLUTION_SLAE_METHOD_CGS_H
