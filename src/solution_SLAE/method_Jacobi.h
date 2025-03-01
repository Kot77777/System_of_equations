#ifndef SOLUTION_SLAE_METHOD_JACOBI_H
#define SOLUTION_SLAE_METHOD_JACOBI_H
#include <primitives/CSR_matrix.h>
#include <primitives/vector_from_vector.h>
#include <algorithms/stop_cond.h>

template<typename T>
Vector<T> method_Jacobi(const CSR_matrix<T> &A, const Vector<T> &b, const T eps) {
    Vector<T> x_i{b.N()};
    Vector<T> x_i_next{b.N()};
    T d_ii{};
    while(!cond_stop(A, b, x_i, eps)) {
        for (std::size_t i = 0; i != A.N(); ++i) {
            const std::size_t k_start = A.rows(i);
            const std::size_t k_end = A.rows(i+1);
            for (std::size_t k = k_start; k < k_end; ++k) {
                if (A.cols(k) != i) {x_i_next(i) += A.values(k) * x_i(A.cols(k));}
                else {d_ii = A.values(k);}
            }
            x_i_next(i) = (b(i) - x_i_next(i))/d_ii;
        }
        x_i = x_i_next;
        x_i_next = Vector<T>{b.N()};
    }
    return x_i;
}

#endif //SOLUTION_SLAE_METHOD_JACOBI_H
