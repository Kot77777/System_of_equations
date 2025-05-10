#ifndef SOLUTION_SLAE_METHOD_GMRES_H
#define SOLUTION_SLAE_METHOD_GMRES_H
#include <primitives/vector_from_vector.h>
#include "primitives/CSR_matrix.h"
#include "primitives/Arnoldi.h"

template<typename T>
Vector<T> method_GMRES(const CSR_matrix<T> &A, const Vector<T> &b, const Vector<T> &x_0,
                       const std::size_t N_iter, const T eps) {
    const Vector<T> r0 = A * x_0 - b;
    Vector<T> t{b.N()};
    Vector<T> x_i{x_0};
    Vector<T> e;
    T h_1{}, h_2{1.}, p{}, cos_{}, sin_{}, rho{r0.norm()}, gamma{r0.norm()};
    std::size_t j{};

    Arnoldi<T> arnoldi;
    Vector<T> v_i = r0 * (1.0 / r0.norm());
    arnoldi.basis_.set_N(b.N());

    while (abs(gamma) > eps && j < N_iter) {
        arnoldi.basis_.resize();
        v_i = v_i * (1 / h_2);
        arnoldi.basis_.col(v_i, j);

        e.push_back(gamma);
        t = A * v_i;

        arnoldi.R_.resize(j + 1);
        for (std::size_t i = 0; i < j + 1; ++i) {
            arnoldi.R_(i, j) = dot(arnoldi.basis_.get_col(i), t);
            t = t - arnoldi.R_(i, j) * arnoldi.basis_.get_col(i);
        }

        for (std::size_t i = 0; i < j; ++i) {
            p = arnoldi.rotations_(i).first * arnoldi.R_(i, j) - arnoldi.rotations_(i).second * arnoldi.R_(i + 1, j);
            arnoldi.R_(i + 1, j) = arnoldi.rotations_(i).second * arnoldi.R_(i, j) + arnoldi.rotations_(i).first *
                                   arnoldi.R_(i + 1, j);
            arnoldi.R_(i, j) = p;
        }

        h_2 = t.norm();
        h_1 = arnoldi.R_(j, j);
        cos_ = h_1 / sqrt(h_1 * h_1 + h_2 * h_2);
        sin_ = -h_2 / sqrt(h_1 * h_1 + h_2 * h_2);
        arnoldi.rotations_.push_back(std::pair<T, T>(cos_, sin_));

        arnoldi.R_(j, j) = arnoldi.rotations_(j).first * arnoldi.R_(j, j) - arnoldi.rotations_(j).second * h_2;
        gamma = e(j) * arnoldi.rotations_(j).second;
        e(j) = arnoldi.rotations_(j).first * e(j);

        v_i = t;

        j += 1;
    }

    const Vector<T> y = solution_SLAE(arnoldi.R_, e);
    x_i = x_0 - arnoldi.basis_ * y;
    return x_i;
}

#endif //SOLUTION_SLAE_METHOD_GMRES_H
