#ifndef ALGORITHMS_ACCEL_FOR_SYMM_METHOD_H
#define ALGORITHMS_ACCEL_FOR_SYMM_METHOD_H
#include <utility>
#include "primitives/CSR_matrix.h"
#include <iostream>

template<typename T, typename Method>
Vector<T> accel(const CSR_matrix<T> &A, const Vector<T> &x_0, const Vector<T> &b, const T rho, const Method &method,
                const std::size_t N_iter, const T eps) {
    std::size_t n_iter{};
    Vector<T> y_0{x_0};
    Vector<T> y_1 = method(b, y_0);
    T w = 2 / (2 - rho * rho);

    while (!cond_stop(A, b, y_1, eps) && n_iter < N_iter) {
        y_0 = std::exchange(y_1, w * (method(b, y_1) - y_0) + y_0);
        w = 1 / (1 - rho * rho * w / 4);
        n_iter += 1;
    }

    return y_1;
}

#endif //ALGORITHMS_ACCEL_FOR_SYMM_METHOD_H
