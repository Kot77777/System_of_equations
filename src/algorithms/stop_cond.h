#ifndef ALGORITHMS_STOP_COND_H
#define ALGORITHMS_STOP_COND_H
#include <primitives/matrix_from_vector.h>

template<typename T, typename T_matrix>
bool cond_stop(const T_matrix &A, const Vector<T> &b, const Vector<T> &x_i, const T eps) {
    if ((A * x_i - b).norm() < eps) {
        return true;
    }
    return false;
}

#endif //ALGORITHMS_STOP_COND_H
