#ifndef ALGORITHMS_LAMDA_MAX_FOR_SOR_H
#define ALGORITHMS_LAMDA_MAX_FOR_SOR_H
#include "algorithms/lamda_max.h"

template<typename T>
T lamda_max_for_SOR(const CSR_matrix<T> &A, const Vector<T> &r_0, const std::size_t N_iter) {
    Matrix<T> res{A.N(), A.N()};
    T d{};
    for (std::size_t i = 0; i < A.N(); ++i) {
        d = 1 / A(i, i);
        for(std::size_t j = 0; j < A.N(); ++j) {
            res(i, j) = -d * A(i, j);
        }
        res(i, i) += 1;
    }
    return lamda_max(res, r_0, N_iter);
}

#endif //ALGORITHMS_LAMDA_MAX_FOR_SOR_H
