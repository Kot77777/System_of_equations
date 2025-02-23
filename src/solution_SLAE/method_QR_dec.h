#ifndef SOLUTION_SLAE_METHOD_QR_DEC_H
#define SOLUTION_SLAE_METHOD_QR_DEC_H
#include "algorithms/QR_dec.h"

template<typename T>
Vector<T> method_QR_dec(const Matrix<T> &A, const Vector<T> &b) {
    Vector<T> x{b.N()};
    const QR_matrix<T> QR{QR_dec(A)};
    const Matrix<T> R{QR.R_};
    const Matrix<T> Q(QR.Q_);
    const Vector<T> RHS = transposed(Q) * b;

    x(b.N() - 1) = RHS(b.N() - 1) / R(b.N() - 1, b.N() - 1);
    for (int i = b.N() - 2; i >= 0; --i) {
        x(i) = RHS(i);
        for (std::size_t j = i + 1; j < b.N(); ++j) {
            x(i) -= R(i, j) * x(j);
        }
        x(i) /= R(i, i);
    }
    return x;
}

#endif //SOLUTION_SLAE_METHOD_QR_DEC_H
