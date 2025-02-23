#ifndef ALGORITHMS_QR_DEC_H
#define ALGORITHMS_QR_DEC_H
#include "primitives/matrix_from_vector.h"

template<typename T>
struct QR_matrix {
    Matrix<T> Q_;
    Matrix<T> R_;
};

template<typename T>
QR_matrix<T> QR_dec(const Matrix<T> &A) {
    Matrix<T> R{A};
    Matrix<T> Q{eye<T>(R.N())};
    for (std::size_t i = 0; i < R.M(); ++i) {
        Vector<T> x_R{R.N() - i};
        Vector<T> x_Q{R.N() - i};
        for (std::size_t j = i; j < R.N(); ++j) {
            x_R(j - i) = R(j, i);
            x_Q(j - i) = Q(0 , j);
        }
        Vector<T> v{x_R};
        v(0) -= x_R.norm();
        T v_dot = dot(v, v);
        if(v_dot != 0){
            for (std::size_t k = i; k < R.M() - 1; ++k) {
                T vx_R_dot = dot(v, x_R);
                for (std::size_t z = i; z < R.N(); ++z) {
                    R(z, k) = R(z, k) - 2 * (vx_R_dot / v_dot) * v(z - i); //применяем преобразование только к значениям ниже или на диагонали
                    x_R(z - i) = R(z, k + 1); //сразу собираем в x_R следующий столбец
                }
            }
            T vx_R_dot = dot(v, x_R);

            for (std::size_t z = i; z < R.N(); ++z) { // вынесена отдельно так как уже не нужно собирать следующий столбец(его нет)
                R(z, R.M() - 1) = R(z, R.M() - 1) - 2 * (vx_R_dot / v_dot) * v(z - i);
            }

            //аналогично для Q только преобразовываем все строки(k с нуля) размером N-1
            for (std::size_t k = 0; k < Q.N() - 1; ++k) {
                T vx_Q_dot = dot(v, x_Q);
                for (std::size_t z = i; z < Q.N(); ++z) {
                    Q(k, z) = Q(k, z) - 2 * (vx_Q_dot / v_dot) * v(z - i);
                    x_Q(z - i) = Q(k + 1, z);
                }
            }
            T vx_Q_dot = dot(v, x_Q);
            for (std::size_t z = i; z < Q.N(); ++z) {
                Q(Q.N() - 1, z) = Q(Q.N() - 1, z) - 2 * (vx_Q_dot / v_dot) * v(z - i);
            }
        }
    }
    QR_matrix<T> res{Q, R};
    return res;
}

#endif //ALGORITHMS_QR_DEC_H
