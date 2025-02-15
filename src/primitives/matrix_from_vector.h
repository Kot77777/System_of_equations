#ifndef PRIMITIVES_MATRIX_FROM_VECTOR_H
#define PRIMITIVES_MATRIX_FROM_VECTOR_H
#include "vector_from_vector.h"

template<typename T>
class Matrix {
    std::size_t N_{};
    std::size_t M_{};
    std::vector<T> data_{};

public:
    Matrix() = default;

    Matrix(const std::vector<T> &matrix, const std::size_t N, const std::size_t M) : data_(matrix), N_(N), M_(M) {
    }

    Matrix(const std::size_t N, const std::size_t M) {
        data_.resize(N * M);
        N_ = N;
        M_ = M;
    }

    T &operator()(const std::size_t i, const std::size_t j) {
        if (i >= N_ || j >= M_) {
            throw std::out_of_range("Index out of range");
        }
        return data_[i * M_ + j];
    }

    T operator()(const std::size_t i, const std::size_t j) const {
        if (i >= N_ || j >= M_) {
            throw std::out_of_range("Index out of range");
        }
        return data_[i * M_ + j];
    }

    const std::size_t N() const { return N_; }

    const std::size_t M() const { return M_; }

    operator T() const {
        if (N_ != 1 || M_ != 1) {
            throw std::out_of_range("This conversion is only valid for 1x1 matrix.");
        }
        return data_[0];
    }
};

template<typename T>
Matrix<T> operator*(const Matrix<T> &matrix1, const Matrix<T> &matrix2) {
    if (matrix1.M() != matrix2.N()) {
        throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
    }
    Matrix<T> matrix{matrix1.N(), matrix2.M()};
    for (std::size_t i = 0; i != matrix.N(); ++i) {
        for (std::size_t j = 0; j != matrix.M(); ++j) {
            for (std::size_t k = 0; k != matrix1.M(); ++k) {
                matrix(i, j) += matrix1(i, k) * matrix2(k, j);
            }
        }
    }
    return matrix;
}

template<typename T>
Vector<T> operator*(const Matrix<T> &matrix, const Vector<T> &vector) {
    if (matrix.M() != vector.N()) {
        throw std::invalid_argument("Incompatible matrix and vector dimensions for multiplication.");
    }
    Vector<T> vector_res{matrix.N()};
    for (std::size_t i = 0; i != matrix.N(); ++i) {
        for (std::size_t k = 0; k != vector.N(); ++k) {
            vector_res(i) += matrix(i, k) * vector(k);
        }
    }
    return vector_res;
}

#endif //PRIMITIVES_MATRIX_FROM_VECTOR_H
