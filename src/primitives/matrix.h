#ifndef PRIMITIVES_MATRIX_H
#define PRIMITIVES_MATRIX_H
#include <array>

template<std::size_t N, std::size_t M, typename T>
class Matrix {
    std::array<T, N * M> data_;

public:
    Matrix() = default;

    Matrix(const std::array<T, N * M> &matrix) : data_(matrix) {
    }

    T &operator()(const std::size_t i, const std::size_t j) {
        if (i >= N || j >= M) {
            throw std::out_of_range("Index out of range");
        }
        return data_[i * M + j];
    }

    T operator()(const std::size_t i, const std::size_t j) const {
        if (i >= N || j >= M) {
            throw std::out_of_range("Index out of range");
        }
        return data_[i * M + j];
    }

    T operator()(const std::size_t i) const {
        static_assert(M == 1, "This operation is only valid for vectors.");
        return data_[i];
    }

    T &operator()(const std::size_t i) {
        static_assert(M == 1, "This operation is only valid for vectors.");
        return data_[i];
    }

    operator T() const {
        static_assert(N == 1 && M == 1, "This conversion is only valid for 1x1 matrix.");
        return data_[0];
    }
};

template<std::size_t N, typename T>
using Vector = Matrix<N, 1, T>;

template<std::size_t N, typename T>
using Vector_transp = Matrix<1, N, T>;

template<std::size_t N, std::size_t M, typename T>
Matrix<N, M, T> operator+(const Matrix<N, M, T> &matrix1, const Matrix<N, M, T> &matrix2) {
    Matrix<N, M, T> matrix;
    for (std::size_t i = 0; i != N; ++i) {
        for (std::size_t j = 0; j != M; ++j) {
            matrix(i, j) = matrix1(i, j) + matrix2(i, j);
        }
    }
    return matrix;
}

template<std::size_t N, std::size_t M, typename T>
Matrix<N, M, T> operator-(const Matrix<N, M, T> &matrix1) {
    Matrix<N, M, T> matrix;
    for (std::size_t i = 0; i != N; ++i) {
        for (std::size_t j = 0; j != M; ++j) {
            matrix(i, j) = -matrix1(i, j);
        }
    }
    return matrix;
}

template<std::size_t N, std::size_t M, typename T>
Matrix<N, M, T> operator-(const Matrix<N, M, T> &matrix1, const Matrix<N, M, T> &matrix2) {
    return matrix1 + (-matrix2);
}

template<std::size_t N, std::size_t M, typename T>
Matrix<N, M, T> operator*(const Matrix<N, M, T> &matrix1, const T c) {
    Matrix<N, M, T> matrix;
    for (std::size_t i = 0; i != N; ++i) {
        for (std::size_t j = 0; j != M; ++j) {
            matrix(i, j) = matrix1(i, j) * c;
        }
    }
    return matrix;
}

template<std::size_t N, std::size_t M, typename T>
Vector<N, T> operator*(const Vector<N, T> &matrix1, const T c) {
    Matrix<N, 1, T> matrix;
    for (std::size_t i = 0; i != N; ++i) {
        matrix(i, 0) = matrix1(i, 0) * c;
    }
    return matrix;
}

template<std::size_t N, std::size_t M, typename T>
Matrix<N, M, T> operator*(const T c, const Matrix<N, M, T> &matrix1) {
    return matrix1 * c;
}

template<std::size_t N, std::size_t M, std::size_t K, typename T>
Matrix<N, M, T> operator*(const Matrix<N, K, T> &matrix1, const Matrix<K, M, T> &matrix2) {
    Matrix<N, M, T> matrix;
    for (std::size_t i = 0; i != N; ++i) {
        for (std::size_t j = 0; j != M; ++j) {
            for (std::size_t k = 0; k != K; ++k) {
                T p = matrix1(i, k) * matrix2(k, j);
                matrix(i, j) += p;
            }
        }
    }
    return matrix;
}

template<std::size_t N, std::size_t M, typename T>
Matrix<M, N, T> transposed(const Matrix<N, M, T> &matrix1) {
    Matrix<M, N, T> matrix;
    for (std::size_t i = 0; i != N; ++i) {
        for (std::size_t j = 0; j != M; ++j) {
            matrix(j, i) = matrix1(i, j);
        }
    }
    return matrix;
}

template<std::size_t N, std::size_t M, typename T>
Matrix<N - 1, M - 1, T> minor(const Matrix<N, M, T> &matrix, const std::size_t i_, const std::size_t j_) {
    std::array<T, (N - 1) * (M - 1)> minor;
    std::size_t count = 0;
    for (std::size_t i = 0; i != N; ++i) {
        for (std::size_t j = 0; j != M; ++j) {
            if (i != i_ and j != j_) { minor[count++] = matrix(i, j); }
        }
    }
    Matrix<N - 1, M - 1, T> minor_matrix{minor};
    return minor_matrix;
}

template<typename T>
T determinant(const Matrix<1, 1, T> &matrix) {
    return matrix(0, 0);
}

template<std::size_t N, typename T>
T determinant(const Matrix<N, N, T> &matrix) {
    T determ = 0;
    for (std::size_t i = 0; i != N; ++i) {
        determ += ((i % 2 == 0) ? 1 : -1) * matrix(0, i) * determinant(minor(matrix, 0, i));
    }
    return determ;
}

template<std::size_t N, std::size_t M, typename T>
bool operator==(const Matrix<N, M, T> &matrix1, const Matrix<N, M, T> &matrix2) {
    for (std::size_t i = 0; i != N; ++i) {
        for (std::size_t j = 0; j != M; ++j) {
            if (matrix1(i, j) != matrix2(i, j)) { return false; }
        }
    }
    return true;
}

#endif //PRIMITIVES_MATRIX_H
