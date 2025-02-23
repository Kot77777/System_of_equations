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
Matrix<T> operator+(const Matrix<T> &matrix1, const Matrix<T> &matrix2) {
    if (matrix1.N() != matrix2.N() || matrix1.M() != matrix2.M()) {
        throw std::invalid_argument("Incompatible matrix dimensions for sum.");
    }
    Matrix<T> result{matrix1.N(), matrix1.M()};
    for (std::size_t i = 0; i != matrix1.N(); ++i) {
        for (std::size_t j = 0; j != matrix1.M(); ++j) {
            result(i, j) = matrix1(i, j) + matrix2(i, j);
        }
    }
    return result;
}

template<typename T>
Matrix<T> operator-(const Matrix<T> &matrix) {
    Matrix<T> result{matrix.N(), matrix.M()};
    for (std::size_t i = 0; i != result.N(); ++i) {
        for (std::size_t j = 0; j != result.M(); ++j) {
            result(i, j) = -matrix(i, j);
        }
    }
    return result;
}

template<typename T>
Matrix<T> operator-(const Matrix<T> &matrix1, const Matrix<T> &matrix2) {
    if (matrix1.N() != matrix2.N() || matrix1.M() != matrix2.M()) {
        throw std::invalid_argument("Incompatible matrix dimensions for sum.");
    }
    return matrix1 + (-matrix2);
}

template<typename T>
Matrix<T> operator*(const Matrix<T> &matrix1, const T c) {
    Matrix<T> result{matrix1.N(), matrix1.M()};
    for (std::size_t i = 0; i != matrix1.N(); ++i) {
        for (std::size_t j = 0; j != matrix1.M(); ++j) {
            result(i, j) = matrix1(i, j) * c;
        }
    }
    return result;
}

template<typename T>
Matrix<T> operator*(const T c, const Matrix<T> &matrix1) {
    return matrix1 * c;
}

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
Matrix<T> diada_for_vector(const Vector<T> &vector) {
    Matrix<T> result{vector.N(), vector.N()};
    for (std::size_t i = 0; i != vector.N(); ++i) {
        for (std::size_t j = i; j != vector.N(); ++j) {
            result(i, j) = vector(i) * vector(j);
            result(j, i) = vector(i) * vector(j);
        }
    }
    return result;
}

template<typename T>
Matrix<T> transposed(const Matrix<T> &matrix) {
    Matrix<T> result{matrix.M(), matrix.N()};
    for (std::size_t i = 0; i != result.N(); ++i) {
        for (std::size_t j = 0; j != result.M(); ++j) {
            result(j, i) = matrix(i, j);
        }
    }
    return result;
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

template<typename T>
bool operator==(const Matrix<T> &matrix1, const Matrix<T> &matrix2) {
    if (matrix1.N() != matrix2.N() || matrix1.M() != matrix2.M()) {
        throw std::invalid_argument("Incompatible matrix dimensions for equal.");
    }
    for (std::size_t i = 0; i != matrix1.N(); ++i) {
        for (std::size_t j = 0; j != matrix2.M(); ++j) {
            if (matrix1(i, j) != matrix2(j, i)) { return false; }
        }
    }
    return true;
}

template<typename T>
Matrix<T> eye(const std::size_t n) {
    Matrix<T> result{n, n};
    for (std::size_t i = 0; i < result.N(); ++i) {
        result(i, i) = 1;
    }
    return result;
}

#endif //PRIMITIVES_MATRIX_FROM_VECTOR_H
