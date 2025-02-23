#ifndef PRIMITIVES_CSR_MATRIX_H
#define PRIMITIVES_CSR_MATRIX_H
#include <map>
#include <vector>
#include "vector_from_vector.h"

template<typename T>
class CSR_matrix {
    std::size_t N_;
    std::size_t M_;
    std::vector<T> values_;
    std::vector<std::size_t> cols_;
    std::vector<std::size_t> rows_;

public:
    CSR_matrix() = default;

    CSR_matrix(const std::map<std::array<std::size_t, 2>, T> &DOK, const std::size_t N, const std::size_t M) {
        N_ = N;
        M_ = M;
        std::size_t count = 0;
        values_.resize(DOK.size());
        cols_.resize(DOK.size());
        rows_.resize(N + 1);
        for (const auto &element: DOK) {
            values_[count] = element.second;
            cols_[count++] = element.first[1];
            rows_[element.first[0] + 1] += rows_[element.first[0]] * (rows_[element.first[0] + 1] == 0) + 1;
        }
    }

    T operator()(const std::size_t i, const std::size_t j) const {
        for (std::size_t k = rows_[i]; k < rows_[i + 1]; ++k) {
            if (cols_[k] == j) { return values_[k]; }
        }
        return 0;
    }

    const std::size_t rows(const std::size_t i) const {
        return rows_[i];
    }

    const T values(const std::size_t i) const {
        return values_[i];
    }

    const std::size_t cols(const std::size_t i) const {
        return cols_[i];
    }

    const std::size_t N() const { return N_; }

    const std::size_t M() const { return M_; }
};

template<typename T>
Vector<T> operator*(const CSR_matrix<T> &matrix, const Vector<T> &vector) {
    Vector<T> vector_res{matrix.N()};
    const std::size_t N = matrix.N();

    for (std::size_t i = 0; i != N; ++i) {
        const std::size_t k_start = matrix.rows(i);
        const std::size_t k_end = matrix.rows(i+1);
        for (std::size_t k = k_start; k < k_end; ++k) {
            vector_res(i) += matrix.values(k) * vector(matrix.cols(k));
        }
    }
    return vector_res;
}

#endif //PRIMITIVES_CSR_MATRIX_H
