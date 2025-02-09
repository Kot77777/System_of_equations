#ifndef PRIMITIVES_CSR_MATRIX_H
#define PRIMITIVES_CSR_MATRIX_H
#include <map>
#include "matrix.h"

template<typename T>
class CSR_matrix {
    std::vector<T> values_;
    std::vector<std::size_t> cols_;
    std::vector<std::size_t> rows_;

public:
    CSR_matrix() = default;

    CSR_matrix(const std::map<std::array<std::size_t, 2>, T> &DOK) {
        std::size_t count = 0;
        const std::size_t size = DOK.size();
        values_.reserve(size);
        cols_.reserve(size);
        rows_.resize(DOK.end()->first[0]);
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

    const std::vector<std::size_t> rows() const {
        return rows_;
    }

    const std::vector<std::size_t> values() const {
        return values_;
    }

    const std::vector<std::size_t> cols() const {
        return cols_;
    }
};

template<std::size_t N, std::size_t M, typename T>
Vector<N, T> operator*(const CSR_matrix<T> &matrix, const Vector<M, T> &vector) {
    Vector<N, T> vector_res;
    for (std::size_t i = 0; i != M; ++i) {
        for (std::size_t k = matrix.rows()[i]; k < matrix.rows()[i + 1]; ++k) {
            vector_res[i] += matrix.values()[k] * vector[matrix.cols()[k]];
        }
    }
    return vector_res;
}

#endif //PRIMITIVES_CSR_MATRIX_H
