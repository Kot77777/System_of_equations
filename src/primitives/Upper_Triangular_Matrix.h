#ifndef UPPER_TRIANGULAR_MATRIX_H
#define UPPER_TRIANGULAR_MATRIX_H
#include <vector>

template<typename T>
class Upper_Triangular_Matrix {
    std::size_t N_{};
    std::vector<T> data_{};

public:
    Upper_Triangular_Matrix() = default;

    Upper_Triangular_Matrix(const std::vector<T> &data, const std::size_t N) : data_(data), N_(N) {
    }

    Upper_Triangular_Matrix(const std::size_t N) {
        data_.resize(N * (N + 1) / 2);
        N_ = N;
    }

    T &operator()(const std::size_t i, const std::size_t j) {
        return data_[i * N_ + j - i * (i+1) / 2];
    }

    T operator()(const std::size_t i, const std::size_t j) const {
        return data_[i * N_ + j - i * (i+1) / 2];
    }
};

#endif //UPPER_TRIANGULAR_MATRIX_H
