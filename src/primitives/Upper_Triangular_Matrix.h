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

    void resize(const std::size_t i) {
        data_.insert(data_.end(), i, 0);
        N_ = i;
    }

    T &operator()(const std::size_t i, const std::size_t j) {
        return data_[j * N_ + i - j * (j+1) / 2];
    }

    T operator()(const std::size_t i, const std::size_t j) const {
        return data_[j * N_ + i - j * (j+1) / 2];
    }
};

template<typename T>
Vector<T> solution_SLAE(const Upper_Triangular_Matrix<T> &R, const Vector<T> &z) {
    Vector<T> res{z.N()};
    for (int i = z.N() - 1; i >= 0; --i) {
        for (int j = z.N() - 1; j > i; --j) {
            res(i) -= res(j) * R(i, j);
        }
        res(i) = (res(i) + z(i)) / R(i, i);
    }
    return res;
}

#endif //UPPER_TRIANGULAR_MATRIX_H
