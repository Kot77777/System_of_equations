#ifndef COL_MAJOR_MATRIX_H
#define COL_MAJOR_MATRIX_H
#include <vector>
#include <span>

template<typename T>
class Col_Major_Matrix {
    std::size_t N_{};
    std::vector<T> data_{};

public:
    Col_Major_Matrix() = default;

    Col_Major_Matrix(const std::vector<T> &data, const std::size_t N) : data_(data), N_(N) {
    }

    void set_N(const std::size_t N) {
        N_ = N;
    }

    void insert(const std::size_t i) {
        data_.insert(data_.end(), i, 0);
    }

    void col(const Vector<T> &v, const std::size_t i) {
        auto column = std::span<T>(data_.data() + N_ * i, N_);
        std::copy(v.data_get().begin(), v.data_get().end(), column.begin());
    }

    Vector<T> get_col(const std::size_t i) const {
        auto column = std::span<const T>(data_.begin() + N_ * i, N_);
        return Vector<T>(column);
    }

    T &operator()(const std::size_t i, const std::size_t j) {
        return data_[j * N_ + i];
    }

    T operator()(const std::size_t i, const std::size_t j) const {
        return data_[j * N_ + i];
    }

};

template<typename T>
Vector<T> operator*(const Col_Major_Matrix<T> &matrix, const Vector<T> &vector) {
    Vector<T> vector_res{vector.N()};
    for (std::size_t i = 0; i != vector.N(); ++i) {
        for (std::size_t k = 0; k != vector.N(); ++k) {
            vector_res(i) += matrix(k, i) * vector(k);
        }
    }
    return vector_res;
}

#endif //COL_MAJOR_MATRIX_H
