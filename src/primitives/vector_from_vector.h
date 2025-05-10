#ifndef PRIMITIVES_VECTOR_FROM_VECTOR_H
#define PRIMITIVES_VECTOR_FROM_VECTOR_H
#include <cmath>
#include <vector>
#include <span>

template<typename T>
class Vector {
    std::size_t N_{};
    std::vector<T> data_{};

public:
    Vector() = default;

    Vector(const std::vector<T> &vector, const std::size_t N) : data_(vector), N_(N) {
    }

    Vector(const std::span<const T> &sp) : N_(sp.size()), data_(sp.begin(), sp.end()) {
    }

    Vector(const std::size_t N) {
        data_.resize(N);
        N_ = N;
    }

    void push_back(const T &i) {
        data_.push_back(i);
        N_ += 1;
    }

    T &operator()(const std::size_t i) {
        if (i >= N_) {
            throw std::out_of_range("Index out of range");
        }
        return data_[i];
    }

    T operator()(const std::size_t i) const {
        if (i >= N_) {
            throw std::out_of_range("Index out of range");
        }
        return data_[i];
    }

    const std::size_t N() const { return N_; }

    const T norm() const {
        T res{};
        for (auto &i: data_) {
            res += i * i;
        }
        return std::sqrt(res);
    }

    const std::vector<T> &data_get() const { return data_; }

    void clean() {
        data_.clear();
        data_.resize(N_);
    }

    operator T() const {
        if (N_ != 1) {
            throw std::out_of_range("This conversion is only valid for vector with 1 element.");
        }
        return data_[0];
    }
};

template<typename T>
T dot(const Vector<T> &vector1, const Vector<T> &vector2) {
    if (vector1.N() != vector2.N()) {
        throw std::invalid_argument("Incompatible vector dimensions for dot.");
    }
    T result{};
    for (std::size_t i = 0; i != vector1.N(); ++i) {
        result += vector1(i) * vector2(i);
    }
    return result;
}

template<typename T>
Vector<T> operator+(const Vector<T> &vector1, const Vector<T> &vector2) {
    if (vector1.N() != vector2.N()) {
        throw std::invalid_argument("Incompatible vectors dimensions for sum.");
    }
    Vector<T> vector_res{vector1.N()};
    for (std::size_t i = 0; i != vector_res.N(); ++i) {
        vector_res(i) = vector1(i) + vector2(i);
    }
    return vector_res;
}

template<typename T>
Vector<T> operator-(const Vector<T> &vector) {
    Vector<T> vector_res{vector.N()};
    for (std::size_t i = 0; i != vector_res.N(); ++i) {
        vector_res(i) = -vector(i);
    }
    return vector_res;
}

template<typename T>
Vector<T> operator-(const Vector<T> &vector1, const Vector<T> &vector2) {
    if (vector1.N() != vector2.N()) {
        throw std::invalid_argument("Incompatible vectors dimensions for sum.");
    }
    return vector1 + (-vector2);
}

template<typename T>
Vector<T> operator*(const Vector<T> &vector, const T c) {
    Vector<T> vector_res{vector.N()};
    for (std::size_t i = 0; i != vector_res.N(); ++i) {
        vector_res(i) = c * vector(i);
    }
    return vector_res;
}

template<typename T>
Vector<T> operator*(const T c, const Vector<T> &vector) {
    return vector * c;
}

template<typename T>
bool operator==(const Vector<T> &vector1, const Vector<T> &vector2) {
    if (vector1.N() != vector2.N()) {
        throw std::invalid_argument("Incompatible vector dimensions for equal.");
    }
    for (std::size_t i = 0; i != vector1.N(); ++i) {
        if (vector1(i) != vector2(i)) { return false; }
    }
    return true;
}

#endif //PRIMITIVES_VECTOR_FROM_VECTOR_H
