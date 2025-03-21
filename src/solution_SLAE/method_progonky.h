#ifndef SOLUTION_SLAE_METHOD_PROGONKY_H
#define SOLUTION_SLAE_METHOD_PROGONKY_H
#include <array>

template<typename T>
struct PQ{
    T p_ = 0;
    T q_ = 0;
};

template<typename T, std::size_t N>
std::array<T, N> method_progonky(const std::array<T, N - 1>& a,
                                 const std::array<T, N>& b,
                                 const std::array<T, N - 1>& c,
                                 const std::array<T, N>& d) {

    std::array<PQ<T>, N-1> pq_arr;
    pq_arr[0] = PQ<T> {-c[0] / b[0], d[0] / b[0]};

    for (std::size_t i = 1; i != N - 1; ++i) {
        T p_i = -c[i] / (a[i-1] * pq_arr[i-1].p_ + b[i]);
        T q_i = (d[i] - a[i-1] * pq_arr[i-1].q_) / (a[i-1] * pq_arr[i-1].p_ + b[i]);
        pq_arr[i] = PQ<T> {p_i, q_i};
    }

    std::array<T, N> solution;
    solution[N-1] = (d[N-1] - a[N-2] * pq_arr[N-2].q_) / (a[N-2] * pq_arr[N-2].p_ + b[N-1]);

    for (std::size_t i = N - 1; i != 0; --i) {
        solution[i - 1] = pq_arr[i-1].p_ * solution[i] + pq_arr[i-1].q_;
    }

    return solution;
}

#endif //SOLUTION_SLAE_METHOD_PROGONKY_H
