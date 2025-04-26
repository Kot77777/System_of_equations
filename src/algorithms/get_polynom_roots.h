#ifndef ALGORITHMS_GET_POLYNOM_ROOTS_H
#define ALGORITHMS_GET_POLYNOM_ROOTS_H

template<typename T>
Vector<T> get_polynom_roots(const std::size_t n, const T k, const T b) {
    Vector<T> t{n};
    Vector<T> res{n};
    const T cos_pi_2n = cos(M_PI / 2 / n);
    //const T sin_pi_2n = sin(M_PI / 2 / n);
    const T cos_pi_n = cos(M_PI / n);
    const T sin_pi_n = sin(M_PI / n);
    t(0) = cos_pi_2n;
    for (std::size_t i = 1; i < n; ++i) {
        t(i) = t(i - 1) * cos_pi_n - sqrt(1 - t(i - 1) * t(i - 1)) * sin_pi_n;
    }
    for (std::size_t i = 0; i < n; ++i) {res(i) = 1 / (k*t(i) + b);}
    return res;
}

#endif //ALGORITHMS_GET_POLYNOM_ROOTS_H
