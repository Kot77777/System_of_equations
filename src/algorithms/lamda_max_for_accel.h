#ifndef ALGORITHMS_LAMDA_MAX_FOR_ACCEL_H
#define ALGORITHMS_LAMDA_MAX_FOR_ACCEL_H
#include "primitives/vector_from_vector.h"

template<typename T, typename Method>
T lamda_max_for_accel(const Method &method, const Vector<T> &r_0, const std::size_t N_iter) {

    const Vector<T> b{r_0.N()};
    Vector<T> Ar = method(b, r_0);

    T mu_i = dot(r_0, Ar) / dot(r_0, r_0);
    std::size_t n_iter{};
    Vector<T> r_i_next{r_0};
    Vector<T> r_i{r_0.N()};

    while (n_iter < N_iter) {
        r_i = r_i_next;

        Ar = method(b, r_i);
        r_i_next = Ar;
        r_i_next = r_i_next * (1 / (Ar).norm());
        mu_i = dot(r_i_next, method(b, r_i_next)) / dot(r_i_next, r_i_next);
        n_iter+=1;
    }
    return mu_i;
}

#endif //ALGORITHMS_LAMDA_MAX_FOR_ACCEL_H
