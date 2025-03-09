#ifndef ALGORITHMS_LAMDA_MAX_H
#define ALGORITHMS_LAMDA_MAX_H

template<typename T>
T lamda_max(const CSR_matrix<T> &A, const Vector<T> &r_0, const std::size_t N_iter) {
    T mu_i = dot(r_0, A * r_0) / dot(r_0, r_0);
    std::size_t n_iter{};
    Vector<T> r_i_next{r_0};
    Vector<T> r_i{r_0.N()};

    while (n_iter < N_iter) {
        r_i = r_i_next;

        r_i_next = A * r_i;
        r_i_next = r_i_next * (1 / (A * r_i).norm());
        mu_i = dot(r_i_next, A * r_i_next) / dot(r_i_next, r_i_next);
        n_iter+=1;
    }
    return mu_i;
}

#endif //ALGORITHMS_LAMDA_MAX_H
