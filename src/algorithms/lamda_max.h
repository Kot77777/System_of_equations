#ifndef ALGORITHMS_LAMDA_MAX_H
#define ALGORITHMS_LAMDA_MAX_H

template<typename T>
T lamda_max(const CSR_matrix<T> &A, const Vector<T> &r_0, const T eps) {
    T mu_i_next = dot(r_0, A * r_0) / dot(r_0, r_0);
    T mu_i{};
    Vector<T> r_i_next{r_0};
    Vector<T> r_i{r_0.N()};

    while (abs(mu_i - mu_i_next) > eps) {
        mu_i = mu_i_next;
        r_i = r_i_next;

        r_i_next = A * r_i;
        r_i_next = r_i_next * (1 / (A * r_i).norm());
        mu_i_next = dot(r_i_next, A * r_i_next) / dot(r_i_next, r_i_next);
    }
    return mu_i_next;
}

#endif //ALGORITHMS_LAMDA_MAX_H
