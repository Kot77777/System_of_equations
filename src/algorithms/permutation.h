#ifndef ALGORITHMS_PERMUTATION_H
#define ALGORITHMS_PERMUTATION_H

Vector<std::size_t> permutation(const Vector<std::size_t> &c, const std::size_t n) {
    if (c.N() == n) return c;
    Vector<std::size_t> d{c.N() * 2};
    for (std::size_t i = 0; i < d.N(); i += 2) {
        d(i) = c(i / 2);
        d(i + 1) = d.N() - d(i) - 1;
    }
    return permutation(d, n);
}


#endif //ALGORITHMS_PERMUTATION_H
