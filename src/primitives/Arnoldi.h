#ifndef PRIMITIVES_ARNOLDI_H
#define PRIMITIVES_ARNOLDI_H
#include "vector_from_vector.h"
#include "primitives/Upper_Triangular_Matrix.h"
#include "primitives/Col_Major_Matrix.h"

template<typename T>
struct Arnoldi {
    using Rotation = std::pair<T, T>;

    std::size_t N_;
    Upper_Triangular_Matrix<T> R_;
    Vector<Rotation> rotations_;
    Col_Major_Matrix<T> basis_;
};

#endif //PRIMITIVES_ARNOLDI_H
