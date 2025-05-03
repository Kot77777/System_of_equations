#ifndef PRIMITIVES_ARNOLDI_H
#define PRIMITIVES_ARNOLDI_H
#include "vector_from_vector.h"
#include "primitives/Upper_Triangular_Matrix.h"
#include "primitives/Col_Major_Matrix.h"

template<typename T, std::size_t N>
struct Arnoldi {
    using Rotation = std::pair<T, T>;
    Upper_Triangular_Matrix<T> R{N};
    Vector<Rotation> rotations{N};
    Col_Major_Matrix<T> basis{N};

};

#endif //PRIMITIVES_ARNOLDI_H
