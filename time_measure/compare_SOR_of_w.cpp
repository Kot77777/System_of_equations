#include <map>
#include <array>
#include <fstream>

#include "solution_SLAE/method_SOR.h"

int main() {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 6.0;
    DOK[{1, 1}] = 5.0;
    DOK[{2, 2}] = 4.0;
    DOK[{3, 3}] = 3.0;
    DOK[{4, 4}] = 4.0;
    DOK[{0, 1}] = 2.0;
    DOK[{1, 0}] = 2.0;
    DOK[{0, 2}] = 1.0;
    DOK[{2, 0}] = 1.0;
    DOK[{1, 2}] = 1.0;
    DOK[{2, 1}] = 1.0;
    DOK[{1, 3}] = 1.0;
    DOK[{3, 1}] = 1.0;
    DOK[{2, 4}] = 1.0;
    DOK[{4, 2}] = 1.0;
    DOK[{3, 4}] = 2.0;
    DOK[{4, 3}] = 2.0;

    const CSR_matrix<double> A{DOK, 5, 5};
    const Vector<double> b{{1., 2., 4., 8., 10.}, 5};
    Vector<double> x_0{5};
    const double eps = 10e-15;
    const double w_start = 0.01;

    std::ofstream data_measure_SOR_of_w("data_measure_SOR_of_w");
    data_measure_SOR_of_w << "w" << "," << "N_iter"<< '\n';
    for (int i = 8; i < 199; ++i) {
        data_measure_SOR_of_w << w_start * i << "," << method_SOR(A, b, x_0, 2000, eps, w_start * i) << '\n';
    }

}