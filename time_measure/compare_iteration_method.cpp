#include <fstream>
#include <map>
#include <array>
#include <solution_SLAE/method_Gauss_Seidel.h>
#include <solution_SLAE/method_Jacobi.h>
#include <solution_SLAE/method_simple_iteration.h>
#include "algorithms/lamda_max.h"
#include "algorithms/get_polynom_roots.h"
#include "algorithms/permutation.h"

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
    const double eps = 10e-15;

    std::ofstream data_measure_Jacobi("data_measure_Jacobi.csv");
    data_measure_Jacobi << "n_iterate" << "," << "nevyzka" << '\n';
    Vector<double> res{5};
    std::size_t count = 1;

    while (!cond_stop(A, b, res, eps)) {
        res = method_Jacobi(A, b, res, 1, eps);
        data_measure_Jacobi << count << "," << (A * res - b).norm() << '\n';
        count += 1;
    }

    std::ofstream data_measure_Gauss_Seidel("data_measure_Gauss_Seidel.csv");
    data_measure_Gauss_Seidel << "n_iterate" << "," << "nevyzka" << '\n';
    res.clean();
    count = 1;

    while (!cond_stop(A, b, res, eps)) {
        res = method_Gauss_Seidel(A, b, res, 1, eps);
        data_measure_Gauss_Seidel << count << "," << (A * res - b).norm() << '\n';
        count += 1;
    }

    std::ofstream data_measure_simple_iteration("data_measure_simple_iteration.csv");
    data_measure_simple_iteration << "n_iterate" << "," << "nevyzka" << '\n';
    res.clean();
    count = 1;
    const double lamd_min = 1.02311881410475;
    const double lamd_max = lamda_max(A, b, 1000);
    const double t = 2. / (lamd_max + lamd_min);

    while (!cond_stop(A, b, res, eps)) {
        res = method_simple_iteration(A, b, res, t, 1, eps);
        data_measure_simple_iteration << count << "," << (A * res - b).norm() << '\n';
        count += 1;
    }

    std::ofstream data_measure_quick_simple_iteration("data_measure_quick_simple_iteration.csv");
    data_measure_quick_simple_iteration << "n_iterate" << "," << "nevyzka" << '\n';
    res.clean();
    const std::size_t n = 256;

    const Vector<double> t_{get_polynom_roots<double>(n, (lamd_max - lamd_min) / 2, (lamd_max + lamd_min) / 2)};
    const Vector<std::size_t> perm = permutation(Vector<size_t>{{0, 1}, 2}, n);

    Vector<double> x_i{res};
    std::size_t n_iter{};
    while (!cond_stop(A, b, x_i, eps) && n_iter < 1000) {
        for (std::size_t i = 0; i < n; ++i) {
            x_i = x_i - t_(perm(i)) * (A * x_i - b);
            n_iter += 1;
            data_measure_quick_simple_iteration << n_iter << "," << (A * x_i - b).norm() << '\n';
            if(cond_stop(A, b, x_i, eps)){n_iter+=1000; break;}
        }
    }
}
