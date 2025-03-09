#include <fstream>
#include <map>
#include <array>
#include <solution_SLAE/method_Gauss_Seidel.h>
#include <solution_SLAE/method_Jacobi.h>
#include <solution_SLAE/method_simple_iteration.h>
#include <solution_SLAE/method_quick_simple_iteration.h>
#include "algorithms/lamda_max.h"
#include "algorithms/get_polynom_roots.h"
#include "algorithms/permutation.h"

int main() {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 28.0;
    DOK[{1, 1}] = 10.0;
    DOK[{2, 2}] = 26.0;
    DOK[{0, 1}] = 12.0;
    DOK[{1, 0}] = 2.0;
    DOK[{2, 0}] = 10.0;
    DOK[{0, 2}] = 2.0;
    DOK[{3, 3}] = 10.0;
    DOK[{3, 0}] = 2.0;
    DOK[{2, 3}] = 8.0;

    const CSR_matrix<double> A{DOK, 4, 4};
    const Vector<double> b{{1., 2., 4., 8.}, 4};
    const double eps = 10e-15;

    std::ofstream data_measure_Jacobi("data_measure_Jacobi.csv");
    data_measure_Jacobi << "n_iterate" << "," << "nevyzka" << '\n';
    Vector<double> res{4};
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
    const double lamd_min = 8.77427527658846;
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
