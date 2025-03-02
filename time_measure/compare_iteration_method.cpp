#include <fstream>
#include <map>
#include <array>
#include <solution_SLAE/method_Gauss_Seidel.h>
#include <solution_SLAE/method_Jacobi.h>
#include <solution_SLAE/method_simple_iteration.h>

int main() {
    std::map<std::array<std::size_t, 2>, double> DOK;

    DOK[{0, 0}] = 10.0;
    DOK[{1, 1}] = 9.0;
    DOK[{1, 2}] = 1.0;
    DOK[{2, 1}] = 2.0;
    DOK[{2, 2}] = 8.0;
    DOK[{2, 9}] = 3.0;
    DOK[{3, 0}] = 1.0;
    DOK[{3, 3}] = 7.0;
    DOK[{4, 2}] = 3.0;
    DOK[{4, 4}] = 12.0;
    DOK[{4, 7}] = 4.0;
    DOK[{5, 1}] = 5.0;
    DOK[{5, 5}] = 11.0;
    DOK[{6, 2}] = 6.0;
    DOK[{6, 6}] = 15.0;
    DOK[{7, 3}] = 1.0;
    DOK[{7, 4}] = 2.0;
    DOK[{7, 7}] = 14.0;
    DOK[{8, 0}] = 4.0;
    DOK[{8, 2}] = 5.0;
    DOK[{8, 3}] = 3.0;
    DOK[{8, 8}] = 13.0;
    DOK[{9, 1}] = 6.0;
    DOK[{9, 4}] = 1.0;
    DOK[{9, 6}] = 2.0;
    DOK[{9, 8}] = 3.0;
    DOK[{9, 9}] = 9.0;

    const CSR_matrix<double> A{DOK, 10, 10};
    const Vector<double> b{{1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}, 10};
    const double eps = 10e-15;

    std::ofstream data_measure_Jacobi("data_measure_Jacobi.csv");
    data_measure_Jacobi << "n_iterate" << "," << "nevyzka" << '\n';
    Vector<double> res{10};
    std::size_t count = 1;

    while (!cond_stop(A, b, res, eps)) {
        res = method_Jacobi(A, b, res, count, eps);
        data_measure_Jacobi << count << "," << (A * res - b).norm() << '\n';
        count += 1;
    }

    std::ofstream data_measure_Gauss_Seidel("data_measure_Gauss_Seidel.csv");
    data_measure_Gauss_Seidel << "n_iterate" << "," << "nevyzka" << '\n';
    res.clean();
    count = 1;

    while (!cond_stop(A, b, res, eps)) {
        res = method_Gauss_Seidel(A, b, res, count, eps);
        data_measure_Gauss_Seidel << count << "," << (A * res - b).norm() << '\n';
        count += 1;
    }

    std::ofstream data_measure_simple_iteration("data_measure_simple_iteration.csv");
    data_measure_simple_iteration << "n_iterate" << "," << "nevyzka" << '\n';
    res.clean();
    count = 1;
    const double t = 2. / 23;

    while (!cond_stop(A, b, res, eps)) {
        res = method_simple_iteration(A, b, res, t, count, eps);
        data_measure_simple_iteration << count << "," << (A * res - b).norm() << '\n';
        count += 1;
    }
}
