#include <fstream>
#include <chrono>
#include <solution_SLAE/method_Gauss_Seidel.h>
#include <solution_SLAE/method_Jacobi.h>
#include "solution_SLAE/method_symmetriz_Gauss_Seidel.h"
#include "algorithms/generate_matrix_of_elliptic_equation.h"
#include "algorithms/lamda_max_for_accel.h"
#include "algorithms/accel_for_symm_method.h"

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
    const double eps = 1e-14;
    const Vector<double> x_0{5};

    std::ofstream data_measure_Jacobi("data_measure_Jacobi_1.csv");
    std::ofstream time_measure_Jacobi("time_measure_Jacobi_1.csv");
    data_measure_Jacobi << "n_iterate" << "," << "nevyzka" << '\n';
    time_measure_Jacobi << "time" << "," << "nevyzka" << '\n';

    Vector<double> res{5};
    std::size_t count = 1;
    double r = 0;
    std::chrono::microseconds accumulated_duration(0);

    while (!cond_stop(A, b, res, eps)) {
        auto start = std::chrono::high_resolution_clock::now();
        res = method_Jacobi(A, b, res, 1, eps);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        accumulated_duration += duration;
        r = (A * res - b).norm();

        data_measure_Jacobi << count << "," << r << '\n';
        time_measure_Jacobi <<  accumulated_duration.count() << "," << r << '\n';
        count += 1;
    }

    std::ofstream data_measure_Gauss_Seidel("data_measure_Gauss_Seidel_1.csv");
    std::ofstream time_measure_Gauss_Seidel("time_measure_Gauss_Seidel_1.csv");
    data_measure_Gauss_Seidel << "n_iterate" << "," << "nevyzka" << '\n';
    time_measure_Gauss_Seidel << "time" << "," << "nevyzka" << '\n';
    res.clean();
    count = 1;
    std::chrono::microseconds accumulated_duration1(0);

    while (!cond_stop(A, b, res, eps)) {
        auto start = std::chrono::high_resolution_clock::now();
        res = method_Gauss_Seidel(A, b, res, 1, eps);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        accumulated_duration1 += duration;
        r = (A * res - b).norm();

        data_measure_Gauss_Seidel << count << "," << r << '\n';
        time_measure_Gauss_Seidel <<  accumulated_duration1.count() << "," << r << '\n';
        count += 1;
    }

    std::ofstream data_measure_accel_symm_Gauss_Seidel("data_measure_accel_symm_Gauss_Seidel.csv");
    std::ofstream time_measure_accel_symm_Gauss_Seidel("time_measure_accel_symm_Gauss_Seidel.csv");
    data_measure_accel_symm_Gauss_Seidel << "n_iterate" << "," << "nevyzka" << '\n';
    time_measure_accel_symm_Gauss_Seidel << "time" << "," << "nevyzka" << '\n';
    res.clean();
    count = 1;
    std::chrono::microseconds accumulated_duration2(0);

    const auto symm_Gauss_Seidel = [&A, eps](const Vector<double> &b, const Vector<double> &x_0) {
        return method_symmetriz_Gauss_Seidel(A, b, x_0, 1, eps);
    };

    const double rho_symm_GS = lamda_max_for_accel(symm_Gauss_Seidel, b, 1000);

    while (!cond_stop(A, b, res, eps)) {
        auto start = std::chrono::high_resolution_clock::now();
        res = accel(A, res, b, rho_symm_GS, symm_Gauss_Seidel, 1, eps);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        accumulated_duration2 += duration;
        r = (A * res - b).norm();

        data_measure_accel_symm_Gauss_Seidel << count << "," << r << '\n';
        time_measure_accel_symm_Gauss_Seidel <<  accumulated_duration2.count() << "," << r << '\n';
        count += 1;
    }
}