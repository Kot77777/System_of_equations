#include <fstream>
#include <algorithms/generate_matrix_of_elliptic_equation.h>
#include <solution_SLAE/method_ConjugateGradients.h>
#include <solution_SLAE/method_Gauss_Seidel.h>
#include <solution_SLAE/method_Jacobi.h>
#include <solution_SLAE/method_simple_iteration.h>
#include <solution_SLAE/method_SOR.h>
#include <solution_SLAE/method_symmetriz_Gauss_Seidel.h>
#include <solution_SLAE/method_steepest_gradient_descent.h>

#include "algorithms/lamda_max_for_SOR.h"
#include "algorithms/get_polynom_roots.h"
#include "algorithms/permutation.h"
#include "algorithms/lamda_max_for_accel.h"
#include "algorithms/accel_for_symm_method.h"

int main() {
    const std::size_t N = 50;
    const CSR_matrix<double> A = generate_matrix<double>(N);
    const Vector<double> b{std::vector<double>(N * N, 1), N * N};
    const double eps = 1e-11;

    Vector<double> res{N * N};
    double r0 = (A * res - b).norm();

    std::ofstream data_measure_Jacobi("data_measure_Jacobi.csv");
    data_measure_Jacobi << "n_iterate" << "," << "nevyzka" << '\n';
    data_measure_Jacobi << 0 << "," << r0 << '\n';
    std::size_t count = 1;

    while (!cond_stop(A, b, res, eps)) {
        res = method_Jacobi(A, b, res, 1, eps);
        data_measure_Jacobi << count << "," << (A * res - b).norm() << '\n';
        count += 1;
    }

    std::ofstream data_measure_Gauss_Seidel("data_measure_Gauss_Seidel.csv");
    data_measure_Gauss_Seidel << "n_iterate" << "," << "nevyzka" << '\n';
    data_measure_Gauss_Seidel << 0 << "," << r0 << '\n';
    res.clean();
    count = 1;

    while (!cond_stop(A, b, res, eps)) {
        res = method_Gauss_Seidel(A, b, res, 1, eps);
        data_measure_Gauss_Seidel << count << "," << (A * res - b).norm() << '\n';
        count += 1;
    }

    std::ofstream data_measure_accel_symm_Gauss_Seidel("data_measure_accel_symm_Gauss_Seidel.csv");
    data_measure_accel_symm_Gauss_Seidel << "n_iterate" << "," << "nevyzka" << '\n';
    data_measure_accel_symm_Gauss_Seidel << 0 << "," << r0 << '\n';
    res.clean();
    count = 2;

    const auto symm_Gauss_Seidel = [&A, eps](const Vector<double> &b, const Vector<double> &x_0) {
        return method_symmetriz_Gauss_Seidel(A, b, x_0, 1, eps);
    };

    const double rho_symm_GS = lamda_max_for_accel(symm_Gauss_Seidel, b, 1000);

    while (!cond_stop(A, b, res, eps)) {
        res = accel(A, res, b, rho_symm_GS, symm_Gauss_Seidel, 1, eps);

        data_measure_accel_symm_Gauss_Seidel << count << "," << (A * res - b).norm() << '\n';
        count += 1;
    }

    std::ofstream data_measure_simple_iteration("data_measure_simple_iteration.csv");
    data_measure_simple_iteration << "n_iterate" << "," << "nevyzka" << '\n';
    data_measure_simple_iteration << 0 << "," << r0 << '\n';
    res.clean();
    count = 1;
    const double lamd_min = 8 * sin(M_PI / (2 * (N + 1))) * sin(M_PI / (2 * (N + 1)));
    const double lamd_max = 8 * sin(N * M_PI / (2 * (N + 1))) * sin(N * M_PI / (2 * (N + 1)));
    const double t = 2. / (lamd_max + lamd_min);

    while (!cond_stop(A, b, res, eps)) {
        res = method_simple_iteration(A, b, res, t, 1, eps);
        data_measure_simple_iteration << count << "," << (A * res - b).norm() << '\n';
        count += 1;
    }

    std::ofstream data_measure_quick_simple_iteration("data_measure_quick_simple_iteration.csv");
    data_measure_quick_simple_iteration << "n_iterate" << "," << "nevyzka" << '\n';
    data_measure_quick_simple_iteration << 0 << "," << r0 << '\n';
    res.clean();
    const std::size_t n = 256;

    const Vector<double> t_{get_polynom_roots<double>(n, (lamd_max - lamd_min) / 2, (lamd_max + lamd_min) / 2)};
    const Vector<std::size_t> perm = permutation(Vector<size_t>{{0, 1}, 2}, n);

    Vector<double> x_i{res};
    std::size_t n_iter{};
    while (!cond_stop(A, b, x_i, eps)) {
        for (std::size_t i = 0; i < n; ++i) {
            x_i = x_i - t_(perm(i)) * (A * x_i - b);
            n_iter += 1;
            data_measure_quick_simple_iteration << n_iter << "," << (A * x_i - b).norm() << '\n';
            if (cond_stop(A, b, x_i, eps)) {
                break;
            }
        }
    }

    std::ofstream data_measure_SOR("data_measure_SOR.csv");
    data_measure_SOR << "n_iterate" << "," << "nevyzka" << '\n';
    data_measure_SOR << 0 << "," << r0 << '\n';
    res.clean();
    count = 1;
    const double w = 1.5;

    while (!cond_stop(A, b, res, eps)) {
        res = method_SOR(A, b, res, 1, eps, w);
        data_measure_SOR << count << "," << (A * res - b).norm() << '\n';
        count += 1;
    }

    std::ofstream data_measure_SOR_fastest("data_measure_SOR_fastest.csv");
    data_measure_SOR_fastest << "n_iterate" << "," << "nevyzka" << '\n';
    data_measure_SOR_fastest << 0 << "," << r0 << '\n';
    res.clean();
    count = 1;
    const double mu = lamda_max_for_SOR(A, b, 1000);
    const double mu_p = mu / (1.0 + sqrt(1.0 - mu * mu));
    const double w_fastest = 1 + mu_p * mu_p;

    while (!cond_stop(A, b, res, eps)) {
        res = method_SOR(A, b, res, 1, eps, w_fastest);
        data_measure_SOR_fastest << count << "," << (A * res - b).norm() << '\n';
        count += 1;
    }

    std::ofstream data_measure_symmetriz_Gauss_Seidel("data_measure_symmetriz_Gauss_Seidel.csv");
    data_measure_symmetriz_Gauss_Seidel << "n_iterate" << "," << "nevyzka" << '\n';
    data_measure_symmetriz_Gauss_Seidel << 0 << "," << r0 << '\n';
    res.clean();
    count = 2;

    while (!cond_stop(A, b, res, eps)) {
        res = method_symmetriz_Gauss_Seidel(A, b, res, 1, eps);
        data_measure_symmetriz_Gauss_Seidel << count << "," << (A * res - b).norm() << '\n';
        count += 2;
    }

    std::ofstream data_measure_steepest_gradient_descent("data_measure_steepest_gradient_descent.csv");
    data_measure_steepest_gradient_descent << "n_iterate" << "," << "nevyzka" << '\n';
    data_measure_steepest_gradient_descent << 0 << "," << r0 << '\n';
    res.clean();
    count = 1;

    while (!cond_stop(A, b, res, eps)) {
        res = method_steepest_gradient_descent(A, b, res, 1, eps);
        data_measure_steepest_gradient_descent << count << "," << (A * res - b).norm() << '\n';
        count += 1;
    }

    std::ofstream data_measure_CG("data_measure_CG.csv");
    data_measure_CG << "n_iterate" << "," << "nevyzka" << '\n';
    data_measure_CG << 0 << "," << r0 << '\n';
    res.clean();
    count = 1;

    Vector<double> r_i = A * res - b;
    Vector<double> r_i_next{r_i};
    Vector<double> d_i{r_i};
    x_i.clean();
    double alpha_i{}, betta_i{};
    n_iter = 0;

    while (r_i.norm() >= eps && n_iter < 1000) {
        alpha_i = dot(r_i, r_i) / dot(d_i, A * d_i);
        x_i = x_i - alpha_i * d_i;
        r_i_next = A * x_i - b;

        betta_i = dot(r_i_next, r_i_next) / dot(r_i, r_i);
        d_i = r_i_next + betta_i * d_i;
        r_i = r_i_next;

        n_iter += 1;

        data_measure_CG << n_iter << "," << r_i.norm() << '\n';
    }
}
