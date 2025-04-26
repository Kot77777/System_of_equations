#include <fstream>
#include <algorithms/generate_matrix_of_elliptic_equation.h>

int main() {
    const std::size_t N = 10;
    const CSR_matrix<double> A = generate_matrix<double>(N);
    const Vector<double> b{std::vector<double>(N * N, 1), N * N};
    constexpr double eps = 1e-11; {
        std::ofstream data_measure_BiCG("data_measure_BiCG.csv");
        data_measure_BiCG << "n_iterate" << "," << "nevyzka" << '\n';
        Vector<double> x_0{N * N};

        Vector<double> x_i{x_0};
        Vector<double> r_i = A * x_0 - b;
        Vector<double> r_i_c{r_i};
        Vector<double> r_i_next{r_i}, r_i_c_next{r_i};
        Vector<double> d_i{r_i}, d_i_c{r_i};
        double alpha_i{}, betta_i{};
        std::size_t n_iter{};

        while (r_i.norm() >= eps) {
            alpha_i = dot(r_i_c, r_i) / dot(d_i_c, A * d_i);
            r_i_next = r_i - alpha_i * (A * d_i);
            r_i_c_next = r_i_c - alpha_i * transpose_multiply(A, d_i_c);

            x_i = x_i - alpha_i * d_i;

            betta_i = dot(r_i_c_next, r_i_next) / dot(r_i_c, r_i);
            d_i = r_i_next + betta_i * d_i;
            d_i_c = r_i_c_next + betta_i * d_i_c;

            r_i = r_i_next;
            r_i_c = r_i_c_next;

            n_iter += 1;
            data_measure_BiCG << n_iter << "," << r_i.norm() << '\n';
        }
    } {
        std::ofstream data_measure_CGS("data_measure_CGS.csv");
        data_measure_CGS << "n_iterate" << "," << "nevyzka" << '\n';
        Vector<double> x_0{N * N};

        Vector<double> x_i{x_0};
        Vector<double> r_0 = A * x_0 - b;
        Vector<double> r_i{r_0}, r_i_next{r_0};
        Vector<double> d_i{r_i}, u_i{r_i}, q_i{x_0.N()};
        double alpha_i{}, betta_i{};
        std::size_t n_iter{};

        while (r_i.norm() >= eps) {
            alpha_i = dot(r_0, r_i) / dot(r_0, A * d_i);
            q_i = u_i - alpha_i * (A * d_i);
            x_i = x_i - alpha_i * (u_i + q_i);
            r_i_next = A * x_i - b;

            betta_i = dot(r_0, r_i_next) / dot(r_0, r_i);
            u_i = r_i_next + betta_i * q_i;
            d_i = u_i + betta_i * (q_i + betta_i * d_i);

            r_i = r_i_next;

            n_iter += 1;
            data_measure_CGS << n_iter << "," << r_i.norm() << '\n';
        }
    }
}
