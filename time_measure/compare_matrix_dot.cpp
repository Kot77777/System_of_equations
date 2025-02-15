#include <fstream>
#include <array>
#include <map>
#include <chrono>
#include "primitives/CSR_matrix.h"
#include "primitives/matrix_from_vector.h"

template<typename T>
std::map<std::array<std::size_t, 2>, T> get_random_CSR_matrix(const std::size_t size, const T densyty) {
    std::map<std::array<std::size_t, 2>, T> DOK;
    for (std::size_t i = 0; i != size; ++i) {
        for (std::size_t j = 0; j != size; ++j) {
            const T rand_num = static_cast<T>(std::rand()) / RAND_MAX;
            if (rand_num > densyty) { DOK[{i, j}] = rand_num; }
        }
    }
    return DOK;
}

int main() {
    std::ofstream data_measure_CSR_matrix_dot("data_measure_CSR_matrix_dot.csv");
    data_measure_CSR_matrix_dot << "size" << "," << "densyty" << "," << "time" << '\n';

    std::srand(static_cast<unsigned int>(std::time(0)));
    for (int i = 0; i < 10; i += 2) {
        for (std::size_t j = 1; j <= 10; ++j) {
            CSR_matrix<double> A{get_random_CSR_matrix<double>(100 * j, 0.1 * i), 100 * j, 100 * j};
            Vector<double> c{100 * j};

            auto start = std::chrono::high_resolution_clock::now();
            Vector<double> res = A * c;
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

            data_measure_CSR_matrix_dot << 100 * j << "," << 0.1 * i << "," << duration.count() << '\n';
        }
    }

    std::ofstream data_measure_Dense_matrix_dot("data_measure_Dense_matrix_dot.csv");
    data_measure_Dense_matrix_dot << "size" << "," << "time" << '\n';

    for(std::size_t j = 1; j <= 10; ++j) {
        std::vector<double> m((10000 * j * j), 1);
        Matrix<double> matrix{m,100 * j, 100 * j};
        Vector<double> v{100 * j};

        auto start = std::chrono::high_resolution_clock::now();
        Vector<double> mult = matrix * v;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        data_measure_Dense_matrix_dot << 100 * j << "," << duration.count() << '\n';
    }

}
