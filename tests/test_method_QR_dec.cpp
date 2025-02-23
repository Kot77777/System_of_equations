#include <gtest/gtest.h>
#include "algorithms/QR_dec.h"
#include "solution_SLAE/method_QR_dec.h"

TEST(QR_dec, QR_dec_1) {
    Matrix<double> A{{4, 9, 3, 4, 0, 6, 1, 2, 3, 3, 4, 5, 5, 2, 7}, 5, 3};
    QR_matrix<double> res{QR_dec(A)};
    Matrix<double> B = res.Q_ * res.R_;

    for (std::size_t i = 0; i < 5; ++i) {
        for(std::size_t j = 0; j < 3; ++j) {
            ASSERT_NEAR(A(i, j), B(i, j), 10e-15);
        }
    }
}

TEST(QR_dec, QR_dec_2) {
    Matrix<double> A{{4, 9, 3, 4, 0, 6, 1, 2, 3, 3, 4, 5, 5, 2, 7, 4, 9, 3, 4, 0, 6, 1, 2, 3, 3, 4, 5, 5, 2, 7}, 6, 5};
    QR_matrix<double> res{QR_dec(A)};
    Matrix<double> B = res.Q_ * res.R_;

    for (std::size_t i = 0; i < 6; ++i) {
        for(std::size_t j = 0; j < 5; ++j) {
            ASSERT_NEAR(A(i, j), B(i, j), 10e-15);
        }
    }
}

TEST(QR_dec, QR_dec_3) {
    Matrix<double> A{
        {
            4, 9, 3, 4, 0, 6, 1, 2, 3, 3, 4, 5, 5, 2, 7, 4, 9, 3, 4, 0, 6, 1, 2, 3, 3, 4, 5, 5, 2, 7,
            4, 9, 3, 4, 0, 6, 1, 2, 3, 3, 4, 5, 5, 2, 7, 4, 9, 3, 4, 0, 6, 1, 2, 3, 3, 4, 5, 5, 2, 7
        }, 12, 5};
    QR_matrix<double> res{QR_dec(A)};
    Matrix<double> B = res.Q_ * res.R_;

    for (std::size_t i = 0; i < 6; ++i) {
        for(std::size_t j = 0; j < 5; ++j) {
            ASSERT_NEAR(A(i, j), B(i, j), 10e-14);
        }
    }
}

TEST(QR_dec, method_QR_dec_1) {
    Matrix<double> A{{4, 9, 3, 4, 0, 6, 1, 2, 3}, 3, 3};
    Vector<double> b{{3, 8, 64}, 3};
    Vector<double> res{method_QR_dec(A, b)};
    Vector<double> exp{{-542./13, 119./13, 1136./39}, 3};

    for (std::size_t i = 0; i < 3; ++i) {
        ASSERT_NEAR(res(i), exp(i), 10e-15);
    }
}

TEST(QR_dec, method_QR_dec_2) {
    Matrix<double> A{{25, 8, 34, 34, 34, 53, 65, 44, 45, 90, 58, 24, 29, 54, 39, 64}, 4, 4};
    Vector<double> b{{9, 95, 30, 92}, 4};
    Vector<double> res{method_QR_dec(A, b)};
    Vector<double> exp{{-4517088./878587, 896333. / 878587, 2021274. / 878587, 2643559. / 1757174}, 4};

    for (std::size_t i = 0; i < 3; ++i) {
        ASSERT_NEAR(res(i), exp(i), 10e-12);
    }
}


