#include <gtest/gtest.h>
#include "solution_SLAE/method_progonky.h"

TEST(method_progonky, equal_of_array_1) {
    std::vector<double> expected = {2.0 / 7, 1.0 / 7, 2.0 / 7};
    std::vector<double> result = method_progonky<double>({1, 1}, {3, 3, 3}, {1, 1}, {1, 1, 1});

    for (size_t i = 0; i != 3; ++i) {
        ASSERT_NEAR(result[i], expected[i], 1e-15);
    }
}

TEST(method_progonky, equal_of_array_2) {
    std::vector<double> expected = {5.0 / 26, 3.0 / 13, 6.0 / 13};
    std::vector<double> result = method_progonky<double>({2, 1}, {4, 5, 6}, {1, 1}, {1, 2, 3});

    for (size_t i = 0; i != 3; ++i) {
        ASSERT_NEAR(result[i], expected[i], 1e-15);
    }
}

TEST(method_progonky, equal_of_array_3) {
    std::vector<float> expected = {3.0 / 14, -1.0 / 7, 15.0 / 28};
    std::vector<float> result = method_progonky<float>({5, 1}, {6, 8, 4}, {2, 2}, {1, 1, 2});

    for (size_t i = 0; i != 3; ++i) {
        ASSERT_NEAR(result[i], expected[i], 1e-7);
    }
}

TEST(method_progonky, equal_of_array_4) {
    std::vector<float> expected = {29.0 / 65, 10.0 / 13, 33.0 / 130, 1.0};
    std::vector<float> result = method_progonky<float>({3, 2, 0}, {5, 8, -10, 1}, {1, 2, 6}, {3, 8, 5, 1});

    for (size_t i = 0; i != 4; ++i) {
        ASSERT_NEAR(result[i], expected[i], 1e-7);
    }
}

TEST(method_progonky, equal_of_array_5) {
    std::vector<float> expected = {6.0 / 11, -1.0 / 11};
    std::vector<float> result = method_progonky<float>({5}, {2, 8}, {1}, {1, 2});

    for (size_t i = 0; i != 2; ++i) {
        ASSERT_NEAR(result[i], expected[i], 1e-7);
    }
}

TEST(method_progonky, equal_of_array_6) {
    std::vector<float> expected = {3.};
    std::vector<float> result = method_progonky<float>({}, {2}, {}, {6});

    for (size_t i = 0; i != 1; ++i) {
        ASSERT_NEAR(result[i], expected[i], 1e-7);
    }
}
