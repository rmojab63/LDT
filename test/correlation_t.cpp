/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the LGPL 3.0.
 See accompanying file LICENSE
*/

#include "correlation.h"
#include "matrix.h"
#include <gtest/gtest.h>
#include <random>

using namespace ldt;

TEST(Correlation_T, cov) {

  Ti m = 4;
  Ti n = 3;
  auto mat =
      Matrix<Tv>(new Tv[m * n]{1, 3, 5, 2, 4, 6, 7, 18, 3, 4, 5, 7}, m, n);

  auto model = Correlation<false, CorrelationType::kCovariance,
                           CorrelationMethod::kPearson>(m, n, true);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(mat, work, storage, true);

  ASSERT_NEAR(2.9166666666666665, model.Result.Get(0, 0), 1e-10);
  ASSERT_NEAR(-1.083333333333333, model.Result.Get(0, 1), 1e-10);
  ASSERT_NEAR(0.5833333333333334, model.Result.Get(0, 2), 1e-10);
  ASSERT_NEAR(10.2500000000000000, model.Result.Get(1, 2), 1e-10);

  delete[] storage;
  delete[] work;
  delete[] mat.Data;
}

TEST(Correlation_T, cov_nan) {

  Ti m = 6;
  Ti n = 3;
  auto mat = Matrix<Tv>(new Tv[m * n]{NAN, 1, 3, 5, 9, 2, 0, 4, 6, 7, NAN, 18,
                                      NAN, 3, 4, 5, NAN, 7},
                        m, n);

  auto model = Correlation<true, CorrelationType::kCovariance,
                           CorrelationMethod::kPearson>(m, n, true);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(mat, work, storage, true);

  ASSERT_NEAR(10.0, model.Result.Get(0, 0), 1e-10);
  ASSERT_NEAR(-1.083333333333333, model.Result.Get(0, 1), 1e-10);
  ASSERT_NEAR(0.5833333333333334, model.Result.Get(0, 2), 1e-10);
  ASSERT_NEAR(10.2500000000000000, model.Result.Get(1, 2), 1e-10);

  delete[] storage;
  delete[] work;
  delete[] mat.Data;
}

TEST(Correlation_T, corr) {

  Ti m = 4;
  Ti n = 2;
  auto mat = Matrix<Tv>(new Tv[m * n]{1, 3, 5, 2, 4, 6, 7, 18}, m, n);

  auto model = Correlation<false, CorrelationType::kCorrelation,
                           CorrelationMethod::kPearson>(m, n, true);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(mat, work, storage);

  ASSERT_EQ(true,
            model.Result.Equals(
                *new Matrix<Tv>(new std::vector<Tv>{1.0, -0.1008236754628326,
                                                    -0.1008236754628326, 1.0},
                                2, 2),
                1e-10));

  delete[] storage;
  delete[] work;
  delete[] mat.Data;
}

TEST(Correlation_T, corr_nan) {

  Ti m = 6;
  Ti n = 2;
  auto mat =
      Matrix<Tv>(new Tv[m * n]{1, 0, 3, 5, 2, NAN, 4, NAN, 6, 7, 18, 2}, m, n);

  auto model = Correlation<true, CorrelationType::kCorrelation,
                           CorrelationMethod::kPearson>(m, n, true);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(mat, work, storage);

  ASSERT_EQ(true,
            model.Result.Equals(
                *new Matrix<Tv>(new std::vector<Tv>{1.0, -0.1008236754628326,
                                                    -0.1008236754628326, 1.0},
                                2, 2),
                1e-10));

  delete[] storage;
  delete[] work;
  delete[] mat.Data;
}

TEST(Correlation_T, corr_spearman) {

  Ti m = 22;
  Ti n = 7;
  auto mat = Matrix<Tv>(
      new Tv[m * n]{
          32.446, 44.145, 17.062, 65.818, 76.19,  40.408, 78.131, 26.695,
          21.992, 68.033, 98.872, 61.154, 71.842, 66.922, 31.142, 58.429,
          45.123, 80.99,  26.345, 50.096, 36.478, 29.377, 27.141, 65.037,
          72.621, 63.391, 68.125, 60.369, 76.384, 5.449,  99.204, 6.87,
          2.514,  98.799, 95.082, 6.048,  59.287, 48.889, 21.464, 54.972,
          14.997, 27.161, 92.522, 65.383, 51.852, 71.011, 55.434, 89.082,
          59.556, 29.524, 21.193, 2.684,  35.457, 69.849, 76.352, 49.455,
          18.762, 8.492,  95.032, 39.042, 32.517, 13.667, 91.408, 23.432,
          56.526, 33.531, 10.67,  72.891, 11.796, 31.202, 96.893, 2.552,
          82.001, 87.786, 96.292, 93.249, 11.688, 19.522, 37.55,  55.967,
          97.026, 14.017, 19.869, 60.988, 91.525, 33.192, 50.666, 97.465,
          58.493, 17.033, 76.138, 3.432,  58.561, 69.172, 56.453, 46.325,
          63.116, 84.577, 12.265, 77.277, 9.141,  69.192, 65.464, 29.827,
          8.261,  26.696, 94.1,   64.958, 68.924, 97.838, 91.389, 76.779,
          56.448, 14.524, 33.549, 39.059, 94.886, 98.52,  80.476, 2.754,
          93.605, 17.733, 37.658, 97.567, 2.705,  74.385, 59.03,  10.732,
          82.043, 92.891, 69.384, 86.848, 40.02,  62.295, 18.609, 61.597,
          22.438, 67.702, 83.393, 96.283, 64.895, 34.39,  42.212, 52.377,
          24.745, 42.534, 64.688, 7.392,  82.462, 22.022, 68.858, 55.901,
          98.156, 96.029},
      m, n);

  auto model = Correlation<false, CorrelationType::kCorrelation,
                           CorrelationMethod::kSpearman>(m, n, true);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(mat, work, storage);

  ASSERT_NEAR(model.Result.Get(1, 0), -0.107848673, 1e-7);
  ASSERT_NEAR(model.Result.Get(4, 1), 0.0750988142, 1e-7);

  delete[] storage;
  delete[] work;
  delete[] mat.Data;
}

TEST(Correlation_T, corr_spearman_nan) {

  Ti m = 22;
  Ti n = 7;
  auto mat = Matrix<Tv>(
      new Tv[m * n]{
          32.446, 44.145, 17.062, 65.818, 76.19,  40.408, 78.131, 26.695,
          21.992, 68.033, 98.872, 61.154, 71.842, 66.922, 31.142, 58.429,
          45.123, 80.99,  26.345, 50.096, 36.478, 29.377, 27.141, 65.037,
          72.621, 63.391, 68.125, 60.369, 76.384, 5.449,  99.204, 6.87,
          2.514,  98.799, 95.082, 6.048,  59.287, 48.889, 21.464, 54.972,
          14.997, 27.161, 92.522, 65.383, 51.852, 71.011, 55.434, 89.082,
          59.556, 29.524, 21.193, 2.684,  35.457, 69.849, 76.352, 49.455,
          18.762, 8.492,  95.032, 39.042, 32.517, 13.667, 91.408, 23.432,
          56.526, 33.531, 10.67,  72.891, 11.796, 31.202, 96.893, 2.552,
          82.001, 87.786, 96.292, 93.249, 11.688, 19.522, 37.55,  55.967,
          97.026, 14.017, 19.869, 60.988, 91.525, 33.192, 50.666, 97.465,
          58.493, 17.033, NAN,    3.432,  58.561, 69.172, 56.453, 46.325,
          63.116, 84.577, 12.265, 77.277, NAN,    69.192, 65.464, 29.827,
          8.261,  26.696, 94.1,   64.958, 68.924, 97.838, 91.389, 76.779,
          56.448, 14.524, 33.549, 39.059, 94.886, 98.52,  80.476, 2.754,
          93.605, 17.733, 37.658, 97.567, 2.705,  74.385, 59.03,  10.732,
          82.043, 92.891, 69.384, 86.848, 40.02,  62.295, 18.609, 61.597,
          22.438, 67.702, 83.393, 96.283, 64.895, 34.39,  NAN,    52.377,
          24.745, 42.534, 64.688, 7.392,  82.462, 22.022, 68.858, 55.901,
          98.156, 96.029},
      m, n);

  auto model = Correlation<true, CorrelationType::kCorrelation,
                           CorrelationMethod::kSpearman>(m, n, true);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(mat, work, storage);

  ASSERT_NEAR(model.Result.Get(1, 0), -0.1078486731, 1e-5);
  ASSERT_NEAR(model.Result.Get(0, 2), -0.1451158, 1e-5);

  delete[] storage;
  delete[] work;
  delete[] mat.Data;
}
