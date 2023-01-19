/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "clustering.h"
#include "matrix.h"
#include <gtest/gtest.h>
#include <random>

using namespace ldt;

TEST(Distance_t, euclidean) {
  Ti n = 9;
  Ti count = 16;
  auto data = Matrix<Tv>(
      new Tv[400]{5.24, 1.29, 1.37, 5.6,  0.76, 4.12, 1.96, 7.49, 7.31, 9.83,
                  5.81, 7.96, 0.18, 4.09, 1.93, 7.99, 4.03, 1.08, 1.95, 0.35,
                  9.62, 7.43, 3.19, 2.99, 1.83, 9.5,  0.69, 3.01, 7.65, 1.76,
                  1.81, 6.85, 8.22, 2.43, 2.22, 5.5,  2.54, 1.15, 1.06, 0.89,
                  2.19, 6.94, 0.08, 6,    9.29, 1.76, 5.73, 6.84, 2.82, 3.5,
                  4.07, 1.92, 9.16, 6.3,  3.45, 9.63, 4.19, 6.49, 8.23, 8.87,
                  8.68, 1.66, 9.82, 9.14, 5.02, 0.06, 0.59, 3.78, 4.64, 5.79,
                  3.45, 1.65, 4.15, 1.68, 3.03, 6.09, 1.27, 1.71, 5.11, 5.66,
                  1.33, 0.51, 8.01, 0.61, 9.78, 5.31, 4.28, 5.68, 4.9,  3.17,
                  3.6,  5.21, 4.41, 1.54, 8.35, 5.09, 5.18, 2.66, 5.5,  3.62,
                  1.41, 5.6,  6.94, 1.63, 6.5,  3.43, 1.89, 0.79, 2.49, 9.17,
                  1.15, 6.7,  4.77, 9.56, 0.81, 9.49, 3.38, 9.7,  8.04, 0.34,
                  0.41, 3.05, 0.89, 4.66, 3.98, 3.62, 9.95, 2.71, 6.66, 7.9,
                  2.37, 5.06, 5.05, 7.44, 3.9,  1.28, 2.87, 2.89, 3.61, 5.39,
                  9.37, 5.72, 4.28, 3.31},
      n, count);

  auto model = Distance<false, DistanceMethod::kEuclidean>(n, count);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(data, storage, work);

  ASSERT_NEAR(11.755096766934757, model.Result.Get(2, 0), 1e-14);

  auto model1 = Distance<false, DistanceMethod::kManhattan>(n, count);
  auto storage1 = new Tv[model1.StorageSize];
  auto work1 = new Tv[model1.WorkSize];
  model1.Calculate(data, storage1, work1);

  ASSERT_NEAR(26.629999999999995, model1.Result.Get(2, 0), 1e-14);

  auto model2 = Distance<false, DistanceMethod::kMaximum>(n, count);
  auto storage2 = new Tv[model2.StorageSize];
  auto work2 = new Tv[model2.WorkSize];
  model2.Calculate(data, storage2, work2);

  ASSERT_NEAR(8.2500000000000000, model2.Result.Get(2, 0), 1e-14);

  delete[] storage;
  delete[] work;
  delete[] storage1;
  delete[] work1;
  delete[] storage2;
  delete[] work2;
}

TEST(Distance_t, euclidean_nan) {
  Ti n = 10;
  Ti count = 16;
  auto data = Matrix<Tv>(
      new Tv[400]{5.24, NAN, 1.29, 1.37, 5.6,  0.76, 4.12, 1.96, 7.49, 7.31,
                  9.83, NAN, 5.81, 7.96, 0.18, 4.09, 1.93, 7.99, 4.03, 1.08,
                  1.95, NAN, 0.35, 9.62, 7.43, 3.19, 2.99, 1.83, 9.5,  0.69,
                  3.01, NAN, 7.65, 1.76, 1.81, 6.85, 8.22, 2.43, 2.22, 5.5,
                  2.54, NAN, 1.15, 1.06, 0.89, 2.19, 6.94, 0.08, 6,    9.29,
                  1.76, NAN, 5.73, 6.84, 2.82, 3.5,  4.07, 1.92, 9.16, 6.3,
                  3.45, NAN, 9.63, 4.19, 6.49, 8.23, 8.87, 8.68, 1.66, 9.82,
                  9.14, NAN, 5.02, 0.06, 0.59, 3.78, 4.64, 5.79, 3.45, 1.65,
                  4.15, NAN, 1.68, 3.03, 6.09, 1.27, 1.71, 5.11, 5.66, 1.33,
                  0.51, NAN, 8.01, 0.61, 9.78, 5.31, 4.28, 5.68, 4.9,  3.17,
                  3.6,  NAN, 5.21, 4.41, 1.54, 8.35, 5.09, 5.18, 2.66, 5.5,
                  3.62, NAN, 1.41, 5.6,  6.94, 1.63, 6.5,  3.43, 1.89, 0.79,
                  2.49, NAN, 9.17, 1.15, 6.7,  4.77, 9.56, 0.81, 9.49, 3.38,
                  9.7,  NAN, 8.04, 0.34, 0.41, 3.05, 0.89, 4.66, 3.98, 3.62,
                  9.95, NAN, 2.71, 6.66, 7.9,  2.37, 5.06, 5.05, 7.44, 3.9,
                  1.28, NAN, 2.87, 2.89, 3.61, 5.39, 9.37, 5.72, 4.28, 3.31},
      n, count);

  auto model = Distance<true, DistanceMethod::kEuclidean>(n, count);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(data, storage, work);

  ASSERT_NEAR(11.755096766934757, model.Result.Get(2, 0), 1e-14);

  auto model1 = Distance<true, DistanceMethod::kManhattan>(n, count);
  auto storage1 = new Tv[model1.StorageSize];
  auto work1 = new Tv[model1.WorkSize];
  model1.Calculate(data, storage1, work1);

  ASSERT_NEAR(26.629999999999995, model1.Result.Get(2, 0), 1e-14);

  auto model2 = Distance<true, DistanceMethod::kMaximum>(n, count);
  auto storage2 = new Tv[model2.StorageSize];
  auto work2 = new Tv[model2.WorkSize];
  model2.Calculate(data, storage2, work2);

  ASSERT_NEAR(8.2500000000000000, model2.Result.Get(2, 0), 1e-14);

  delete[] storage;
  delete[] work;
  delete[] storage1;
  delete[] work1;
  delete[] storage2;
  delete[] work2;
}

TEST(Distance_t, correlation_t) {
  Ti n = 9;
  Ti count = 16;
  auto data = Matrix<Tv>(
      new Tv[400]{5.24, 1.29, 1.37, 5.6,  0.76, 4.12, 1.96, 7.49, 7.31, 9.83,
                  5.81, 7.96, 0.18, 4.09, 1.93, 7.99, 4.03, 1.08, 1.95, 0.35,
                  9.62, 7.43, 3.19, 2.99, 1.83, 9.5,  0.69, 3.01, 7.65, 1.76,
                  1.81, 6.85, 8.22, 2.43, 2.22, 5.5,  2.54, 1.15, 1.06, 0.89,
                  2.19, 6.94, 0.08, 6,    9.29, 1.76, 5.73, 6.84, 2.82, 3.5,
                  4.07, 1.92, 9.16, 6.3,  3.45, 9.63, 4.19, 6.49, 8.23, 8.87,
                  8.68, 1.66, 9.82, 9.14, 5.02, 0.06, 0.59, 3.78, 4.64, 5.79,
                  3.45, 1.65, 4.15, 1.68, 3.03, 6.09, 1.27, 1.71, 5.11, 5.66,
                  1.33, 0.51, 8.01, 0.61, 9.78, 5.31, 4.28, 5.68, 4.9,  3.17,
                  3.6,  5.21, 4.41, 1.54, 8.35, 5.09, 5.18, 2.66, 5.5,  3.62,
                  1.41, 5.6,  6.94, 1.63, 6.5,  3.43, 1.89, 0.79, 2.49, 9.17,
                  1.15, 6.7,  4.77, 9.56, 0.81, 9.49, 3.38, 9.7,  8.04, 0.34,
                  0.41, 3.05, 0.89, 4.66, 3.98, 3.62, 9.95, 2.71, 6.66, 7.9,
                  2.37, 5.06, 5.05, 7.44, 3.9,  1.28, 2.87, 2.89, 3.61, 5.39,
                  9.37, 5.72, 4.28, 3.31},
      n, count);

  auto model = Distance<false, DistanceMethod::kCorrelation,
                        CorrelationMethod::kPearson>(n, count);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(data, storage, work);

  ASSERT_NEAR(0.84997521712902702, model.Result.Get(1, 0), 1e-14);
  ASSERT_NEAR(0.73301473408209417, model.Result.Get(2, 1), 1e-14);

  auto model1 = Distance<false, DistanceMethod::kAbsCorrelation,
                         CorrelationMethod::kPearson>(n, count);
  auto storage1 = new Tv[model1.StorageSize];
  auto work1 = new Tv[model1.WorkSize];
  model1.Calculate(data, storage1, work1);

  ASSERT_NEAR(0.89557243413052168, model1.Result.Get(1, 0), 1e-14);
  ASSERT_NEAR(0.99721195159138809, model1.Result.Get(2, 1), 1e-14);

  delete[] storage;
  delete[] work;
  delete[] storage1;
  delete[] work1;
}

TEST(Distance_t, correlation_nan_t) {
  Ti n = 10;
  Ti count = 16;
  auto data = Matrix<Tv>(
      new Tv[400]{5.24, NAN, 1.29, 1.37, 5.6,  0.76, 4.12, 1.96, 7.49, 7.31,
                  9.83, NAN, 5.81, 7.96, 0.18, 4.09, 1.93, 7.99, 4.03, 1.08,
                  1.95, NAN, 0.35, 9.62, 7.43, 3.19, 2.99, 1.83, 9.5,  0.69,
                  3.01, NAN, 7.65, 1.76, 1.81, 6.85, 8.22, 2.43, 2.22, 5.5,
                  2.54, NAN, 1.15, 1.06, 0.89, 2.19, 6.94, 0.08, 6,    9.29,
                  1.76, NAN, 5.73, 6.84, 2.82, 3.5,  4.07, 1.92, 9.16, 6.3,
                  3.45, NAN, 9.63, 4.19, 6.49, 8.23, 8.87, 8.68, 1.66, 9.82,
                  9.14, NAN, 5.02, 0.06, 0.59, 3.78, 4.64, 5.79, 3.45, 1.65,
                  4.15, NAN, 1.68, 3.03, 6.09, 1.27, 1.71, 5.11, 5.66, 1.33,
                  0.51, NAN, 8.01, 0.61, 9.78, 5.31, 4.28, 5.68, 4.9,  3.17,
                  3.6,  NAN, 5.21, 4.41, 1.54, 8.35, 5.09, 5.18, 2.66, 5.5,
                  3.62, NAN, 1.41, 5.6,  6.94, 1.63, 6.5,  3.43, 1.89, 0.79,
                  2.49, NAN, 9.17, 1.15, 6.7,  4.77, 9.56, 0.81, 9.49, 3.38,
                  9.7,  NAN, 8.04, 0.34, 0.41, 3.05, 0.89, 4.66, 3.98, 3.62,
                  9.95, NAN, 2.71, 6.66, 7.9,  2.37, 5.06, 5.05, 7.44, 3.9,
                  1.28, NAN, 2.87, 2.89, 3.61, 5.39, 9.37, 5.72, 4.28, 3.31},
      n, count);

  auto model =
      Distance<true, DistanceMethod::kCorrelation, CorrelationMethod::kPearson>(
          n, count);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(data, storage, work);

  ASSERT_NEAR(0.84997521712902702, model.Result.Get(1, 0), 1e-14);
  ASSERT_NEAR(0.73301473408209417, model.Result.Get(2, 1), 1e-14);

  auto model1 = Distance<true, DistanceMethod::kAbsCorrelation,
                         CorrelationMethod::kPearson>(n, count);
  auto storage1 = new Tv[model1.StorageSize];
  auto work1 = new Tv[model1.WorkSize];
  model1.Calculate(data, storage1, work1);

  ASSERT_NEAR(0.89557243413052168, model1.Result.Get(1, 0), 1e-14);
  ASSERT_NEAR(0.99721195159138809, model1.Result.Get(2, 1), 1e-14);

  delete[] storage;
  delete[] work;
  delete[] storage1;
  delete[] work1;
}
