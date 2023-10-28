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

TEST(Cluster_t, single_t) {
  Ti n = 16;
  auto dist = MatrixSym<false>(
      new Tv[120]{
          7, 10, 8, 8, 7,  8, 8,  8,  10, 10, 7, 3, 8, 8, 4, 8, 10, 7, 6, 7,
          8, 5,  7, 5, 8,  4, 10, 7,  7,  7,  4, 6, 5, 6, 8, 8, 9,  6, 8, 7,
          7, 10, 7, 6, 8,  6, 7,  8,  10, 8,  6, 5, 9, 6, 5, 7, 6,  7, 8, 8,
          6, 6,  7, 5, 8,  8, 8,  3,  5,  6,  6, 6, 6, 7, 7, 8, 8,  9, 7, 7,
          6, 8,  8, 8, 8,  8, 7,  8,  6,  6,  8, 7, 8, 8, 8, 6, 7,  8, 8, 7,
          9, 8,  8, 7, 10, 9, 8,  10, 8,  10, 7, 8, 9, 9, 6, 6, 5,  9, 8, 8,
      },
      n);

  auto model = HCluster<HClusterLinkage::kSingle>(n);
  model.Calculate(dist);

  /*Assert::IsTrue(model.height->Equals(Matrix<Tv>(new Tv[] {6, 6, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 8, 8, 8}, (Ti)15)));

  Assert::IsTrue(model.merge->Equals(Matrix<Ti, Ti>(new Ti[] {
-1,-14,-16,-10,-2,-9,-12,-7,-8,-11,-15,-6,-4,-5,-13,
-3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
}, (Ti)15, (Ti)2)));*/

  std::vector<std::unique_ptr<std::vector<Ti>>> map;
  for (int i = 0; i < 4; i++)
    map.push_back(std::make_unique<std::vector<Ti>>());
  model.Group(map);

  ASSERT_EQ(map.at(0)->size(), (Ti)12);
  ASSERT_EQ(map.at(1)->at(0), (Ti)3);
  ASSERT_EQ(map.at(1)->at(1), (Ti)13);
  ASSERT_EQ(map.at(2)->at(0), (Ti)11);
  ASSERT_EQ(map.at(3)->at(0), (Ti)7);

  Matrix<Tv> heights = Matrix<Tv>(new Tv[n - 1], n - 1);
  Matrix<Ti> merge = Matrix<Ti>(new Ti[2 * (n - 1)], n - 1, (Ti)2);
  std::vector<Ti> order = std::vector<Ti>();
  model.MergeR(merge, heights, order);

  ASSERT_NEAR(3, heights.Data[0], 1e-10);
  ASSERT_NEAR(4, heights.Data[2], 1e-10);
  ASSERT_EQ(14, merge.Get(14, 0));
  ASSERT_EQ(12, merge.Get(14, 1));
}

TEST(Cluster_t, average_t) {
  Ti n = 16;
  auto dist = MatrixSym<false>(
      new Tv[120]{
          5.79, 6.68, 4.72, 8.45, 7.07, 7.16, 7.33, 6.03, 5.47, 8.14, 6.05,
          5.86, 6.46, 7.02, 6.32, 8.75, 4.31, 9.27, 6.35, 4.28, 9.22, 9.15,
          8.54, 8.67, 6.67, 6.86, 6.51, 7.84, 6.85, 8.68, 7.67, 8.89, 7.23,
          6.12, 5.94, 8.46, 4.44, 6.59, 7.07, 6.47, 9.14, 6.62, 9.27, 7.08,
          4.44, 9.15, 9.08, 9.15, 8.6,  6.95, 7.3,  5.87, 7.9,  7.22, 6.32,
          6.43, 6.73, 7.79, 9.07, 8.93, 7.2,  6.75, 8.24, 7.81, 7.23, 4.24,
          9.36, 9.29, 6.65, 8.81, 5.04, 8.14, 6.08, 5.62, 6.99, 7.7,  7.63,
          7.87, 7.15, 6.15, 8.23, 4.42, 6.37, 6.03, 5.44, 6.51, 2.33, 5.11,
          8.4,  7.97, 9.61, 5.95, 7.67, 4.61, 5.71, 7.13, 3.36, 9.54, 5.02,
          8.81, 6.49, 9.65, 7.74, 7.9,  3.98, 5.92, 9.21, 7.28, 9.06, 6.76,
          7.78, 8.38, 6.03, 5.91, 7.53, 5.97, 7.81, 8.16, 7.48, 7.24,
      },
      n);

  auto model = HCluster<HClusterLinkage::kAverageU>(n);
  model.Calculate(dist);

  /*Assert::IsTrue(model.height->Equals(Matrix<Tv>(new Tv[] {6, 6, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 8, 8, 8}, (Ti)15)));

  Assert::IsTrue(model.merge->Equals(Matrix<Ti, Ti>(new Ti[] {
-1,-14,-16,-10,-2,-9,-12,-7,-8,-11,-15,-6,-4,-5,-13,
-3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
}, (Ti)15, (Ti)2)));*/

  std::vector<std::unique_ptr<std::vector<Ti>>> map;
  for (int i = 0; i < 4; i++)
    map.push_back(std::make_unique<std::vector<Ti>>());
  model.Group(map);

  ASSERT_EQ(map.at(0)->size(), (Ti)8);
  ASSERT_EQ(map.at(1)->at(0), (Ti)12);
  ASSERT_EQ(map.at(1)->at(1), (Ti)14);

  Matrix<Tv> heights = Matrix<Tv>(new Tv[n - 1], n - 1);
  Matrix<Ti> merge = Matrix<Ti>(new Ti[2 * (n - 1)], n - 1, (Ti)2);
  std::vector<Ti> order = std::vector<Ti>();
  model.MergeR(merge, heights, order);

  ASSERT_NEAR(2.33, heights.Data[0], 1e-10);
  ASSERT_NEAR(3.36, heights.Data[1], 1e-10);
  ASSERT_EQ(14, merge.Get(14, 0));
  ASSERT_EQ(13, merge.Get(14, 1));
}

TEST(Cluster_t, group_t) {
  Tv *D = new Tv[154]{
      32.446, 44.145, 17.062, 65.818, 76.19,  40.408, 78.131, 26.695, 21.992,
      68.033, 98.872, 61.154, 71.842, 66.922, 31.142, 58.429, 45.123, 80.99,
      26.345, 50.096, 36.478, 29.377, 27.141, 65.037, 72.621, 63.391, 68.125,
      60.369, 76.384, 5.449,  99.204, 6.87,   2.514,  98.799, 95.082, 6.048,
      59.287, 48.889, 21.464, 54.972, 14.997, 27.161, 92.522, 65.383, 51.852,
      71.011, 55.434, 89.082, 59.556, 29.524, 21.193, 2.684,  35.457, 69.849,
      76.352, 49.455, 18.762, 8.492,  95.032, 39.042, 32.517, 13.667, 91.408,
      23.432, 56.526, 33.531, 10.67,  72.891, 11.796, 31.202, 96.893, 2.552,
      82.001, 87.786, 96.292, 93.249, 11.688, 19.522, 37.55,  55.967, 97.026,
      14.017, 19.869, 60.988, 91.525, 33.192, 50.666, 97.465, 58.493, 17.033,
      76.138, 3.432,  58.561, 69.172, 56.453, 46.325, 63.116, 84.577, 12.265,
      77.277, 9.141,  69.192, 65.464, 29.827, 8.261,  26.696, 94.1,   64.958,
      68.924, 97.838, 91.389, 76.779, 56.448, 14.524, 33.549, 39.059, 94.886,
      98.52,  80.476, 2.754,  93.605, 17.733, 37.658, 97.567, 2.705,  74.385,
      59.03,  10.732, 82.043, 92.891, 69.384, 86.848, 40.02,  62.295, 18.609,
      61.597, 22.438, 67.702, 83.393, 96.283, 64.895, 34.39,  42.212, 52.377,
      24.745, 42.534, 64.688, 7.392,  82.462, 22.022, 68.858, 55.901, 98.156,
      96.029};
  Matrix<Tv> data = Matrix<Tv>(D, 22, 7);

  auto model = GroupDataBase::GetFromType(
      HClusterLinkage::kSingle, DistanceMethod::kEuclidean,
      CorrelationMethod::kPearson, data.RowsCount, data.ColsCount);

  auto W = new Tv[model->WorkSize];
  model->Calculate(data, W, 4, 0.0);

  ASSERT_EQ(4, model->Groups.at(0)->at(0));
  ASSERT_EQ(5, model->Groups.at(0)->at(1));
  ASSERT_EQ(3, model->Groups.at(0)->size());

  delete[] W;
  delete[] D;
}
