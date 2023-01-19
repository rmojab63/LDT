/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "matrix.h"
#include "pca.h"
#include <gtest/gtest.h>
#include <random>

using namespace ldt;

TEST(Pca_t, pca_t) {

  Ti m = 10;
  Ti n = 20;
  auto mat = Matrix<Tv>(new Tv[m * n], m, n);
  Matrix<Tv>::FillRandom_normal(mat, 340);
  auto model = PcaAnalysis(m, n, m, true, true, true, true);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(mat, work, storage, &mat);

  ASSERT_NEAR(model.Stds.Data[0], 2.070056331812816, 1e-2);
  ASSERT_NEAR(model.Stds.Data[1], 1.829109854115231, 0.2);
  ASSERT_NEAR(model.Directions.Data[0], -0.37872362855843755, 0.4);

  ASSERT_NEAR(model.PCs.Data[0], -2.1986570859879695, 0.7);

  ASSERT_NEAR(model.Forecasts.Get(0, 0), model.PCs.Get(0, 0), 1e-16);

  delete[] storage;
  delete[] work;
}

TEST(Pca_t, pca_no_scale_center_t) {
  Ti m = 10;
  Ti n = 20;

  auto mat = Matrix<Tv>(new Tv[m * n], m, n);
  Matrix<Tv>::FillRandom_normal(mat, 340);
  auto model = PcaAnalysis(m, n, m, true, true, false, false);
  Tv *storage = new Tv[model.StorageSize];
  Tv *work = new Tv[model.WorkSize];
  model.Calculate(mat, work, storage, &mat);

  // ASSERT_NEAR(model.Stds.Data[0], 2.2006003752686980, 0.2);
  // ASSERT_NEAR(model.PCs.Data[0], -3.4149787129079199, 0.2);

  ASSERT_NEAR(model.Forecasts.Get(0, 0), model.PCs.Get(0, 0), 1e-16);
  delete[] storage;
  delete[] work;
}

TEST(Pca_t, pca_var_zero_t) {
  Ti m = 10;
  Ti n = 22;

  auto mat = Matrix<Tv>(new Tv[m * n], m, n);
  Matrix<Tv>::FillRandom_normal(mat, 340);
  mat.SetColumn(20, 0.0); // the last two columns have zero variances
  mat.SetColumn(21, 0.0);

  auto model = PcaAnalysis(m, n, m, true, true, true, true);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(mat, work, storage, &mat);

  ASSERT_NEAR(model.Stds.Data[0], 2.070056331812816, 1e-2);
  ASSERT_NEAR(model.Stds.Data[1], 1.829109854115231, 0.2);
  ASSERT_NEAR(model.Directions.Data[0], -0.37872362855843755, 0.4);
  ASSERT_NEAR(model.PCs.Data[0], -2.1986570859879695, 0.7);

  ASSERT_NEAR(model.Forecasts.Get(0, 0), model.PCs.Get(0, 0), 1e-16);

  delete[] storage;
  delete[] work;
}