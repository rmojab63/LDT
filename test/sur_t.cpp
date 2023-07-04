/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include <gtest/gtest.h>

#include <random>

#include "matrix.h"
#include "sur.h"

using namespace ldt;

TEST(Sur_t, unresricted_t) {
  auto Y = Matrix<Tv>(new Tv[40]{1,  9, 2, 3, 7, 2, 9, 2, 6, 10, 7,  6, 4, 3,
                                 5,  7, 9, 1, 4, 4, 3, 7, 2, 9,  0,  5, 1, 5,
                                 10, 2, 4, 7, 3, 9, 8, 2, 3, 8,  10, 10},
                      20, 2);
  auto X = Matrix<Tv>(new Tv[60]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 8, 4, 4, 6, 5, 3, 9, 4, 2, 6,
                                 3, 2, 3, 8, 3, 7, 3, 6, 7, 5, 7, 8, 2, 4, 10,
                                 2, 4, 2, 6, 2, 5, 3, 6, 5, 8, 3, 9, 3, 7, 7},
                      20, 3);

  auto N = Y.RowsCount;
  auto m = Y.ColsCount;
  auto k = X.ColsCount;

  auto model = Sur(N, m, k, false, true);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(Y, X, storage, work);

  ASSERT_NEAR(model.beta.Get(1, 0), -0.0764430204198711, 1e-10);
  ASSERT_NEAR(model.gamma.Data[1], -0.0764430204198711, 1e-10);
  ASSERT_NEAR(model.resid_var.Get(0, 0), 7.2380672851609207,
              1e-10); // dof adjusted: 8.515373276659906
  ASSERT_NEAR(model.resid_var.Get(1, 0), -3.3455358423796087,
              1e-10); //-3.935924520446599
  ASSERT_NEAR(model.e_beta_std.Get(1, 0), 0.29126987042974234,
              1e-10); // 0.3159265317055236
  ASSERT_NEAR(model.e_beta_prob.Get(0, 1), 0.028084973367330734,
              1e-10); // 0.043347800043200824

  ASSERT_NEAR(model.r2, 0.045871940324214266, 1e-10);
  ASSERT_NEAR(model.logLikelihood, -98.252596228066253, 1e-10);

  delete[] work;

  // prediction
  auto forc = SurProjection(1, m, k, false, true);
  auto newx = Matrix<Tv>(new Tv[3]{5, 7, 5}, 1, 3);
  auto fstorage = new Tv[forc.StorageSize];
  auto fwork = new Tv[forc.WorkSize];
  forc.Calculate(model, newx, fstorage, fwork);

  ASSERT_NEAR(forc.Means.Data[0], 20.27870498780536, 1e-10);
  ASSERT_NEAR(forc.Covariance.Get(0, 0), 71.712767221393293, 1e-10);
  ASSERT_NEAR(forc.Covariance.Get(1, 0), -33.146643108342261, 1e-10);

  delete[] fstorage;
  delete[] fwork;
  delete[] storage;
}

TEST(Sur_t, restricted_t) {
  auto Y = Matrix<Tv>(new Tv[40]{1,  9, 2, 3, 7, 2, 9, 2, 6, 10, 7,  6, 4, 3,
                                 5,  7, 9, 1, 4, 4, 3, 7, 2, 9,  0,  5, 1, 5,
                                 10, 2, 4, 7, 3, 9, 8, 2, 3, 8,  10, 10},
                      20, 2);
  auto X = Matrix<Tv>(new Tv[60]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 8, 4, 4, 6, 5, 3, 9, 4, 2, 6,
                                 3, 2, 3, 8, 3, 7, 3, 6, 7, 5, 7, 8, 2, 4, 10,
                                 2, 4, 2, 6, 2, 5, 3, 6, 5, 8, 3, 9, 3, 7, 7},
                      20, 3);

  auto N = Y.RowsCount;
  auto m = Y.ColsCount;
  auto k = X.ColsCount;

  auto R = Matrix<Tv>(new Tv[6 * 5], 6, 5);
  R.SetValue(0); // an identity without the second column
  R.Set(0, 0, 1);
  R.Set(2, 1, 1);
  R.Set(3, 2, 1);
  R.Set(4, 3, 1);
  R.Set(5, 4, 1);

  auto model = Sur(N, m, k, true, true);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(Y, X, storage, work, &R);

  ASSERT_NEAR(model.beta.Get(2, 0), 0.3088535291717698, 1e-10);
  ASSERT_NEAR(model.gamma.Data[1], 0.3088535291717698, 1e-10);
  ASSERT_NEAR(model.beta.Get(0, 1), 6.0046405569804353, 1e-10);
  ASSERT_NEAR(model.gamma.Data[2], 6.0046405569804353, 1e-10);
  ASSERT_NEAR(model.resid_var.Get(1, 1), 10.311287947999730, 1e-10);
  ASSERT_NEAR(model.e_beta_std.Get(0, 1), 2.3552559004872551, 1e-10);

  ASSERT_NEAR(model.r2, 0.044226640704708720, 1e-10);

  delete[] work;

  // prediction
  auto forc = SurProjection(1, m, k, true, true);
  auto newx = Matrix<Tv>(new Tv[3]{5, 7, 5}, 1, 3);
  auto fstorage = new Tv[forc.StorageSize];
  auto fwork = new Tv[forc.WorkSize];
  forc.Calculate(model, newx, fstorage, fwork);

  ASSERT_NEAR(forc.Means.Data[0], 18.84128926968576, 1e-10);
  ASSERT_NEAR(forc.Covariance.Get(0, 0), 41.740528581906879, 1e-10);
  ASSERT_NEAR(forc.Covariance.Get(1, 0), -19.293055583626728, 1e-10);

  delete[] fstorage;
  delete[] fwork;
  delete[] storage;
}

TEST(Sur_t, sig_search_t) {
  auto Y = Matrix<Tv>(new Tv[40]{1,  9, 2, 3, 7, 2, 9, 2, 6, 10, 7,  6, 4, 3,
                                 5,  7, 9, 1, 4, 4, 3, 7, 2, 9,  0,  5, 1, 5,
                                 10, 2, 4, 7, 3, 9, 8, 2, 3, 8,  10, 10},
                      20, 2);
  auto X = Matrix<Tv>(
      new Tv[100]{1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1,  1,  1,  1,
                  1,  1,  1,  1,  1,  8,  4,  4,  6,  5,   3,  9,  4,  2,  6,
                  3,  2,  3,  8,  3,  7,  3,  6,  7,  5,   7,  8,  2,  4,  10,
                  2,  4,  2,  6,  2,  5,  3,  6,  5,  8,   3,  9,  3,  7,  7,
                  3,  5,  3,  8,  20, 3,  5,  6,  8,  2,   6,  6,  7,  8,  4,
                  2,  6,  8,  2,  4,  13, 52, 32, 18, 220, 35, 53, 62, 48, 52,
                  36, 56, 77, 98, 34, 42, 66, 78, 32, 54},
      20, 5);

  auto N = Y.RowsCount;
  auto m = Y.ColsCount;
  auto k = X.ColsCount;
  auto km = k * m;
  auto R = Matrix<Tv>(new Tv[km * km], km, km);

  auto model = Sur(N, m, k, true, true, 10);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(Y, X, storage, work, &R, 0.3);

  ASSERT_EQ(model.mSigSearchIter, 2);
  ASSERT_NEAR(model.beta.Get(2, 0), 0.38220977181777094, 1e-10);
  ASSERT_NEAR(model.gamma.Data[1], 0.38220977181777094, 1e-10);

  delete[] work;

  // prediction
  auto forc = SurProjection(5, m, k, true, true);
  auto newx = Matrix<Tv>(
      new Tv[20]{1, 9, 2, 3, 7, 2, 9, 2, 6, 10, 7, 6, 4, 3, 5, 7, 9, 1, 4, 4},
      4, 5);
  auto fstorage = new Tv[forc.StorageSize];
  auto fwork = new Tv[forc.WorkSize];
  forc.Calculate(model, newx, fstorage, fwork);

  ASSERT_NEAR(forc.Variances.Data[0], 7.6717282164384795, 1e-10); // self

  delete[] fstorage;
  delete[] fwork;
  delete[] R.Data;
  delete[] storage;
}

TEST(Sur_t, extended_t) {
  // using 'sig_search_t' data and results
  auto data = Matrix<Tv>(
      new Tv[147]{1,  NAN, 9,  2,  3,  7,   2,  9,  2,   6,  10, 7,  6,  4,
                  3,  5,   7,  9,  1,  4,   4,  3,  100, 7,  2,  9,  0,  5,
                  1,  5,   10, 2,  4,  7,   3,  9,  8,   2,  3,  8,  10, 10,
                  1,  100, 1,  1,  1,  1,   1,  1,  1,   1,  1,  1,  1,  1,
                  1,  1,   1,  1,  1,  1,   1,  8,  100, 4,  4,  6,  5,  3,
                  9,  4,   2,  6,  3,  2,   3,  8,  3,   7,  3,  6,  7,  5,
                  7,  100, 8,  2,  4,  10,  2,  4,  2,   6,  2,  5,  3,  6,
                  5,  8,   3,  9,  3,  7,   7,  3,  100, 5,  3,  8,  20, 3,
                  5,  6,   8,  2,  6,  6,   7,  8,  4,   2,  6,  8,  2,  4,
                  13, 100, 52, 32, 18, 220, 35, 53, 62,  48, 52, 36, 56, 77,
                  98, 34,  42, 66, 78, 32,  54},
      21, 7);
  auto newX = Matrix<Tv>(
      new Tv[20]{1, 9, 2, 3, 7, 2, 9, 2, 6, 10, 7, 6, 4, 3, 5, 7, 9, 1, 4, 4},
      4, 5);

  auto N = data.RowsCount;
  auto m = 2;
  auto k = 5;
  auto km = k * m;
  auto R = Matrix<Tv>(new Tv[km * km], km, km);

  auto model =
      SurExtended(N, m, k, true, true, true, 4, 10, true, nullptr, nullptr);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(data, m, storage, work, &R, 0.3, &newX, nullptr);

  ASSERT_EQ(model.Model.mSigSearchIter, 2);
  ASSERT_NEAR(model.Model.beta.Get(2, 0), 0.38220977181777094, 1e-10);
  ASSERT_NEAR(model.Model.gamma.Data[1], 0.38220977181777094, 1e-10);

  ASSERT_NEAR(model.Projections.Variances.Data[0], 7.6717282164384795,
              1e-10); // from 'sig_search_t' (self)

  delete[] work;
  delete[] R.Data;
  delete[] storage;
}

TEST(Sur_t, extended_pca_t) {
  // using 'sig_search_t' data and results
  auto data = Matrix<Tv>(
      new Tv[147]{1,  NAN, 9,  2,  3,  7,  2,  9,  2,   6,  10, 7,  6,   4,
                  3,  5,   7,  9,  1,  4,  4,  3,  100, 7,  2,  9,  0,   5,
                  1,  5,   10, 2,  4,  7,  3,  9,  8,   2,  3,  8,  10,  10,
                  8,  100, 4,  4,  6,  5,  3,  9,  4,   2,  6,  3,  2,   3,
                  8,  3,   7,  3,  6,  7,  5,  7,  100, 8,  2,  4,  10,  2,
                  4,  2,   6,  2,  5,  3,  6,  5,  8,   3,  9,  3,  7,   7,
                  3,  100, 5,  3,  8,  20, 3,  5,  6,   8,  2,  6,  6,   7,
                  8,  4,   2,  6,  8,  2,  4,  13, 100, 52, 32, 18, 220, 35,
                  53, 62,  48, 52, 36, 56, 77, 98, 34,  42, 66, 78, 32,  54},
      21, 6);
  auto newX = Matrix<Tv>(
      new Tv[20]{1, 9, 2, 3, 7, 2, 9, 2, 6, 10, 7, 6, 4, 3, 5, 7, 9, 1, 4, 4},
      4, 3);

  auto N = data.RowsCount;
  auto m = 3;
  auto k = 3;
  auto km = k * m;
  auto R = Matrix<Tv>(new Tv[km * km], km, km);

  auto pca_y = PcaAnalysisOptions();
  pca_y.ExactCount = 1;
  pca_y.IgnoreFirstCount = 1;
  pca_y.CutoffCountMax = 10000;

  auto pca_x = PcaAnalysisOptions();
  pca_x.ExactCount = 1;
  pca_x.IgnoreFirstCount = 1;
  pca_x.CutoffCountMax = 10000;

  auto model =
      SurExtended(N, m, k, true, true, true, 4, 10, true, &pca_y, &pca_x);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(data, m, storage, work, &R, 0.3, &newX);

  // test values ?!

  delete[] work;
  delete[] R.Data;
  delete[] storage;
}

TEST(Sur_t, sim_t) {
  auto data = Matrix<Tv>(
      new Tv[140]{1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,  14,
                  15, 16, 17, 18, 19, 20, 3,  7,  2,  9,  0,  5,  1,   5,
                  10, 2,  4,  7,  3,  9,  8,  2,  3,  8,  10, 10, 1,   1,
                  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,
                  1,  1,  1,  1,  8,  4,  4,  6,  5,  3,  9,  4,  2,   6,
                  3,  2,  3,  8,  3,  7,  3,  6,  7,  5,  7,  8,  2,   4,
                  10, 2,  4,  2,  6,  2,  5,  3,  6,  5,  8,  3,  9,   3,
                  7,  7,  3,  5,  3,  8,  20, 3,  5,  6,  8,  2,  6,   6,
                  7,  8,  4,  2,  6,  8,  2,  4,  13, 52, 32, 18, 220, 35,
                  53, 62, 48, 52, 36, 56, 77, 98, 34, 42, 66, 78, 32,  54},
      20, 7);

  auto N = data.RowsCount;
  auto m = 2;
  auto k = 5;
  auto km = k * m;
  auto R = Matrix<Tv>(new Tv[km * km], km, km);

  auto model =
      SurExtended(N, m, k, true, true, true, 4, 10, true, nullptr, nullptr);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(data, m, storage, work, &R, 0.3, nullptr, nullptr);

  auto meangamma = model.Model.gamma.Mean();

  auto metrics = std::vector<ScoringType>(
      {ScoringType::kRmse, ScoringType::kMae, ScoringType::kRmse});

  auto simmodel =
      SurSimulation(N, m, k, 0.8, 0, metrics, true, 10, nullptr, nullptr);
  auto storage0 = new Tv[simmodel.StorageSize];
  auto work0 = new Tv[simmodel.WorkSize];
  bool cancel = false;
  simmodel.Calculate(data, m, storage0, work0, &R, cancel, 10, 340, 0.3,
                     INFINITY);

  ASSERT_NEAR(simmodel.Results.Get(0, 0), simmodel.Results.Get(2, 0), 1e-15);

  ASSERT_NEAR(meangamma, model.Model.gamma.Mean(), 1e-15);

  delete[] work;
  delete[] R.Data;
  delete[] storage;
}

TEST(Sur_t, search_t) {
  auto data = Matrix<Tv>(
      new Tv[140]{1,  9,  2,  3,  7,  2,  9,  2,  6,  10, 7,  6,  4,   3,
                  5,  7,  9,  1,  4,  4,  3,  7,  2,  9,  0,  5,  1,   5,
                  10, 2,  4,  7,  3,  9,  8,  2,  3,  8,  10, 10, 1,   1,
                  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,
                  1,  1,  1,  1,  8,  4,  4,  6,  5,  3,  9,  4,  2,   6,
                  3,  2,  3,  8,  3,  7,  3,  6,  7,  5,  7,  8,  2,   4,
                  10, 2,  4,  2,  6,  2,  5,  3,  6,  5,  8,  3,  9,   3,
                  7,  7,  3,  5,  3,  8,  20, 3,  5,  6,  8,  2,  6,   6,
                  7,  8,  4,  2,  6,  8,  2,  4,  13, 52, 32, 18, 220, 35,
                  53, 62, 48, 52, 36, 56, 77, 98, 34, 42, 66, 78, 32,  54},
      20, 7);

  auto N = data.RowsCount;
  auto m = 2;
  auto k = 5;

  auto parallel = false;
  auto simfixsize = 4;
  auto sigseaItr = 0;

  auto items = SearchItems();
  auto searchOptions = SearchOptions();
  auto metrics = SearchMetricOptions();
  auto checks = SearchModelChecks();

  items.KeepBestCount = 2;
  items.KeepAll = false;
  items.LengthTargets = 2;

  items.KeepModelEvaluations = true;
  items.Length1 = 0;
  items.LengthDependents = m;
  items.LengthExogenouses = k;
  items.ExtremeBoundsMultiplier = 2.0;

  metrics.SimFixSize = simfixsize;
  searchOptions.Parallel = parallel;
  metrics.Seed = 340;
  metrics.TrainRatio = 0.8;
  metrics.MetricsIn = std::vector<GoodnessOfFitType>({GoodnessOfFitType::kAic});
  metrics.MetricsOut =
      std::vector<ScoringType>({ScoringType::kMae, ScoringType::kCrps}); //

  auto exosizes = std::vector<Ti>({1, 2, 3});
  auto endoIndexes = std::vector<std::vector<Ti>>({{0}, {0, 1}, {0, 1, 2}});
  auto exogroups = std::vector<std::vector<Ti>>({{3}, {4}, {5}, {6}});

  auto model = SurModelset(searchOptions, items, metrics, checks, exosizes,
                           exogroups, 3, data, endoIndexes, sigseaItr, 0.5);
  auto W = new Tv[model.Modelset.WorkSize];
  model.Modelset.Start(W, nullptr);

  SearcherModelingInfo info;
  auto list1 = std::vector<SearcherSummary *>();
  auto list2 = std::vector<SearcherSummary *>();
  auto list3 = std::vector<SearcherSummary *>();
  model.Modelset.CombineInfo(info, list1, list2, list3);

  auto res = std::vector<EstimationKeep *>();
  model.Modelset.CombineAll(0, 0, 0, list1, res);

  delete[] W;
}
