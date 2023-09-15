/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include <gtest/gtest.h>

#include "matrix.h"
#include "varma.h"

using namespace ldt;

TEST(Varma_T, sizes) {
  auto sizes = VarmaSizes(100, 2, 3, 1, 0, 2, 2, 0, 1,
                          3); // y(-1),y(-3),y(-6), e(-1), e(-2), e(-3)
  auto W = new Ti[sizes.WorkSizeI];
  sizes.Calculate(W);

  ASSERT_EQ(true, std::equal(sizes.ArLags.begin(), sizes.ArLags.end(),
                             (std::vector<Ti>{1, 3, 6}).begin()));
  ASSERT_EQ(true, std::equal(sizes.MaLags.begin(), sizes.MaLags.end(),
                             (std::vector<Ti>{1, 2, 3}).begin()));

  delete[] W;
}

TEST(Varma_T, sizes_diff) {
  auto sizes = VarmaSizes(100, 2, 3, 1, 1, 2, 2, 0, 1, 3);

  ASSERT_EQ(true, std::equal(sizes.DiffPoly.begin(), sizes.DiffPoly.end(),
                             (std::vector<Ti>{1, -1}).begin()));

  auto sizes1 = VarmaSizes(100, 2, 3, 1, 0, 2, 2, 1, 1, 3);

  ASSERT_EQ(true, std::equal(sizes1.DiffPoly.begin(), sizes1.DiffPoly.end(),
                             (std::vector<Ti>{1, 0, 0, -1}).begin()));

  auto sizes2 = VarmaSizes(100, 2, 3, 1, 2, 2, 2, 3, 1, 3);

  ASSERT_EQ(
      true,
      std::equal(
          sizes2.DiffPoly.begin(), sizes2.DiffPoly.end(),
          (std::vector<Ti>{1, -2, 1, -3, 6, -3, 3, -6, 3, -1, 2, -1}).begin()));
}

TEST(Varma_T, getARMAs) {
  Ti g = 2;
  Ti exo = 3;

  auto sizes = VarmaSizes(100, g, exo, 3, 0, 3, 0, 0, 0, 0);

  auto V = new Ti[sizes.WorkSizeI];
  sizes.Calculate(V);
  delete[] V;

  auto arma = VarmaArma(sizes, 0);
  auto S = std::unique_ptr<Tv[]>(new Tv[arma.StorageSize]);
  auto W = std::unique_ptr<Tv[]>(new Tv[arma.WorkSize]);

  auto Pi =
      Matrix<Tv>(new Tv[26]{1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
                            21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33},
                 g, g * 3 + exo + g * 3);
  arma.Calculate(Pi, S.get(), W.get());

  ASSERT_EQ(true, arma.Ar.Coefficients.at(0)->Equals(
                      Matrix<Tv>(new Tv[4]{1, 0, 0, 1}, 2, 2)));
  ASSERT_EQ(true, arma.Ar.Coefficients.at(1)->Equals(
                      Matrix<Tv>(new Tv[4]{-1, -2, -3, -4}, 2, 2)));

  // MA
  ASSERT_EQ(true, arma.Ma.Coefficients.at(0)->Equals(
                      Matrix<Tv>(new Tv[4]{1, 0, 0, 1}, 2, 2)));
  ASSERT_EQ(true, arma.Ma.Coefficients.at(1)->Equals(
                      Matrix<Tv>(new Tv[4]{26, 27, 28, 29}, 2, 2)));
}

TEST(Varma_T, getMAinf) {
  // Lutkpohl, 2.1.14
  Ti g = 3;

  auto sizes = VarmaSizes(100, g, 0, 1, 0, 0, 0, 0, 0, 0);
  auto V = new Ti[sizes.WorkSizeI];
  sizes.Calculate(V);
  delete[] V;

  auto arma = VarmaArma(sizes, 4);
  auto S = std::unique_ptr<Tv[]>(new Tv[arma.StorageSize]);
  auto W = std::unique_ptr<Tv[]>(new Tv[arma.WorkSize]);

  auto Pi =
      Matrix<Tv>(new Tv[9]{0.5, 0.1, 0.0, 0.0, 0.1, 0.2, 0.0, 0.3, 0.3}, g, g);
  arma.Calculate(Pi, S.get(), W.get());

  // AR
  ASSERT_EQ(true, arma.Ar.Coefficients.at(1)->Equals(
                      Matrix<Tv>(new Tv[9]{-0.5, -0.1, -0.0, -0.0, -0.1, -0.2,
                                           -0.0, -0.3, -0.3},
                                 3, 3),
                      1e-8));

  // MAinf

  ASSERT_EQ(true,
            arma.MaInf.Coefficients.at(0)->Equals(
                Matrix<Tv>(new Tv[9]{1, 0, 0, 0, 1, 0, 0, 0, 1}, 3, 3), 1e-8));
  ASSERT_EQ(true, arma.MaInf.Coefficients.at(1)->Equals(
                      Matrix<Tv>(new Tv[9]{0.5, 0.1, 0.0, 0.0, 0.1, 0.2, 0.0,
                                           0.3, 0.3},
                                 3, 3),
                      1e-8));
  ASSERT_EQ(true, arma.MaInf.Coefficients.at(2)->Equals(
                      Matrix<Tv>(new Tv[9]{0.25, 0.06, 0.02, 0.0, 0.07, 0.08,
                                           0.0, 0.12, 0.15},
                                 3, 3),
                      1e-8));

  delete[] Pi.Data;
}

TEST(Varma_T, difference) {
  auto sizes1 = VarmaSizes(15, 3, 0, 1, 1, 2, 2, 0, 1, 3);

  const Matrix<Tv> y = Matrix<Tv>(
      new std::vector<Tv>{1, 2, 3, 4, 5, 3, 4, 6, 7, 5,  4, 5, 6, 7, 5,
                          4, 6, 5, 4, 5, 7, 5, 8, 9, 7,  3, 7, 8, 5, 3,
                          5, 3, 3, 4, 6, 5, 8, 3, 6, 10, 2, 4, 7, 4, 6},
      3, 15);

  auto diffy1 = Matrix<Tv>(1000, new Tv[42], 3, 14);
  Varma::Difference(sizes1.DiffPoly, y, diffy1);
  ASSERT_EQ(diffy1.Get(1, 2), -2.0);
  ASSERT_EQ(diffy1.Get(2, 3), 0.0);

  // d=2
  auto sizes2 = VarmaSizes(15, 3, 0, 1, 2, 2, 2, 0, 1, 3);
  auto diffy2 = Matrix<Tv>(1000, new Tv[39], 3, 15 - (3 - 1));
  Varma::Difference(sizes2.DiffPoly, y, diffy2);
  ASSERT_EQ(diffy2.Get(1, 2), 5.0);

  delete y.Data;
}

TEST(Varma_T, ols) {
  Tv *D = new Tv[300];
  Matrix<Tv> data = Matrix<Tv>(D, 3, 100);
  Matrix<Tv>::FillRandom_uniform(data, 3243, 1, 100);
  Tv *E = new Tv[200];
  Matrix<Tv> exo = Matrix<Tv>(E, 2, 100);
  Matrix<Tv>::FillRandom_uniform(exo, 87, 1, 100);
  Ti info = 0;
  auto sizes = VarmaSizes(100, 3, 2, 2, 0, 2, 0, 0, 0, 0);
  auto var = Varma(sizes, false, true);
  auto S = std::unique_ptr<Tv[]>(new Tv[var.Result.StorageSize]);
  auto W = std::unique_ptr<Tv[]>(new Tv[var.Result.WorkSize]);
  var.EstimateOls(data, &exo, nullptr, nullptr, W.get(), S.get(), 0, false,
                  nullptr);

  ASSERT_NEAR(-0.44874739420702131, var.Result.gamma.Data[0], 1e-6);
  ASSERT_NEAR(0.65034107480637338, var.Result.gamma.Data[3 * 14 - 1], 1e-6);

  delete[] D;
  delete[] E;
}

TEST(Varma_T, finalFormMA) {
  Ti equCount = 3;
  Ti exoCount = 2;

  auto sizes = VarmaSizes(100, 3, 2, 2, 0, 2, 0, 0, 0, 0);
  auto rest = VarmaRestriction(sizes, VarmaRestrictionType::kMaFinal);
  Tv *F = new Tv[rest.StorageSize];
  rest.Calculate(F);
  Tv *G = new Tv[rest.R.ColsCount];
  Matrix<Tv> gamma = Matrix<Tv>(1.0, G, rest.R.ColsCount, 1);
  Tv *D = new Tv[rest.R.RowsCount];
  Matrix<Tv> coef = Matrix<Tv>(D, equCount, rest.R.RowsCount / equCount);
  rest.R.Dot0(gamma, coef);

  ASSERT_NEAR(30, coef.Sum(), 1e-16);

  delete[] F;
  delete[] G;
  delete[] D;
}

TEST(Varma_T, olsRest) {
  Tv *D = new Tv[300];
  Matrix<Tv> data = Matrix<Tv>(D, 3, 100);
  Matrix<Tv>::FillRandom_uniform(data, 3243, 1, 100);
  Tv *E = new Tv[200];
  Matrix<Tv> exo = Matrix<Tv>(E, 2, 100);
  Matrix<Tv>::FillRandom_uniform(exo, 87, 1, 100);

  // define model
  auto sizes = VarmaSizes(100, 3, 2, 2, 0, 2, 0, 0, 0, 0);
  auto var = Varma(sizes, true, true);

  // create restriction Matrix
  auto rest = VarmaRestriction(sizes, VarmaRestrictionType::kMaFinal);
  Tv *F = new Tv[rest.StorageSize];
  rest.Calculate(F);

  // set storage
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateOls(data, &exo, &rest.R, nullptr, W, S, 0, false, nullptr);

  ASSERT_NEAR(0.2411267168945912, var.Result.gamma.Data[0],
              1e-6); // self-evaluated

  delete[] D;
  delete[] E;
  delete[] F;
  delete[] S;
  delete[] W;
}

// there are two types of tests for ARMA
//     m: just gamma and mean of forecasts
//     v: variances, too

TEST(Varma_T, arma_m_ar2) {
  const Matrix<Tv> y = Matrix<Tv>(
      new Tv[16]{1, 2, 3, 4, 5, 3, 4, 6, 7, 5, 4, 6, 5, 4, 5, 7}, 1, 16);
  const Matrix<Tv> exo =
      Matrix<Tv>(new Tv[25]{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                 1, 25);

  auto sizes = VarmaSizes(16, 1, 1, 0, 0, 0, 1, 0, 0, 2);

  // estimate
  auto var = Varma(sizes, false, true, false);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S);

  ASSERT_NEAR(var.Result.gamma.Get(1), 0.311079, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, false, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, nullptr, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(2), 7.016449, 1e-5);
  ASSERT_NEAR(forc.Forecast.Get(4), 8.460990, 1e-5);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, arma_m_ar2_diff) {
  const Matrix<Tv> y = Matrix<Tv>(
      new Tv[16]{1, 2, 3, 4, 5, 3, 4, 6, 7, 5, 4, 6, 5, 4, 5, 7}, 1, 16);
  const Matrix<Tv> exo =
      Matrix<Tv>(new Tv[25]{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                 1, 25);

  auto sizes = VarmaSizes(16, 1, 1, 0, 1, 0, 1, 0, 0, 2);

  // estimate
  auto var = Varma(sizes, false, true, false);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S);

  ASSERT_NEAR(var.Result.gamma.Get(0), -0.60386833, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, false, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(2), 7, 1e-5);
  ASSERT_NEAR(forc.Forecast.Get(4), 6.2521868, 1e-5);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, arma_m_ar2_diff2) {
  const Matrix<Tv> y = Matrix<Tv>(
      new Tv[16]{1, 2, 3, 4, 5, 3, 4, 6, 7, 5, 4, 6, 5, 4, 5, 7}, 1, 16);
  const Matrix<Tv> exo =
      Matrix<Tv>(new Tv[25]{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                 1, 25);

  auto sizes = VarmaSizes(16, 1, 1, 0, 2, 0, 1, 0, 0, 2);

  // estimate
  auto var = Varma(sizes, false, true, false);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S);

  ASSERT_NEAR(var.Result.gamma.Get(0), -0.57874717, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, false, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(2), 5, 1e-5);
  ASSERT_NEAR(forc.Forecast.Get(5), 8.596566377, 1e-5);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, arma_m_ma1) {
  const Matrix<Tv> y = Matrix<Tv>(
      new Tv[16]{1, 2, 3, 4, 5, 3, 4, 6, 7, 5, 4, 6, 5, 4, 5, 7}, 1, 16);

  auto sizes = VarmaSizes(16, 1, 0, 0, 0, 1, 0, 0, 0, 2);

  // estimate
  auto var = Varma(sizes, false, true, false);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, nullptr, W, S);

  ASSERT_NEAR(var.Result.gamma.Get(0), 0.96108589, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, false, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, nullptr, &y, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(0), 3.43409279, 1e-4);
  ASSERT_NEAR(forc.Forecast.Get(2), 0, 1e-4);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, arma_m_ma1_season) {
  const Matrix<Tv> y = Matrix<Tv>(
      new Tv[16]{1, 2, 3, 4, 5, 3, 4, 6, 7, 5, 4, 6, 5, 4, 5, 7}, 1, 16);
  const Matrix<Tv> exo =
      Matrix<Tv>(new Tv[25]{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                 1, 25);

  auto sizes = VarmaSizes(16, 1, 1, 0, 0, 0, 0, 0, 1, 2);

  // estimate
  auto var = Varma(sizes, false, true, false);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S);

  ASSERT_NEAR(var.Result.gamma.Get(0), 0.49071285, 1e-5);
  ASSERT_NEAR(var.Result.gamma.Get(1), 0.2599969, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, false, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(0), 7.3957475067466092, 1e-5);
  ASSERT_NEAR(forc.Forecast.Get(4), 9.8142570659198078, 1e-5);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, arma_m_arma_season) {
  const Matrix<Tv> y = Matrix<Tv>(
      new Tv[16]{1, 2, 3, 4, 5, 3, 4, 6, 7, 5, 4, 6, 5, 4, 5, 7}, 1, 16);
  const Matrix<Tv> exo =
      Matrix<Tv>(new Tv[25]{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                 1, 25);

  auto sizes = VarmaSizes(16, 1, 1, 0, 0, 0, 1, 0, 1, 2);

  // estimate
  auto var = Varma(sizes, false, true, false);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S);

  ASSERT_NEAR(var.Result.gamma.Get(0), 1.5188137, 1e-5);
  ASSERT_NEAR(var.Result.gamma.Get(2), -1.0092161, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, false, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(2), 5.0130877, 1e-4);
  ASSERT_NEAR(forc.Forecast.Get(3), 6.0992993, 1e-4);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, arma_m_arma_diff) {
  const Matrix<Tv> y = Matrix<Tv>(
      new Tv[16]{1, 2, 3, 4, 5, 3, 4, 6, 7, 5, 4, 6, 5, 4, 5, 7}, 1, 16);
  const Matrix<Tv> exo =
      Matrix<Tv>(new Tv[25]{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                 1, 25);

  auto sizes = VarmaSizes(16, 1, 1, 0, 2, 0, 1, 0, 1, 2);

  // estimate
  auto var = Varma(sizes, false, true, false);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S);

  ASSERT_NEAR(var.Result.gamma.Get(0), -1.0016817, 1e-5);
  ASSERT_NEAR(var.Result.gamma.Get(2), 0.5973201, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, false, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(2), 5, 1e-5);
  ASSERT_NEAR(forc.Forecast.Get(4), 6.9558527, 1e-5);
  ASSERT_NEAR(forc.Forecast.Get(5), 5.5386389, 1e-4);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, arma_v_ar2) {
  const Matrix<Tv> y = Matrix<Tv>(
      new Tv[16]{1, 2, 3, 4, 5, 3, 4, 6, 7, 5, 4, 6, 5, 4, 5, 7}, 1, 16);
  const Matrix<Tv> exo =
      Matrix<Tv>(new Tv[25]{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                 1, 25);

  auto sizes = VarmaSizes(16, 1, 1, 0, 0, 0, 1, 0, 0, 2);

  // estimate
  auto var = Varma(sizes, false, true, true); // new.var=true
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S);

  ASSERT_NEAR(var.Result.gamma.Get(1), 0.311079, 1e-5);
  ASSERT_NEAR(var.Result.sigma2.Get(0, 0), 3.103117, 1e-5); // new

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, true, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, nullptr, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(2), 7.016449, 1e-5);
  ASSERT_NEAR(forc.Forecast.Get(4), 8.460990, 1e-5);
  ASSERT_NEAR(forc.Variance.Get(4), 3.619262, 1e-5); // new
  ASSERT_NEAR(forc.Variance.Get(6), 3.705114, 1e-5); // new

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, arma_v_ma1) {
  const Matrix<Tv> y = Matrix<Tv>(
      new Tv[16]{1, 2, 3, 4, 5, 3, 4, 6, 7, 5, 4, 6, 5, 4, 5, 7}, 1, 16);

  auto sizes = VarmaSizes(16, 1, 0, 0, 0, 1, 0, 0, 0, 2);

  // estimate
  auto var = Varma(sizes, false, true, true);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, nullptr, W, S);

  ASSERT_NEAR(var.Result.gamma.Get(0), 0.96108589, 1e-5);
  ASSERT_NEAR(var.Result.sigma2.Get(0, 0), 7.3243989, 1e-5);
  ASSERT_NEAR(var.Result.gammavar.Get(0, 0), 0.0079551, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, true, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, nullptr, &y, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(0), 3.43409279, 1e-4);
  ASSERT_NEAR(forc.Forecast.Get(2), 0, 1e-4);
  ASSERT_NEAR(forc.Variance.Get(0), 7.32439892, 1e-5); // new
  ASSERT_NEAR(forc.Variance.Get(3), 14.0898442, 1e-5); // new

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, arma_v_ma1_season) {
  const Matrix<Tv> y = Matrix<Tv>(
      new Tv[16]{1, 2, 3, 4, 5, 3, 4, 6, 7, 5, 4, 6, 5, 4, 5, 7}, 1, 16);
  const Matrix<Tv> exo =
      Matrix<Tv>(new Tv[25]{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                 1, 25);

  auto sizes = VarmaSizes(16, 1, 1, 0, 0, 0, 0, 0, 1, 2);

  // estimate
  auto var = Varma(sizes, false, true, true);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S);

  ASSERT_NEAR(var.Result.gamma.Get(0), 0.49071285, 1e-5);
  ASSERT_NEAR(var.Result.gamma.Get(1), 0.2599969, 1e-5);
  ASSERT_NEAR(var.Result.sigma2.Get(0), 2.888547044, 1e-5);
  ASSERT_NEAR(var.Result.gammavar.Get(2), -0.0003021813, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, true, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(0), 7.3957475067466092, 1e-5);
  ASSERT_NEAR(forc.Forecast.Get(4), 9.8142570659198078, 1e-5);
  ASSERT_NEAR(forc.Variance.Get(0), 2.88854704, 1e-5);
  ASSERT_NEAR(forc.Variance.Get(3), 3.0838082, 1e-5);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, arma_v_arma_season) {
  const Matrix<Tv> y = Matrix<Tv>(
      new Tv[16]{1, 2, 3, 4, 5, 3, 4, 6, 7, 5, 4, 6, 5, 4, 5, 7}, 1, 16);
  const Matrix<Tv> exo =
      Matrix<Tv>(new Tv[25]{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                 1, 25);

  auto sizes = VarmaSizes(16, 1, 1, 0, 0, 0, 1, 0, 1, 2);

  // estimate
  auto var = Varma(sizes, false, true, true);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S);

  ASSERT_NEAR(var.Result.gamma.Get(0), 1.5188137, 1e-5);
  ASSERT_NEAR(var.Result.gamma.Get(2), -1.0092161, 1e-5);
  ASSERT_NEAR(var.Result.sigma2.Get(0), 2.922668803, 1e-5);
  ASSERT_NEAR(var.Result.gammavar.Get(2), -0.1546095288, 1e-4);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, true, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(2), 5.0130877, 1e-4);
  ASSERT_NEAR(forc.Forecast.Get(3), 6.0992993, 1e-4);
  ASSERT_NEAR(forc.Variance.Get(4), 3.68165561, 1e-5);
  ASSERT_NEAR(forc.Variance.Get(6), 5.43248269, 1e-4);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, arma_v_ar1_diff) {
  const Matrix<Tv> y = Matrix<Tv>(
      new Tv[16]{1, 2, 3, 4, 5, 3, 4, 6, 7, 5, 4, 6, 5, 4, 5, 7}, 1, 16);
  const Matrix<Tv> exo =
      Matrix<Tv>(new Tv[25]{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                 1, 25);

  auto sizes = VarmaSizes(16, 1, 1, 1, 1, 0, 0, 0, 0, 0);

  // estimate
  auto var = Varma(sizes, false, true, true);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S);

  // ASSERT_NEAR(var.Result.gamma.Get(0), 1.5188137, 1e-5);
  // ASSERT_NEAR(var.Result.gamma.Get(2), -1.0092161, 1e-5);
  // ASSERT_NEAR(var.Result.sigma2.Get(0), 2.922668803, 1e-5);
  // ASSERT_NEAR(var.Result.gammavar.Get(2), -0.1546095288, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, true, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  // ASSERT_NEAR(forc.Forecast.Get(2), 5.0130877, 1e-4);
  // ASSERT_NEAR(forc.Forecast.Get(3), 6.0992993, 1e-4);
  ASSERT_NEAR(forc.Variance.Get(3), 4.036736, 1e-5);
  ASSERT_NEAR(forc.Variance.Get(6), 10.21648, 1e-4);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, arma_v_arma_diff) {
  const Matrix<Tv> y = Matrix<Tv>(
      new Tv[16]{1, 2, 3, 4, 5, 3, 4, 6, 7, 5, 4, 6, 5, 4, 5, 7}, 1, 16);
  const Matrix<Tv> exo =
      Matrix<Tv>(new Tv[25]{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                 1, 25);

  auto sizes = VarmaSizes(16, 1, 1, 0, 2, 0, 1, 0, 1, 2);

  // estimate
  auto var = Varma(sizes, false, true, true);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S);

  // ASSERT_NEAR(var.Result.gamma.Get(0), 1.5188137, 1e-5);
  // ASSERT_NEAR(var.Result.gamma.Get(2), -1.0092161, 1e-5);
  // ASSERT_NEAR(var.Result.sigma2.Get(0), 2.922668803, 1e-5);
  // ASSERT_NEAR(var.Result.gammavar.Get(2), -0.1546095288, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, true, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  // ASSERT_NEAR(forc.Forecast.Get(2), 5.0130877, 1e-4);
  // ASSERT_NEAR(forc.Forecast.Get(3), 6.0992993, 1e-4);
  ASSERT_NEAR(forc.Variance.Get(4), 2.676136, 1e-5);
  ASSERT_NEAR(forc.Variance.Get(5), 13.38068, 1e-4);

  delete[] S;
  delete[] W;
  delete[] S1;
}

// Test VAR

TEST(Varma_T, var) {
  Matrix<Tv> y =
      Matrix<Tv>(new Tv[45]{2, 3, 4, 6, 8, 4, 3, 7, 6, 5, 0, 0, 6, 8, 8,
                            4, 7, 8, 7, 4, 9, 7, 4, 2, 0, 4, 4, 3, 6, 0,
                            3, 9, 0, 4, 4, 4, 8, 6, 0, 9, 9, 2, 5, 0, 6},
                 3, 15);
  Matrix<Tv> exo = Matrix<Tv>(
      new Tv[40]{8, 4, 9, 1, 1, 3, 6, 4, 1, 7, 9, 2, 8, 4, 2, 4, 3, 4, 3, 2,
                 4, 9, 4, 6, 0, 9, 8, 9, 9, 7, 6, 1, 6, 7, 0, 3, 2, 2, 1, 2},
      2, 20);

  auto sizes = VarmaSizes(15, 3, 2, 2, 0, 0, 0, 0, 0, 0);

  // estimate
  auto var = Varma(sizes, false, true, true);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S);

  ASSERT_NEAR(var.Result.gamma.Data[0], 0.2353662, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, true, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(0, 1), 5, 1e-5);
  ASSERT_NEAR(forc.Forecast.Get(0, 2), 7.377207, 1e-5);
  ASSERT_NEAR(forc.Forecast.Get(0, 4), -8.65123656, 1e-5);
  ASSERT_NEAR(forc.Forecast.Get(1, 2), -8.273805, 1e-5);
  ASSERT_NEAR(forc.Forecast.Get(1, 5), -0.1645592, 1e-5);
  ASSERT_NEAR(forc.Forecast.Get(2, 5), -1.442606, 1e-5);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, var_diff) {
  Matrix<Tv> y =
      Matrix<Tv>(new Tv[45]{2, 3, 4, 6, 8, 4, 3, 7, 6, 5, 0, 0, 6, 8, 8,
                            4, 7, 8, 7, 4, 9, 7, 4, 2, 0, 4, 4, 3, 6, 0,
                            3, 9, 0, 4, 4, 4, 8, 6, 0, 9, 9, 2, 5, 0, 6},
                 3, 15);
  Matrix<Tv> exo = Matrix<Tv>(
      new Tv[40]{8, 4, 9, 1, 1, 3, 6, 4, 1, 7, 9, 2, 8, 4, 2, 4, 3, 4, 3, 2,
                 4, 9, 4, 6, 0, 9, 8, 9, 9, 7, 6, 1, 6, 7, 0, 3, 2, 2, 1, 2},
      2, 20);

  auto sizes = VarmaSizes(15, 3, 2, 0, 1, 0, 1, 0, 0, 2);

  // estimate
  auto var = Varma(sizes, false, true, true);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S);

  ASSERT_NEAR(var.Result.gamma.Data[0], -0.3231539, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, true, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(1, 5), 6.58414214, 1e-5);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, var_sample) {
  Matrix<Tv> y =
      Matrix<Tv>(new Tv[45]{2, 3, 4, 6, 8, 4, 3, 7, 6, 5, 0, 0, 6, 8, 8,
                            4, 7, 8, 7, 4, 9, 7, 4, 2, 0, 4, 4, 3, 6, 0,
                            3, 9, 0, 4, 4, 4, 8, 6, 0, 9, 9, 2, 5, 0, 6},
                 3, 15);
  Matrix<Tv> exo = Matrix<Tv>(
      new Tv[40]{8, 4, 9, 1, 1, 3, 6, 4, 1, 7, 9, 2, 8, 4, 2, 4, 3, 4, 3, 2,
                 4, 9, 4, 6, 0, 9, 8, 9, 9, 7, 6, 1, 6, 7, 0, 3, 2, 2, 1, 2},
      2, 20);

  auto sizes = VarmaSizes(15, 3, 2, 2, 0, 0, 0, 0, 0, 0);

  // estimate
  auto var = Varma(sizes, false, true, true);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S, nullptr, nullptr, 1);

  ASSERT_NEAR(var.Result.gamma.Data[0], 0.5008786, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, true, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(0, 3), 10.1541927, 1e-5);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, var_diff_sample) {
  Matrix<Tv> y =
      Matrix<Tv>(new Tv[45]{2, 3, 4, 6, 8, 4, 3, 7, 6, 5, 0, 0, 6, 8, 8,
                            4, 7, 8, 7, 4, 9, 7, 4, 2, 0, 4, 4, 3, 6, 0,
                            3, 9, 0, 4, 4, 4, 8, 6, 0, 9, 9, 2, 5, 0, 6},
                 3, 15);
  Matrix<Tv> exo = Matrix<Tv>(
      new Tv[40]{8, 4, 9, 1, 1, 3, 6, 4, 1, 7, 9, 2, 8, 4, 2, 4, 3, 4, 3, 2,
                 4, 9, 4, 6, 0, 9, 8, 9, 9, 7, 6, 1, 6, 7, 0, 3, 2, 2, 1, 2},
      2, 20);

  auto sizes = VarmaSizes(15, 3, 2, 0, 1, 0, 1, 0, 0, 2);

  // estimate
  auto var = Varma(sizes, false, true, true);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S, nullptr, nullptr, 1);

  ASSERT_NEAR(var.Result.gamma.Data[0], -0.2397216, 1e-5);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, true, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  ASSERT_NEAR(forc.Forecast.Get(0, 3), 11.086097, 1e-5);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, var_extended) {
  Matrix<Tv> data = Matrix<Tv>(
      new Tv[200]{1.1,  0.4,  1.8,  -1.6, 2.6,  0.7,  -1.8, 2.8,  0.3,  -0.5,
                  2.3,  4.1,  0.9,  3.5,  -1.1, -0.3, 0.6,  2,    0.1,  1.1,
                  0.3,  3.3,  4.2,  -1.1, 1.1,  2.7,  2.5,  -0.5, 2.6,  0.3,
                  3.2,  1.9,  0.4,  4.4,  -1.6, 2.9,  2,    3,    2.3,  0.7,
                  1.3,  0.4,  -2.6, 3.1,  0.8,  1.1,  3.5,  4.9,  -1.7, -2.8,
                  1.2,  1.7,  1.3,  0.9,  2,    0.1,  -1.3, 1.9,  0.3,  0.1,
                  0.3,  0.7,  1.1,  0.2,  1.1,  1.5,  -0.4, 0,    -0.8, 1.2,
                  1.1,  -2.1, 0,    0,    -0.5, -0.3, -0.3, 0.6,  0.2,  0.6,
                  -0.5, -1.2, 1.1,  1.3,  0.3,  -0.2, -0.5, 0.1,  1.1,  -1.1,
                  0.3,  0.8,  1.5,  0.2,  -0.3, 0.3,  -0.7, 0.6,  0.9,  0.9,
                  1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
                  1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
                  1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
                  1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
                  1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
                  -1.6, -1.5, -1.6, -0.5, -1.5, 0.7,  2.1,  -1.3, 0.8,  0.8,
                  0.3,  -1,   -0.1, -0.3, 0.6,  -0.4, 1,    -0.4, 1.1,  -1,
                  -1.3, 3.2,  -0.4, 0.3,  0.6,  -0.5, 0.5,  0.4,  -0.2, 0.1,
                  0,    2.1,  -0.7, -1.1, 0,    0.3,  0.4,  -0.5, -1.1, 1.3,
                  -0.3, -0.9, -0.2, -0.2, 1.1,  0.1,  0.8,  -0.5, 0.2,  -0.3},
      50, 4);
  auto sizes = VarmaSizes(50, 2, 2, 1, 0, 0, 0, 0, 0, 0);

  // estimate
  auto var = VarmaExtended(sizes, VarmaRestrictionType::kMaFinal, true, true,
                           true, 0, nullptr, nullptr, nullptr);
  Tv *S = new Tv[var.StorageSize];
  Tv *W = new Tv[var.WorkSize];
  var.Calculate(data, S, W, false, 0, 0);

  ASSERT_NEAR(var.Model.Result.gamma.Get(1, 3), -0.76724074489110317, 1e-10);

  delete[] S;
  delete[] W;
}

TEST(Varma_T, var_extended2) {
  Matrix<Tv> data = Matrix<Tv>(
      new Tv[75]{
          32.446, 44.145, 17.062, 65.818, 76.19,  40.408, 78.131, 26.695,
          21.992, 68.033, 98.872, 61.154, 71.842, 66.922, 31.142, 58.429,
          45.123, 80.99,  26.345, 50.096, 36.478, 29.377, NAN,    NAN,
          NAN,    27.141, 65.037, 72.621, 63.391, 68.125, 60.369, 76.384,
          5.449,  99.204, 6.87,   2.514,  98.799, 95.082, 6.048,  59.287,
          48.889, 21.464, 54.972, 14.997, 27.161, 92.522, 65.383, NAN,
          NAN,    NAN,    1,      1,      1,      1,      1,      1,
          1,      1,      1,      1,      1,      1,      1,      1,
          1,      1,      1,      1,      1,      1,      1,      1,
          1,      1,      1,
      },
      25, 3);

  auto sizes = VarmaSizes(25, 2, 1, 1, 0, 1, 0, 0, 0, 0);

  auto pcaX = PcaAnalysisOptions();
  pcaX.ExactCount = 1;
  pcaX.IgnoreFirstCount = 1;
  auto pcaY = PcaAnalysisOptions();
  pcaY.ExactCount = 1;
  pcaY.IgnoreFirstCount = 1;

  // estimate
  auto var = VarmaExtended(sizes, VarmaRestrictionType::kMaFinal, true, true,
                           true, 3, nullptr, nullptr, nullptr);
  Tv *S = new Tv[var.StorageSize];
  Tv *W = new Tv[var.WorkSize];
  var.Calculate(data, S, W, true, 3, 0);
  delete[] S;
  delete[] W;

  // simulation
  auto horizon = 3;
  auto metrics = std::vector<ScoringType>(
      {ScoringType::kDirection, ScoringType::kRmspe, ScoringType::kCrps});
  auto horizons = std::vector<Ti>(horizon);
  std::iota(horizons.begin(), horizons.end(), 1);
  auto sim = VarmaSimulation(sizes, 1, horizons, metrics, nullptr, true);
  sim.KeepDetails = true;
  S = new Tv[sim.StorageSize];
  W = new Tv[sim.WorkSize];
  sim.CalculateE(S, W, data);

  delete[] S;
  delete[] W;
}

TEST(Varma_T, varma_sim) {
  auto ars = std::vector<Matrix<Tv> *>(
      {new Matrix<Tv>(new Tv[4]{0.2, -0.6, 0.3, 1.1}, 2, 2)});
  auto mas = std::vector<Matrix<Tv> *>(
      {new Matrix<Tv>(new Tv[4]{-0.5, 0, 0, -0.5}, 2, 2)}); // MA final form
  auto sigma = Matrix<Tv>(new Tv[4]{0.2, 0.0, 0.0, 0.4}, 2, 2);
  auto data = Varma::Simulate(&ars, &mas, {}, {}, &sigma, 1000, 100, 4888);
  const Matrix<Tv> y = std::get<0>(data);
  const Matrix<Tv> exo = std::get<1>(data);

  auto sizes = VarmaSizes(1000, 2, exo.ColsCount, 1, 0, 1, 0, 0, 0, 0);

  // restriction
  auto rest = VarmaRestriction(sizes, VarmaRestrictionType::kMaFinal);
  Tv *F = new Tv[rest.StorageSize];
  rest.Calculate(F);

  // estimate
  auto var = Varma(sizes, true, true, true);
  Tv *S = new Tv[var.Result.StorageSize];
  Tv *W = new Tv[var.Result.WorkSize];
  var.EstimateMl(y, &exo, W, S, &rest.R, nullptr, 0);

  ASSERT_NEAR(var.Result.gamma.Get(0), 0.2, 1e-1);
  ASSERT_NEAR(var.Result.gamma.Get(2), 0.4, 1e-1);

  delete[] W;

  auto forc = VarmaForecast(sizes, 5, true, false);
  Tv *S1 = new Tv[forc.StorageSize];
  W = new Tv[forc.WorkSize];
  forc.Calculate(var, &exo, &y, S1, W);

  // TODO:

  // ASSERT_NEAR(forc.Forecast.Get(0, 3), 10.1541927, 1e-5);

  delete[] S;
  delete[] W;
  delete[] S1;
}

TEST(Varma_T, varma_changeIndex) {
  Matrix<Tv> y0 =
      Matrix<Tv>(new Tv[50]{0.4, 0.9, 0.9, 0.6, 0.1, 0,   0.2, 1,   0.3, 0.4,
                            0.7, 0.7, 0.6, 0.9, 0.1, 0.8, 0.1, 0.2, 0.5, 0.4,
                            0.3, 0.7, 0.8, 0.6, 0.6, 0.8, 0.9, 0.8, 0.4, 0.3,
                            0.6, 0.3, 0.6, 0.5, 0.2, 0.6, 0.2, 0.7, 0.4, 0,
                            0.3, 0.4, 0.8, 0.4, 0.6, 0.5, 0.2, 0.7, 0.7, 0.9},
                 25, 2);
  Matrix<Tv> x0 =
      Matrix<Tv>(new Tv[50]{1,   0.4, 0.3, 0.7, 0.6, 0.6, 0.3, 0.9, 0.1, 0.1,
                            0.8, 0.7, 0.4, 0.8, 0.6, 0,   0.3, 0.9, 0.7, 0.6,
                            0.1, 0,   0.8, 0.1, 0.1, 0.7, 0.3, 0.1, 0.1, 0.9,
                            0.6, 0.8, 0.4, 0.8, 0.3, 0.5, 0.8, 0.6, 0.4, 0.8,
                            0.3, 0.7, 0.4, 0.2, 1,   0,   0,   0.5, 1,   0.4},
                 25, 2);
  Matrix<Tv> y1 =
      Matrix<Tv>(new Tv[50]{0.8, 0.9, 0.8, 0.4, 0.3, 0.6, 0.3, 0.6, 0.5, 0.2,
                            0.6, 0.2, 0.7, 0.4, 0,   0.3, 0.4, 0.8, 0.4, 0.6,
                            0.5, 0.2, 0.7, 0.7, 0.9, 0.4, 0.9, 0.9, 0.6, 0.1,
                            0,   0.2, 1,   0.3, 0.4, 0.7, 0.7, 0.6, 0.9, 0.1,
                            0.8, 0.1, 0.2, 0.5, 0.4, 0.3, 0.7, 0.8, 0.6, 0.6},
                 25, 2);
  Matrix<Tv> x1 =
      Matrix<Tv>(new Tv[50]{0.7, 0.3, 0.1, 0.1, 0.9, 0.6, 0.8, 0.4, 0.8, 0.3,
                            0.5, 0.8, 0.6, 0.4, 0.8, 0.3, 0.7, 0.4, 0.2, 1,
                            0,   0,   0.5, 1,   0.4, 1,   0.4, 0.3, 0.7, 0.6,
                            0.6, 0.3, 0.9, 0.1, 0.1, 0.8, 0.7, 0.4, 0.8, 0.6,
                            0,   0.3, 0.9, 0.7, 0.6, 0.1, 0,   0.8, 0.1, 0.1},
                 25, 2);
  y0.Transpose();
  x0.Transpose();
  y1.Transpose();
  x1.Transpose();
  auto sizes = VarmaSizes(25, 2, 2, 1, 0, 1, 0, 0, 0, 0);

  // restriction
  auto rest = VarmaRestriction(sizes, VarmaRestrictionType::kMaFinal);
  Tv *F = new Tv[rest.StorageSize];
  rest.Calculate(F);

  // estimate
  auto var0 = Varma(sizes, true, true, true);
  Tv *S0 = new Tv[var0.Result.StorageSize];
  Tv *W0 = new Tv[var0.Result.WorkSize];
  var0.EstimateMl(y0, &x0, W0, S0, &rest.R, nullptr, 0);

  auto var1 = Varma(sizes, true, true, true);
  Tv *S1 = new Tv[var1.Result.StorageSize];
  Tv *W1 = new Tv[var1.Result.WorkSize];
  var1.EstimateMl(y1, &x1, W1, S1, &rest.R, nullptr, 0);

  ASSERT_NEAR(var0.Result.LogLikelihood, var1.Result.LogLikelihood,
              1e-4); // !!!!?? 1e-4

  delete[] W0;
  delete[] S0;
  delete[] W1;
  delete[] S1;
}

TEST(Varma_T, outofsample) {
  auto y0 = Matrix<Tv>(new Tv[40]{1, 8, 7, 6, 5, 4, 4, 4, 1, 2, 3, 2, 2, 7,
                                  4, 3, 4, 3, 2, 4, 4, 7, 4, 2, 6, 0, 5, 9,
                                  3, 8, 6, 3, 1, 5, 3, 4, 6, 2, 2, 5},
                       20, 2);
  auto x0 =
      Matrix<Tv>(new Tv[50]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 3, 4, 5, 9, 4, 6, 6, 1, 5,
                            8, 7, 2, 2, 3, 4, 5, 4, 1, 3, 6, 4, 2, 1, 2, 3},
                 25, 2);
  Matrix<Tv> y = Matrix<Tv>(new Tv[40], 2, 20);
  Matrix<Tv> x = Matrix<Tv>(new Tv[50], 2, 25);
  y0.Transpose(y);
  x0.Transpose(x);

  auto sizes = VarmaSizes(20, 2, 2, 1, 0, 0, 0, 0, 0, 0);

  auto horizons = std::vector<Ti>({1, 3});
  Ti outcount = 4;
  // bool keepdetails = true;
  auto metrics = std::vector<ScoringType>({ScoringType::kRmse});
  auto sim = VarmaSimulation(sizes, outcount, horizons, metrics, nullptr);
  sim.KeepDetails = true;
  auto W = new Tv[sim.WorkSize];
  auto S = new Tv[sim.StorageSize];
  bool cancel = false;
  sim.Calculate(S, W, y, cancel, &x);

  /*
  // 16
  ASSERT_NEAR(result.details.at(0).sampleEnd, (Ti)4);
  ASSERT_NEAR(result.details.at(0).forecast.Get(0, 1), 4.1201140, 1e-5);
  ASSERT_NEAR(result.details.at(0).actual.Get(1, 0), 6, 1e-6);
  ASSERT_NEAR(result.details.at(0).metrics.at(0).Get(1, 1), 5.9865124, 1e-5);

  // 17
  ASSERT_NEAR(result.details.at(1).sampleEnd, (Ti)3);
  ASSERT_NEAR(result.details.at(1).forecast.Get(1, 1), 4.35570892, 1e-5);
  ASSERT_NEAR(result.details.at(1).actual.Get(0, 1), 4, 1e-6);
  ASSERT_NEAR(result.details.at(1).metrics.at(0).Get(0, 1), 0.00543269,
  1e-5);

  // 19
  ASSERT_NEAR(result.details.at(3).sampleEnd, (Ti)1);
  ASSERT_NEAR(result.details.at(3).forecast.Get(0, 0), 4.03157163, 1e-5);
  ASSERT_NEAR(result.details.at(3).actual.Get(1, 0), 5, 1e-6);
  ASSERT_NEAR(result.details.at(3).metrics.at(0).Get(1, 0), 0.13348641,
  1e-5);
  */

  // average (second variable,1 period ahead)
  auto sum1 = 1.295428409821386 + 5.869232939915652 + 9.480065543353192 +
              0.1334864116986188;
  // ASSERT_NEAR(std::sqrt(sum1 / 4), sim.Results.at(0).Get(1, 0), 1e-13);
  // average (second variable, second h, which is 3 periods ahead)
  // auto sum2 = 5.986512374677352 + 0.4151109869193009;
  // ASSERT_NEAR(std::sqrt(sum2 / 2), sim.Results.at(0).Get(1, 1), 1e-13);
}

TEST(Varma_T, outofsample2) { // similar to the previous test, just change the
                              // metrics & horizons length

  auto y0 = Matrix<Tv>(new Tv[40]{1, 8, 7, 6, 5, 4, 4, 4, 1, 2, 3, 2, 2, 7,
                                  4, 3, 4, 3, 2, 4, 4, 7, 4, 2, 6, 0, 5, 9,
                                  3, 8, 6, 3, 1, 5, 3, 4, 6, 2, 2, 5},
                       20, 2);
  auto x0 =
      Matrix<Tv>(new Tv[50]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 3, 4, 5, 9, 4, 6, 6, 1, 5,
                            8, 7, 2, 2, 3, 4, 5, 4, 1, 3, 6, 4, 2, 1, 2, 3},
                 25, 2);
  Matrix<Tv> y = Matrix<Tv>(new Tv[40], 2, 20);
  Matrix<Tv> x = Matrix<Tv>(new Tv[50], 2, 25);
  y0.Transpose(y);
  x0.Transpose(x);

  auto sizes = VarmaSizes(20, 2, 2, 1, 0, 0, 0, 0, 0, 0);

  auto horizons = std::vector<Ti>({1, 2, 3});
  Ti outcount = 4;
  // bool keepdetails = true;
  auto metrics =
      std::vector<ScoringType>({ScoringType::kCrps, ScoringType::kRmse});
  auto sim = VarmaSimulation(sizes, outcount, horizons, metrics);
  auto W = new Tv[sim.WorkSize];
  auto S = new Tv[sim.StorageSize];
  bool cancel = false;
  sim.Calculate(S, W, y, cancel, &x, nullptr, nullptr);

  // average (second variable,1 period ahead)
  auto sum1 = 1.295428409821386 + 5.869232939915652 + 9.480065543353192 +
              0.1334864116986188;
  // ASSERT_NEAR(std::sqrt(sum1 / 4), sim.Results.at(1).Get(1, 0), 1e-13);
  // average (second variable, second h, which is 3 periods ahead)
  // auto sum2 = 5.986512374677352 + 0.4151109869193009;
  // ASSERT_NEAR(std::sqrt(sum2 / 2), sim.Results.at(1).Get(1, 2), 1e-13);
}

TEST(Varma_T, outofsample_compare) { // row-wise and extended

  auto data =
      Matrix<Tv>(new Tv[100]{1, 8, 7, 6, 5, 4, 4,   4,   1,   2,   3,   2,   2,
                             7, 4, 3, 4, 3, 2, 4,   NAN, NAN, NAN, NAN, NAN, 4,
                             7, 4, 2, 6, 0, 5, 9,   3,   8,   6,   3,   1,   5,
                             3, 4, 6, 2, 2, 5, NAN, NAN, NAN, NAN, NAN, 1,   1,
                             1, 1, 1, 1, 1, 1, 1,   1,   1,   1,   1,   1,   1,
                             1, 1, 1, 1, 1, 1, 1,   1,   1,   1,   3,   4,   5,
                             9, 4, 6, 6, 1, 5, 8,   7,   2,   2,   3,   4,   5,
                             4, 1, 3, 6, 4, 2, 1,   2,   3},
                 25, 4);
  auto y = Matrix<Tv>(new Tv[40]{1, 8, 7, 6, 5, 4, 4, 4, 1, 2, 3, 2, 2, 7,
                                 4, 3, 4, 3, 2, 4, 4, 7, 4, 2, 6, 0, 5, 9,
                                 3, 8, 6, 3, 1, 5, 3, 4, 6, 2, 2, 5},
                      20, 2);
  auto x =
      Matrix<Tv>(new Tv[50]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 3, 4, 5, 9, 4, 6, 6, 1, 5,
                            8, 7, 2, 2, 3, 4, 5, 4, 1, 3, 6, 4, 2, 1, 2, 3},
                 25, 2);
  y.Transpose();
  x.Transpose();

  auto sizes = VarmaSizes(20, 2, 2, 1, 0, 0, 0, 0, 0, 0);

  auto horizons = std::vector<Ti>({1});
  Ti outcount = 1;
  // bool keepdetails = true;
  auto metrics =
      std::vector<ScoringType>({ScoringType::kCrps, ScoringType::kRmse});

  auto sim_ex =
      VarmaSimulation(sizes, outcount, horizons, metrics, nullptr, true);
  sim_ex.KeepDetails = true;
  auto W = new Tv[sim_ex.WorkSize];
  auto S = new Tv[sim_ex.StorageSize];
  sim_ex.CalculateE(S, W, data, 1e12, 2, false, true);

  auto sim_r =
      VarmaSimulation(sizes, outcount, horizons, metrics, nullptr, false);
  auto W_r = new Tv[sim_r.WorkSize];
  auto S_r = new Tv[sim_r.StorageSize];
  auto restriction = VarmaRestriction(sizes, VarmaRestrictionType::kMaFinal);
  restriction.Calculate(new Tv[restriction.StorageSize]);
  bool cancel = false;
  sim_r.Calculate(S_r, W_r, y, cancel, &x, &restriction.R, nullptr, true, 1e12,
                  2, false);

  ASSERT_EQ(true, sim_ex.ResultAggs.Equals(sim_r.ResultAggs, 1e-14));
}

TEST(Varma_T, outofsample_dir) {
  auto y0 = Matrix<Tv>(new Tv[54]{7, 6, 5, 4, 4, 4, 1, 2, 3, 2, 2, 7, 4, 3,
                                  4, 3, 2, 4, 4, 2, 6, 0, 5, 9, 3, 8, 6, 3,
                                  1, 5, 3, 4, 6, 2, 2, 5, 5, 9, 4, 6, 6, 1,
                                  5, 8, 7, 2, 2, 3, 4, 5, 4, 1, 3, 6},
                       18, 3);
  auto x0 = Matrix<Tv>(new Tv[46]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 7, 7, 3, 9,
                                  3, 8, 4, 3, 5, 4, 2, 3, 2, 4, 7, 2, 1, 2, 6},
                       23, 2);
  Matrix<Tv> y = Matrix<Tv>(new Tv[54], 3, 18);
  Matrix<Tv> x = Matrix<Tv>(new Tv[46], 2, 23);
  y0.Transpose(y);
  x0.Transpose(x);

  auto sizes = VarmaSizes(18, 3, 2, 1, 0, 0, 0, 0, 0, 0);

  auto horizons = std::vector<Ti>({1, 3});
  Ti outcount = 4;
  // bool keepdetails = true;
  auto metrics = std::vector<ScoringType>(
      {ScoringType::kCrps, ScoringType::kDirection, ScoringType::kRmse});
  auto sim = VarmaSimulation(sizes, outcount, horizons, metrics);
  sim.KeepDetails = true;
  auto W = new Tv[sim.WorkSize];
  auto S = new Tv[sim.StorageSize];
  bool cancel = false;
  sim.Calculate(S, W, y, cancel, &x, nullptr, nullptr);

  /*
   auto y1_15 = y.Get(0, y.ColsCount - 5);
   auto y1_16 = y.Get(0, y.ColsCount - 4);
   auto y1_17 = y.Get(0, y.ColsCount - 3);
   auto y1_18 = y.Get(0, y.ColsCount - 2);
   auto y1_19 = y.Get(0, y.ColsCount - 1);

   // first variable
   auto f1_h1_2016 = result.details.at(0).forecast.Get(0, 0);
   auto f1_h1_2017 = result.details.at(1).forecast.Get(0, 0);
   auto f1_h1_2018 = result.details.at(2).forecast.Get(0, 0);
   auto f1_h1_2019 = result.details.at(3).forecast.Get(0, 0);

   Tv v1 = ((f1_h1_2016 < y1_15&& y1_16 < y1_15) ? 1 : 0) +
           ((f1_h1_2016 > y1_15 && y1_16 > y1_15) ? 1 : 0);
   Tv v2 = ((f1_h1_2017 < y1_16&& y1_17 < y1_16) ? 1 : 0) +
           ((f1_h1_2017 > y1_16 && y1_17 > y1_16) ? 1 : 0);
   Tv v3 = ((f1_h1_2018 < y1_17&& y1_18 < y1_17) ? 1 : 0) +
           ((f1_h1_2018 > y1_17 && y1_18 > y1_17) ? 1 : 0);
   Tv v4 = ((f1_h1_2019 < y1_18&& y1_19 < y1_18) ? 1 : 0) +
           ((f1_h1_2019 > y1_18 && y1_19 > y1_18) ? 1 : 0);

   ASSERT_NEAR(result.details.at(0).metrics.at(1).Get(0, 0), v1);
   ASSERT_NEAR(result.details.at(1).metrics.at(1).Get(0, 0), v2);
   ASSERT_NEAR(result.details.at(2).metrics.at(1).Get(0, 0), v3);
   ASSERT_NEAR(result.details.at(3).metrics.at(1).Get(0, 0), v4);
   */
  // ASSERT_NEAR(sim.Results.at(1).Get(0, 0), 0.750, 1e-16);

  // second variable
  /*
  auto y2_15 = y.Get(1, y.ColsCount - 5);
  auto y2_16 = y.Get(1, y.ColsCount - 4);
  auto y2_17 = y.Get(1, y.ColsCount - 3);
  auto y2_18 = y.Get(1, y.ColsCount - 2);
  auto y2_19 = y.Get(1, y.ColsCount - 1);

  auto f2_h1_2016 = result.details.at(0).forecast.Get(1, 0);
  auto f2_h1_2017 = result.details.at(1).forecast.Get(1, 0);
  auto f2_h1_2018 = result.details.at(2).forecast.Get(1, 0);
  auto f2_h1_2019 = result.details.at(3).forecast.Get(1, 0);


  v1 = ((f2_h1_2016 < y2_15&& y2_16 < y2_15) ? 1 : 0) +
          ((f2_h1_2016 > y2_15 && y2_16 > y2_15) ? 1 : 0);
  v2 = ((f2_h1_2017 < y2_16&& y2_17 < y2_16) ? 1 : 0) +
          ((f2_h1_2017 > y2_16 && y2_17 > y2_16) ? 1 : 0);
  v3 = ((f2_h1_2018 < y2_17&& y2_18 < y2_17) ? 1 : 0) +
          ((f2_h1_2018 > y2_17 && y2_18 > y2_17) ? 1 : 0);
  v4 = ((f2_h1_2019 < y2_18&& y2_19 < y2_18) ? 1 : 0) +
          ((f2_h1_2019 > y2_18 && y2_19 > y2_18) ? 1 : 0);

  ASSERT_NEAR(result.details.at(0).metrics.at(1).Get(1, 0), v1);
  ASSERT_NEAR(result.details.at(1).metrics.at(1).Get(1, 0), v2);
  ASSERT_NEAR(result.details.at(2).metrics.at(1).Get(1, 0), v3);
  ASSERT_NEAR(result.details.at(3).metrics.at(1).Get(1, 0), v4);
  */
  // ASSERT_NEAR(sim.Results.at(1).Get(1, 0), 0.750, 1e-16);

  // Horizon 3
  // third variable
  /*
  auto y3_15 = y.Get(2, y.ColsCount - 5);
  auto y3_16 = y.Get(2, y.ColsCount - 4);
  auto y3_17 = y.Get(2, y.ColsCount - 3);
  auto y3_18 = y.Get(2, y.ColsCount - 2);
  auto y3_19 = y.Get(2, y.ColsCount - 1);

  auto f3_h3_2018 = result.details.at(0).forecast.Get(2, 1);
  auto f3_h3_2019 = result.details.at(1).forecast.Get(2, 1);

  v3 = ((f3_h3_2018 < y3_17&& y3_18 < y3_17) ? 1 : 0) +
          ((f3_h3_2018 > y3_17 && y3_18 > y3_17) ? 1 : 0);
  v4 = ((f3_h3_2019 < y3_18&& y3_19 < y3_18) ? 1 : 0) +
          ((f3_h3_2019 > y3_18 && y3_19 > y3_18) ? 1 : 0);

  ASSERT_NEAR(result.details.at(0).metrics.at(1).Get(2, 1), v1);
  ASSERT_NEAR(result.details.at(1).metrics.at(1).Get(2, 1), v2);
  */
  // ASSERT_NEAR(sim.Results.at(1).Get(2, 1), 1.0, 1e-16);
}

TEST(Varma_T, var_search) {
  Matrix<Tv> data0 = Matrix<Tv>(
      new Tv[175]{
          99.72, 97.92, 56.54, 83.42, 92.26, 70.87, 76.07, 4.4,   80.39, 17.36,
          15.97, 8.2,   63.25, 22.21, 73.89, 97.61, 25.91, 77.64, 1.99,  47.91,
          NAN,   NAN,   NAN,   NAN,   NAN,   22.69, 85.4,  86.52, 12.38, 30.98,
          55.28, 82.9,  91.49, 94.46, 25.23, 17.81, 68.72, 51.19, 58.44, 87.2,
          78.37, 46.01, 78.36, 23.55, 19.98, NAN,   NAN,   NAN,   NAN,   NAN,
          29.87, 87.62, 93.55, 22.62, 56.17, 42.02, 17.68, 7.23,  44.26, 12.48,
          59.15, 97.62, 19.1,  97.52, 91.92, 9.64,  50.57, 69.87, 34.69, 41.93,
          NAN,   NAN,   NAN,   NAN,   NAN,   66.65, 54.83, 78.89, 48.38, 90.27,
          1.06,  5.73,  69.28, 43.3,  59.16, 93.13, 18.33, 67.54, 39.53, 59.2,
          81.9,  58.92, 61.82, 47.84, 84.39, 2,     3,     4,     5,     6,
          62.17, 27.3,  82.17, 27.67, 38.42, 62.38, 50.52, 82.43, 36.42, 90.27,
          92.14, 92.19, 83.92, 54.37, 70.26, 39.16, 28.02, 2.9,   32.99, 51.88,
          7,     8,     6,     3,     5,     79.18, 99.12, 97.16, 44.54, 10.64,
          95.93, 57.89, 17.41, 34.78, 58.1,  99.85, 64.26, 35.68, 33.85, 28.57,
          39.84, 15.76, 94.6,  48.69, 5.85,  11,    12,    13,    14,    15,
          1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
          1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
          1,     1,     1,     1,     1},
      25, 7);
  data0.Transpose();
  auto data = DatasetTs<true>(data0.RowsCount, data0.ColsCount, true, true);
  data.Data(data0);

  auto items = SearchItems();
  auto searchOptions = SearchOptions();
  auto metrics = SearchMetricOptions();
  auto checks = SearchModelChecks();

  items.KeepBestCount = 2;
  items.KeepAll = true;
  items.LengthTargets = 2;
  items.LengthDependents = 3;
  items.KeepMixture = true;

  items.Length1 = 2;
  checks.Prediction = true;

  auto sizes = std::vector<Ti>({1});
  auto exo = std::vector<std::vector<Ti>>({{6}});
  auto endogroups = std::vector<std::vector<Ti>>({{0}, {1}, {2}});
  metrics.SimFixSize = 2;
  metrics.Horizons = std::vector<Ti>({1});
  metrics.MetricsIn = std::vector<GoodnessOfFitType>(
      {GoodnessOfFitType::kAic, GoodnessOfFitType::kSic});
  metrics.MetricsOut =
      std::vector<ScoringType>({ScoringType::kDirection, ScoringType::kRmse});
  metrics.TrainFixSize = 4;
  auto parm = std::vector<Ti>({2, 1, 2, 0, 0, 0});

  auto modelset =
      VarmaModelset(searchOptions, items, metrics, checks, sizes, endogroups,
                    data, parm, 0, exo, true, nullptr, 2, items.Length1);
  auto W = new Tv[modelset.Modelset.WorkSize];
  modelset.Modelset.Start(W, nullptr);

  auto result = SearcherModelingInfo();
  auto list0 = std::vector<SearcherSummary *>();
  auto list1 = std::vector<SearcherSummary *>();
  auto list2 = std::vector<SearcherSummary *>();
  modelset.Modelset.CombineInfo(result, list0, list1, list2);
  auto mixture = RunningMoments<4, true, true, Tv>();
  modelset.Modelset.CombineMixture(0, 0, 1, list1, mixture);

  delete[] W;
}