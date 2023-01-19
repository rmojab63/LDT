/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "functionx.h"

#include "matrix.h"
#include <gtest/gtest.h>
#include <random>

using namespace ldt;

TEST(Functionx_T, simple) {
  std::function<Tv(const Matrix<Tv> &)> f = [](const Matrix<Tv> &x) -> Tv {
    auto storage = Matrix<Tv>(new std::vector<Tv>(x.length()), x.length());
    x.Apply([](Tv a) -> Tv { return (a - 1) * (a - 1); }, storage);
    return storage.Sum();
  };
  std::function<void(Matrix<Tv> &, Matrix<Tv> &)> g =
      [](Matrix<Tv> &x, Matrix<Tv> &storage) -> void {
    x.Apply([](Tv a) -> Tv { return 2 * (a - 1); }, storage);
  };

  Ti m = 3;
  Matrix<Tv> x0 = Matrix<Tv>(new std::vector<Tv>{1.0, 3.0, 4.0}, m);
  Matrix<Tv> actual = Matrix<Tv>(new std::vector<Tv>(m), m);
  g(x0, actual);

  auto derf = Derivative(m, true, false);

  Tv *W0 = new Tv[derf.WorkSize];
  Tv *S0 = new Tv[derf.StorageSize1];

  derf.CalculateFirst(f, x0, S0, W0);
  ASSERT_EQ(true, actual.Equals(derf.Result1, 1e-7));

  delete[] W0;
  delete[] S0;
}

TEST(Functionx_T, rosenbrock_der) {
  std::function<Tv(const Matrix<Tv> &)> f = [](const Matrix<Tv> &x) -> Tv {
    auto x0 = x.Data[0];
    auto x1 = x.Data[1];
    return (1 - x0) * (1 - x0) + 100 * (x1 - x0 * x0) * (x1 - x0 * x0);
  };
  std::function<void(Matrix<Tv> &, Matrix<Tv> &)> g =
      [](Matrix<Tv> &x, Matrix<Tv> &storage) -> void {
    auto x0 = x.Data[0];
    auto x1 = x.Data[1];
    storage.Data[0] = -2 * (1 - x0) - 400 * x0 * (x1 - x0 * x0);
    storage.Data[1] = 200 * (x1 - x0 * x0);
  };

  Ti m = 2;
  Matrix<Tv> x0 = Matrix<Tv>(new std::vector<Tv>({-1.0, 1.0}), m);
  Matrix<Tv> actual = Matrix<Tv>(new std::vector<Tv>(m), m);
  g(x0, actual);

  auto derf = Derivative(m, true, false);

  Tv *W0 = new Tv[derf.WorkSize];
  Tv *S0 = new Tv[derf.StorageSize1];

  derf.CalculateFirst(f, x0, S0, W0);
  ASSERT_EQ(true, actual.Equals(derf.Result1, 1e-7));

  delete[] W0;
  delete[] S0;
}

TEST(Functionx_T, simple_hes) {
  std::function<Tv(const Matrix<Tv> &)> f = [](const Matrix<Tv> &x) -> Tv {
    auto storage = Matrix<Tv>(new std::vector<Tv>(x.length()), x.length());
    x.Apply([](Tv a) -> Tv { return (a - 1) * (a - 1); }, storage);
    return storage.Sum();
  };

  Ti m = 3;
  Matrix<Tv> x0 = Matrix<Tv>(new std::vector<Tv>({1.0, 3.0, 4.0}), m);
  Matrix<Tv> actual = Matrix<Tv>(0.0, new Tv[m * m], m, m);
  Matrix<Tv>::Diagonal(actual, 2.0);

  auto derf = Derivative(m, false, true);

  Tv *W0 = new Tv[derf.WorkSize];
  Tv *S0 = new Tv[derf.StorageSize2];

  derf.CalculateSecond(f, x0, S0, W0);
  ASSERT_EQ(true, actual.Equals(derf.Result2, 1e-5));

  delete[] W0;
  delete[] S0;
}

TEST(Functionx_T, simple_hes2) {
  std::function<Tv(const Matrix<Tv> &)> f = [](const Matrix<Tv> &a) -> Tv {
    auto x = a.Data[0];
    auto y = a.Data[1];
    return x * x + 2 * x * x * y * y + 2 * y * y * y;
  };

  std::function<void(Matrix<Tv> &, Matrix<Tv> &)> hes =
      [](Matrix<Tv> &a, Matrix<Tv> &storage) -> void {
    auto x = a.Data[0];
    auto y = a.Data[1];
    auto fxx = 2 + 4 * y * y;
    auto fyx = 8 * x * y;
    auto fyy = 12 * y + 4 * x * x;
    storage.Set0(0, 0, fxx);
    storage.Set0(1, 0, fyx);
    storage.Set0(0, 1, fyx);
    storage.Set0(1, 1, fyy);
  };

  Ti m = 2;
  auto z = Matrix<Tv>(new std::vector<Tv>{1.0, 2.0}, m);
  Matrix<Tv> h = Matrix<Tv>(0.0, new Tv[4]{0, 0, 0, 0}, m, m);
  hes(z, h);

  auto derf = Derivative(m, false, true);

  Tv *W0 = new Tv[derf.WorkSize];
  Tv *S0 = new Tv[derf.StorageSize2];

  derf.CalculateSecond(f, z, S0, W0);

  ASSERT_EQ(true, h.Equals(derf.Result2, 1e-4));

  delete[] W0;
  delete[] S0;
}

TEST(Functionx_T, mle_hes) {
  auto n = 10000;
  auto y = Matrix<Tv>(new std::vector<Tv>(n), n);
  Matrix<Tv>::FillRandom_normal(y, 3403, 0.0, 2.0);

  std::function<Tv(const Matrix<Tv> &)> f = [&y](const Matrix<Tv> &co) -> Tv {
    auto s2 = co.Data[1];
    auto m = co.Data[0];
    Tv sum = 0.0;
    for (Ti i = 0; i < y.length(); i++) {
      auto pdf = (1 / (sqrt(2 * c_pi * s2))) *
                 exp(-0.5 * std::pow((y.Data[i] - m), 2.0) / s2);
      sum += std::log(pdf);
    }
    return -sum;
  };

  Ti m = 2;
  auto derf = Derivative(m, false, true);

  Tv *W0 = new Tv[derf.WorkSize];
  Tv *S0 = new Tv[derf.StorageSize2];

  auto z = Matrix<Tv>(new std::vector<Tv>{0.1, 2.2}, m);

  derf.CalculateSecond(f, z, S0, W0);

  // ASSERT_EQ(true,derf.Result2.Equals(Matrix<Tv>(new std::vector<Tv>{
  // 456,-14,-14,249 }, 2, 2), 2));
}