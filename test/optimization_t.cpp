/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "functionx.h"
#include "matrix.h"
#include "optimization.h"
#include <algorithm>
#include <gtest/gtest.h>
#include <iostream>
// #include <ppl.h>

using namespace ldt;

TEST(LMBFGSB_T, simple) {
  std::function<Tv(const Matrix<Tv> &)> f = [](const Matrix<Tv> &x) -> Tv {
    auto storage = Matrix<Tv>(new std::vector<Tv>(x.length()), x.length());
    std::function<Tv(Tv)> f = [](Tv a) -> Tv { return (a - 1) * (a - 1); };
    x.Apply(f, storage);
    return storage.Sum();
  };
  std::function<void(const Matrix<Tv> &, Matrix<Tv> &)> g =
      [](const Matrix<Tv> &x, Matrix<Tv> &storage) -> void {
    std::function<Tv(Tv)> f = [](Tv a) -> Tv { return 2 * (a - 1); };
    x.Apply(f, storage);
  };

  Ti m = 3;
  Matrix<Tv> x0 = Matrix<Tv>(new Tv[3]{2, 2, 2}, m);
  auto prob = LimitedMemoryBFGSB(m);

  Tv *W = new Tv[prob.WorkSize];
  Tv *S = new Tv[prob.StorageSize];

  prob.Minimize(f, g, x0, S, W);
  ASSERT_EQ(true, Matrix<Tv>(new Tv[3]{1.0, 1.0, 1.0}, 3).Equals(x0, 1e-8));
  ASSERT_NEAR(0, prob.FunctionValue, 1e-8);

  auto lower = Matrix<Tv>(new Tv[3]{6.0, -2, -2}, 3);
  auto upper = Matrix<Tv>(new Tv[3]{6.0, INFINITY, INFINITY}, 3);

  prob.Minimize(f, g, x0, S, W, &lower, &upper);
  ASSERT_EQ(true, Matrix<Tv>(new Tv[3]{6.0, 1.0, 1.0}, 3).Equals(x0, 1e-8));

  delete[] W;
  delete[] S;
  delete[] lower.Data;
  delete[] upper.Data;
  delete[] x0.Data;
}

TEST(LMBFGSB_T, rosenbrock_der) {
  std::function<Tv(const Matrix<Tv> &)> f = [](const Matrix<Tv> &x) -> Tv {
    auto x0 = x.Data[0];
    auto x1 = x.Data[1];
    return (1 - x0) * (1 - x0) + 100 * (x1 - x0 * x0) * (x1 - x0 * x0);
  };
  std::function<void(const Matrix<Tv> &, Matrix<Tv> &)> g =
      [](const Matrix<Tv> &x, Matrix<Tv> &storage) -> void {
    auto x0 = x.Data[0];
    auto x1 = x.Data[1];
    storage.Data[0] = -2 * (1 - x0) - 400 * x0 * (x1 - x0 * x0);
    storage.Data[1] = 200 * (x1 - x0 * x0);
  };
  Ti m = 2;
  Matrix<Tv> x0 = Matrix<Tv>(new Tv[2]{-3, 4}, m);

  auto prob = LimitedMemoryBFGSB(m);
  prob.Options.IterationMax = 1000000;
  prob.Options.Factor = 1e-10;
  prob.Options.ProjectedGradientTol = 1e-10;

  Tv *W = new Tv[prob.WorkSize];
  Tv *S = new Tv[prob.StorageSize];

  prob.Minimize(f, g, x0, S, W);
  ASSERT_EQ(true, Matrix<Tv>(new Tv[2]{1.0, 1.0}, 2).Equals(x0, 1e-8));
  ASSERT_NEAR(0, prob.FunctionValue, 1e-8);

  delete[] W;
  delete[] S;
  delete[] x0.Data;
}

TEST(LMBFGSB_T, rosenbrock) {
  std::function<Tv(const Matrix<Tv> &)> f = [](const Matrix<Tv> &x) -> Tv {
    auto x0 = x.Data[0];
    auto x1 = x.Data[1];
    return (1 - x0) * (1 - x0) + 100 * (x1 - x0 * x0) * (x1 - x0 * x0);
  };
  Ti m = 2;

  // numerical derivative
  auto derf = Derivative(m, true, true);
  auto W = new Tv[derf.WorkSize];
  auto S = new Tv[derf.StorageSize1];
  std::function<void(const Matrix<Tv> &, Matrix<Tv> &)> g =
      [&derf, &f, W](const Matrix<Tv> &x, Matrix<Tv> &storage) -> void {
    derf.CalculateFirst(f, x, storage.Data, W);
  };

  Matrix<Tv> x0 = Matrix<Tv>(new Tv[2]{-3, 4}, m);
  auto prob = LimitedMemoryBFGSB(m);
  prob.Options.IterationMax = 1000000;
  prob.Options.Factor = 1e-10;
  prob.Options.ProjectedGradientTol = 1e-10;

  Tv *W1 = new Tv[prob.WorkSize];
  Tv *S1 = new Tv[prob.StorageSize];

  prob.Minimize(f, g, x0, S1, W1);
  ASSERT_EQ(true, Matrix<Tv>(new Tv[2]{1.0, 1.0}, 2).Equals(x0, 1e-8));
  ASSERT_NEAR(0, prob.FunctionValue, 1e-8);

  delete[] W1;
  delete[] W;
  delete[] S1;
  delete[] S;
}

/*
TEST(LMBFGSB_T, rosenbrock_der_parallel)
{
        concurrency::parallel_for(Ti(0), Ti(100), [&](Ti i)
                {

                        std::function<Tv(const Matrix<Tv>&)> f = [](const
Matrix<Tv>& x) -> Tv { auto x0 = x.Data[0]; auto x1 = x.Data[1]; return	(1 - x0)
* (1 - x0) + 100 * (x1 - x0 * x0) * (x1 - x0 * x0);
                        };
                        std::function<void(const Matrix<Tv>&, Matrix<Tv>&)> g =
[](const Matrix<Tv>& x, Matrix<Tv>& storage) -> void { auto x0 = x.Data[0]; auto
x1 = x.Data[1]; storage.Data[0] = -2 * (1 - x0) - 400 * x0 * (x1 - x0 * x0);
                                storage.Data[1] = 200 * (x1 - x0 * x0);
                        };


                        Ti m = 2;
                        Matrix<Tv> x0 = Matrix<Tv>(new Tv[2] {-3, 4}, m);
                        auto prob = LimitedMemoryBFGSB(m);
                        prob.Options.IterationMax = 1000000;
                        prob.Options.Factor = 1e-10;
                        prob.Options.ProjectedGradientTol = 1e-10;

                        Tv* W = new Tv[prob.WorkSize];
                        Tv* S = new Tv[prob.StorageSize];

                        prob.Minimize(f, g, x0, S, W);
                        ASSERT_EQ(true,Matrix<Tv>(new Tv[2] {1.0, 1.0},
2).Equals(x0, 1e-8)); ASSERT_NEAR(0, prob.FunctionValue, 1e-8);


                });
}
*/

TEST(Newton_T, simple) {
  std::function<Tv(const Matrix<Tv> &)> f = [](const Matrix<Tv> &x) -> Tv {
    auto storage = Matrix<Tv>(new std::vector<Tv>(x.length()), x.length());
    std::function<Tv(Tv)> f = [](Tv a) -> Tv { return (a - 1) * (a - 1); };
    x.Apply(f, storage);
    return storage.Sum();
  };
  std::function<void(const Matrix<Tv> &, Matrix<Tv> &)> g =
      [](const Matrix<Tv> &x, Matrix<Tv> &storage) -> void {
    std::function<Tv(Tv)> f = [](Tv a) -> Tv { return 2 * (a - 1); };
    x.Apply(f, storage);
  };
  std::function<void(const Matrix<Tv> &, Matrix<Tv> &)> h =
      [](const Matrix<Tv> &x, Matrix<Tv> &storage) -> void {
    Matrix<Tv>::Diagonal(storage, 2.0, 0.0);
  };
  Ti m = 3;
  Matrix<Tv> x0 = Matrix<Tv>(new Tv[3]{2, 2, 2}, m);
  Newton prob = Newton(m);
  prob.IterationMax = 1000;
  prob.TolFunction = 1e-8;

  auto W = new Tv[prob.WorkSize];
  auto S = new Tv[prob.StorageSize];

  prob.Minimize(f, g, h, x0, S, W);
  ASSERT_EQ(true, Matrix<Tv>(new Tv[3]{1.0, 1.0, 1.0}, 3).Equals(x0, 1e-8));
  ASSERT_NEAR(0.0, prob.FunctionValue, 1e-8);

  delete[] W;
  delete[] S;
}

TEST(Newton_T, rosenbrock) {
  std::function<Tv(const Matrix<Tv> &)> f = [](const Matrix<Tv> &x) -> Tv {
    auto x0 = x.Data[0];
    auto x1 = x.Data[1];
    return (1 - x0) * (1 - x0) + 100 * (x1 - x0 * x0) * (x1 - x0 * x0);
  };
  std::function<void(const Matrix<Tv> &, Matrix<Tv> &)> g =
      [](const Matrix<Tv> &x, Matrix<Tv> &storage) -> void {
    auto x0 = x.Data[0];
    auto x1 = x.Data[1];
    storage.Data[0] = -2 * (1 - x0) - 400 * x0 * (x1 - x0 * x0);
    storage.Data[1] = 200 * (x1 - x0 * x0);
  };

  std::function<void(const Matrix<Tv> &, Matrix<Tv> &)> h =
      [](const Matrix<Tv> &x, Matrix<Tv> &storage) -> void {
    auto x0 = x.Data[0];
    auto x1 = x.Data[1];
    storage.Set(0, 0, 2 - 400 * x1 + 1200 * x0 * x0);
    storage.Set(1, 1, 200);
    storage.Set(0, 1, -400 * x0);
    storage.Set(1, 0, -400 * x0);
  };
  Ti m = 2;
  Matrix<Tv> x0 = Matrix<Tv>(new Tv[2]{2, 2}, m);
  Newton prob = Newton(m);
  prob.IterationMax = 1000;
  prob.TolFunction = 1e-8;

  auto W = new Tv[prob.WorkSize];
  auto S = new Tv[prob.StorageSize];

  prob.Minimize(f, g, h, x0, S, W);

  ASSERT_EQ(true, Matrix<Tv>(new Tv[2]{1.0, 1.0}, 2).Equals(x0, 1e-7));
  ASSERT_NEAR(0.0, prob.FunctionValue, 1e-8);

  delete[] W;
  delete[] S;
}

TEST(NelderMead_T, univariate) {
  auto quandratic = [](Tv x) { return std::pow(x - 2, 2); };
  auto res = NelderMead::Minimize1(quandratic, 0.0);
  ASSERT_NEAR(2.0, std::get<0>(res), 1e-6);

  auto absf = [](Tv x) { return std::abs(x - 4.5); };
  res = NelderMead::Minimize1(absf, 0.0);
  ASSERT_NEAR(4.5, std::get<0>(res), 1e-6);

  auto absf_sin = [](Tv x) { return std::abs(std::sin(x)); };

  res = NelderMead::Minimize1(absf_sin, M_PI / 2);
  ASSERT_NEAR(M_PI, std::get<0>(res), 1e-4);

  res =
      NelderMead::Minimize1(absf_sin, 0.0, 0.1, 100, 1e-7, M_PI, 3 * M_PI / 2);
  ASSERT_NEAR(M_PI, std::get<0>(res), 1e-4);
}

TEST(NelderMead_T, multivar) {

  auto optim = NelderMead(3);
  auto start = Matrix<Tv>(new double[3]{2.0, 2.0, 2.0}, 3, 1);
  std::function<double(const Matrix<Tv> &)> sphere =
      [](const Matrix<Tv> &x) -> double {
    double sum = 0.0;
    for (size_t i = 0; i < x.length(); ++i)
      sum += x.Data[i] * x.Data[i];
    return sum;
  };
  auto work = new double[optim.WorkSize];
  auto storage = new double[optim.StorageSize];
  optim.Minimize(sphere, start, work, storage);

  ASSERT_TRUE(std::abs(optim.Result.Data[0]) < 1e-3);
  ASSERT_TRUE(std::abs(optim.Result.Data[1]) < 1e-3);
  ASSERT_TRUE(std::abs(optim.Result.Data[2]) < 1e-3);

  // restrict
  auto lower = Matrix<Tv>(new double[3]{-4, 2.0, -INFINITY}, 3, 1);
  auto upper = Matrix<Tv>(new double[3]{4, 5.0, INFINITY}, 3, 1);
  optim.Minimize(sphere, start, work, storage, &lower, &upper);
}