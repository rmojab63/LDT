/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "matrix.h"
#include "polynomial.h"
#include <gtest/gtest.h>

using namespace ldt;

TEST(Polynomial_T, polynomial_simple) {
  auto a = Matrix<Tv>(new Tv[5]{2, 4, 5, 0, 0}, 5);
  auto pol = Polynomial<Tv>();
  pol.Data(a);
  ASSERT_EQ(pol.GetDegree(), 2);

  // multiply
  auto pol1 = Polynomial<Tv>();
  pol1.Data(*new Matrix<Tv>(new Tv[2]{2, 3}, 2));
  auto pol2 = Polynomial<Tv>();
  pol2.Data(*new Matrix<Tv>(new Tv[2]{4, 5}, 2));
  auto mul = PolynomialMultiply<Tv>(2, 2);
  auto S = new Tv[mul.StorageSize];
  mul.Calculate(pol1, pol2, S);
  ASSERT_EQ(true, mul.Result.Coefficients.Equals(
                      Matrix<Tv>(new Tv[3]{8, 22, 15}, 3), 1e-10));

  mul.Calculate(pol1, pol2, S, 2);
  ASSERT_EQ(true, mul.Result.Coefficients.Equals(
                      Matrix<Tv>(new Tv[2]{8, 22}, 2), 1e-10));

  delete[] S;

  // power
  auto pol3 = Polynomial<Tv>();
  pol3.Data(*new Matrix<Tv>(new Tv[2]{1, -1}, 2));

  for (auto i = 1; i < 4; i++) {
    auto pow = PolynomialPower<Tv>(i, 1);
    auto S1 = new Tv[pow.StorageSize];
    auto W1 = new Tv[pow.WorkSize];
    pow.Calculate(pol3, i, S1, W1);

    if (i == 1)
      ASSERT_EQ(pow.Result.Coefficients.Equals(Matrix<Tv>(new Tv[2]{1, -1}, 2),
                                               1e-10),
                true);
    else if (i == 2)
      ASSERT_EQ(pow.Result.Coefficients.Equals(
                    Matrix<Tv>(new Tv[3]{1, -2, 1}, 3), 1e-10),
                true);
    else if (i == 3)
      ASSERT_EQ(pow.Result.Coefficients.Equals(
                    Matrix<Tv>(new Tv[4]{1, -3, 3, -1}, 4), 1e-10),
                true);

    delete[] S1;
    delete[] W1;
  }
}

TEST(Polynomial_T, matrix_polynomial_simple) {
  auto a1 = Matrix<Tv>(new Tv[4]{1, 2, 3, 4}, 2, 2);
  auto a2 = Matrix<Tv>(new Tv[4]{5, 6, 7, 8}, 2, 2);
  auto a3 = Matrix<Tv>(new Tv[4]{9, 10, 11, 12}, 2, 2);
  auto a4 = Matrix<Tv>(new Tv[4]{0, 0, 0, 0}, 2, 2);

  auto pol = PolynomialM();
  auto vec = std::vector<Matrix<Tv> *>({&a1, &a2, &a3, &a4});
  pol.Data(vec, true);
  ASSERT_EQ(pol.GetDegree(), 2);

  a2.Set(1, 5.0);
  ASSERT_EQ(pol.Coefficients[1]->Get(1, 0), 5.0); // not a copy

  ASSERT_EQ(pol.IsMonic(), false);
}

TEST(Polynomial_T, matrix_polynomial_inv1) {
  // 1x1 simple
  // f(x) = 1-x
  auto b0 = Matrix<Tv>(new Tv[1]{1}, 1, 1);
  auto b1 = Matrix<Tv>(new Tv[1]{-1}, 1, 1);
  auto fx = PolynomialM();
  auto vec = std::vector<Matrix<Tv> *>({&b0, &b1});
  fx.Data(vec, true);

  auto inv1 = PolynomialMInvert(fx.GetSize(), fx.GetDegree(), 10);
  auto S1 = new Tv[inv1.StorageSize];
  auto W1 = new Tv[inv1.WorkSize];
  inv1.Calculate(fx, S1, W1, 10);

  ASSERT_EQ(inv1.Result.Coefficients.at(0)->Data[0], 1.0);
  ASSERT_EQ(inv1.Result.Coefficients.at(0)->Data[0], 1.0); //?
  ASSERT_EQ(inv1.Result.Coefficients.at(0)->Data[0], 1.0); //?

  // 1x1
  auto a1 = Matrix<Tv>(new Tv[1]{1}, 1, 1);
  auto a2 = Matrix<Tv>(new Tv[1]{0}, 1, 1);
  auto a3 = Matrix<Tv>(new Tv[1]{-0.4}, 1, 1);
  auto mp = PolynomialM();
  auto vec1 = std::vector<Matrix<Tv> *>({&a1, &a2, &a3});
  mp.Data(vec1, true);

  auto inv2 = PolynomialMInvert(mp.GetSize(), mp.GetDegree(), 20);
  auto S2 = new Tv[inv2.StorageSize];
  auto W2 = new Tv[inv2.WorkSize];
  inv2.Calculate(mp, S2, W2, 20);

  auto mul = PolynomialMMultiply(mp.GetSize(), mp.GetDegree(),
                                 inv2.Result.GetDegree(), 20);
  auto S3 = new Tv[mul.StorageSize];
  mul.Calculate(mp, inv2.Result, S3, 20);

  ASSERT_EQ(mul.Result.Coefficients.at(0)->Data[0], 1.0);
  ASSERT_EQ(mul.Result.Coefficients.at(1)->Data[0], 0.0);
  ASSERT_EQ(mul.Result.Coefficients.at(2)->Data[0], 0.0);
}

TEST(Polynomial_T, matrix_polynomial_inv3) {
  auto a1 = Matrix<Tv>(new Tv[9]{1, 2, 6, 3, 3, 7, 2, 5, 8}, 3, 3);
  auto a2 = Matrix<Tv>(new Tv[9]{4, 6, 6, 2, 4, 1, 3, 5, 5}, 3, 3);
  auto a3 = Matrix<Tv>(new Tv[9]{3, 5, 6, 2, 7, 3, 4, 5, 3}, 3, 3);
  auto mp = PolynomialM();
  auto vec1 = std::vector<Matrix<Tv> *>({&a1, &a2, &a3});
  mp.Data(vec1, true);

  auto inv2 = PolynomialMInvert(mp.GetSize(), mp.GetDegree(), 20);
  auto S2 = new Tv[inv2.StorageSize];
  auto W2 = new Tv[inv2.WorkSize];
  inv2.Calculate(mp, S2, W2, 20);

  auto a1inv = Matrix<Tv>(0.0, new Tv[9]{0, 0, 0, 0, 0, 0, 0, 0, 0}, 3, 3);
  a1.Inv(a1inv);
  ASSERT_EQ(true, a1inv.Equals(*inv2.Result.Coefficients.at(0)));

  // check with multiply
  auto mul = PolynomialMMultiply(mp.GetSize(), mp.GetDegree(),
                                 inv2.Result.GetDegree(), 20);
  auto S3 = new Tv[mul.StorageSize];
  mul.Calculate(mp, inv2.Result, S3, 20);

  auto I = Matrix<Tv>(new Tv[9]{1, 0, 0, 0, 1, 0, 0, 0, 1}, 3, 3);
  auto Z = Matrix<Tv>(new Tv[9]{0, 0, 0, 0, 0, 0, 0, 0, 0}, 3, 3);
  ASSERT_EQ(true, mul.Result.Coefficients.at(0)->Equals(I, 1e-10));
  ASSERT_EQ(true, mul.Result.Coefficients.at(1)->Equals(Z, 1e-10));
  ASSERT_EQ(true, mul.Result.Coefficients.at(2)->Equals(Z, 1e-10));
  ASSERT_EQ(true, mul.Result.Coefficients.at(9)->Equals(Z, 1e-10));
  ASSERT_EQ(true, mul.Result.Coefficients.at(10)->Equals(Z, 1e-10));
  ASSERT_EQ(true, mul.Result.Coefficients.at(11)->Equals(Z, 1e-10));
}

TEST(Polynomial_T, matrix_polynomial_mul) {
  auto a1 = Matrix<Tv>(new Tv[9]{1, 2, 6, 3, 3, 7, 2, 5, 8}, 3, 3);
  auto a2 = Matrix<Tv>(new Tv[9]{4, 6, 6, 2, 4, 1, 3, 5, 5}, 3, 3);
  auto a3 = Matrix<Tv>(new Tv[9]{3, 5, 6, 2, 7, 3, 4, 5, 3}, 3, 3);
  auto pol1 = PolynomialM();
  auto vec = std::vector<Matrix<Tv> *>({&a1, &a2, &a3});
  pol1.Data(vec, true);

  auto b1 = Matrix<Tv>(new Tv[9]{1, 2, 6, 3, 3, 7, 2, 5, 8}, 3, 3);
  auto b2 = Matrix<Tv>(new Tv[9]{4, 6, 6, 2, 4, 1, 3, 5, 5}, 3, 3);
  auto pol2 = PolynomialM();
  auto vec1 = std::vector<Matrix<Tv> *>({&b1, &b2});
  pol2.Data(vec1);

  auto mul =
      PolynomialMMultiply(pol1.GetSize(), pol1.GetDegree(), pol2.GetDegree());
  auto S = new Tv[mul.StorageSize];
  mul.Calculate(pol1, pol2, S);

  ASSERT_EQ(true, mul.Result.Coefficients.at(0)->Equals(
                      Matrix<Tv>(new Tv[9]{19, 38, 68, 26, 50, 95, 33, 59, 111},
                                 3, 3),
                      1e-8));
  ASSERT_EQ(true,
            mul.Result.Coefficients.at(3)->Equals(
                Matrix<Tv>(new Tv[9]{48, 92, 60, 18, 43, 27, 39, 75, 48}, 3, 3),
                1e-8));
}

TEST(Polynomial_T, matrix_polynomial_mul_1x1) {
  auto a1 = Matrix<Tv>(new Tv[1]{1}, 1, 1);
  auto a2 = Matrix<Tv>(new Tv[1]{0}, 1, 1);
  auto a3 = Matrix<Tv>(new Tv[1]{2}, 1, 1);
  auto a4 = Matrix<Tv>(new Tv[1]{0}, 1, 1);
  auto a5 = Matrix<Tv>(new Tv[1]{4}, 1, 1);
  auto pol1 = PolynomialM();
  auto vec = std::vector<Matrix<Tv> *>({&a1, &a2, &a3, &a4, &a5});
  pol1.Data(vec);

  auto b1 = Matrix<Tv>(new Tv[1]{1}, 1, 1);
  auto b2 = Matrix<Tv>(new Tv[1]{0}, 1, 1);
  auto b3 = Matrix<Tv>(new Tv[1]{-1}, 1, 1);
  auto pol2 = PolynomialM();
  auto vec1 = std::vector<Matrix<Tv> *>({&b1, &b2, &b3});
  pol2.Data(vec1);

  auto mul =
      PolynomialMMultiply(pol1.GetSize(), pol1.GetDegree(), pol2.GetDegree());
  auto S = new Tv[mul.StorageSize];
  mul.Calculate(pol1, pol2, S);

  ASSERT_EQ(mul.Result.Coefficients.at(0)->Get(0), 1.0);
  ASSERT_EQ(mul.Result.Coefficients.at(4)->Get(0), 2.0);
}
