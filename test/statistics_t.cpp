/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the LGPL 3.0.
 See accompanying file LICENSE
*/

#include "matrix.h"
#include "statistics.h"
#include <gtest/gtest.h>
#include <random>

using namespace ldt;

TEST(Statistics_T, rank) {

  auto mat = Matrix<Tv>(new Tv[10 * 2]{1.1, 2,   4, 0,    2.2, 2.1, 2.03,
                                       3.2, 0.2, 3, -2.4, 1,   1,   1.1,
                                       1,   3,   2, 2,    2.1, -2.4},
                        10, 2);

  auto model = Rank(mat.RowsCount, mat.ColsCount);
  auto storage = new Tv[model.StorageSize];
  auto work = new Tv[model.WorkSize];
  model.Calculate(mat, work, storage, true);

  auto eq = Matrix<Tv>(new Tv[20]{2, 3, 9, 0, 6, 5, 4, 8, 1, 7,
                                  0, 2, 2, 5, 2, 9, 6, 6, 8, 0},
                       10, 2)
                .Equals(model.Result);
  ASSERT_EQ(true, eq);

  delete[] storage;
  delete[] work;
  delete[] mat.Data;
}

TEST(Statistics_T, ols) {

  // univariate
  Tv yy[9]{1, 1, 0, 0, 1, 0, 0, 1, 0};
  auto myy = Matrix<Tv>(yy, 9, 1);
  Tv xx[18]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 3, 6, 7, 3, 3, 5};

  Matrix<Tv> mxx = Matrix<Tv>(xx, 9, 2);

  auto ols = Ols(myy.RowsCount, myy.ColsCount, mxx.ColsCount, false, false);
  auto W = new Tv[ols.WorkSize];
  auto S = new Tv[ols.StorageSize];
  ols.Calculate(myy, mxx, S, W);

  ASSERT_NEAR(0.8842794, ols.Beta.Get0(0, 0), 1e-6);
  ASSERT_NEAR(-0.106986899, ols.Beta.Get0(1, 0), 1e-6);

  delete[] W;
  delete[] S;

  // multi
  Tv y[20]{1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
           20, 19, 18, 17, 16, 15, 14, 13, 12, 11};
  Matrix<Tv> my = Matrix<Tv>(y, 10, 2);

  Tv x[30]{1, 1, 1, 1, 1, 1,  1,  1,  1,  1,  1,  4,  2,  4, 6,
           6, 8, 8, 3, 3, 22, 10, 11, 12, 16, 12, 12, 11, 1, 11};
  Matrix<Tv> mx = Matrix<Tv>(x, 10, 3);

  auto ols1 = Ols(my.RowsCount, my.ColsCount, mx.ColsCount, false, false);
  auto W1 = new Tv[ols1.WorkSize];
  auto S1 = new Tv[ols1.StorageSize];
  ols1.Calculate(my, mx, S1, W1);

  ASSERT_NEAR(7.31859804231, ols1.Beta.Get0(0, 0), 1e-6);
  ASSERT_NEAR(13.68140195, ols1.Beta.Get0(0, 1), 1e-6);

  delete[] W1;
  delete[] S1;
}

TEST(Statistics_T, gls) {

  // univariate
  Tv yy[9]{1, 1, 0, 0, 1, 0, 0, 1, 0};
  Matrix<Tv> myy = Matrix<Tv>(yy, 9, 1);
  Tv xx[18]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 3, 6, 7, 3, 3, 5};
  Matrix<Tv> mxx = Matrix<Tv>(xx, 9, 2);
  Tv omega2[81];
  Matrix<Tv> momega2 = Matrix<Tv>(omega2, 9, 9);
  Matrix<Tv>::Diagonal(momega2, 1.0, 0.0);

  auto gls = Gls(myy.RowsCount, myy.ColsCount, mxx.ColsCount, false, false);
  auto W = new Tv[gls.WorkSize];
  auto S = new Tv[gls.StorageSize];
  gls.Calculate(myy, mxx, momega2, S, W);

  ASSERT_NEAR(0.8842794, gls.Beta.Get0(0, 0), 1e-5);
  ASSERT_NEAR(-0.106986899, gls.Beta.Get0(1, 0), 1e-5);

  delete[] W;
  delete[] S;

  // multi
  Tv y[20]{1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
           20, 19, 18, 17, 16, 15, 14, 13, 12, 11};
  Matrix<Tv> my = Matrix<Tv>(y, 10, 2);
  Tv x[30]{1, 1, 1, 1, 1, 1,  1,  1,  1,  1,  1,  4,  2,  4, 6,
           6, 8, 8, 3, 3, 22, 10, 11, 12, 16, 12, 12, 11, 1, 11};
  Matrix<Tv> mx = Matrix<Tv>(x, 10, 3);
  Tv omega[100];
  Matrix<Tv> momega = Matrix<Tv>(omega, 10, 10);
  Matrix<Tv>::Diagonal(momega, 1.0, 0.0);

  auto gls1 = Gls(my.RowsCount, my.ColsCount, mx.ColsCount, false, false);
  auto W1 = new Tv[gls1.WorkSize];
  auto S1 = new Tv[gls1.StorageSize];
  gls1.Calculate(my, mx, momega, S1, W1);

  ASSERT_NEAR(7.31859804231, gls1.Beta.Get0(0, 0), 1e-6);
  ASSERT_NEAR(13.68140195, gls1.Beta.Get0(0, 1), 1e-6);

  delete[] W1;
  delete[] S1;
}
