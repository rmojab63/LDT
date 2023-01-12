/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the LGPL 3.0.
 See accompanying file LICENSE
*/

#include "running.h"
#include <gtest/gtest.h>

using namespace ldt;

TEST(Running_T, runningweighted_mean) {

  auto r = RunningWeightedMean();
  r.PushNew(2, 1);
  r.PushNew(3, 1);
  r.PushNew(5, 1);
  r.PushNew(4, 1);
  r.PushNew(2, 1);
  r.PushNew(-3, 1);
  ASSERT_NEAR(r.GetMean(), 2.1666666666666666666, 1e-14);
}

TEST(Running_T, runningweighted_mean_combine) {

  auto r1 = RunningWeightedMean();
  r1.PushNew(2, 2);
  r1.PushNew(3, 1);
  r1.PushNew(5, 3);
  r1.PushNew(4, 0.4);
  auto r2 = RunningWeightedMean();
  r2.PushNew(2, 0.3);
  r2.PushNew(-3, 3);
  r2.PushNew(3, 2);
  r2.PushNew(-3, 3);
  r2.PushNew(4, 2);
  r2.PushNew(-2, 1);
  auto r = RunningWeightedMean();
  r.PushNew(2, 2);
  r.PushNew(3, 1);
  r.PushNew(5, 3);
  r.PushNew(4, 0.4);
  r.PushNew(2, 0.3);
  r.PushNew(-3, 3);
  r.PushNew(3, 2);
  r.PushNew(-3, 3);
  r.PushNew(4, 2);
  r.PushNew(-2, 1);

  r1.Combine(r2);
  ASSERT_NEAR(r.GetMean(), r1.GetMean(), 1e-14);
}

TEST(Running_T, runningweighted_variance) {

  auto r = RunningWeightedVariance();
  r.PushNew(2, 1);
  r.PushNew(3, 1);
  r.PushNew(5, 1);
  r.PushNew(4, 1);
  r.PushNew(2, 1);
  r.PushNew(-3, 1);
  ASSERT_NEAR(r.GetMean(), 2.1666666666666666666, 1e-14);
  ASSERT_NEAR(r.GetVariancePopulation(), 6.47222222222222, 1e-14);
}

TEST(Running_T, runningweighted_4) {

  RunningWeighted4 r = RunningWeighted4();
  r.PushNew(1, 1);
  r.PushNew(2, 1);
  r.PushNew(3, 1);
  r.PushNew(4, 1);
  r.PushNew(5, 1);
  r.PushNew(6, 1);
  r.PushNew(1, 1);
  ASSERT_NEAR(r.GetMean(), 3.14285714285714000, 1e-14);
  ASSERT_NEAR(r.GetVariancePopulation(), 3.26530612244898000, 1e-14);
  ASSERT_NEAR(r.GetSkewnessPopulation(), 0.22234764798058900, 1e-14);
  ASSERT_NEAR(r.GetKurtosisPopulation(), -1.3526562500, 1e-6);

  // weight
  RunningWeighted4 r1 = RunningWeighted4();
  r1.PushNew(1, 2);
  r1.PushNew(2, 1);
  r1.PushNew(3, 1);
  r1.PushNew(4, 1);
  r1.PushNew(5, 1);
  r1.PushNew(6, 1);
  ASSERT_NEAR(r1.GetMean(), 3.14285714285714000, 1e-14);
  ASSERT_NEAR(r1.GetVariancePopulation(), 3.26530612244898000, 1e-14);
  ASSERT_NEAR(r1.GetSkewnessPopulation(), 0.22234764798058900, 1e-14);
  ASSERT_NEAR(r1.GetKurtosisPopulation(), -1.3526562500, 1e-6);
}

TEST(Running_T, runningweighted_4_combine) {

  RunningWeighted4 r1 = RunningWeighted4();
  r1.PushNew(1, 1);
  r1.PushNew(2, 1);
  r1.PushNew(3, 1);
  r1.PushNew(4, 1);

  RunningWeighted4 r2 = RunningWeighted4();
  r2.PushNew(5, 1);
  r2.PushNew(6, 1);
  r2.PushNew(1, 1);

  r1.Combine(r2);

  ASSERT_NEAR(r1.GetMean(), 3.14285714285714000, 1e-14);
  ASSERT_NEAR(r1.GetVariancePopulation(), 3.26530612244898000, 1e-14);
  ASSERT_NEAR(r1.GetSkewnessPopulation(), 0.22234764798058900, 1e-14);
  ASSERT_NEAR(r1.GetKurtosisPopulation(), -1.3526562500, 1e-6);
}
