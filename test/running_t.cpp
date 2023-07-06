/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "running.h"
#include <gtest/gtest.h>
#include <random>

using namespace ldt;

TEST(Running_T, mean) {

  auto r = RunningMoments<1>();
  r.PushNew(2, 1);
  r.PushNew(3, 1);
  r.PushNew(5, 1);
  r.PushNew(4, 1);
  r.PushNew(2, 1);
  r.PushNew(-3, 1);
  ASSERT_NEAR(r.GetMean(), 2.1666666666666666666, 1e-14);
}

TEST(Running_T, mean_combine) {

  auto r1 = RunningMoments<1, true, true, Tv>();
  r1.PushNew(2, 2);
  r1.PushNew(3, 1);
  r1.PushNew(5, 3);
  r1.PushNew(4, 0.4);
  auto r2 = RunningMoments<1, true, true, Tv>();
  r2.PushNew(2, 0.3);
  r2.PushNew(-3, 3);
  r2.PushNew(3, 2);
  r2.PushNew(-3, 3);
  r2.PushNew(4, 2);
  r2.PushNew(-2, 1);
  auto r = RunningMoments<1, true, true, Tv>();
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

TEST(Running_T, variance) {

  auto r = RunningMoments<2, true, true, Tv>();
  r.PushNew(2, 1);
  r.PushNew(3, 1);
  r.PushNew(5, 1);
  r.PushNew(4, 1);
  r.PushNew(2, 1);
  r.PushNew(-3, 1);
  ASSERT_NEAR(r.GetMean(), 2.1666666666666666666, 1e-14);
  ASSERT_NEAR(r.GetVariance(), 6.47222222222222, 1e-14);
}

TEST(Running_T, kurtosis) {

  auto r = RunningMoments<4, true, true, Tv>();
  r.PushNew(1, 1);
  r.PushNew(2, 1);
  r.PushNew(3, 1);
  r.PushNew(4, 1);
  r.PushNew(5, 1);
  r.PushNew(6, 1);
  r.PushNew(1, 1);
  ASSERT_NEAR(r.GetMean(), 3.14285714285714000, 1e-14);
  ASSERT_NEAR(r.GetVariance(), 3.26530612244898000, 1e-14);
  ASSERT_NEAR(r.GetSkewness(), 0.22234764798058900, 1e-14);
  ASSERT_NEAR(r.GetKurtosis(), -1.3526562500, 1e-6);

  // weight
  auto r1 = RunningMoments<4, true, true, Tv>();
  r1.PushNew(1, 2);
  r1.PushNew(2, 1);
  r1.PushNew(3, 1);
  r1.PushNew(4, 1);
  r1.PushNew(5, 1);
  r1.PushNew(6, 1);
  ASSERT_NEAR(r1.GetMean(), 3.14285714285714000, 1e-14);
  ASSERT_NEAR(r1.GetVariance(), 3.26530612244898000, 1e-14);
  ASSERT_NEAR(r1.GetSkewness(), 0.22234764798058900, 1e-14);
  ASSERT_NEAR(r1.GetKurtosis(), -1.3526562500, 1e-6);
}

TEST(Running_T, combine) {

  auto r1 = RunningMoments<4, true, true, Tv>();
  r1.PushNew(1, 1);
  r1.PushNew(2, 1);
  r1.PushNew(3, 1);
  r1.PushNew(4, 1);

  auto r2 = RunningMoments<4, true, true, Tv>();
  r2.PushNew(5, 1);
  r2.PushNew(6, 1);
  r2.PushNew(1, 1);

  r1.Combine(r2);

  ASSERT_NEAR(r1.GetMean(), 3.14285714285714000, 1e-14);
  ASSERT_NEAR(r1.GetVariance(), 3.26530612244898000, 1e-14);
  ASSERT_NEAR(r1.GetSkewness(), 0.22234764798058900, 1e-14);
  ASSERT_NEAR(r1.GetKurtosis(), -1.3526562500, 1e-6);
}

TEST(Running_T, Combine_sample) {

  // Sample 1
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(5.0, 2.0);
  int n1 = 100000;
  std::vector<double> sample1(n1);
  for (int i = 0; i < n1; i++) {
    sample1[i] = distribution(generator);
  }

  // Sample 2
  int n2 = 200000;
  std::vector<double> sample2(n2);
  for (int i = 0; i < n2; i++) {
    sample2[i] = distribution(generator);
  }

  // Combine
  std::vector<double> combined_sample(n1 + n2);
  copy(sample1.begin(), sample1.end(), combined_sample.begin());
  copy(sample2.begin(), sample2.end(), combined_sample.begin() + n1);

  Tv mean1 = 0, var1 = 0, skew1 = 0, kurt1 = 0;
  Array<Tv>::template Moments<false, false, 4>(
      &sample1[0], sample1.size(), nullptr, mean1, &var1, &skew1, &kurt1);

  Tv mean2 = 0, var2 = 0, skew2 = 0, kurt2 = 0;
  Array<Tv>::template Moments<false, false, 4>(
      &sample2[0], sample2.size(), nullptr, mean2, &var2, &skew2, &kurt2);

  Tv mean = 0, var = 0, skew = 0, kurt = 0;
  Array<Tv>::template Moments<false, false, 4>(&combined_sample[0],
                                               combined_sample.size(), nullptr,
                                               mean, &var, &skew, &kurt);

  // Combine moments
  auto result = RunningMoments<4, true, false>();
  result.Combine(mean1, var1, skew1, kurt1, n1);
  result.Combine(mean2, var2, skew2, kurt2, n2);

  // Compare combined moments
  ASSERT_NEAR(result.GetMean(), mean, 1e-4);
  ASSERT_NEAR(result.GetVariance(), var, 1e-4);
  ASSERT_NEAR(result.GetSkewness(), skew, 1e-4);
  ASSERT_NEAR(result.GetKurtosis(), kurt, 1e-4);
}