/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "array.h"
#include <gtest/gtest.h>
#include <random>

using namespace ldt;

TEST(Array_t, range) {

  Tv data1[] = {NAN, NAN, 2, 3, NAN, 5, NAN, NAN, NAN};
  bool hasmissing;
  auto range = Array<Tv>::GetRange(data1, 9, hasmissing);
  ASSERT_EQ(range.StartIndex, (Ti)2);
  ASSERT_EQ(range.EndIndex, (Ti)5);
  ASSERT_EQ(hasmissing, true);

  Tv data2[] = {NAN, 1, 2, 3, 4, 5, 6, NAN, NAN};
  range = Array<Tv>::GetRange(data2, 9, hasmissing);
  ASSERT_EQ(range.StartIndex, (Ti)1);
  ASSERT_EQ(range.EndIndex, (Ti)6);
  ASSERT_EQ(hasmissing, false);
}

TEST(Array_t, interpolate) {

  // by column

  Tv data[] = {NAN, NAN, 2, 3, NAN, 4, 5, NAN, NAN, NAN, 9, 7, NAN};
  Ti count = 0;
  auto range = Array<Tv>::Interpolate(data, 13, count);
  ASSERT_NEAR(data[7], 6.0, 1e-16);
  ASSERT_NEAR(data[8], 7.0, 1e-16);
  ASSERT_NEAR(data[9], 8.0, 1e-16);
  ASSERT_EQ(count, 4);

  Tv data1[] = {NAN, NAN, 3, NAN, 2, 1, 2, 3, 2, 4, 5, 1, 7};
  range = Array<Tv>::Interpolate(data1, 13, count);
  ASSERT_EQ(count, 1);
}

TEST(Array_t, partition) {

  auto v = std::vector<Tv>({1, 2, 3, 4, 5, 6, 7});
  std::vector<std::vector<Tv>> partitions;
  Array<Tv>::PartitionEqual(v, partitions, 2, false);
  ASSERT_EQ(4, partitions.size());
  ASSERT_EQ(std::vector<Tv>({1, 2}), partitions.at(0));
  ASSERT_EQ(std::vector<Tv>({7}), partitions.at(3));

  // backward
  Array<Tv>::PartitionEqual(v, partitions, 2, true);
  ASSERT_EQ(4, partitions.size());
  ASSERT_EQ(std::vector<Tv>({1}), partitions.at(0));
  ASSERT_EQ(std::vector<Tv>({6, 7}), partitions.at(3));
}

TEST(Array_t, stat_moments) {

  // Unweighted (population)
  const double arr[] = {1, 2, 5, 5, 10};
  auto v = Array<Tv>::Mean<false, false>(arr, 5, nullptr);
  ASSERT_NEAR(4.6, v, 1e-6);
  v = Array<Tv>::Variance<false, false>(arr, 5, nullptr);
  ASSERT_NEAR(9.84, v, 1e-5);
  v = Array<Tv>::Skewness<false, false>(arr, 5, nullptr);
  ASSERT_NEAR(0.60492, v, 1e-5);
  v = Array<Tv>::Kurtosis<false, false>(arr, 5, nullptr);
  ASSERT_NEAR(-0.802201, v, 1e-5);

  // Weighted (population)
  const double arr1[] = {1, 2, 5, 10};
  const double weig[] = {1, 1, 2, 1};
  v = Array<Tv>::Mean<true, false>(arr1, 4, weig);
  ASSERT_NEAR(4.6, v, 1e-6);
  v = Array<Tv>::Variance<true, false>(arr1, 4, weig);
  ASSERT_NEAR(9.84, v, 1e-6);
  v = Array<Tv>::Skewness<true, false>(arr1, 4, weig);
  ASSERT_NEAR(0.60492, v, 1e-3);
  v = Array<Tv>::Kurtosis<true, false>(arr1, 4, weig);
  ASSERT_NEAR(-0.802201, v, 1e-3);

  // NAN
}