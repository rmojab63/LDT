/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "data_split.h"
#include "matrix.h"
#include <gtest/gtest.h>
#include <random>

using namespace ldt;

TEST(Data_Split_t, discrete) {
  const Matrix<Tv> data = Matrix<Tv>(
      new std::vector<Tv>{1,  2,  3,  0,  2,  1,  3,  0,  0,  0,  1,  2,
                          0,  1,  0,  2,  0,  1,  3,  0,  1,  1,  1,  1,
                          1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                          1,  1,  1,  1,  10, 11, 12, 13, 14, 15, 16, 17,
                          18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29},
      20, 3);

  auto model = DataSplitDiscrete(data.RowsCount, data.ColsCount, 4);
  auto S = new Tv[model.StorageSize];
  model.Calculate(data, S, 0.8, 0);

  ASSERT_EQ(model.Rows.at(0)->size(), 8);
  ASSERT_EQ(model.Rows.at(1)->size(), 5);

  std::random_device rdev{};
  auto eng = std::mt19937(rdev());
  auto Wi = new Ti[model.WorkSizeI];
  auto row = Matrix<Tv>(new Tv[data.ColsCount], data.ColsCount, 1);
  auto row0 = Matrix<Tv>(new Tv[data.ColsCount], data.ColsCount, 1);

  auto vec = std::set<Tv>();
  for (auto k = 0; k < 10; k++) {

    model.Shuffle(data, Wi, eng);
    vec.insert(model.Sample0.Get(0, 2));
    // have we distributed all observations?

    bool found = false;
    for (auto i = 0; i < data.RowsCount; i++) {
      data.GetRow(i, row);
      for (auto j = 0; j < model.Sample0.RowsCount; j++) {
        model.Sample0.GetRow(j, row0);
        if (row.Equals(row0, 0)) {
          found = true;
          break;
        }
      }
      if (found == false) {
        for (auto j = 0; j < model.Sample1.RowsCount; j++) {
          model.Sample1.GetRow(j, row0);
          if (row.Equals(row0, 0)) {
            found = true;
            break;
          }
        }
      }
      ASSERT_EQ(true, found);
    }
  }

  ASSERT_NE(vec.size(), 1);
}