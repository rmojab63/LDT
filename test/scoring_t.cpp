/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "matrix.h"
#include "scoring.h"
#include <gtest/gtest.h>

using namespace ldt;

TEST(Scoring_T, normal) {

  ASSERT_EQ(0, Scoring::GetScoreCrpsNormal(0, 0, 0));
  ASSERT_EQ(5, Scoring::GetScoreCrpsNormal(13, 18, 0));
  ASSERT_NEAR(0.2336949772551092, Scoring::GetScoreCrpsNormal(0, 0, 1), 1e-14);
  ASSERT_NEAR(0.2336949772551092, Scoring::GetScoreCrpsNormal(1, 1, 1), 1e-14);
  ASSERT_NEAR(0.4673899545102183, Scoring::GetScoreCrpsNormal(0, 0, 2), 1e-14);
  ASSERT_NEAR(0.4673899545102183, Scoring::GetScoreCrpsNormal(1, 1, 2), 1e-14);
  ASSERT_NEAR(0.6628070625097118, Scoring::GetScoreCrpsNormal(0, 1, 2), 1e-14);
  ASSERT_NEAR(1.988848007954906, Scoring::GetScoreCrpsNormal(-2, 1, 2), 1e-14);
  ASSERT_NEAR(0.8328479351511628, Scoring::GetScoreCrpsNormal(1, 2, 3), 1e-14);
}

TEST(Scoring_T, logNormal) {
  ASSERT_EQ(1, Scoring::GetScoreCrpsLogNormal(0, 0, 0));
  ASSERT_NEAR(65659956.13733051, Scoring::GetScoreCrpsLogNormal(13, 18, 0), 1);
  ASSERT_NEAR(0.7905620507529407, Scoring::GetScoreCrpsLogNormal(0, 0, 1),
              1e-14);
  ASSERT_NEAR(1.262362929292143, Scoring::GetScoreCrpsLogNormal(1, 1, 1),
              1e-14);
  ASSERT_NEAR(1.162292665211865, Scoring::GetScoreCrpsLogNormal(0, 0, 2),
              1e-14);
  ASSERT_NEAR(3.159439031196646, Scoring::GetScoreCrpsLogNormal(0, 1, 2),
              1e-14);
  ASSERT_NEAR(2.427583576155807, Scoring::GetScoreCrpsLogNormal(-2, -1, 2),
              1e-14);
  ASSERT_NEAR(21.88641619776969, Scoring::GetScoreCrpsLogNormal(1, 2, 3),
              1e-14);
}

TEST(Scoring_T, auc_binary) {
  auto y = Matrix<Tv>(new Tv[6]{1, 0, 1, 0, 1, 1}, 6, 1);
  auto scores = Matrix<Tv>(new Tv[6]{0.5, 0.5, 0.5, 0.5, 0.5, 0.5}, 6, 1);
  auto auc = ROC<false, true>(6);
  auc.Calculate(y, scores, nullptr, true);
  ASSERT_NEAR(auc.Result, 0.5, 1e-15);

  scores = Matrix<Tv>(new Tv[6]{0.2, 0.2, 0.2, 0.2, 0.2, 0.2}, 6, 1);
  auc.Calculate(y, scores, nullptr, true);
  ASSERT_NEAR(auc.Result, 0.5, 1e-15);

  scores = Matrix<Tv>(new Tv[6]{0.1, 0.9, 0.1, 0.9, 0.1, 0.9}, 6, 1);
  auc.Calculate(y, scores, nullptr, true);
  ASSERT_NEAR(auc.Result, 0.875, 1e-15);
  ASSERT_NEAR(std::get<1>(auc.Points.at(1)), 0.75, 1e-15);

  scores = Matrix<Tv>(new Tv[6]{0.9, 0.1, 0.9, 0.1, 0.9, 0.1}, 6, 1);
  auc.Calculate(y, scores, nullptr, true);
  ASSERT_NEAR(auc.Result, 0.125, 1e-15);
  ASSERT_NEAR(std::get<1>(auc.Points.at(1)), 0.25, 1e-15);

  // partial
  scores = Matrix<Tv>(new Tv[6]{0.5, 0.5, 0.5, 0.5, 0.5, 0.5}, 6, 1);
  auc.Calculate(y, scores, nullptr, true, 0.2, 0.8);
  ASSERT_NEAR(auc.Result, 0.3, 1e-15);

  scores = Matrix<Tv>(new Tv[6]{0.1, 0.9, 0.1, 0.9, 0.1, 0.9}, 6, 1);
  auc.Calculate(y, scores, nullptr, true, 0.2, 0.8);
  ASSERT_NEAR(auc.Result, 0.525, 1e-15);
}

TEST(Scoring_T, costMatrix) {

  auto y = Matrix<Tv>(new Tv[8]{0, 1, 0, 1, 0, 1, 0, 1}, 8, 1);
  auto scores = Matrix<Tv>(new Tv[16]{0.4, 0.6, 0.45, 0.9, 0.6, 0.3, 0.5, 0.1,
                                      0.6, 0.4, 0.55, 0.1, 0.4, 0.7, 0.5, 0.9},
                           8, 2);

  Matrix<Tv> cost_table =
      Matrix<Tv>(new Tv[12]{0.3, 0.5, 0.7, 1.0, 0, 1, 2, 3, 9, 8, 7, 0}, 4, 3);
  CostMatrix<true>::Check(cost_table, 2);
  auto weights =
      Matrix<Tv>(new Tv[8]{0.4, 0.6, 0.45, 0.9, 0.6, 0.3, 0.5, 0.1}, 8, 1);

  auto cm = CostMatrix<false>(1);
  auto S = new Tv[cm.StorageSize];
  cm.Calculate(std::vector<Matrix<Tv>>({cost_table, cost_table}), y, scores,
               &weights, S);

  //??!!!
}
