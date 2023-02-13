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
  RocOptions rocOptions;
  auto auc = ROC<false, false>(6);
  auc.Calculate(y, scores, nullptr, rocOptions);
  ASSERT_NEAR(auc.Result, 0.5, 1e-15);

  scores = Matrix<Tv>(new Tv[6]{0.2, 0.2, 0.2, 0.2, 0.2, 0.2}, 6, 1);
  auc.Calculate(y, scores, nullptr, rocOptions);
  ASSERT_NEAR(auc.Result, 0.5, 1e-15);

  scores = Matrix<Tv>(new Tv[6]{0.1, 0.9, 0.1, 0.9, 0.1, 0.9}, 6, 1);
  auc.Calculate(y, scores, nullptr, rocOptions);
  ASSERT_NEAR(auc.Result, 0.875, 1e-15);
  ASSERT_NEAR(std::get<1>(auc.Points.at(1)), 0.75, 1e-15);

  scores = Matrix<Tv>(new Tv[6]{0.9, 0.1, 0.9, 0.1, 0.9, 0.1}, 6, 1);
  auc.Calculate(y, scores, nullptr, rocOptions);
  ASSERT_NEAR(auc.Result, 0.125, 1e-15);
  ASSERT_NEAR(std::get<1>(auc.Points.at(1)), 0.25, 1e-15);

  // partial
  scores = Matrix<Tv>(new Tv[6]{0.5, 0.5, 0.5, 0.5, 0.5, 0.5}, 6, 1);
  rocOptions.LowerThreshold = 0.2;
  rocOptions.UpperThreshold = 0.8;
  auc.Calculate(y, scores, nullptr, rocOptions);
  ASSERT_NEAR(auc.Result, 0.3 / 0.6, 1e-15);

  scores = Matrix<Tv>(new Tv[6]{0.1, 0.9, 0.1, 0.9, 0.1, 0.9}, 6, 1);
  rocOptions.LowerThreshold = 0.2;
  rocOptions.UpperThreshold = 0.8;
  auc.Calculate(y, scores, nullptr, rocOptions);
  ASSERT_NEAR(auc.Result, 0.525 / 0.6, 1e-15);

  // another y (more complex)
  y = Matrix<Tv>(new Tv[10]{1, 0, 1, 0, 1, 1, 0, 0, 1, 0}, 10, 1);
  scores = Matrix<Tv>(
      new Tv[10]{0.1, 0.2, 0.3, 0.5, 0.5, 0.5, 0.7, 0.8, 0.9, 1}, 10, 1);
  auc = ROC<false, false>(10);
  rocOptions.LowerThreshold = 0.2;
  rocOptions.UpperThreshold = 0.8;
  auc.Calculate(y, scores, nullptr, rocOptions);
  ASSERT_NEAR(auc.Result, 0.44 / 0.6, 1e-15);
  rocOptions.LowerThreshold = 0.2;
  rocOptions.UpperThreshold = 0.4;
  auc.Calculate(y, scores, nullptr, rocOptions);
  ASSERT_NEAR(auc.Result, 0.12 / 0.2, 1e-15);
  rocOptions.LowerThreshold = 0.3;
  rocOptions.UpperThreshold = 0.4;
  auc.Calculate(y, scores, nullptr, rocOptions);
  ASSERT_NEAR(auc.Result, 0.07 / 0.1, 1e-15);
  rocOptions.LowerThreshold = 0.3;
  rocOptions.UpperThreshold = 0.8;
  auc.Calculate(y, scores, nullptr, rocOptions);
  ASSERT_NEAR(auc.Result, 0.39 / 0.5, 1e-15);
  rocOptions.LowerThreshold = 0.4;
  rocOptions.UpperThreshold = 0.6;
  auc.Calculate(y, scores, nullptr, rocOptions);
  ASSERT_NEAR(auc.Result, 0.16 / 0.2, 1e-15);
  rocOptions.LowerThreshold = 0;
  rocOptions.UpperThreshold = 1;
  auc.Calculate(y, scores, nullptr, rocOptions);
  ASSERT_NEAR(auc.Result, 0.68 / 1, 1e-15);
  auc.Calculate(y, scores, nullptr, rocOptions);
  ASSERT_NEAR(auc.Result, 0.68, 1e-15);
}
