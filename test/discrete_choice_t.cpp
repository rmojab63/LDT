/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "discrete_choice.h"
#include "matrix.h"
#include "optimization.h"
#include <gtest/gtest.h>
#include <random>

using namespace ldt;

TEST(DiscreteChoice_t, logit) {
  const Matrix<Tv> y = Matrix<Tv>(new std::vector<Tv>{
      1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0});
  const Matrix<Tv> x =
      Matrix<Tv>(new std::vector<Tv>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 1, 1, 7, 5, 8, 9, 7, 3, 7, 8,
                                     5, 3, 5, 3, 3, 4, 6, 5, 8, 3, 6, 10},
                 20, 2);

  auto model = DiscreteChoice<DiscreteChoiceModelType::kBinary,
                              DiscreteChoiceDistType::kLogit>(
      y.RowsCount, x.ColsCount, 2, true);
  auto W = new Tv[model.WorkSize];
  auto S = new Tv[model.StorageSize];
  model.Calculate(y, x, nullptr, S, W, 2);

  // ASSERT_NEAR(model.resid.Data[0], 0.56312260, 1e-4);
  // ASSERT_NEAR(model.resid.Data[19], -0.297344, 1e-4);

  ASSERT_NEAR(model.Beta.Data[0], 1.1604, 1e-4);
  ASSERT_NEAR(model.Beta.Data[1], -0.20204, 1e-4);
  ASSERT_NEAR(model.BetaStd.Data[0], 1.33173, 1e-4);
  ASSERT_NEAR(model.BetaStd.Data[1], 0.21821, 1e-4);
  ASSERT_NEAR(model.BetaZ.Data[0], 0.8713, 1e-4);
  ASSERT_NEAR(model.BetaProb.Data[1], 0.3545, 1e-4);
  ASSERT_NEAR(model.LogL, -13.41503, 1e-4);
  ASSERT_NEAR(model.Aic, 1.54150, 1e-4);
  ASSERT_NEAR(model.Sic, 1.64107, 1e-4);

  // test probabilities calculations
  Tv *Q = new Tv[x.RowsCount + model.NumChoices - 2];
  Tv *Z = new Tv[(model.NumChoices) * x.RowsCount];
  Matrix<Tv> stor = Matrix<Tv>(Z, x.RowsCount, model.NumChoices);
  model.GetProbabilities(x, stor, Q);

  ASSERT_NEAR(stor.Get(0, 0), 0.4368774, 1e-4);
  ASSERT_NEAR(stor.Get(1, 1), 1.0 - 0.53749, 1e-4);

  delete[] Q;
  delete[] Z;
  delete[] W;
  delete[] S;
}

TEST(DiscreteChoice_t, probit) {
  const Matrix<Tv> y = Matrix<Tv>(new std::vector<Tv>{
      1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0});
  const Matrix<Tv> x =
      Matrix<Tv>(new std::vector<Tv>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 1, 1, 7, 5, 8, 9, 7, 3, 7, 8,
                                     5, 3, 5, 3, 3, 4, 6, 5, 8, 3, 6, 10},
                 20, 2);

  auto model = DiscreteChoice<DiscreteChoiceModelType::kBinary,
                              DiscreteChoiceDistType::kProbit>(
      y.RowsCount, x.ColsCount, 2, true);
  auto W = new Tv[model.WorkSize];
  auto S = new Tv[model.StorageSize];
  model.Calculate(y, x, nullptr, S, W, 2);

  // ASSERT_NEAR(model.resid.Data[0], 0.564637, 1e-4);
  // ASSERT_NEAR(model.resid.Data[19], -0.29243322, 1e-4);

  ASSERT_NEAR(model.Beta.Data[0], 0.732, 1e-3);
  ASSERT_NEAR(model.Beta.Data[1], -0.127, 1e-3);
  ASSERT_NEAR(model.BetaStd.Data[0], 0.828, 1e-3);
  ASSERT_NEAR(model.BetaStd.Data[1], 0.136, 1e-3);
  ASSERT_NEAR(model.BetaZ.Data[0], 0.883, 1e-3);
  ASSERT_NEAR(model.BetaProb.Data[1], 0.349, 1e-3);
  ASSERT_NEAR(model.LogL, -13.410, 1e-3);
  ASSERT_NEAR(model.Aic, 1.54150, 1e-3);
  ASSERT_NEAR(model.Sic, 1.64107, 1e-3);

  // test probabilities calculations

  Tv *Q = new Tv[x.RowsCount + model.NumChoices - 2];
  Tv *Z = new Tv[(model.NumChoices) * x.RowsCount];
  Matrix<Tv> stor = Matrix<Tv>(Z, x.RowsCount, model.NumChoices);
  model.GetProbabilities(x, stor, Q);

  ASSERT_NEAR(stor.Get(0, 0), 0.4353621, 1e-4);
  ASSERT_NEAR(stor.Get(1, 1), 1.0 - 0.5370332, 1e-4);

  delete[] Q;
  delete[] Z;
  delete[] W;
  delete[] S;
}

TEST(DiscreteChoice_t, logit_ordered) {
  const Matrix<Tv> y = Matrix<Tv>(new std::vector<Tv>{
      1, 2, 3, 0, 2, 1, 3, 0, 0, 0, 1, 2, 0, 1, 0, 2, 0, 1, 3, 0});
  const Matrix<Tv> x =
      Matrix<Tv>(new std::vector<Tv>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 1, 1, 7, 5, 8, 9, 7, 3, 7, 8,
                                     5, 3, 5, 3, 3, 4, 6, 5, 8, 3, 6, 10},
                 20, 2);

  auto model = DiscreteChoice<DiscreteChoiceModelType::kOrdered,
                              DiscreteChoiceDistType::kLogit>(
      y.RowsCount, x.ColsCount, 4, true);
  auto W = new Tv[model.WorkSize];
  auto S = new Tv[model.StorageSize];
  model.Calculate(y, x, nullptr, S, W);

  ASSERT_NEAR(model.Counts.Data[0], 8, 1e-4);
  ASSERT_NEAR(model.Counts.Data[1], 5, 1e-4);
  ASSERT_NEAR(model.Counts.Data[2], 4, 1e-4);
  ASSERT_NEAR(model.Counts.Data[3], 3, 1e-4);
  // ASSERT_NEAR(model.resid.Data[0], -0.116003507, 1e-4);
  // ASSERT_NEAR(model.resid.Data[19], 0.50331641, 1e-4);

  ASSERT_NEAR(model.Beta.Data[1], -0.0862511988674, 1e-4);
  ASSERT_NEAR(model.Beta.Data[0], 0.875777841676, 1e-4);
  ASSERT_NEAR(model.Beta.Data[2], 1.042348717, 1e-4);
  ASSERT_NEAR(model.Beta.Data[3], 2.147200467, 1e-4);

  ASSERT_NEAR(model.BetaStd.Data[0], 1.1639, 1e-4);
  ASSERT_NEAR(model.BetaStd.Data[1], 0.19577, 1e-4);
  ASSERT_NEAR(model.BetaZ.Data[0], 0.7524, 1e-4);
  ASSERT_NEAR(model.BetaProb.Data[1], 0.6595, 1e-4);
  ASSERT_NEAR(model.LogL, -26.2937, 1e-4);
  ASSERT_NEAR(model.Aic, 3.0293, 1e-4);
  ASSERT_NEAR(model.Sic, 3.2285, 1e-4);

  // test probabilities calculations

  Tv *Q = new Tv[x.RowsCount + model.NumChoices - 2];
  Tv *Z = new Tv[(model.NumChoices) * x.RowsCount];
  Matrix<Tv> stor = Matrix<Tv>(Z, x.RowsCount, model.NumChoices);
  model.GetProbabilities(x, stor, Q);

  ASSERT_NEAR(stor.Get(1, 0), 0.390664, 1e-4);
  ASSERT_NEAR(stor.Get(2, 1), 0.2482547, 1e-4);
  ASSERT_NEAR(stor.Get(3, 2), 0.166042, 1e-4);
  ASSERT_NEAR(stor.Get(4, 3), 0.1329434, 1e-4);

  delete[] Q;
  delete[] Z;

  delete[] W;
  delete[] S;
}

TEST(DiscreteChoice_t, logit_ordered2) {
  const Matrix<Tv> y =
      Matrix<Tv>(new std::vector<Tv>{1, 0, 1, 2, 3, 3, 1, 3, 0});
  const Matrix<Tv> x = Matrix<Tv>(
      new std::vector<Tv>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 9, 3, 5, 5, 4, 7},
      9, 2);

  auto model = DiscreteChoice<DiscreteChoiceModelType::kOrdered,
                              DiscreteChoiceDistType::kLogit>(
      y.RowsCount, x.ColsCount, 4, true);
  auto W = new Tv[model.WorkSize];
  auto S = new Tv[model.StorageSize];
  model.Calculate(y, x, nullptr, S, W);

  ASSERT_NEAR(model.Beta.Data[1], 0.09236752, 1e-4);
  ASSERT_NEAR(model.Beta.Data[0], 0.88666679, 1e-4);

  delete[] W;
  delete[] S;
}

TEST(DiscreteChoice_t, probit_ordered) {
  const Matrix<Tv> y = Matrix<Tv>(new std::vector<Tv>{
      1, 2, 3, 0, 2, 1, 3, 0, 0, 0, 1, 2, 0, 1, 0, 2, 0, 1, 3, 0});
  const Matrix<Tv> x =
      Matrix<Tv>(new std::vector<Tv>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 1, 1, 7, 5, 8, 9, 7, 3, 7, 8,
                                     5, 3, 5, 3, 3, 4, 6, 5, 8, 3, 6, 10},
                 20, 2);

  auto model = DiscreteChoice<DiscreteChoiceModelType::kOrdered,
                              DiscreteChoiceDistType::kProbit>(
      y.RowsCount, x.ColsCount, 4, true);
  auto W = new Tv[model.WorkSize];
  auto S = new Tv[model.StorageSize];
  model.Calculate(y, x, nullptr, S, W);

  ASSERT_NEAR(model.Counts.Data[0], 8, 1e-4);
  ASSERT_NEAR(model.Counts.Data[1], 5, 1e-4);
  ASSERT_NEAR(model.Counts.Data[2], 4, 1e-4);
  ASSERT_NEAR(model.Counts.Data[3], 3, 1e-4);
  // ASSERT_NEAR(model.resid.Data[0], -0.1121832, 1e-4);
  // ASSERT_NEAR(model.resid.Data[1], -0.667030, 1e-4);
  // ASSERT_NEAR(model.resid.Data[2], -1.6187411, 1e-4);
  // ASSERT_NEAR(model.resid.Data[19], 0.8638701, 1e-4);

  ASSERT_NEAR(model.Beta.Data[1], -0.0346, 1e-4);
  ASSERT_NEAR(model.Beta.Data[0], 0.448015, 1e-4);
  ASSERT_NEAR(model.Beta.Data[2], 0.643751189, 1e-4);
  ASSERT_NEAR(model.Beta.Data[3], 1.28704564, 1e-4);

  ASSERT_NEAR(model.BetaStd.Data[0], 0.7261112, 1e-4);
  ASSERT_NEAR(model.BetaStd.Data[1], 0.119042, 1e-4);
  ASSERT_NEAR(model.BetaZ.Data[0], 0.6170063, 1e-4);
  ASSERT_NEAR(model.BetaProb.Data[1], 0.7712450, 1e-4);

  ASSERT_NEAR(model.LogL, -26.34849, 1e-4);
  ASSERT_NEAR(model.Aic, 3.034849, 1e-4);
  ASSERT_NEAR(model.Sic, 3.23399, 1e-4);

  // test probabilities calculations
  Tv *Q = new Tv[x.RowsCount + model.NumChoices - 2];
  Tv *Z = new Tv[(model.NumChoices) * x.RowsCount];
  Matrix<Tv> stor = Matrix<Tv>(Z, x.RowsCount, model.NumChoices);
  model.GetProbabilities(x, stor, Q);

  ASSERT_NEAR(stor.Get(1, 0), 0.3916737, 1e-4);
  ASSERT_NEAR(stor.Get(2, 1), 0.2496973, 1e-4);
  ASSERT_NEAR(stor.Get(3, 2), 0.1810318, 1e-4);
  ASSERT_NEAR(stor.Get(4, 3), 0.1397799, 1e-4);

  delete[] Q;
  delete[] Z;

  delete[] W;
  delete[] S;
}

TEST(DiscreteChoice_t, probit_ordered_large) {
  Ti n = 10000;
  auto NumCutoff = 20;
  Matrix<Tv> y = Matrix<Tv>(new std::vector<Tv>(n), n, 1);
  Matrix<Tv> x = Matrix<Tv>(new std::vector<Tv>(3 * n), n, 3);

  std::mt19937 reg = std::mt19937(340);
  auto dis = Distribution<DistributionType::kUniformDis>(0, NumCutoff);
  dis.GetSample(y.Data, n, 340);
  for (Ti i = 0; i < n; i++) {
    x.Set(i, 0, 1.0);
    x.Set(i, 1, sqrt(dis.GetSample1(reg) * dis.GetSample1(reg)));
    x.Set(i, 2, sqrt(dis.GetSample1(reg) * dis.GetSample1(reg)));
  }
  // auto ystr = y.ToString0();
  // auto xstr = x.ToString0();

  auto model = DiscreteChoice<DiscreteChoiceModelType::kOrdered,
                              DiscreteChoiceDistType::kProbit>(
      y.RowsCount, x.ColsCount, NumCutoff + 1, true);
  auto W = new Tv[model.WorkSize];
  auto S = new Tv[model.StorageSize];
  model.Calculate(y, x, nullptr, S, W);

  ASSERT_NEAR(model.Beta.Data[0], 1.711136501632428, 1e-5);

  delete[] W;
  delete[] S;
}

template <DiscreteChoiceModelType type, DiscreteChoiceDistType dist>
void test_weight() {

  // first three observations are the same
  const Matrix<Tv> y =
      Matrix<Tv>(type == DiscreteChoiceModelType::kBinary
                     ? new std::vector<Tv>{1, 1, 1, 0, 1, 1, 1, 0, 0, 0,
                                           1, 1, 0, 1, 0, 0, 0, 1, 1, 0}
                     : new std::vector<Tv>{2, 2, 2, 0, 1, 2, 1, 0, 2, 0,
                                           1, 2, 0, 1, 0, 2, 0, 1, 1, 0},
                 20, 1);
  const Matrix<Tv> x = Matrix<Tv>(
      new std::vector<Tv>{1, 1, 1, 1, 1, 1, 1, 1,  1, 1,  1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 7, 7, 7,  9, 7,  3, 7, 8, 5, 3,
                          5, 3, 3, 4, 6, 5, 8, 3,  6, 10, 2, 2, 2, 6, 4,
                          3, 5, 6, 7, 4, 5, 6, 63, 0, 3,  2, 4, 3, 4, 5},
      20, 3);
  const Matrix<Tv> w =
      Matrix<Tv>(new std::vector<Tv>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                 20, 1);

  auto numChoices = type == DiscreteChoiceModelType::kBinary ? 2 : 3;
  auto model1 =
      DiscreteChoice<type, dist>(y.RowsCount, x.ColsCount, numChoices, true);
  auto W = new Tv[model1.WorkSize];
  auto S = new Tv[model1.StorageSize];
  model1.Calculate(y, x, nullptr, S, W, numChoices);

  // estimate with weight 1
  auto model2 =
      DiscreteChoice<type, dist>(y.RowsCount, x.ColsCount, numChoices, true);
  auto W1 = new Tv[model2.WorkSize];
  auto S1 = new Tv[model2.StorageSize];
  model2.Calculate(y, x, &w, S1, W1, numChoices);

  ASSERT_EQ(true, model1.Beta.Equals(model2.Beta));
  ASSERT_EQ(true, model1.BetaVar.Equals(model2.BetaVar));

  // combine first 3 observations and give a weight
  const Matrix<Tv> y2 =
      Matrix<Tv>(type == DiscreteChoiceModelType::kBinary
                     ? new std::vector<Tv>{1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1,
                                           0, 0, 0, 1, 1, 0}
                     : new std::vector<Tv>{2, 0, 1, 2, 1, 0, 2, 0, 1, 2, 0, 1,
                                           0, 2, 0, 1, 1, 0},
                 18, 1);
  const Matrix<Tv> x2 =
      Matrix<Tv>(new std::vector<Tv>{1, 1, 1, 1, 1,  1, 1, 1,  1, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 7,  9, 7, 3,  7, 8, 5, 3, 5, 3,
                                     3, 4, 6, 5, 8,  3, 6, 10, 2, 6, 4, 3, 5, 6,
                                     7, 4, 5, 6, 63, 0, 3, 2,  4, 3, 4, 5},
                 18, 3);
  const Matrix<Tv> w2 = Matrix<Tv>(new std::vector<Tv>{
      3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});

  auto model3 =
      DiscreteChoice<type, dist>(y.RowsCount, x.ColsCount, numChoices, true);
  auto W2 = new Tv[model3.WorkSize];
  auto S2 = new Tv[model3.StorageSize];
  model3.Calculate(y, x, &w, S2, W2, numChoices);

  ASSERT_EQ(true, model1.Beta.Equals(model3.Beta, 1e-14));
  ASSERT_EQ(true, model1.BetaVar.Equals(model3.BetaVar, 1e-12));
  ASSERT_NEAR(model1.LogL, model3.LogL, 1e-14);

  delete[] W;
  delete[] S;
  delete[] W1;
  delete[] S1;
  delete[] W2;
  delete[] S2;
}

TEST(DiscreteChoice_t, with_weight) {
  test_weight<DiscreteChoiceModelType::kBinary,
              DiscreteChoiceDistType::kLogit>();
  test_weight<DiscreteChoiceModelType::kBinary,
              DiscreteChoiceDistType::kProbit>();
  test_weight<DiscreteChoiceModelType::kOrdered,
              DiscreteChoiceDistType::kLogit>();
  test_weight<DiscreteChoiceModelType::kOrdered,
              DiscreteChoiceDistType::kProbit>();
}

TEST(DiscreteChoice_t, CrossValidate) {
  Matrix<Tv> source = Matrix<Tv>(
      new std::vector<Tv>{
          0,   1,   1,   1,   0,   1,   0,   0,   0,   0,   0,   1,   1,   1,
          0,   1,   1,   1,   1,   1,   0,   1,   1,   1,   0,   1,   1,   0,
          1,   0,   0,   0,   0,   0,   1,   1,   0,   1,   0,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0.1, 0.9, 0.2, 0.4,
          0.7, 0.2, 0.7, 0.9, 0.4, 0.1, 0.6, 0.1, 0.2, 0.6, 0.8, 0.6, 0.7, 0.7,
          0.9, 0.5, 0.3, 0.6, 0.2, 0.4, 0.8, 0.6, 0.5, 0.1, 0.6, 0.9, 0.5, 0.8,
          0.5, 0.8, 0,   0.2, 0.2, 0.3, 0.4, 0.2, 0.5, 0.1, 0.8, 0.3, 1,   0.4,
          0.6, 0.4, 1,   0.3, 0.9, 0.1, 0,   0.8, 0.3, 0.1, 0.7, 0.1, 0.9, 0.7,
          0.5, 0,   0.3, 0.1, 0.9, 1,   0.6, 0.9, 0.3, 0.8, 0.5, 0.6, 0.5, 0.1,
          0.9, 0,   0.7, 0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.3, 0.8, 0.4, 0.9, 0.8,
          0.6, 0.2, 0.3, 0.2, 0.1, 0.9, 0.1, 0.1, 0,   0.3, 0.8, 0.7, 0.4, 0.1,
          0.8, 0.5, 0.8, 1,   0.7, 0.5, 0.9, 0.4, 1,   0.6, 0.2, 0.2, 0.2, 0.4,
          0.1, 0.9, 0.1, 0.4, 0.7, 0.8, 0.8, 0.4, 0.6, 0.3, 0.8, 0.3, 0.7, 0.1,
          0.1, 0.2, 0.2, 0.5, 0.8, 0.9, 0.3, 0.9, 0.9, 0.9, 0.5, 0.4, 0.7, 0.6,
          0.2, 0.9, 0.4, 0.3, 0.7, 0.8, 0.1, 0.5, 0,   0.1, 0.9, 0,   0.2, 0.1,
          0.9, 0.7,
      },
      40, 6);

  auto trainRatio = 0.5;
  Matrix<Tv> cost_table =
      Matrix<Tv>(new Tv[10]{0.2, 1.0, 1, 0, 1, 0, 1, 0, 1, 0}, 2, 4);
  auto ct_vec = std::vector<Matrix<Tv>>({cost_table, cost_table, cost_table});

  auto model = DiscreteChoiceSim<true, DiscreteChoiceModelType::kBinary,
                                 DiscreteChoiceDistType::kLogit>(
      source.RowsCount, source.ColsCount, 2, trainRatio, 0, ct_vec.size(), true,
      true, nullptr);
  model.Seed = 340;
  model.SimulationMax = 100;
  auto W = new Tv[model.WorkSize];
  auto S = new Tv[model.StorageSize];
  auto Wi = new Ti[model.WorkSizeI];
  bool cancel = false;
  model.Calculate(source, &ct_vec, S, W, Wi, cancel);

  ASSERT_NEAR(model.CostRatios.Data[0], 0.3320, 1e-14); // self
  ASSERT_NEAR(model.CostRatios.Data[0], model.CostRatios.Data[1], 1e-16);
  ASSERT_NEAR(model.CostRatios.Data[0], model.CostRatios.Data[2], 1e-16);

  // how can I test the simulation?!

  delete[] S;
  delete[] Wi;
  delete[] W;
}

TEST(DiscreteChoice_t, CrossValidate_pca) {
  Matrix<Tv> source = Matrix<Tv>(
      new std::vector<Tv>{
          1,   0,   1,   0,   1,   1,   0,   0,   1,   1,   0,   1,   0,   0,
          1,   0,   1,   0,   0,   0,   0,   1,   0,   1,   0,   1,   1,   0,
          1,   1,   1,   0,   1,   1,   0,   1,   1,   1,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2,   4.7, 9.9, 3,
          9.9, 7.5, 2.9, 2.6, 2.4, 2.8, 3.7, 0.6, 9.8, 0.7, 5.8, 6.3, 6.2, 5.1,
          8.1, 4,   3.2, 1.6, 7.9, 8.9, 4.6, 2.2, 7.3, 8.5, 3.1, 4.6, 4.7, 4.2,
          9.6, 3,   3.1, 9.6, 9,   2.6, 1.5, 1.9, 4.7, 6.3, 6.1, 7.2, 1,   6.4,
          7.1, 1.6, 2.9, 8.1, 2,   1,   9,   2.4, 4.9, 0.4, 6.8, 1,   8.4, 7.8,
          0.2, 4.5, 9.5, 7.8, 8,   9.5, 0.8, 6.9, 7.6, 6.2, 7.3, 3.7, 3.7, 6.2,
          9.1, 6.9, 5.5, 0.7, 7.9, 8.1, 5.1, 5.7, 5.6, 5,   2.6, 3.7, 5.2, 7.3,
          8.1, 2.7, 2.8, 8,   4.2, 3.7, 1.7, 0.3, 6.8, 8.3, 9.1, 1.5, 5.6, 0.2,
          4,   6.7, 0.8, 8.3, 7.4, 2.8, 3.2, 7.9, 4.9, 6.1, 3.2, 0.2, 2.7, 4.5,
          8.9, 4.5, 5.2, 4.3, 4.3, 3.6, 8.8, 4.9, 3.4, 9.4, 9.5, 9.9, 4.5, 9,
          0.2, 7.1, 1.6, 6.3, 2.6, 4.8, 1.7, 8.9, 4.2, 4.3, 6.3, 2.5, 4.1, 0.2,
          7.8, 6.7, 4.8, 3.9, 1,   7.2, 1.4, 2.2, 3.1, 4.9, 7.3, 2,   1.5, 9.1,
          3.3, 1.4, 8.8, 7.4, 8.3, 4.2, 7.7, 7.8, 4.5, 6.9, 7.6, 7.8, 2.8, 0,
          0.3, 1,   7.9, 1.8, 9.5, 8.6, 7.1, 6.7, 4.6, 5.6, 8.9, 7.3, 7.6, 8,
          4.1, 6.7, 9.1, 0.6, 3.4, 2.6, 5.6, 2.8, 5.5, 8.9, 4.1, 9.8, 2.4, 10,
          2.2, 7.9, 6.3, 3.2, 4.4, 1.3, 1.9, 2.1, 6.9, 5.5, 2.5, 2.9, 2,   6.6,
          6.9, 3.3, 8.5, 1.1, 6.6, 9.8, 3.7, 7.8, 7.8, 7,   4.5, 0.2, 4.4, 7.5,
          5.5, 9,   7.7, 9.6, 6.2, 7.1, 4.5, 6.7, 9.3, 5.7, 2.8, 6.3, 6,   2,
          9.6, 3.6, 8.1, 0.5, 0.6, 9,   5.9, 2.2, 0.4, 3,   3.2, 2.6, 6.3, 4.4,
          0.1, 9.3, 8.6, 4.2, 6.9, 4.4, 6.6, 7.9, 7.8, 2.8, 6.7, 4.8, 2.1, 0.2,
          4,   2.7, 8.1, 0.3, 3.8, 1.7, 8.1, 4.2, 8.9, 4,   9.6, 5.4, 1.6, 3.3,
          8.6, 4.7, 6.4, 3.6, 0.8, 0.1, 6.6, 10,  7,   8.6, 9.4, 1.6, 5.3, 1.8,
          5.8, 5.3, 3.9, 5.8, 6.9, 5.7, 8.1, 7.7, 4.8, 2.9, 8.9, 4.6, 6.3, 5.9,
          0.3, 4.7, 4.6, 1.4, 2.8, 3.4, 9.3, 3.1, 6.4, 3.1, 4,   2.1, 3.6, 4,
          5.7, 9.6, 0.9, 2.2, 6.4, 6.3, 8.4, 3.6, 0.1, 3.8, 8.6, 5.9, 5.1, 0.9,
          9.9, 3.5, 4.7, 6.3, 3.8, 0.7, 9.6, 2.3, 4,   8,   8.8, 1.9, 2.4, 2.7,
          8.9, 0.8, 1.7, 4.8, 2.1, 4.2, 2,   5.2, 2,   6.7, 8.1, 9.9, 7.4, 6.7,
          4.6, 0.8, 2.9, 6.8, 8.3, 8.4, 3.3, 2.8, 8.4, 7.7, 8.1, 2.9, 4.4, 3.4,
          5.5, 3,   3.3, 4.5, 3.9, 5.8, 4.7, 9.4, 2.1, 3.8, 3.8, 0.9, 8.9, 4.7,
          1.6, 8,   8.1, 3.9, 3.3, 2.9, 7.7, 3.3, 5.1, 2,   3.6, 4,   9.4, 6.1,
          2.3, 4,   6.6, 9.1, 9.8, 8.6, 9.6, 9,   3.4, 3.6, 6,   1.7, 2,   5.3,
          8.6, 9.1, 3,   7.4, 9.5, 8.6, 1.3, 0.4, 9.6, 4.3, 9.6, 1.5, 5.7, 2.3,
          7.9, 1.3, 0.8, 9.9, 6.4, 8.6, 7,   8.1, 7.2, 4.6, 1.1, 0.3, 5,   2.5,
          4.6, 4.5, 7.8, 1.6, 1,   2.4, 4.6, 6,   9,   6.2, 4.5, 0.8, 4.8, 1.9,
          4.4, 2.2, 9.6, 8.3, 2.9, 7.4, 1.2, 2,   2.8, 6.2, 5.3, 3.3, 8.3, 0.6,
          1.9, 8.4, 6.3, 9.1, 3.5, 8.8, 5.7, 5.1, 6.5, 8,   0.5, 3.3, 7.7, 1.8,
          2.3, 4.6, 5.3, 5.5, 3.9, 5.1, 3.6, 5,   5,   2.8, 7.5, 7.7, 1.4, 2.6,
          1.7, 1.7, 6.3, 7.7, 4.2, 4,   4.2, 0.6, 9.9, 5.4, 7.3, 6,   6.2, 9.4,
          5.8, 1.6, 7.8, 9.6, 9.6, 9.4, 9.6, 8.1, 0.8, 9.8, 9.3, 3.2, 7,   6.9,
          3.4, 3.5, 1.7, 5.9, 7.4, 3.1, 3.4, 4.8, 1.7, 2.8, 9.1, 3.7, 7.6, 1.8,
          4.7, 10,  8.1, 7.9, 8.9, 4.3, 9.1, 3.5, 5,   7.8, 8,   2,   2.9, 3.5,
          0.8, 6.6, 8.4, 1,   7,   7.3, 7.2, 6.4, 3.7, 1.5, 8.6, 3.1, 6,   0.9,
          9.9, 5,   0.5, 6,   1.7, 7.9, 3.5, 0.9, 1.6, 1.4, 4.3, 2.7, 7.7, 9.9,
          8.5, 3.9, 7.9, 5.4, 0.2, 9.8, 4,   5.1, 4.4, 8.1, 6.6, 0.3, 4.5, 8.3,
          6.3, 6.7, 9.3, 5.8, 8.3, 8.3, 3.5, 10,  7.4, 1.9, 0.3, 5,   0,   10,
          5,   6.9, 6.6, 4.3, 9.6, 7.3, 8,   9.8, 3.6, 3.2, 2.9, 6.2, 6.9, 4.8,
          0.6, 2.3, 9.4, 1.4, 3.4, 5.2, 4.6, 4.9, 2.1, 8.5, 4.4, 9.4, 6.4, 8.5,
          6.5, 6,   5.1, 5.2, 9.8, 4,   0.1, 7.5, 1.5, 2.3, 0.8, 7.2, 0.3, 6.5,
          3.2, 6.1, 2.9, 2.1, 1.3, 8.6, 9.4, 5,   2.5, 8.9, 1.6, 4.8, 5.4, 0.1,
          4.9, 8.1, 9.9, 0.1, 1,   9.5, 3.3, 7.6, 8.1, 7.6, 6.3, 2.9, 4.7, 3.1,
          0.8, 7.1, 0.6, 6.2, 4.9, 8,   1.7, 6,   5.6, 2.3, 8.1, 7.7, 9.1, 1.5,
          3.5, 6.2, 3.6, 5,   8.6, 10,  2.6, 3.4, 2.4, 0.6, 6.3, 0.8, 3.2, 5.6,
          8,   5.3, 5.8, 8.4, 4,   6.3, 4.2, 5.4, 8.3, 5.4, 9.9, 7.4, 8.4, 8.9,
          2.5, 7.6, 2,   3.3, 0.4, 1.8, 4.1, 1.8, 4,   1.3, 4.3, 4.7, 0.5, 7.4,
          4.7, 10,  9.9, 4.1, 4.8, 6.1, 7.9, 0.3, 6.6, 9.6, 3.7, 2.8, 7.2, 8.1,
          8,   8.1, 5.1, 1.3, 1.6, 3.7, 0.8, 4.8, 8.6, 3,   0.4, 1.4, 1.5, 6.2,
          1.3, 4.4, 4.2, 9.3, 5.8, 6.8, 0.6, 8.3, 5.6, 5.8, 2.5, 2.3, 1,   2.3,
          1,   2.3, 9.6, 7.6, 4.6, 6,   3.9, 5.9, 6.2, 3,   3.1, 6.8, 6.6, 9.9,
          1.4, 6.9, 3.7, 3.5, 4.6, 1.7, 9,   2.1, 5.6, 4.4, 7,   0.5, 7.5, 9.6,
          2.6, 9.3, 8.8, 9.3, 2.8, 8.6, 3.8, 6.4, 3.3, 2.8, 6.8, 5.5, 2.3, 5.4,
          8.1, 1.2, 6.1, 9,   0.2, 0.7, 6.9, 8.6, 2.1, 9,   7.7, 9.1, 3.2, 8.9,
          6.7, 5.3, 0.6, 4.2, 8.5, 4.1, 3.3, 6.9, 4.2, 3.9, 4,   8,   6.2, 5,
          3.6, 0.9, 6.7, 9.1, 9.4, 3.6, 3.9, 4.4, 8.3, 8,   0.2, 6.8, 2.8, 4,
          7.2, 9.6, 1.5, 3.5, 5.3, 4.8, 9.6, 8.1, 2.2, 2.6, 7.1, 8.1, 2.1, 1.3,
          7,   2.2, 5,   2.2, 4.5, 1.9, 5,   1.3, 2.3, 8.6, 9.8, 9.6, 5,   2.3,
          7,   2.6, 3.9, 1.6, 1.3, 9.3, 9.6, 4.5, 2.4, 0.5, 2.7, 1.4, 5.3, 5.5,
          8.6, 4.1, 2,   1.2, 1.5, 2.6, 7.4, 6,   3.1, 8,   2.3, 0.9, 8.9, 5.9,
          7.4, 2.9, 3.7, 5.7, 9.8, 7.4, 1.1, 9.6, 8.6, 5.4, 6.9, 2.8, 1.6, 2.5,
          9.9, 4.5, 9.5, 3.1, 6.6, 9.5, 2.6, 8.1, 6.7, 8.6, 8.1, 2.5, 4.9, 1,
          9.5, 5.9, 9.3, 0.4, 0.4, 3.8, 7,   9.5, 6.8, 7,   1.6, 7.9, 4.4, 9,
          5.4, 6.5, 4.3, 8.6, 0.9, 4.8, 6.2, 8.6, 8.5, 9.1, 4.3, 9.4, 3.3, 9.1,
          6.2, 5,   7.6, 7.3, 2.4, 0.1, 0.2, 2.8, 6.1, 6.5, 2.5, 7.8, 5.5, 0.4,
          7.7, 2.6, 1.6, 8.6, 5.1, 3.6, 2.8, 0.1, 4.4, 9.5, 1.5, 0.1, 4.8, 0.9,
          7.3, 0.1, 6.2, 1,   0.6, 7.8, 3.2, 1.8, 9.4, 4.2, 5.5, 4.3, 8.8, 2.6,
          0.4, 7.8, 6.2, 3.2, 9.3, 5.7, 0.9, 6,   7.6, 5.5, 4.3, 9.6, 6.4, 6.1,
          0.7, 8.1, 5.8, 2.3, 5.3, 9.3, 5.5, 6,   0.3, 1.5, 7.2, 7.7, 7.2, 5.9,
          3.2, 7.3, 3.3, 0.8, 4.8, 0.2, 9.6, 6.6, 3.3, 8.4, 0.8, 9.1, 9.9, 1.4,
          2.5, 0.6, 6.7, 1,   1.4, 7.1, 1.1, 1.3, 8.5, 3.1, 9,   9.5, 7.1, 4.5,
          7.4, 4.7, 5.9, 1.6, 8.6, 5,   2.8, 0.3, 2.2, 3.4, 4.8, 6.2, 6.4, 4.1,
          4.8, 7.4, 1.5, 6.2, 9.4, 4.1, 3.5, 9.6, 4.8, 6.3, 6.3, 0,   4.8, 1.8,
          8.9, 0.1, 1.8, 5.3, 3.1, 5.2, 9.1, 8.6, 2.4, 9.1, 0.5, 4.9, 8.6, 2.2,
          8.4, 3.3, 6.2, 5.1, 0.7, 7,   6.6, 1.1, 1.6, 2.7, 0.2, 3,   1,   7.7,
          9.9, 9.6, 5.3, 8.6, 0.2, 1,   1.4, 2.1, 1.6, 9.4, 3.1, 1.2, 5,   6.1,
          5,   2.4, 5.8, 4.5, 5.5, 1.2, 1.3, 1.5, 0.1, 1.2, 9.3, 7.9, 5.3, 7.8,
          8.5, 9.1, 5.2, 0.1, 4.2, 1.4, 8.6, 2.7, 3.8, 1.9, 4.7, 7.6, 0.8, 6.3,
          1.7, 6.7, 2,   9.5, 2,   0.4, 6.1, 9.1, 7,   5.7, 7.5, 9.8, 6.4, 8.5,
          9.6, 4.2, 9.3, 2.3, 5.7, 3.7},
      40, 32);

  auto projX = Matrix<Tv>(&source.Data[80], 2, 31);

  auto pcaoptions = PcaAnalysisOptions();
  pcaoptions.IgnoreFirstCount = 1;
  pcaoptions.ExactCount = 3;

  auto trainRatio = 0.5;
  Matrix<Tv> cost_table = Matrix<Tv>(new Tv[6]{0.2, 1.0, 1, 0, 1, 0}, 2, 3);
  auto ct_vec = std::vector<Matrix<Tv>>({cost_table, cost_table, cost_table});

  auto model = DiscreteChoiceSim<false, DiscreteChoiceModelType::kBinary,
                                 DiscreteChoiceDistType::kLogit>(
      source.RowsCount, source.ColsCount, 2, trainRatio, 0, ct_vec.size(), true,
      true, &pcaoptions);
  model.Seed = 340;
  model.SimulationMax = 100;
  auto W = new Tv[model.WorkSize];
  auto S = new Tv[model.StorageSize];
  auto Wi = new Ti[model.WorkSizeI];
  bool cancel = false;
  model.Calculate(source, &ct_vec, S, W, Wi, cancel);

  delete[] S;
  delete[] Wi;
  delete[] W;
}

TEST(DiscreteChoice_t, searcherSmall) {
  Matrix<Tv> source = Matrix<Tv>(
      new std::vector<Tv>{
          1,  2,  3,  0,  2,  1,  3,  0,  0,  0,  1,  2,  0,  1,  0,  2,
          0,  1,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
          1,  1,  1,  1,  1,  1,  1,  1,  7,  5,  8,  9,  7,  3,  7,  8,
          5,  3,  5,  3,  3,  4,  6,  5,  8,  3,  6,  10, 17, 15, 18, 19,
          17, 13, 17, 18, 15, 13, 15, 13, 13, 14, 16, 15, 18, 13, 16, 110},
      20, 4);

  auto items = SearchItems();
  auto searchOptions = SearchOptions();
  auto measures = SearchMeasureOptions();
  auto checks = SearchModelChecks();

  items.KeepBestCount = 2;
  items.KeepAll = true;
  items.KeepInclusionWeights = true;
  items.ExtremeBoundsMultiplier = 2;
  items.CdfsAt = std::vector<Tv>({0});
  items.KeepMixture = true;
  items.Length1 = 2 + 1 + 4 - 2;

  measures.TrainRatio = 0.8;
  measures.SimFixSize = 10;
  measures.MeasuresOut.push_back(ScoringType::kFrequencyCost);

  checks.Estimation = true;

  auto sizes = std::vector<Ti>({1, 2});

  // groups
  auto gr1 = std::vector<Ti>({0});
  auto gr2 = std::vector<Ti>({1});

  // cost tables
  auto cost1 = Matrix<Tv>(new Tv[10]{0.15, 1, 0, 1, 0, 1, 0, 1, 0, 1}, 2,
                          5); // new? the owner will be the model set
  auto cost2 = Matrix<Tv>(new Tv[10]{0.15, 1, 0, 1, 0, 1, 0, 1, 0, 1}, 2, 5);
  auto costs = std::vector<Matrix<Tv>>({cost1, cost2});
  auto groups = std::vector<std::vector<Ti>>({gr1, gr2});
  auto newton = Newton();
  auto modelset =
      DiscreteChoiceModelset<false, DiscreteChoiceModelType::kOrdered>(
          searchOptions, items, measures, checks, sizes, source, costs, groups,
          newton, true, true);

  auto W = new Tv[modelset.Modelset.WorkSize];
  auto Wi = new Ti[modelset.Modelset.WorkSizeI];

  modelset.Modelset.Start(W, Wi);

  auto summaries = std::vector<SearcherSummary *>();
  auto list1 = std::vector<SearcherSummary *>();
  auto list2 = std::vector<SearcherSummary *>();
  auto result = SearcherModelingInfo();
  modelset.Modelset.CombineInfo(result, summaries, list1, list2);

  ASSERT_EQ(6, result.SearchedCount);
  ASSERT_EQ(6, result.ExpectedCount);

  auto all = std::vector<EstimationKeep *>();
  modelset.Modelset.CombineAll(0, 0, 0, summaries, all);

  auto bests = std::vector<EstimationKeep *>();
  modelset.Modelset.CombineBests(0, 0, 0, summaries, bests);
}

TEST(DiscreteChoice_t, searcherLarager) {
  Matrix<Tv> source = Matrix<Tv>(
      new std::vector<Tv>{
          1,   0,   1,   0,   1,   1,   0,   0,   1,   1,   0,   1,   0,   0,
          1,   0,   1,   0,   0,   0,   0,   1,   0,   1,   0,   1,   1,   0,
          1,   1,   1,   0,   1,   1,   0,   1,   1,   1,   1,   1,   2,   3,
          4,   2,   3,   4,   5,   6,   7,   8,   9,   4,   5,   3,   4,   7,
          8,   4,   8,   9,   3,   5,   7,   8,   4,   5,   6,   7,   3,   4,
          6,   7,   3,   5,   6,   7,   8,   2,   4,   6,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   2,   4.7, 9.9, 3,   9.9, 7.5,
          2.9, 2.6, 2.4, 2.8, 3.7, 0.6, 9.8, 0.7, 5.8, 6.3, 6.2, 5.1, 8.1, 4,
          3.2, 1.6, 7.9, 8.9, 4.6, 2.2, 7.3, 8.5, 3.1, 4.6, 4.7, 4.2, 9.6, 3,
          3.1, 9.6, 9,   2.6, 1.5, 1.9, 4.7, 6.3, 6.1, 7.2, 1,   6.4, 7.1, 1.6,
          2.9, 8.1, 2,   1,   9,   2.4, 4.9, 0.4, 6.8, 1,   8.4, 7.8, 0.2, 4.5,
          9.5, 7.8, 8,   9.5, 0.8, 6.9, 7.6, 6.2, 7.3, 3.7, 3.7, 6.2, 9.1, 6.9,
          5.5, 0.7, 7.9, 8.1, 5.1, 5.7, 5.6, 5,   2.6, 3.7, 5.2, 7.3, 8.1, 2.7,
          2.8, 8,   4.2, 3.7, 1.7, 0.3, 6.8, 8.3, 9.1, 1.5, 5.6, 0.2, 4,   6.7,
          0.8, 8.3, 7.4, 2.8, 3.2, 7.9, 4.9, 6.1, 3.2, 0.2, 2.7, 4.5, 8.9, 4.5,
          5.2, 4.3, 4.3, 3.6, 8.8, 4.9, 3.4, 9.4, 9.5, 9.9, 4.5, 9,   0.2, 7.1,
          1.6, 6.3, 2.6, 4.8, 1.7, 8.9, 4.2, 4.3, 6.3, 2.5, 4.1, 0.2, 7.8, 6.7,
          4.8, 3.9, 1,   7.2, 1.4, 2.2, 3.1, 4.9, 7.3, 2,   1.5, 9.1, 3.3, 1.4,
          8.8, 7.4, 8.3, 4.2, 7.7, 7.8, 4.5, 6.9, 7.6, 7.8, 2.8, 0,   0.3, 1,
          7.9, 1.8, 9.5, 8.6, 7.1, 6.7, 4.6, 5.6, 8.9, 7.3, 7.6, 8,   4.1, 6.7,
          9.1, 0.6, 3.4, 2.6, 5.6, 2.8, 5.5, 8.9, 4.1, 9.8, 2.4, 10,  2.2, 7.9,
          6.3, 3.2, 4.4, 1.3, 1.9, 2.1, 6.9, 5.5, 2.5, 2.9, 2,   6.6, 6.9, 3.3,
          8.5, 1.1, 6.6, 9.8, 3.7, 7.8, 7.8, 7,   4.5, 0.2, 4.4, 7.5, 5.5, 9,
          7.7, 9.6, 6.2, 7.1, 4.5, 6.7, 9.3, 5.7, 2.8, 6.3, 6,   2,   9.6, 3.6,
          8.1, 0.5, 0.6, 9,   5.9, 2.2, 0.4, 3,   3.2, 2.6, 6.3, 4.4, 0.1, 9.3,
          8.6, 4.2, 6.9, 4.4, 6.6, 7.9, 7.8, 2.8, 6.7, 4.8, 2.1, 0.2, 4,   2.7,
          8.1, 0.3, 3.8, 1.7, 8.1, 4.2, 8.9, 4,   9.6, 5.4, 1.6, 3.3, 8.6, 4.7,
          6.4, 3.6, 0.8, 0.1, 6.6, 10,  7,   8.6, 9.4, 1.6, 5.3, 1.8, 5.8, 5.3,
          3.9, 5.8, 6.9, 5.7, 8.1, 7.7, 4.8, 2.9, 8.9, 4.6, 6.3, 5.9, 0.3, 4.7,
          4.6, 1.4, 2.8, 3.4, 9.3, 3.1, 6.4, 3.1, 4,   2.1, 3.6, 4,   5.7, 9.6,
          0.9, 2.2, 6.4, 6.3, 8.4, 3.6, 0.1, 3.8, 8.6, 5.9, 5.1, 0.9, 9.9, 3.5,
          4.7, 6.3, 3.8, 0.7, 9.6, 2.3, 4,   8,   8.8, 1.9, 2.4, 2.7, 8.9, 0.8,
          1.7, 4.8, 2.1, 4.2, 2,   5.2, 2,   6.7, 8.1, 9.9, 7.4, 6.7, 4.6, 0.8,
          2.9, 6.8, 8.3, 8.4, 3.3, 2.8, 8.4, 7.7, 8.1, 2.9, 4.4, 3.4, 5.5, 3,
          3.3, 4.5, 3.9, 5.8, 4.7, 9.4, 2.1, 3.8, 3.8, 0.9, 8.9, 4.7, 1.6, 8,
          8.1, 3.9, 3.3, 2.9, 7.7, 3.3, 5.1, 2,   3.6, 4,   9.4, 6.1, 2.3, 4,
          6.6, 9.1, 9.8, 8.6, 9.6, 9,   3.4, 3.6, 6,   1.7, 2,   5.3, 8.6, 9.1,
          3,   7.4, 9.5, 8.6, 1.3, 0.4, 9.6, 4.3, 9.6, 1.5, 5.7, 2.3, 7.9, 1.3,
          0.8, 9.9, 6.4, 8.6, 7,   8.1, 7.2, 4.6, 1.1, 0.3, 5,   2.5, 4.6, 4.5,
          7.8, 1.6, 1,   2.4, 4.6, 6,   9,   6.2, 4.5, 0.8, 4.8, 1.9, 4.4, 2.2,
          9.6, 8.3, 2.9, 7.4, 1.2, 2,   2.8, 6.2, 5.3, 3.3, 8.3, 0.6, 1.9, 8.4,
          6.3, 9.1, 3.5, 8.8, 5.7, 5.1, 6.5, 8,   0.5, 3.3, 7.7, 1.8, 2.3, 4.6,
          5.3, 5.5, 3.9, 5.1, 3.6, 5,   5,   2.8, 7.5, 7.7, 1.4, 2.6, 1.7, 1.7,
          6.3, 7.7, 4.2, 4,   4.2, 0.6, 9.9, 5.4, 7.3, 6,   6.2, 9.4, 5.8, 1.6,
          7.8, 9.6, 9.6, 9.4, 9.6, 8.1, 0.8, 9.8, 9.3, 3.2, 7,   6.9, 3.4, 3.5,
          1.7, 5.9, 7.4, 3.1, 3.4, 4.8, 1.7, 2.8, 9.1, 3.7, 7.6, 1.8, 4.7, 10,
          8.1, 7.9, 8.9, 4.3, 9.1, 3.5, 5,   7.8, 8,   2,   2.9, 3.5, 0.8, 6.6,
          8.4, 1,   7,   7.3, 7.2, 6.4, 3.7, 1.5, 8.6, 3.1, 6,   0.9, 9.9, 5,
          0.5, 6,   1.7, 7.9, 3.5, 0.9, 1.6, 1.4, 4.3, 2.7, 7.7, 9.9, 8.5, 3.9,
          7.9, 5.4, 0.2, 9.8, 4,   5.1, 4.4, 8.1, 6.6, 0.3, 4.5, 8.3, 6.3, 6.7,
          9.3, 5.8, 8.3, 8.3, 3.5, 10,  7.4, 1.9, 0.3, 5,   0,   10,  5,   6.9,
          6.6, 4.3, 9.6, 7.3, 8,   9.8, 3.6, 3.2, 2.9, 6.2, 6.9, 4.8, 0.6, 2.3,
          9.4, 1.4, 3.4, 5.2, 4.6, 4.9, 2.1, 8.5, 4.4, 9.4, 6.4, 8.5, 6.5, 6,
          5.1, 5.2, 9.8, 4,   0.1, 7.5, 1.5, 2.3, 0.8, 7.2, 0.3, 6.5, 3.2, 6.1,
          2.9, 2.1, 1.3, 8.6, 9.4, 5,   2.5, 8.9, 1.6, 4.8, 5.4, 0.1, 4.9, 8.1,
          9.9, 0.1, 1,   9.5, 3.3, 7.6, 8.1, 7.6, 6.3, 2.9, 4.7, 3.1, 0.8, 7.1,
          0.6, 6.2, 4.9, 8,   1.7, 6,   5.6, 2.3, 8.1, 7.7, 9.1, 1.5, 3.5, 6.2,
          3.6, 5,   8.6, 10,  2.6, 3.4, 2.4, 0.6, 6.3, 0.8, 3.2, 5.6, 8,   5.3,
          5.8, 8.4, 4,   6.3, 4.2, 5.4, 8.3, 5.4, 9.9, 7.4, 8.4, 8.9, 2.5, 7.6,
          2,   3.3, 0.4, 1.8, 4.1, 1.8, 4,   1.3, 4.3, 4.7, 0.5, 7.4, 4.7, 10,
          9.9, 4.1, 4.8, 6.1, 7.9, 0.3, 6.6, 9.6, 3.7, 2.8, 7.2, 8.1, 8,   8.1,
          5.1, 1.3, 1.6, 3.7, 0.8, 4.8, 8.6, 3,   0.4, 1.4, 1.5, 6.2, 1.3, 4.4,
          4.2, 9.3, 5.8, 6.8, 0.6, 8.3, 5.6, 5.8, 2.5, 2.3, 1,   2.3, 1,   2.3,
          9.6, 7.6, 4.6, 6,   3.9, 5.9, 6.2, 3,   3.1, 6.8, 6.6, 9.9, 1.4, 6.9,
          3.7, 3.5, 4.6, 1.7, 9,   2.1, 5.6, 4.4, 7,   0.5, 7.5, 9.6, 2.6, 9.3,
          8.8, 9.3, 2.8, 8.6, 3.8, 6.4, 3.3, 2.8, 6.8, 5.5, 2.3, 5.4, 8.1, 1.2,
          6.1, 9,   0.2, 0.7, 6.9, 8.6, 2.1, 9,   7.7, 9.1, 3.2, 8.9, 6.7, 5.3,
          0.6, 4.2, 8.5, 4.1, 3.3, 6.9, 4.2, 3.9, 4,   8,   6.2, 5,   3.6, 0.9,
          6.7, 9.1, 9.4, 3.6, 3.9, 4.4, 8.3, 8,   0.2, 6.8, 2.8, 4,   7.2, 9.6,
          1.5, 3.5, 5.3, 4.8, 9.6, 8.1, 2.2, 2.6, 7.1, 8.1, 2.1, 1.3, 7,   2.2,
          5,   2.2, 4.5, 1.9, 5,   1.3, 2.3, 8.6, 9.8, 9.6, 5,   2.3, 7,   2.6,
          3.9, 1.6, 1.3, 9.3, 9.6, 4.5, 2.4, 0.5, 2.7, 1.4, 5.3, 5.5, 8.6, 4.1,
          2,   1.2, 1.5, 2.6, 7.4, 6,   3.1, 8,   2.3, 0.9, 8.9, 5.9, 7.4, 2.9,
          3.7, 5.7, 9.8, 7.4, 1.1, 9.6, 8.6, 5.4, 6.9, 2.8, 1.6, 2.5, 9.9, 4.5,
          9.5, 3.1, 6.6, 9.5, 2.6, 8.1, 6.7, 8.6, 8.1, 2.5, 4.9, 1,   9.5, 5.9,
          9.3, 0.4, 0.4, 3.8, 7,   9.5, 6.8, 7,   1.6, 7.9, 4.4, 9,   5.4, 6.5,
          4.3, 8.6, 0.9, 4.8, 6.2, 8.6, 8.5, 9.1, 4.3, 9.4, 3.3, 9.1, 6.2, 5,
          7.6, 7.3, 2.4, 0.1, 0.2, 2.8, 6.1, 6.5, 2.5, 7.8, 5.5, 0.4, 7.7, 2.6,
          1.6, 8.6, 5.1, 3.6, 2.8, 0.1, 4.4, 9.5, 1.5, 0.1, 4.8, 0.9, 7.3, 0.1,
          6.2, 1,   0.6, 7.8, 3.2, 1.8, 9.4, 4.2, 5.5, 4.3, 8.8, 2.6, 0.4, 7.8,
          6.2, 3.2, 9.3, 5.7, 0.9, 6,   7.6, 5.5, 4.3, 9.6, 6.4, 6.1, 0.7, 8.1,
          5.8, 2.3, 5.3, 9.3, 5.5, 6,   0.3, 1.5, 7.2, 7.7, 7.2, 5.9, 3.2, 7.3,
          3.3, 0.8, 4.8, 0.2, 9.6, 6.6, 3.3, 8.4, 0.8, 9.1, 9.9, 1.4, 2.5, 0.6,
          6.7, 1,   1.4, 7.1, 1.1, 1.3, 8.5, 3.1, 9,   9.5, 7.1, 4.5, 7.4, 4.7,
          5.9, 1.6, 8.6, 5,   2.8, 0.3, 2.2, 3.4, 4.8, 6.2, 6.4, 4.1, 4.8, 7.4,
          1.5, 6.2, 9.4, 4.1, 3.5, 9.6, 4.8, 6.3, 6.3, 0,   4.8, 1.8, 8.9, 0.1,
          1.8, 5.3, 3.1, 5.2, 9.1, 8.6, 2.4, 9.1, 0.5, 4.9, 8.6, 2.2, 8.4, 3.3,
          6.2, 5.1, 0.7, 7,   6.6, 1.1, 1.6, 2.7, 0.2, 3,   1,   7.7, 9.9, 9.6,
          5.3, 8.6, 0.2, 1,   1.4, 2.1, 1.6, 9.4, 3.1, 1.2, 5,   6.1, 5,   2.4,
          5.8, 4.5, 5.5, 1.2, 1.3, 1.5, 0.1, 1.2, 9.3, 7.9, 5.3, 7.8, 8.5, 9.1,
          5.2, 0.1, 4.2, 1.4, 8.6, 2.7, 3.8, 1.9, 4.7, 7.6, 0.8, 6.3, 1.7, 6.7,
          2,   9.5, 2,   0.4, 6.1, 9.1, 7,   5.7, 7.5, 9.8, 6.4, 8.5, 9.6, 4.2,
          9.3, 2.3, 5.7, 3.7},
      40, 32);

  auto items = SearchItems();
  auto searchOptions = SearchOptions();
  auto measures = SearchMeasureOptions();
  auto checks = SearchModelChecks();

  items.KeepBestCount = 4;
  items.KeepAll = true;
  items.KeepInclusionWeights = false;
  items.ExtremeBoundsMultiplier = 0;
  items.CdfsAt = std::vector<Tv>({});
  items.KeepMixture = false;

  auto sizes = std::vector<Ti>({2});

  measures.TrainRatio = 0.8;
  measures.Seed = -340;
  measures.SimFixSize = 200;
  measures.MeasuresIn.push_back(GoodnessOfFitType::kAuc);
  measures.MeasuresIn.push_back(GoodnessOfFitType::kFrequencyCost);
  measures.MeasuresOut.push_back(ScoringType::kFrequencyCost);
  measures.MeasuresOut.push_back(ScoringType::kAuc);

  auto gr1 = std::vector<Ti>({0});
  auto gr2 = std::vector<Ti>({1});

  auto cost1 = *new Matrix<Tv>(new Tv[6]{0.5, 1, 1, 0, 1, 0}, 2, 3);
  auto costs = std::vector<Matrix<Tv>>({cost1});
  auto groups = std::vector<std::vector<Ti>>({gr1, gr2});
  auto newton = Newton();
  auto modelset = *DiscreteChoiceModelsetBase::GetFromTypes(
      true, true, searchOptions, items, measures, checks, sizes, source, costs,
      groups, true, false, newton);

  // add groups

  auto W = new Tv[modelset.Modelset.WorkSize];
  auto Wi = new Ti[modelset.Modelset.WorkSizeI];

  modelset.Modelset.Start(W, Wi);

  auto summaries = std::vector<SearcherSummary *>();
  auto list11 = std::vector<SearcherSummary *>();
  auto list12 = std::vector<SearcherSummary *>();
  auto result = SearcherModelingInfo();
  modelset.Modelset.CombineInfo(result, summaries, list11, list12);

  auto all = std::vector<EstimationKeep *>();
  modelset.Modelset.CombineAll(1, 0, 0, summaries, all);

  auto bests = std::vector<EstimationKeep *>();
  modelset.Modelset.CombineBests(0, 0, 0, summaries, bests);
}

TEST(DiscreteChoice_t, searcherNan) {
  Matrix<Tv> source1 = Matrix<Tv>(
      new std::vector<Tv>{
          1,   0,   1,   0,   1,   1,   0,   0,   1,   1,   0,   1,   0,   0,
          1,   0,   1,   0,   0,   0,   0,   1,   0,   1,   0,   1,   1,   0,
          1,   1,   1,   0,   1,   1,   0,   1,   1,   1,   1,   1,   2,   3,
          4,   2,   3,   4,   5,   6,   7,   8,   9,   4,   5,   3,   4,   7,
          8,   4,   8,   9,   3,   5,   7,   8,   4,   5,   6,   7,   3,   4,
          6,   7,   3,   5,   6,   7,   8,   2,   4,   6,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   9,   9,   6.6, 8.2, 3.5, 2.4,
          0.7, 0.1, 8.7, 2.9, 1.1, 9,   6.8, 1.6, 9.9, 0.2, 6.1, 9.9, 10,  2.1,
          8.1, 7.9, 2.4, 1.6, 7.2, 0.8, 9.4, 8.2, 9.4, 3.2, 4.9, 2.1, 8.7, 7.4,
          9.7, 0.4, 8.3, 0.2, 7.5, 0.1, 6.7, 9.4, 2.2, 5.7, 3.6, 1.1, 6.8, 6.2,
          9.5, 6.3, 5.4, 7.7, 8.1, 7.3, 2.8, 9.9, 3.8, 8.7, 0.3, 7.2, 8.4, 7.7,
          6.7, 1.8, 7.2, 2.8, 4.6, 2.4, 2.5, 7.5, 8,   2.3, 4.3, 1,   1.4, 5.2,
          8.8, 4.8, 5.2, 4.7},
      40, 5);
  Matrix<Tv> source2 = Matrix<Tv>(
      new std::vector<Tv>{
          1,   NAN, 0,   1,    0,   1,   1,   0,    0,   1,    1,   0,   1,
          0,   0,   1,   0,    1,   0,   0,   0,    0,   1,    0,   1,   0,
          1,   1,   0,   1,    1,   1,   0,   1,    1,   0,    1,   1,   1,
          1,   1,   2,   1000, 3,   4,   2,   3,    4,   5,    6,   7,   8,
          9,   4,   5,   3,    4,   7,   8,   4,    8,   9,    3,   5,   7,
          8,   4,   5,   6,    7,   3,   4,   6,    7,   3,    5,   6,   7,
          8,   2,   4,   6,    1,   1,   1,   1,    1,   1,    1,   1,   1,
          1,   1,   1,   1,    1,   1,   1,   1,    1,   1,    1,   1,   1,
          1,   1,   1,   1,    1,   1,   1,   1,    1,   1,    1,   1,   1,
          1,   1,   1,   1,    1,   1,   9,   1000, 9,   6.6,  8.2, 3.5, 2.4,
          0.7, 0.1, 8.7, 2.9,  1.1, 9,   6.8, 1.6,  9.9, 0.2,  6.1, 9.9, 10,
          2.1, 8.1, 7.9, 2.4,  1.6, 7.2, 0.8, 9.4,  8.2, 9.4,  3.2, 4.9, 2.1,
          8.7, 7.4, 9.7, 0.4,  8.3, 0.2, 7.5, 0.1,  6.7, 1000, 9.4, 2.2, 5.7,
          3.6, 1.1, 6.8, 6.2,  9.5, 6.3, 5.4, 7.7,  8.1, 7.3,  2.8, 9.9, 3.8,
          8.7, 0.3, 7.2, 8.4,  7.7, 6.7, 1.8, 7.2,  2.8, 4.6,  2.4, 2.5, 7.5,
          8,   2.3, 4.3, 1,    1.4, 5.2, 8.8, 4.8,  5.2, 4.7},
      41, 5);

  auto items = SearchItems();
  auto searchOptions = SearchOptions();
  auto measures = SearchMeasureOptions();
  auto checks = SearchModelChecks();

  items.KeepBestCount = 4;
  items.KeepAll = true;
  items.KeepInclusionWeights = false;
  items.ExtremeBoundsMultiplier = 0;
  items.CdfsAt = std::vector<Tv>({});
  items.KeepMixture = false;

  auto sizes = std::vector<Ti>({2});

  measures.TrainRatio = 0.8;
  measures.Seed = 340;
  measures.SimFixSize = 0;
  measures.MeasuresIn.push_back(GoodnessOfFitType::kAic);

  auto gr1 = std::vector<Ti>({0});
  auto gr2 = std::vector<Ti>({1});
  auto groups = std::vector<std::vector<Ti>>({gr1, gr2});
  auto newton = Newton();
  auto vms = std::vector<Matrix<Tv>>();
  auto modelset1 = *DiscreteChoiceModelsetBase::GetFromTypes(
      true, true, searchOptions, items, measures, checks, sizes, source1, vms,
      groups, true, false, newton);
  auto modelset2 = *DiscreteChoiceModelsetBase::GetFromTypes(
      true, true, searchOptions, items, measures, checks, sizes, source2, vms,
      groups, true, false, newton);

  // add groups

  auto W1 = new Tv[modelset1.Modelset.WorkSize];
  auto Wi1 = new Ti[modelset1.Modelset.WorkSizeI];

  auto W2 = new Tv[modelset2.Modelset.WorkSize];
  auto Wi2 = new Ti[modelset2.Modelset.WorkSizeI];

  modelset1.Modelset.Start(W1, Wi1);
  modelset2.Modelset.Start(W2, Wi2);

  auto summaries1 = std::vector<SearcherSummary *>();
  auto list11 = std::vector<SearcherSummary *>();
  auto list12 = std::vector<SearcherSummary *>();
  auto result1 = SearcherModelingInfo();
  modelset1.Modelset.CombineInfo(result1, summaries1, list11, list12);

  auto summaries2 = std::vector<SearcherSummary *>();
  auto result2 = SearcherModelingInfo();
  modelset2.Modelset.CombineInfo(result2, summaries2, list11, list12);

  auto bests1 = std::vector<EstimationKeep *>();
  modelset1.Modelset.CombineBests(0, 0, 0, summaries1, bests1);
  auto bests2 = std::vector<EstimationKeep *>();
  modelset2.Modelset.CombineBests(0, 0, 0, summaries2, bests2);

  ASSERT_NEAR(bests1.at(0)->Weight, bests2.at(0)->Weight, 1e-15);

  delete[] W1;
  delete[] Wi1;
  delete[] W2;
  delete[] Wi2;
}

TEST(DiscreteChoice_t, extended) {
  Matrix<Tv> source1 = Matrix<Tv>(
      new std::vector<Tv>{
          1,   0,   1,   0,   1,   1,   0,   0,   1,   1,   0,   1,   0,   0,
          1,   0,   1,   0,   0,   0,   0,   1,   0,   1,   0,   1,   1,   0,
          1,   1,   1,   0,   1,   1,   0,   1,   1,   1,   1,   1,   2,   3,
          4,   2,   3,   4,   5,   6,   7,   8,   9,   4,   5,   3,   4,   7,
          8,   4,   8,   9,   3,   5,   7,   8,   4,   5,   6,   7,   3,   4,
          6,   7,   3,   5,   6,   7,   8,   2,   4,   6,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
          1,   1,   1,   1,   1,   1,   1,   1,   9,   9,   6.6, 8.2, 3.5, 2.4,
          0.7, 0.1, 8.7, 2.9, 1.1, 9,   6.8, 1.6, 9.9, 0.2, 6.1, 9.9, 10,  2.1,
          8.1, 7.9, 2.4, 1.6, 7.2, 0.8, 9.4, 8.2, 9.4, 3.2, 4.9, 2.1, 8.7, 7.4,
          9.7, 0.4, 8.3, 0.2, 7.5, 0.1, 6.7, 9.4, 2.2, 5.7, 3.6, 1.1, 6.8, 6.2,
          9.5, 6.3, 5.4, 7.7, 8.1, 7.3, 2.8, 9.9, 3.8, 8.7, 0.3, 7.2, 8.4, 7.7,
          6.7, 1.8, 7.2, 2.8, 4.6, 2.4, 2.5, 7.5, 8,   2.3, 4.3, 1,   1.4, 5.2,
          8.8, 4.8, 5.2, 4.7},
      40, 5);
  Matrix<Tv> source2 = Matrix<Tv>(
      new std::vector<Tv>{
          1,   NAN, 0,   1,    0,   1,   1,   0,    0,   1,    1,   0,   1,
          0,   0,   1,   0,    1,   0,   0,   0,    0,   1,    0,   1,   0,
          1,   1,   0,   1,    1,   1,   0,   1,    1,   0,    1,   1,   1,
          1,   1,   2,   1000, 3,   4,   2,   3,    4,   5,    6,   7,   8,
          9,   4,   5,   3,    4,   7,   8,   4,    8,   9,    3,   5,   7,
          8,   4,   5,   6,    7,   3,   4,   6,    7,   3,    5,   6,   7,
          8,   2,   4,   6,    1,   1,   1,   1,    1,   1,    1,   1,   1,
          1,   1,   1,   1,    1,   1,   1,   1,    1,   1,    1,   1,   1,
          1,   1,   1,   1,    1,   1,   1,   1,    1,   1,    1,   1,   1,
          1,   1,   1,   1,    1,   1,   9,   1000, 9,   6.6,  8.2, 3.5, 2.4,
          0.7, 0.1, 8.7, 2.9,  1.1, 9,   6.8, 1.6,  9.9, 0.2,  6.1, 9.9, 10,
          2.1, 8.1, 7.9, 2.4,  1.6, 7.2, 0.8, 9.4,  8.2, 9.4,  3.2, 4.9, 2.1,
          8.7, 7.4, 9.7, 0.4,  8.3, 0.2, 7.5, 0.1,  6.7, 1000, 9.4, 2.2, 5.7,
          3.6, 1.1, 6.8, 6.2,  9.5, 6.3, 5.4, 7.7,  8.1, 7.3,  2.8, 9.9, 3.8,
          8.7, 0.3, 7.2, 8.4,  7.7, 6.7, 1.8, 7.2,  2.8, 4.6,  2.4, 2.5, 7.5,
          8,   2.3, 4.3, 1,    1.4, 5.2, 8.8, 4.8,  5.2, 4.7},
      41, 5);

  Matrix<Tv> xforecast = Matrix<Tv>(
      new std::vector<Tv>{1, 1, 1, 1, 1, 1, 2, 4, 5, 3, 6, 7, 6, 8, 3, 4, 9, 0},
      6, 3);

  auto model1 = DiscreteChoiceExtended(DiscreteChoiceModelType::kBinary,
                                       DiscreteChoiceDistType::kLogit, 40, 5,
                                       true, false, 2, true, 6, nullptr);
  auto model2 = DiscreteChoiceExtended(DiscreteChoiceModelType::kBinary,
                                       DiscreteChoiceDistType::kLogit, 41, 5,
                                       true, true, 2, true, 6, nullptr);

  auto W1 = new Tv[model1.WorkSize];
  auto S1 = new Tv[model1.StorageSize];
  auto W2 = new Tv[model2.WorkSize];
  auto S2 = new Tv[model2.StorageSize];

  model1.Calculate(source1, S1, W1, true, &xforecast);
  model2.Calculate(source2, S2, W2, true, &xforecast);

  ASSERT_EQ(true, model1.Model->Beta.Equals(model2.Model->Beta, 1e-15));

  delete[] W1;
  delete[] S1;
  delete[] W2;
  delete[] S2;
}

TEST(DiscreteChoice_t, extended_pca) {
  Matrix<Tv> source1 = Matrix<Tv>(new std::vector<Tv>{1,
                                                      0,
                                                      1,
                                                      0,
                                                      1,
                                                      1,
                                                      0,
                                                      0,
                                                      1,
                                                      1,
                                                      0,
                                                      1,
                                                      0,
                                                      0,
                                                      1,
                                                      0,
                                                      1,
                                                      0,
                                                      0,
                                                      0,
                                                      0,
                                                      1,
                                                      0,
                                                      1,
                                                      0,
                                                      1,
                                                      1,
                                                      0,
                                                      1,
                                                      1,
                                                      1,
                                                      0,
                                                      1,
                                                      1,
                                                      0,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      2,
                                                      3,
                                                      4,
                                                      2,
                                                      3,
                                                      4,
                                                      5,
                                                      6,
                                                      7,
                                                      8,
                                                      9,
                                                      4,
                                                      5,
                                                      3,
                                                      4,
                                                      7,
                                                      8,
                                                      4,
                                                      8,
                                                      9,
                                                      3,
                                                      5,
                                                      7,
                                                      8,
                                                      4,
                                                      5,
                                                      6,
                                                      7,
                                                      3,
                                                      4,
                                                      6,
                                                      7,
                                                      3,
                                                      5,
                                                      6,
                                                      7,
                                                      8,
                                                      2,
                                                      4,
                                                      6,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      -0.515747485750720,
                                                      -1.518528612793070,
                                                      1.177399496576260,
                                                      -0.137065367019942,
                                                      0.685655504751216,
                                                      1.624168887607020,
                                                      -0.477339876425763,
                                                      -0.249038349991145,
                                                      -1.552938016544990,
                                                      -0.311664346310363,
                                                      0.038979857312718,
                                                      -0.887147903173811,
                                                      -1.015683391546910,
                                                      -0.671231999108461,
                                                      0.924522228228185,
                                                      -1.624130107120040,
                                                      0.587709892016612,
                                                      -1.266740234568050,
                                                      1.852113059122460,
                                                      -0.638643020683457,
                                                      -1.138936281398830,
                                                      -0.877135563875742,
                                                      -0.455673449962293,
                                                      1.371470296718540,
                                                      -0.685063866519969,
                                                      1.007351580603140,
                                                      0.260552540183927,
                                                      1.088556010476260,
                                                      1.040493416772420,
                                                      -0.760075485208456,
                                                      -0.961249309198837,
                                                      1.181219024689690,
                                                      0.378344154055079,
                                                      1.615798296176280,
                                                      1.446303237947440,
                                                      0.119631429441564,
                                                      -1.289316873694990,
                                                      0.270012021737722,
                                                      0.055006330335831,
                                                      0.308062276143493,
                                                      -0.994461326881295,
                                                      -0.976109240381340,
                                                      -0.359418824645812,
                                                      -0.779381958265770,
                                                      0.509868268104989,
                                                      0.797955697243019,
                                                      1.308186421055500,
                                                      1.470515507878250,
                                                      -0.892225869340340,
                                                      0.694627682872129,
                                                      1.187732305507030,
                                                      -0.987664257807236,
                                                      -0.374785226531268,
                                                      1.061973963191750,
                                                      -1.270580888670890,
                                                      1.467930108741070,
                                                      -0.209870740571333,
                                                      -1.230478181133960,
                                                      -1.315308116067240,
                                                      0.922621482728360,
                                                      -0.733295317054620,
                                                      -0.682584155984064,
                                                      0.836019284057737,
                                                      1.024590083284440,
                                                      -0.491840807542710,
                                                      1.253263590048070,
                                                      -1.119673390781600,
                                                      -0.801812286210160,
                                                      -1.133947235837120,
                                                      0.619580501627405,
                                                      0.151491606074077,
                                                      0.889315844265479,
                                                      -0.927570628525437,
                                                      -0.589451745224259,
                                                      -1.224627675952180,
                                                      1.380514774670600,
                                                      -0.786045598847391,
                                                      1.433265056463370,
                                                      -0.588638609824417,
                                                      1.460319904267170},
                                  40, 5);
  Matrix<Tv> source2 = Matrix<Tv>(new std::vector<Tv>{1,
                                                      NAN,
                                                      0,
                                                      1,
                                                      0,
                                                      1,
                                                      1,
                                                      0,
                                                      0,
                                                      1,
                                                      1,
                                                      0,
                                                      1,
                                                      0,
                                                      0,
                                                      1,
                                                      0,
                                                      1,
                                                      0,
                                                      0,
                                                      0,
                                                      0,
                                                      1,
                                                      0,
                                                      1,
                                                      0,
                                                      1,
                                                      1,
                                                      0,
                                                      1,
                                                      1,
                                                      1,
                                                      0,
                                                      1,
                                                      1,
                                                      0,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      2,
                                                      1000,
                                                      3,
                                                      4,
                                                      2,
                                                      3,
                                                      4,
                                                      5,
                                                      6,
                                                      7,
                                                      8,
                                                      9,
                                                      4,
                                                      5,
                                                      3,
                                                      4,
                                                      7,
                                                      8,
                                                      4,
                                                      8,
                                                      9,
                                                      3,
                                                      5,
                                                      7,
                                                      8,
                                                      4,
                                                      5,
                                                      6,
                                                      7,
                                                      3,
                                                      4,
                                                      6,
                                                      7,
                                                      3,
                                                      5,
                                                      6,
                                                      7,
                                                      8,
                                                      2,
                                                      4,
                                                      6,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      1,
                                                      -0.515747485750720,
                                                      1000,
                                                      -1.518528612793070,
                                                      1.177399496576260,
                                                      -0.137065367019942,
                                                      0.685655504751216,
                                                      1.624168887607020,
                                                      -0.477339876425763,
                                                      -0.249038349991145,
                                                      -1.552938016544990,
                                                      -0.311664346310363,
                                                      0.038979857312718,
                                                      -0.887147903173811,
                                                      -1.015683391546910,
                                                      -0.671231999108461,
                                                      0.924522228228185,
                                                      -1.624130107120040,
                                                      0.587709892016612,
                                                      -1.266740234568050,
                                                      1.852113059122460,
                                                      -0.638643020683457,
                                                      -1.138936281398830,
                                                      -0.877135563875742,
                                                      -0.455673449962293,
                                                      1.371470296718540,
                                                      -0.685063866519969,
                                                      1.007351580603140,
                                                      0.260552540183927,
                                                      1.088556010476260,
                                                      1.040493416772420,
                                                      -0.760075485208456,
                                                      -0.961249309198837,
                                                      1.181219024689690,
                                                      0.378344154055079,
                                                      1.615798296176280,
                                                      1.446303237947440,
                                                      0.119631429441564,
                                                      -1.289316873694990,
                                                      0.270012021737722,
                                                      0.055006330335831,
                                                      0.308062276143493,
                                                      -0.994461326881295,
                                                      1000,
                                                      -0.976109240381340,
                                                      -0.359418824645812,
                                                      -0.779381958265770,
                                                      0.509868268104989,
                                                      0.797955697243019,
                                                      1.308186421055500,
                                                      1.470515507878250,
                                                      -0.892225869340340,
                                                      0.694627682872129,
                                                      1.187732305507030,
                                                      -0.987664257807236,
                                                      -0.374785226531268,
                                                      1.061973963191750,
                                                      -1.270580888670890,
                                                      1.467930108741070,
                                                      -0.209870740571333,
                                                      -1.230478181133960,
                                                      -1.315308116067240,
                                                      0.922621482728360,
                                                      -0.733295317054620,
                                                      -0.682584155984064,
                                                      0.836019284057737,
                                                      1.024590083284440,
                                                      -0.491840807542710,
                                                      1.253263590048070,
                                                      -1.119673390781600,
                                                      -0.801812286210160,
                                                      -1.133947235837120,
                                                      0.619580501627405,
                                                      0.151491606074077,
                                                      0.889315844265479,
                                                      -0.927570628525437,
                                                      -0.589451745224259,
                                                      -1.224627675952180,
                                                      1.380514774670600,
                                                      -0.786045598847391,
                                                      1.433265056463370,
                                                      -0.588638609824417,
                                                      1.460319904267170},
                                  41, 5);
  Matrix<Tv> source3 = Matrix<Tv>(
      new std::vector<Tv>{
          1,
          0,
          1,
          0,
          1,
          1,
          0,
          0,
          1,
          1,
          0,
          1,
          0,
          0,
          1,
          0,
          1,
          0,
          0,
          0,
          0,
          1,
          0,
          1,
          0,
          1,
          1,
          0,
          1,
          1,
          1,
          0,
          1,
          1,
          0,
          1,
          1,
          1,
          1,
          1,
          2,
          3,
          4,
          2,
          3,
          4,
          5,
          6,
          7,
          8,
          9,
          4,
          5,
          3,
          4,
          7,
          8,
          4,
          8,
          9,
          3,
          5,
          7,
          8,
          4,
          5,
          6,
          7,
          3,
          4,
          6,
          7,
          3,
          5,
          6,
          7,
          8,
          2,
          4,
          6,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          -0.515747485750720,
          -1.518528612793070,
          1.177399496576260,
          -0.137065367019942,
          0.685655504751216,
          1.624168887607020,
          -0.477339876425763,
          -0.249038349991145,
          -1.552938016544990,
          -0.311664346310363,
          0.038979857312718,
          -0.887147903173811,
          -1.015683391546910,
          -0.671231999108461,
          0.924522228228185,
          -1.624130107120040,
          0.587709892016612,
          -1.266740234568050,
          1.852113059122460,
          -0.638643020683457,
          -1.138936281398830,
          -0.877135563875742,
          -0.455673449962293,
          1.371470296718540,
          -0.685063866519969,
          1.007351580603140,
          0.260552540183927,
          1.088556010476260,
          1.040493416772420,
          -0.760075485208456,
          -0.961249309198837,
          1.181219024689690,
          0.378344154055079,
          1.615798296176280,
          1.446303237947440,
          0.119631429441564,
          -1.289316873694990,
          0.270012021737722,
          0.055006330335831,
          0.308062276143493,
      },
      40, 4);

  Matrix<Tv> xforecast = Matrix<Tv>(new std::vector<Tv>{1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        1,
                                                        -0.515747485750720,
                                                        -1.518528612793070,
                                                        1.177399496576260,
                                                        -0.137065367019942,
                                                        0.685655504751216,
                                                        1.624168887607020,
                                                        -0.477339876425763,
                                                        -0.249038349991145,
                                                        -1.552938016544990,
                                                        -0.311664346310363,
                                                        0.038979857312718,
                                                        -0.887147903173811,
                                                        -1.015683391546910,
                                                        -0.671231999108461,
                                                        0.924522228228185,
                                                        -1.624130107120040,
                                                        0.587709892016612,
                                                        -1.266740234568050,
                                                        1.852113059122460,
                                                        -0.638643020683457,
                                                        -1.138936281398830,
                                                        -0.877135563875742,
                                                        -0.455673449962293,
                                                        1.371470296718540,
                                                        -0.685063866519969,
                                                        1.007351580603140,
                                                        0.260552540183927,
                                                        1.088556010476260,
                                                        1.040493416772420,
                                                        -0.760075485208456,
                                                        -0.961249309198837,
                                                        1.181219024689690,
                                                        0.378344154055079,
                                                        1.615798296176280,
                                                        1.446303237947440,
                                                        0.119631429441564,
                                                        -1.289316873694990,
                                                        0.270012021737722,
                                                        0.055006330335831,
                                                        0.308062276143493,
                                                        -0.994461326881295,
                                                        -0.976109240381340,
                                                        -0.359418824645812,
                                                        -0.779381958265770,
                                                        0.509868268104989,
                                                        0.797955697243019,
                                                        1.308186421055500,
                                                        1.470515507878250,
                                                        -0.892225869340340,
                                                        0.694627682872129,
                                                        1.187732305507030,
                                                        -0.987664257807236,
                                                        -0.374785226531268,
                                                        1.061973963191750,
                                                        -1.270580888670890,
                                                        1.467930108741070,
                                                        -0.209870740571333,
                                                        -1.230478181133960,
                                                        -1.315308116067240,
                                                        0.922621482728360,
                                                        -0.733295317054620,
                                                        -0.682584155984064,
                                                        0.836019284057737,
                                                        1.024590083284440,
                                                        -0.491840807542710,
                                                        1.253263590048070,
                                                        -1.119673390781600,
                                                        -0.801812286210160,
                                                        -1.133947235837120,
                                                        0.619580501627405,
                                                        0.151491606074077,
                                                        0.889315844265479,
                                                        -0.927570628525437,
                                                        -0.589451745224259,
                                                        -1.224627675952180,
                                                        1.380514774670600,
                                                        -0.786045598847391,
                                                        1.433265056463370,
                                                        -0.588638609824417,
                                                        1.460319904267170},
                                    40, 3);
  Matrix<Tv> xforecast3 = Matrix<Tv>(
      new std::vector<Tv>{
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          -0.515747485750720,
          -1.518528612793070,
          1.177399496576260,
          -0.137065367019942,
          0.685655504751216,
          1.624168887607020,
          -0.477339876425763,
          -0.249038349991145,
          -1.552938016544990,
          -0.311664346310363,
          0.038979857312718,
          -0.887147903173811,
          -1.015683391546910,
          -0.671231999108461,
          0.924522228228185,
          -1.624130107120040,
          0.587709892016612,
          -1.266740234568050,
          1.852113059122460,
          -0.638643020683457,
          -1.138936281398830,
          -0.877135563875742,
          -0.455673449962293,
          1.371470296718540,
          -0.685063866519969,
          1.007351580603140,
          0.260552540183927,
          1.088556010476260,
          1.040493416772420,
          -0.760075485208456,
          -0.961249309198837,
          1.181219024689690,
          0.378344154055079,
          1.615798296176280,
          1.446303237947440,
          0.119631429441564,
          -1.289316873694990,
          0.270012021737722,
          0.055006330335831,
          0.308062276143493,
      },
      40, 2);

  auto pcaoptions = PcaAnalysisOptions();
  pcaoptions.IgnoreFirstCount = 1;
  pcaoptions.ExactCount = 1;

  auto model1 = DiscreteChoiceExtended(
      DiscreteChoiceModelType::kBinary, DiscreteChoiceDistType::kLogit, 40, 5,
      true, false, 2, true, xforecast.RowsCount, &pcaoptions);
  auto model2 = DiscreteChoiceExtended(
      DiscreteChoiceModelType::kBinary, DiscreteChoiceDistType::kLogit, 41, 5,
      true, true, 2, true, xforecast.RowsCount, &pcaoptions);
  auto model3 = DiscreteChoiceExtended(
      DiscreteChoiceModelType::kBinary, DiscreteChoiceDistType::kLogit, 40, 5,
      true, true, 2, true, xforecast.RowsCount,
      nullptr); // they are orthogonal, must be the same

  auto W1 = new Tv[model1.WorkSize];
  auto S1 = new Tv[model1.StorageSize];
  auto W2 = new Tv[model2.WorkSize];
  auto S2 = new Tv[model2.StorageSize];
  auto W3 = new Tv[model3.WorkSize];
  auto S3 = new Tv[model3.StorageSize];

  model1.Calculate(source1, S1, W1, true, &xforecast);
  model2.Calculate(source2, S2, W2, true, &xforecast);
  model3.Calculate(source3, S3, W3, true, &xforecast3);

  ASSERT_EQ(true, model1.Model->Beta.Equals(model2.Model->Beta, 1e-10));
  ASSERT_EQ(true, model1.Model->Beta.Equals(model3.Model->Beta, 1e-10));

  // forecasts are orthogonal too
  ASSERT_EQ(true, model1.PredProbs.Equals(model2.PredProbs, 1e-10));
  ASSERT_EQ(true, model1.PredProbs.Equals(model3.PredProbs, 1e-10));

  delete[] W1;
  delete[] S1;
  delete[] W2;
  delete[] S2;
}
