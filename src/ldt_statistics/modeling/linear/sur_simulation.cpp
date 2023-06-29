/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "sur.h"

using namespace ldt;

SurSimulation::SurSimulation(Ti N, Ti m, Ti k, Tv trainRatio, Ti trainFixSize,
                             const std::vector<ScoringType> &measuresOut,
                             bool isRestricted, Ti maxSigSearchIter,
                             PcaAnalysisOptions *pcaOptionsY,
                             PcaAnalysisOptions *pcaOptionsX) {
  pMeasuresOut = &measuresOut;

  mTrainPerc = trainRatio;
  mTrainFixSize = trainFixSize;

  Ti N0 =
      mTrainFixSize > 0
          ? mTrainFixSize
          : static_cast<Ti>(std::round(N * trainRatio)); // TODO: use data_split
  Ti N1 = N - N0;

  Split = DataSplit(N, m + k);
  Split_d = std::unique_ptr<Ti[]>(new Ti[Split.WorkSizeI]);

  auto doForcVariance = false;
  for (auto &s : measuresOut)
    if (Scoring::RequiresVariance(s)) {
      doForcVariance = true;
      break;
    }

  Model = SurExtended(
      N0, m, k, isRestricted, false, // NANs must be removed before
      false, N1, maxSigSearchIter, doForcVariance, pcaOptionsY, pcaOptionsX);

  WorkSize = Split.StorageSize + Model.StorageSize + Model.WorkSize;
  // I have generated Split.WorkSizeI here
  WorkSize += N1 * m; // for errors matrix

  StorageSize = m * (Ti)measuresOut.size();
}

void SurSimulation::AddError(std::string state) {
  if (KeepErrors) {
    if (state.empty())
      return;
    if (Errors.find(state) != Errors.end())
      Errors.at(state)++;
    else {
      Errors.insert(std::pair<std::string, Ti>(state, 1));
    }
  }
}

Tv sumScores(const ScoringType &e, const Ti &length, const Tv *actuals,
             const Tv *errors, const Tv *means, const Tv *variances) {
  Tv sum = 0;
  switch (e) {
  case ldt::ScoringType::kDirection:
    throw std::logic_error("not implemented (direction)");
  case ldt::ScoringType::kSign:
    for (Ti i = 0; i < length; i++)
      sum += actuals[i] * means[i] > 0 ? 1 : 0;
    break;
  case ldt::ScoringType::kMae:
    for (Ti i = 0; i < length; i++)
      sum += std::abs(errors[i]);
    break;
  case ldt::ScoringType::kMape:
    for (Ti i = 0; i < length; i++)
      sum += std::abs(errors[i] / actuals[i]);
    break;
  case ldt::ScoringType::kRmse:
    for (Ti i = 0; i < length; i++)
      sum += std::pow(errors[i], 2.0);
    break;
  case ldt::ScoringType::kRmspe:
    for (Ti i = 0; i < length; i++)
      sum += std::pow(errors[i] / actuals[i], 2.0);
    break;
  case ldt::ScoringType::kCrps:
    for (Ti i = 0; i < length; i++)
      sum += Scoring::GetScoreCrpsNormal(errors[i], 0, std::sqrt(variances[i]));
    break;
  default:
    throw std::logic_error("not implemented (averaging scores)");
  }

  return sum;
}

void SurSimulation::Calculate(Matrix<Tv> &data, Ti m, Tv *storage, Tv *work,
                              Matrix<Tv> *R, bool &cancel, Ti maxIteration,
                              unsigned int seed, Tv sigSearchMaxProb,
                              Tv maxCondNum, Ti maxInvalidSim) {
  if (cancel)
    return;
  if (maxIteration <= 0)
    throw std::logic_error("Number of iterations must be positive.");

  // Ti N = data.RowsCount;
  Ti k = data.ColsCount - m;
  // Ti km = k * m;

  Results = Matrix<Tv>(storage, (Ti)pMeasuresOut->size(), m);
  Results.SetValue(0.0);

  std::mt19937 eng;
  if (seed == 0) {
    std::random_device rdev{};
    eng = std::mt19937(rdev());
  } else
    eng = std::mt19937(seed);

  Ti p = 0;
  Split.Calculate(data, &work[p], mTrainPerc, mTrainFixSize);
  p += Split.StorageSize;
  auto model_storage = &work[p];
  p += Model.StorageSize;
  auto model_work = &work[p];
  p += Model.WorkSize;
  auto error_work = &work[p];
  p += Split.Sample1.RowsCount * m;

  auto newY = Matrix<Tv>(Split.Sample1.Data, Split.Sample1.RowsCount, m);
  auto newX = Matrix<Tv>(&Split.Sample1.Data[m * Split.Sample1.RowsCount],
                         Split.Sample1.RowsCount, k);
  auto errors = Matrix<Tv>(error_work, Split.Sample1.RowsCount, m);

  Ti i, j, s, invalidCounts = 0;
  ValidCounts = 0;
  ValidIters = 0;
  for (Iteration = 0; Iteration < maxIteration; Iteration++) {
    invalidCounts++;
    if (cancel)
      return;

    Split.Shuffle(data, Split_d.get(), eng);

    try {
      Model.Calculate(Split.Sample0, m, model_storage, model_work, R,
                      sigSearchMaxProb, &newX, nullptr);
      if (Model.Model.condition_number > maxCondNum)
        throw std::logic_error("Model check failed: Maximum CN");
    } catch (std::exception &ex) {
      AddError(ex.what());
      continue;
    } catch (...) {
      AddError("unknown error!");
      continue;
    }

    if (cancel)
      return;

    newY.Subtract(Model.Projections.Means, errors);

    if (errors.Any(NAN)) { // should I check variance too?
      AddError("NAN found");
      continue;
    }

    ValidCounts += newY.RowsCount;
    ValidIters++;
    invalidCounts--;
    if (invalidCounts > maxInvalidSim)
      throw std::logic_error("Model check failed: Minimum Valid Simulations");

    i = -1;
    for (auto &e : *pMeasuresOut) {
      i++;
      for (j = 0; j < m; j++) {
        if (cancel)
          return;
        s = j * errors.RowsCount;
        Results.Set_Plus0(i, j,
                          sumScores(e, errors.RowsCount, &newY.Data[s],
                                    &errors.Data[s],
                                    &Model.Projections.Means.Data[s],
                                    &Model.Projections.Variances.Data[s]));
      }
    }
  }

  if (cancel)
    return;

  if (invalidCounts > maxInvalidSim)
    throw std::logic_error("Model check failed: Minimum Valid Simulations");

  Results.Divide_in((Tv)ValidCounts);

  // sqrt for RMSE or RMSPE
  i = -1;
  for (auto &e : *pMeasuresOut) {
    i++;
    if (e == ScoringType::kRmse || e == ScoringType::kRmspe)
      for (j = 0; j < m; j++) {
        s = j * errors.RowsCount;
        Results.Set0(i, j, std::sqrt(Results.Get0(i, j)));
      }
  }
}
