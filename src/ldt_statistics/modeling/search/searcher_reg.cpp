/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "searchers.h"

using namespace ldt;

SearcherReg::SearcherReg(const SearchData &data,
                         const SearchCombinations &combinations,
                         SearchOptions &options, const SearchItems &items,
                         const SearchMetricOptions &metrics,
                         const SearchModelChecks &checks,
                         const Ti &numPartitions, const bool &isInnerExogenous,
                         const std::vector<Ti> &zInnerIndices,
                         const Ti extraLength, bool checkForEmpty)
    : Searcher::Searcher(data, combinations, options, items, metrics, checks,
                         numPartitions, checkForEmpty) {
  IsInnerExogenous = isInnerExogenous;
  mExtraLength = extraLength;
  auto w = data.HasWeight ? 1 : 0;
  ColIndices = std::vector<Ti>(
      numPartitions + w +
      zInnerIndices.size()); // it gets updated in EstimateOne method

  if (data.HasWeight)
    ColIndices.at(numPartitions) = data.NumEndo;

  if (isInnerExogenous) {
    // zInnerIndices is zero-based. Move InnerIndices
    // Also note that this set is unchanged in the main loop
    for (auto &i : zInnerIndices)
      InnerIndices.push_back(i + numPartitions + w);

    // the second part of ColIndices is constant:
    for (Ti i = 0; i < zInnerIndices.size(); i++)
      ColIndices.at(i + numPartitions + w) = zInnerIndices.at(i);

    // targets are in current indices and should be updated in the main loop
  } else {
    // inner is endogenous, we should move outer which is given by
    // CurrentIndices. Note that this set changes in the main loop.
    this->IndicesOffset = zInnerIndices.size();
    InnerIndices = zInnerIndices;

    // put it in the first part of column indices
    for (Ti i = 0; i < zInnerIndices.size(); i++)
      ColIndices.at(i) = zInnerIndices.at(i);

    // we can set targets here
    for (auto &a : zInnerIndices) {
      if (a < items.LengthTargets)
        TargetsPositions.push_back(a);
    }
    if (TargetsPositions.size() == 0)
      throw LdtException(ErrorType::kLogic, "sur-modelset",
                         "a searcher with no target is not valid");
  }
}

std::string SearcherReg::EstimateOne(Tv *work, Ti *workI) {
  // update column and inner indices based on current indices

  if (IsInnerExogenous) {
    // update the first part of column indices
    for (Ti i = 0; i < CurrentIndices.Vec.size(); i++)
      ColIndices.at(i) = CurrentIndices.Vec.at(i);

    // update target indices
    for (auto &a : CurrentIndices.Vec) {
      if (a < pItems->LengthTargets)
        TargetsPositions.push_back(a);
    }
  } else {
    auto w = pData->HasWeight ? 1 : 0;
    // update the second part of column indices
    for (Ti i = 0; i < CurrentIndices.Vec.size(); i++)
      ColIndices.at(i + NumPartitions + w) = CurrentIndices.Vec.at(i);
  }

  auto numMeas = (Ti)(this->pMetrics->MetricsIn.size() +
                      this->pMetrics->MetricsOut.size());

  auto metrics = VMatrix<Tv>(numMeas, (Ti)TargetsPositions.size());

  auto extra = VMatrix<Ti>(mExtraLength);

  auto type1_mean =
      VMatrix<Tv>(this->pItems->Length1, (Ti)TargetsPositions.size());
  auto type1_var =
      VMatrix<Tv>(this->pItems->Length1, (Ti)TargetsPositions.size());

  auto res = EstimateOneReg(work, workI, metrics, type1_mean, type1_var, extra);
  if (res.length() != 0)
    return res;

  // push:
  Ti i, j;
  bool allNan = true;
  for (i = 0; i < metrics.Mat.RowsCount; i++) { // metric index
    for (j = 0; j < metrics.Mat.ColsCount; j++) {

      // j is the target index (use TargetsPositions to interpret)
      Ti t_pos = TargetsPositions.at(j);

      auto metric = metrics.Mat.Get0(i, j);
      if (std::isnan(metric))
        continue;
      allNan = false;

      Tv weight = NAN;
      if (i < this->pMetrics->MetricsIn.size())
        weight = GoodnessOfFit::ToWeight(
            pMetrics->MetricsIn.at(i), metric,
            this->pMetrics->MinMetrics.Mat.Get(i, t_pos));
      else
        weight = Scoring::ToWeight(
            pMetrics->MetricsOut.at(i - this->pMetrics->MetricsIn.size()),
            metric, this->pMetrics->MinMetrics.Mat.Get(i, t_pos));

      // add model information:
      if (this->pItems->KeepModelEvaluations) {
        EstimationKeep *ek = nullptr;
        if (IsInnerExogenous)
          ek = new EstimationKeep(metric, weight, InnerIndices, extra.Vec,
                                  this->CurrentIndices.Vec);
        else
          ek = new EstimationKeep(metric, weight, this->CurrentIndices.Vec,
                                  extra.Vec, InnerIndices);
        this->Push0(*ek, i, t_pos);
      }

      if (this->pItems->Length1 > 0) { // Add coefficients
        for (Ti i1 = 0; i1 < type1_mean.Mat.RowsCount; i1++) {

          EstimationKeep *ek = nullptr;

          if (IsInnerExogenous)
            ek = new EstimationKeep(metric, weight, InnerIndices, extra.Vec,
                                    this->CurrentIndices.Vec,
                                    type1_mean.Mat.Get(i1, j),
                                    type1_var.Mat.Get(i1, j));
          else
            ek = new EstimationKeep(metric, weight, this->CurrentIndices.Vec,
                                    extra.Vec, InnerIndices,
                                    type1_mean.Mat.Get(i1, j),
                                    type1_var.Mat.Get(i1, j));

          this->Push1(*ek, i, t_pos, i1);
        }
      }

      if (this->pItems->Length2 > 0) { // TODO: for saving projections
      }
    }
  }

  if (allNan) {
    throw LdtException(
        ErrorType::kLogic, "searcher-reg",
        "all metrics are NaN"); // +
                                // Sur::ModelToString(DModel.Sizes).c_str();
                                // // weights are all nan
  }

  return "";
}