/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "rprog.h"

using namespace ldt;

// #pragma region Searcher

RFuncSearcher::RFuncSearcher(
    const SearchData &data, const SearchCombinations &combinations,
    SearchOptions &options, const SearchItems &items,
    const SearchMetricOptions &metrics, const SearchModelChecks &checks,
    const Ti &numPartitions, const std::vector<Ti> &innerIndices,
    const bool &isInnerExogenous,
    std::function<std::string(const std::vector<Ti> &columnIndices,
                              const int &numEndo, VMatrix<Tv> &metrics,
                              VMatrix<Tv> &type1Mean, VMatrix<Tv> &type1Var,
                              VMatrix<Ti> &extra)>
        func)
    : SearcherReg::SearcherReg(data, combinations, options, items, metrics,
                               checks, numPartitions, isInnerExogenous,
                               innerIndices, 0, isInnerExogenous) {

  mFunc = func;
}

std::string RFuncSearcher::EstimateOneReg(Tv *work, Ti *workI,
                                          VMatrix<Tv> &metrics,
                                          VMatrix<Tv> &type1Mean,
                                          VMatrix<Tv> &type1Var,
                                          VMatrix<Ti> &extra) {
  Ti num_endo;
  if (this->IsInnerExogenous)
    num_endo = this->CurrentIndices.Vec.size();
  else
    num_endo = this->InnerIndices.size();
  return mFunc(this->ColIndices, num_endo, metrics, type1Mean, type1Var, extra);
}

// #pragma endregion

// #pragma region Modelset

RFuncModelset::RFuncModelset(
    const SearchData &data, const SearchCombinations &combinations,
    SearchOptions &options, SearchItems &items, SearchMetricOptions &metrics,
    SearchModelChecks &checks, bool isTimeSeries, bool isOutOfSampleRandom,
    const bool &isInnerExogenous,
    std::function<std::string(const std::vector<Ti> &columnIndices,
                              const int &numEndo, VMatrix<Tv> &metrics,
                              VMatrix<Tv> &type1Mean, VMatrix<Tv> &type1Var,
                              VMatrix<Ti> &extra)>
        func) {

  metrics.Update(isOutOfSampleRandom, isTimeSeries);
  checks.Update(metrics);
  items.Update(metrics, items.LengthTargets, items.LengthDependents,
               items.LengthExogenouses);

  unsigned int co = 0;
  for (auto const &size : combinations.Sizes) {
    if (size <= 0)
      throw LdtException(
          ErrorType::kLogic, "sur-modelset",
          "invalid exogenous size (zero or negative). Make sure array is "
          "initialized properly");

    if (combinations.NumFixPartitions > size)
      continue;

    for (auto &inner_group : combinations.InnerGroups) {
      if (inner_group.size() == 0)
        throw LdtException(ErrorType::kLogic, "sur-modelset",
                           "empty endogenous indexes");

      if (isInnerExogenous == false) {
        if (inner_group.at(0) > items.LengthTargets)
          continue; // no target here
      }
      /*auto seed = metrics.Seed == 0
                      ? 0
                      : (unsigned int)(metrics.Seed < 0 ? -metrics.Seed
                                                        : (metrics.Seed +
         co));*/
      co++;
      auto se =
          new RFuncSearcher(data, combinations, options, items, metrics, checks,
                            size, inner_group, isInnerExogenous, func);
      Searchers.push_back(se);
    }
  }
  this->Modelset =
      ModelSet(Searchers, data, combinations, options, items, metrics, checks);
}

// #pragma endregion
