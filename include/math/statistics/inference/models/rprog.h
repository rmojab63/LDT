#pragma once

#include <map>
#include <string>

#include "ldt_base.h"
#include "matrix.h"
#include "searchers.h"

namespace ldt {

/// @brief A searcher class for general R functions
class LDT_EXPORT RFuncSearcher : public SearcherReg {

  std::function<std::string(const std::vector<Ti> &columnIndices,
                            const int &numEndo, VMatrix<Tv> &metrics,
                            VMatrix<Tv> &type1Mean, VMatrix<Tv> &type1Var,
                            VMatrix<Ti> &extra)>
      mFunc;

  std::string EstimateOneReg(Tv *work, Ti *workI, VMatrix<Tv> &metrics,
                             VMatrix<Tv> &type1Mean, VMatrix<Tv> &type1Var,
                             VMatrix<Ti> &extra) override;

public:
  RFuncSearcher(
      const SearchData &data, const SearchCombinations &combinations,
      SearchOptions &options, const SearchItems &items,
      const SearchMetricOptions &metrics, const SearchModelChecks &checks,
      const Ti &numPartitions, const std::vector<Ti> &innerIndices,
      const bool &isInnerExogenous,
      std::function<std::string(const std::vector<Ti> &columnIndices,
                                const int &numEndo, VMatrix<Tv> &metrics,
                                VMatrix<Tv> &type1Mean, VMatrix<Tv> &type1Var,
                                VMatrix<Ti> &extra)>
          func);
};

/// @brief A model set with no estimation (for R)
class LDT_EXPORT RFuncModelset {
public:
  /// @brief The inner model set
  ModelSet Modelset;

  /// @brief List of searchers
  std::vector<Searcher *> Searchers;

  /// @brief Initializes a new instance of this class
  RFuncModelset(){};

  RFuncModelset(
      const SearchData &data, const SearchCombinations &combinations,
      SearchOptions &options, SearchItems &items, SearchMetricOptions &metrics,
      SearchModelChecks &checks, bool isTimeSeries, bool isOutOfSampleRandom,
      const bool &isInnerExogenous,
      std::function<std::string(const std::vector<Ti> &columnIndices,
                                const int &numEndo, VMatrix<Tv> &metrics,
                                VMatrix<Tv> &type1Mean, VMatrix<Tv> &type1Var,
                                VMatrix<Ti> &extra)>
          func);

  ~RFuncModelset() {
    for (auto s : Searchers)
      delete s;
  }
};

} // namespace ldt
