#include "r_ldt.h"
#include "rprog.h"

using namespace Rcpp;
using namespace ldt;

// [[Rcpp::export(.SearchRFunc)]]
SEXP SearchRFunc(List data, List combinations, List metrics, List modelChecks,
                 List items, List options, Function funcR, int length1,
                 bool isInnerExogenous) {

  bool isTimeSeries = false;
  bool isOutOfSampleRamdom = false;

  auto options_ = SearchOptions();
  UpdateSearchOptions(options, options_);

  auto data_ = SearchData();
  UpdateSearchData(data, data_);

  auto exoStart = data_.NumEndo;
  auto colNames = as<std::vector<std::string>>(colnames(data["data"]));
  auto exoNames =
      std::vector<std::string>(colNames.begin() + exoStart, colNames.end());

  auto combinations_ = SearchCombinations();
  UpdateSearchCombinations(combinations, combinations_);

  auto metrics_ = SearchMetricOptions();
  auto metricsNames = std::vector<std::string>();
  auto items_ = SearchItems();
  auto checks_ = SearchModelChecks();

  UpdateOptions(items, metrics, modelChecks, metrics_, items_, checks_,
                metricsNames, length1, data_.NumExo, combinations_.NumTargets,
                data_.NumEndo, false, true, "Coefficients", false);

  auto targetNames = std::vector<std::string>(
      colNames.begin(), colNames.begin() + items_.LengthTargets);

  std::function<std::string(const std::vector<Ti> &columnIndices,
                            const int &numEndo, VMatrix<Tv> &metrics,
                            VMatrix<Tv> &type1Mean, VMatrix<Tv> &type1Var,
                            VMatrix<Ti> &extra)>
      func = [&](const std::vector<Ti> &columnIndices, const int &numEndo,
                 VMatrix<Tv> &metrics, VMatrix<Tv> &type1Mean,
                 VMatrix<Tv> &type1Var, VMatrix<Ti> &extra) -> std::string {
    SEXP columnIndicesR = PROTECT(Rf_allocVector(INTSXP, columnIndices.size()));
    SEXP numEndoR = PROTECT(Rf_ScalarInteger(numEndo));

    try {
      for (int i = 0; i < (Ti)columnIndices.size(); i++)
        INTEGER(columnIndicesR)[i] = columnIndices[i];

      if (options_.RequestCancel)
        throw std::logic_error("Cancel");

      Rcpp::Environment G = Rcpp::Environment::global_env();
      Rcpp::Function F = G["func"];
      F(columnIndicesR, numEndoR);

      // List result =
      //     funcR(_["columnIndices"] = columnIndicesR, _["numEndo"] = numEndoR
      //,
      //_["metrics"] = metrics, _["type1Mean"] = type1Mean,
      //_["type1Var"] = type1Var, _["extra"] = extra
      //   );

      // TODO:

    } catch (...) {
      UNPROTECT(2);
      throw;
    }

    UNPROTECT(2);
    return "TODO";
  };

  // Modelset
  auto model =
      RFuncModelset(data_, combinations_, options_, items_, metrics_, checks_,
                    isTimeSeries, isOutOfSampleRamdom, isInnerExogenous, func);
  model.Modelset.ShuffleSearchers = false;

  bool estimating = true;

  std::unique_ptr<double[]> W;
  try {
    W = std::unique_ptr<double[]>(new double[model.Modelset.WorkSize]);
  } catch (...) {
    throw LdtException(ErrorType::kLogic, "R-sur",
                       "more memory is required for running the project");
  }

  auto alli = model.Modelset.GetExpectedNumberOfModels();

  auto f = std::async(std::launch::async, [&] {
    model.Modelset.Start(W.get(), nullptr);
    estimating = false;
  });

  ReportProgress(model.Modelset, estimating, options_, alli);

  if (options_.RequestCancel)
    return R_NilValue;

  auto extraNames = std::vector<std::string>({"extra"});
  List L = GetModelSetResults(
      model.Modelset, items_, metricsNames, colNames, targetNames, extraNames,
      exoNames, colNames, std::string("coefs"), options_.ReportInterval > 0);

  return L;
}
