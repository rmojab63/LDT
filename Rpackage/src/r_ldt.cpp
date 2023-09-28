#include "r_ldt.h"

using namespace Rcpp;
using namespace ldt;

// [[Rcpp::export(.SupportsParallel)]]
bool SupportsParallel() {
#ifndef _OPENMP
  return false;
#else
  return true;
#endif
}

std::vector<std::vector<int>> listToVectorOfVectors(List list) {
  std::vector<std::vector<int>> result;
  for(int i = 0; i < list.size(); i++) {
    std::vector<int> vec = as<std::vector<int>>(list[i]);
    result.push_back(vec);
  }
  return result;
}

void UpdateSearchData(List& dataR, SearchData& data) {

  auto mat = as<NumericMatrix>(list["data"]);
  data.Data.SetData(&mat[0], mat.nrow(), mat.ncol());

  data.NumEndo = as<int>(list["numEndo"]);
  data.NumExo = as<int>(list["numExo"]);
  data.ObsCount = as<int>(list["obsCount"]);
  data.NewObsCount = as<int>(list["newObsCount"]);
  if (data.NewObsCount > 0){
    auto mat1 = as<NumericMatrix>(list["newX"]);
    data.newX.SetData(&mat1[0], mat1.nrow(), mat1.ncol());
  }
  data.Lambdas = as<std::vector<int>>(list["lambdas"]);
  data.HasIntercept = as<bool>(list["hasIntercept"]);
  data.HasWeight = as<bool>(list["hasWeight"]);

}

void UpdateSearchCombinations(List combinationsR, SearchCombinations& combinations){

  combinations.Sizes = as<std::vector<int>>(list["sizes"]);
  combinations.Partitions = listToVectorOfVectors(list["partitions"]);
  combinations.NumFixPartitions = as<int>(list["numFixPartitions"]);
  combinations.InnerGroups = listToVectorOfVectors(list["innerGroups"]);
  combinations.NumTargets = as<int>(list["numTargets"]);

}

void UpdateSearchOptions(List &optionsR, SearchOptions &options) {

  options.Parallel = as<bool>(optionsR["parallel"]);
  options.ReportInterval = as<int>(optionsR["reportInterval"]);
}

void UpdateModelCheckItems(List &checksR,
                           SearchModelChecks &checks,
                           const SearchMetricOptions &metrics,
                           const SearchItems &items) {

  checks.Estimation = as<bool>(checksR["estimation"]);
  checks.MinObsCount = as<int>(checksR["minObsCount"]);
  checks.MinDof = as<int>(checksR["minDof"]);
  checks.MinOutSim = as<int>(checksR["minOutSim"]);
  checks.PredictionBoundMultiplier =
    as<double>(checksR["predictionBoundMultiplier"]);

  checks.MinR2 = as<double>(checksR["minR2"]);
  checks.MaxAic = as<double>(checksR["maxAic"]);
  checks.MaxSic = as<double>(checksR["maxSic"]);
  checks.MaxConditionNumber = as<double>(checksR["maxConditionNumber"]);
  checks.Prediction = as<bool>(checksR["prediction"]);

  checks.Update(metrics);
}


void UpdateSearchItems(List &itemsR, SearchItems &items,
                       int length1, int length2, const char *length1Informtion,
                       const char *length2Informtion, bool type1NeedsModelEstim,
                       bool type2NeedsModelEstim) {

  items.KeepModelEvaluations = as<bool>(itemsR["model"]);
  items.KeepAll = as<bool>(itemsR["all"]);
  items.KeepMixture = as<bool>(itemsR["mixture4"]);
  items.KeepInclusionWeights = as<bool>(itemsR["inclusion"]);
  items.KeepBestCount = as<int>(itemsR["bestK"]);
  items.ExtremeBoundsMultiplier = as<double>(itemsR["extremeMultiplier"]);

  items.CdfsAt = as<std::vector<double>>(itemsR["cdfs"]);

  // update length1 and 2
  bool type1 = as<bool>(itemsR["type1"]);
  bool type2 = as<bool>(itemsR["type2"]);
  items.Length1 = type1 ? length1 : 0;
  items.Length2 = type2 ? length2 : 0;

  if (type1NeedsModelEstim && items.Length1 > 0)
    items.KeepModelEvaluations = true;
  if (type2NeedsModelEstim && items.Length2 > 0)
    items.KeepModelEvaluations = true;

  if (items.KeepInclusionWeights)
    items.KeepModelEvaluations = true;
  if (items.KeepModelEvaluations == false && items.Length1 == 0 &&
      items.Length2 == 0)
    throw LdtException(ErrorType::kLogic, "R-ldt",
                       "no evaluation data is saved");
  if (hasGoal == false)
    throw LdtException(ErrorType::kLogic, "R-ldt", "no goal is set");
}

static std::vector<double> helper_convertMinMetrics(NumericVector lst, int k) {

  std::vector<double> vec(k);
  if (lst.size() == 1) {
    std::fill(vec.begin(), vec.end(), lst[0]);
  } else if (lst.size() == k) {
    for (int i = 0; i < k; i++) {
      vec[i] = lst[i];
    }
  } else {
    throw std::logic_error("Invalid length for a member of 'minMetrics'.");
  }
  return vec;
}

void UpdatemetricOptions(List &metricsR,
                         SearchMetricOptions &metrics,
                         std::vector<std::string> &metricNames,
                         bool isTimeSeries, bool isDc, int numTargets) {

  bool isOutOfSampleRandom = isTimeSeries == false;
  // bool supportsSimRatio = isTimeSeries;

  auto metricsOut0 = as<StringVector>(metricsR["typesOut"]);
  auto metricsIn0 = as<StringVector>(metricsR["typesIn"]);
  auto lmetricOut = metricsOut0.length();
  auto lmetricIn = metricsIn0.length();

  if (lmetricIn == 0 && lmetricOut == 0)
    throw LdtException(
        ErrorType::kLogic, "R-ldt",
        "No metric is specified. Check the inputs (also, check the number "
        "of simulations)");
  if (lmetricIn > 0) {
    for (auto i = 0; i < lmetricIn; i++) {
      auto a = as<std::string>(metricsIn0[i]);
      boost::algorithm::to_lower(a);
      auto eval = FromString_GoodnessOfFitType(a.c_str());
      metrics.MetricsIn.push_back(eval);
    }
  }
  if (lmetricOut > 0) {
    for (auto i = 0; i < lmetricOut; i++) {
      auto a = as<std::string>(metricsOut0[i]);
      boost::algorithm::to_lower(a);
      auto eval = FromString_ScoringType(a.c_str());
      metrics.MetricsOut.push_back(eval);
    }
  }

  metrics.SimFixSize = as<int>(metricsR["simFixSize"]);
  // metrics.SimRatio = Rf_asInteger(GetListElement(metrics,
  // "simratio"));
  metrics.Seed = as<int>(metricsR["seed"]);

  if (isTimeSeries && lmetricOut > 0) {

    IntegerVector hors = metricsR["horizons"];
    for (auto i = 0; i < hors.length(); i++)
      metrics.Horizons.push_back(hors[i]);

    metrics.TrainFixSize = 0;
    metrics.TrainRatio = 0;
  } else {
    metrics.TrainFixSize = as<int>(metricsR["trainFixSize"]);
    metrics.TrainRatio = as<double>(metricsR["trainRatio"]);
  }

  // set minimum value for metrics with AIC weight formula
  auto minMetrics = as<List>(metricsR["minMetrics"]);

  metrics.minAic = helper_convertMinMetrics(minMetrics["aic"], numTargets);
  metrics.minSic = helper_convertMinMetrics(minMetrics["sic"], numTargets);
  metrics.minBrierIn =
    helper_convertMinMetrics(minMetrics["brierIn"], numTargets);
  metrics.minBrierOut =
    helper_convertMinMetrics(minMetrics["brierOut"], numTargets);
  metrics.minRmse = helper_convertMinMetrics(minMetrics["rmse"], numTargets);
  metrics.minMae = helper_convertMinMetrics(minMetrics["mae"], numTargets);
  metrics.minRmspe = helper_convertMinMetrics(minMetrics["rmspe"], numTargets);
  metrics.minMape = helper_convertMinMetrics(minMetrics["mape"], numTargets);
  metrics.minCrps = helper_convertMinMetrics(minMetrics["crps"], numTargets);

  metrics.Update(isOutOfSampleRandom,
                 isTimeSeries); // update after filling metric vectors

  if (metrics.MetricsIn.size() > 0) {
    for (auto i = 0; i < lmetricIn; i++) {
      auto str = ToString(metrics.MetricsIn.at(i), true);
      metricNames.push_back(str);
      if (printMsg) {
        Rprintf(str);
        if (i != lmetricIn - 1)
          Rprintf(", ");
      }
    }
  }

  if (metrics.MetricsOut.size() > 0) {

    for (auto i = 0; i < lmetricOut; i++) {
      auto str = ToString(metrics.MetricsOut.at(i), true);
      metricNames.push_back(str);
      if (printMsg) {
        Rprintf(str);
        if (i != lmetricOut - 1)
          Rprintf(", ");
      }
    }
  }

  if (isDc)
    metrics.WeightedEval = as<bool>(metrics["weightedEval"]);
}

void UpdateOptions(List &itemsR,
                   List &metricsR,
                   List &modelChecksR,
                   SearchMetricOptions &res_metric,
                   SearchItems &res_items,
                   SearchModelChecks &res_checks,
                   std::vector<std::string> &metricsNames,
                   int length1,
                   int exoCount,
                   int numTargets,
                   int numDependents,
                   bool isTimeSeries,
                   bool type1NeedsModelEstim,
                   const char *length1Informtion,
                   bool isDc) {


  UpdatemetricOptions(metrics, res_metric, metricsNames, isTimeSeries,
                      isDc, numTargets);

  UpdateSearchItems(items, res_items, length1, 0,
                    length1Informtion, nullptr, type1NeedsModelEstim, false);

  UpdateModelCheckItems(modelChecks, res_checks, res_metric,
                        res_items);

  res_items.LengthTargets = numTargets; // Modelset will use it
  res_items.LengthDependents = numDependents;
  res_items.LengthExogenouses = exoCount;
}








std::unique_ptr<double[]>
CombineEndoExo(bool printMsg, ldt::Matrix<double> &result,
               std::vector<std::string> &colNames, ldt::Matrix<double> &my,
               ldt::Matrix<double> &mx, ldt::Matrix<double> &mw,
               ldt::Matrix<double> &mnewX, SEXP &y, SEXP &x, SEXP &w,
               SEXP &newX, bool removeNan, bool addIntercept, int minExpectedY,
               int minExpectedX, int minExpectedW, int minExpectedNewX,
               bool appendNewX) {
  if (minExpectedY == 0 && minExpectedX == 0)
    throw LdtException(ErrorType::kLogic, "R-ldt",
                       "combining two null matrices is not implemented");

  double *y_ = nullptr;
  double *x_ = nullptr;
  double *w_ = nullptr;
  double *newX_ = nullptr;
  int yCols = 0, xCols = 0, yRows = 0, xRows = 0, wRows = 0, newxRows = 0,
      newxCols = 0;

  if (minExpectedY > 0 && y == R_NilValue)
    throw LdtException(ErrorType::kLogic, "R-ldt", "'y' is empty");
  if (minExpectedX > 0 && x == R_NilValue)
    throw LdtException(ErrorType::kLogic, "R-ldt", "'x' is empty");
  if (minExpectedW > 0 && w == R_NilValue)
    throw LdtException(ErrorType::kLogic, "R-ldt", "'w' is empty");

  NumericMatrix y0;
  if (y != R_NilValue && is<NumericMatrix>(y) == false)
    throw LdtException(ErrorType::kLogic, "R-ldt",
                       "'y' must be a 'numeric matrix'");
  else
    y0 = as<NumericMatrix>(y);

  if (y0.nrow() == 0 && y0.ncol() == 0)
    throw LdtException(ErrorType::kLogic, "R-ldt",
                       "invalid data: 'y' is empty");
  y_ = &y0[0];
  yCols = y0.ncol();
  yRows = y0.nrow();
  // names
  SEXP ynames = colnames(y0);

  if (ynames != R_NilValue) {
    StringVector ynames0(ynames);
    for (int i = 0; i < ynames0.length(); i++)
      colNames.push_back(as<std::string>(ynames0[i]));
  } else
    for (int i = 0; i < yCols; i++)
      colNames.push_back(std::string("Y") + std::to_string(i + 1));

  if (w != R_NilValue) {
    if (is<NumericMatrix>(w) == false)
      throw LdtException(ErrorType::kLogic, "R-ldt",
                         "'w' must be a 'numeric matrix'");
    auto w0 = as<NumericMatrix>(w);
    w_ = &w0[0];
    wRows = w0.nrow();
    colNames.push_back("Weight");
  }

  if (addIntercept) {
    colNames.push_back("Intercept");
    if (x == R_NilValue)
      x = NumericMatrix(yRows, 0);
  }

  if (x != R_NilValue) {
    if (is<NumericMatrix>(x) == false)
      throw LdtException(ErrorType::kLogic, "R-ldt",
                         "'x' must be a 'numeric matrix'");
    NumericMatrix x0 = as<NumericMatrix>(x);
    if (x0.nrow() == 0 && x0.ncol() == 0)
      throw LdtException(ErrorType::kLogic, "R-ldt",
                         "invalid data: 'x' is empty");
    x_ = &x0[0];
    xCols = x0.ncol();
    xRows = x0.nrow();
    // names
    if (xCols > 0) {
      SEXP xnames = colnames(x0);
      if (xnames != R_NilValue) {
        StringVector xnames0(xnames);
        for (int i = 0; i < xnames0.length(); i++)
          colNames.push_back(as<std::string>(xnames0[i]));
      } else
        for (int i = 0; i < xCols; i++)
          colNames.push_back(std::string("X") + std::to_string(i + 1));
    }

    // handle new data if x exists
    if (newX != R_NilValue) {
      if (is<NumericMatrix>(newX) == false)
        throw LdtException(ErrorType::kLogic, "R-ldt",
                           "'newX' must be a 'numeric matrix'");
      NumericMatrix newX0 = as<NumericMatrix>(newX);

      newX_ = &newX0[0];
      newxCols = newX0.ncol();
      newxRows = newX0.nrow();

      if (newxCols != xCols + (int)addIntercept) {

        Rcout << "Number of columns in:\n    - newX = " << newxCols
              << "\n    - x = " << (xCols + (int)addIntercept) << "\n";
        throw LdtException(ErrorType::kLogic, "R-ldt",
                           "invalid number of columns in 'newX'. It must "
                           "be equal to the number of columns in 'x'");
      }
    }

    if (newxRows < minExpectedNewX) {
      Rcout << "Required number of new exogenous data = " << minExpectedNewX;
      throw LdtException(ErrorType::kLogic, "R-ldt",
                         "there is not enough new exogenous data");
    }
  }

  if (y != R_NilValue && x != R_NilValue) {
    if (yRows != xRows)
      throw LdtException(
          ErrorType::kLogic, "R-ldt",
          "invalid data. Different number of observations in 'x' and 'y'");
  }
  if (y != R_NilValue && w != R_NilValue) {
    if (yRows != wRows)
      throw LdtException(
          ErrorType::kLogic, "R-ldt",
          "invalid data. Different number of observations in 'w' and 'y'");
  }

  mx.SetData(x_, xRows, xCols);
  my.SetData(y_, yRows, yCols);
  if (w_)
    mw.SetData(w_, wRows, 1);
  if (newX_)
    mnewX.SetData(newX_, newxRows, newxCols);

  if (printMsg) {
    Rprintf("Data:\n");
    Rprintf("    - Number of 'Observations' = %i\n", my.RowsCount);
    Rprintf("    - Number of 'Endogenous' = %i\n", my.ColsCount);
    Rprintf("    - Number of 'Weight' = %i\n", mw.ColsCount);
    Rprintf("    - Number of 'Exogenous' = %i\n", mx.ColsCount);
    Rprintf("    - Adds 'Intercept' = %s\n", addIntercept ? "true" : "false");
  }

  int mat_rows = my.RowsCount + (appendNewX ? newxRows : 0);
  int mat_cols = mx.ColsCount + my.ColsCount + mw.ColsCount + (int)addIntercept;
  auto mat_data = std::unique_ptr<double[]>(new double[mat_rows * mat_cols]);
  result.SetData(NAN, mat_data.get(), mat_rows, mat_cols);
  int j = 0;
  if (y_) {
    result.SetSub(0, j, my, 0, 0, my.RowsCount, my.ColsCount);
    j += my.ColsCount;
  }

  // Rcout<< "\n\nResult after y:" << result.ToString('\t','\n',2) << "\n\n";

  if (w_) {
    result.SetSub(0, j, mw, 0, 0, mw.RowsCount, mw.ColsCount);
    j += mw.ColsCount;
  }
  if (addIntercept) {
    result.SetColumn(j, 1);
    j++;
  }

  // Rcout<< "\n\nResult after intercept:" << result.ToString('\t','\n',2) <<
  // "\n\n";

  if (x_) {
    result.SetSub(0, j, mx, 0, 0, mx.RowsCount, mx.ColsCount);

    // Rcout<< "\n\nResult after x:" << result.ToString('\t','\n',2) << "\n\n";

    if (appendNewX) {
      result.SetSub(mx.RowsCount,
                    addIntercept ? j - 1 : j, // newX has intercept
                    mnewX, 0, 0, mnewX.RowsCount, mnewX.ColsCount);
    }
  }

  if (removeNan) {
    auto dataset =
        ldt::Dataset<double>(result.RowsCount, result.ColsCount, true, false);
    auto mat_data0 = std::unique_ptr<double[]>(new double[dataset.StorageSize]);
    dataset.Calculate(result, nullptr, mat_data0.get());

    if (dataset.Result.RowsCount != result.RowsCount) {
      result.SetData(dataset.Result.Data, dataset.Result.RowsCount,
                     dataset.Result.ColsCount);
      if (printMsg)
        Rprintf("    - Number of 'Observations' (NAN Removed) = %i\n",
                dataset.Result.RowsCount);
      mat_data = std::move(mat_data0);
    }
  }

  if (printMsg) { // TODO: print a summary
    // Rcout<< "\n---\n" << result.ToString_R_Matrix(2) << "\n---\n";
  }

  return mat_data;
}

void GetSizes(bool printMsg, std::vector<int> &result, SEXP &sizes,
              int variableCount, bool isX) {

  if (sizes == R_NilValue)
    result.push_back(1);
  else {
    if (is<IntegerVector>(sizes) == false)
      throw LdtException(ErrorType::kLogic, "R-ldt",
                         "'sizes' must be an 'integer vector'");
    IntegerVector sizes_ = as<IntegerVector>(sizes);
    for (int i = 0; i < sizes_.length(); i++)
      result.push_back((sizes_[i]));
  }
  if (result.size() == 0 ||
      *std::min_element(result.begin(), result.end()) < 1 ||
      *std::max_element(result.begin(), result.end()) > variableCount)
    throw LdtException(
        ErrorType::kLogic, "R-ldt",
        "invalid sizes array. It cannot be empty and elements must larger than "
        "1 and less than the number of variables");

  if (printMsg) {
    if (isX)
      Rprintf("Exogenous Sizes=%s\n", VectorToCsv<int>(result, ',').c_str());
    else
      Rprintf("Endogenous Sizes=%s\n", VectorToCsv<int>(result, ',').c_str());
  }
}

void GetPartitions(bool printMsg, std::vector<std::vector<int>> &result,
                   SEXP &partitions, int variableCount, int adjustPos,
                   bool isX) {

  if (partitions != R_NilValue) {
    if (is<List>(partitions) == false)
      throw LdtException(ErrorType::kLogic, "R-ldt",
                         "'partitions' must be a 'List'");
    List partitions0 = as<List>(partitions);
    for (int i = 0; i < partitions0.length(); i++) {
      if (is<IntegerVector>(partitions0[i]) == false)
        throw LdtException(ErrorType::kLogic, "R-ldt",
                           "'partitions[i]' must be an 'integer vector'");
      IntegerVector par_i = as<IntegerVector>(partitions0[i]);
      auto ar = std::vector<int>();
      for (int j = 0; j < par_i.length(); j++) {
        ar.push_back(par_i[j] + adjustPos -
                     1); // adjust elements for zero-based indexing and their
                         // position in bind data
      }
      result.push_back(ar);
    }
  } else {
    for (int i = adjustPos; i < adjustPos + variableCount; i++)
      result.push_back(std::vector<int>({i}));
  }

  for (auto &p : result) {
    for (auto &pi : p) {
      if (pi > adjustPos + variableCount) {
        Rcout << "Position Adjustment =" << adjustPos
              << "\nNumber of Variables = " << variableCount
              << "\nIndex of Element = " << pi << "\n";
        throw LdtException(
            ErrorType::kLogic, "R-ldt",
            "invalid element in a partition. Elements cannot be larger than "
            "the number of variables");
      }
    }
  }

  if (printMsg) {
    if (isX)
      Rprintf("Number of Exogenous Partitons=%i\n", (int)result.size());
    else
      Rprintf("Number of Endogenous Partitons=%i\n", (int)result.size());
    for (int i = 0; i < (int)result.size(); i++) {
      if (i < 10)
        Rprintf(" %i. Partition:%s\n", i,
                VectorToCsv<int>(result.at(i), ',').c_str());
      else if (i == 10)
        Rprintf("     . . .\n");
    }
  }
}

void GetGroups(bool printMsg, std::vector<std::vector<int>> &result,
               SEXP &groups, int variableCount, int adjustPos, bool isX) {

  if (groups != R_NilValue) {

    if (is<List>(groups) == false)
      throw LdtException(ErrorType::kLogic, "R-ldt",
                         "'groups' must be a 'List'");
    List groups0 = as<List>(groups);
    for (int i = 0; i < groups0.length(); i++) {
      if (is<IntegerVector>(groups0[i]) == false)
        throw LdtException(ErrorType::kLogic, "R-ldt",
                           "'groups[i]' must be an 'integer vector'");
      IntegerVector gro_i = as<IntegerVector>(groups0[i]);
      auto ar = std::vector<int>();
      for (int j = 0; j < gro_i.length(); j++) {
        ar.push_back(gro_i[j] + adjustPos -
                     1); // adjust elements for zero-based indexing and their
                         // position in bind data
      }
      result.push_back(ar);
    }
  } else if (variableCount >
             0) { // just one element (TODO: combination as an option)
    auto os = std::vector<int>(variableCount);
    std::iota(os.begin(), os.end(), adjustPos);
    result.push_back(os);
  }

  if (printMsg) {
    if (isX)
      Rprintf("Number of Exogenous Groups=%i\n", (int)result.size());
    else
      Rprintf("Number of Endogenous Groups=%i\n", (int)result.size());
    for (int i = 0; i < (int)result.size(); i++) {
      if (i < 10)
        Rprintf(" %i. Group:%s\n", i,
                VectorToCsv<int>(result.at(i), ',')
                    .c_str()); // TODO: print the names, indices are misleading
      else if (i == 10)
        Rprintf("     . . .\n");
    }
  }

  for (auto &a : result) {
    for (auto &b : a)
      if (b > variableCount + adjustPos || b + adjustPos < 0) {
        Rcout << "---------------\n";
        Rcout << "Position Adjustment = " << adjustPos << "\n";
        Rcout << "Element of a Groups = " << b << "\n";
        throw LdtException(ErrorType::kLogic, "R-ldt",
                           "invalid variable group. An element is negative "
                           "or larger than the number of variables");
      }
  }
}



std::vector<std::string> GetDefaultColNames(std::string pre, int length) {
  std::vector<std::string> nms;
  for (int i = 1; i <= length; i++)
    nms.push_back(pre + std::to_string(i));
  return nms;
}

NumericMatrix as_matrix(ldt::Matrix<double> &mat,
                        std::vector<std::string> *rowNames,
                        std::vector<std::string> *colNames) {
  NumericMatrix res = NumericMatrix(mat.RowsCount, mat.ColsCount, mat.Data);
  if (rowNames) {
    if ((int)rowNames->size() != mat.RowsCount) {
      Rcout << "---------------\n";
      Rcout << "Number of Rows: " << mat.RowsCount << "\n";
      Rcout << "Row Names:" << VectorToCsv(*rowNames) << "\n";
      throw LdtException(ErrorType::kLogic, "R-ldt",
                         std::string("invalid number of rows/row_names."));
    }
    rownames(res) = wrap(*rowNames);
  }
  if (colNames) {
    if ((int)colNames->size() != mat.ColsCount) {
      Rcout << "---------------\n";
      Rcout << "Number of Columns: " << mat.ColsCount << "\n";
      Rcout << "Column Names:" << VectorToCsv(*colNames) << "\n";
      throw LdtException(ErrorType::kLogic, "R-ldt",
                         std::string("invalid number of columns/col_names."));
    }
    colnames(res) = wrap(*colNames);
  }
  return res;
}

NumericVector as_vector(ldt::Matrix<double> &vec,
                        std::vector<std::string> *names) {
  auto res = NumericVector(vec.Data, vec.Data + vec.length());
  if (names) {
    if ((int)names->size() != vec.length()) {
      Rcout << "names:" << VectorToCsv(*names);
      throw LdtException(ErrorType::kLogic, "R-ldt",
                         std::string("invalid number of elements/names."));
    }
    res.names() = wrap(*names);
  }
  return res;
}

IntegerMatrix as_imatrix(ldt::Matrix<int> &mat,
                         std::vector<std::string> *rowNames,
                         std::vector<std::string> *colNames) {
  auto res = IntegerMatrix(mat.RowsCount, mat.ColsCount, mat.Data);
  if (rowNames) {
    if ((int)rowNames->size() != mat.RowsCount) {
      Rcout << "Row names:" << VectorToCsv(*rowNames);
      throw LdtException(ErrorType::kLogic, "R-ldt",
                         std::string("invalid number of rows/row_names."));
    }
    rownames(res) = wrap(*rowNames);
  }
  if (colNames) {
    if ((int)colNames->size() != mat.ColsCount) {
      Rcout << "Column names:" << VectorToCsv(*colNames);
      throw LdtException(ErrorType::kLogic, "R-ldt",
                         std::string("invalid number of columns/col_names."));
    }
    colnames(res) = wrap(*colNames);
  }
  return res;
}

IntegerVector as_ivector(ldt::Matrix<int> &vec,
                         std::vector<std::string> *names) {
  auto res = IntegerVector(vec.Data, vec.Data + vec.length());
  if (names) {
    if ((int)names->size() != vec.length()) {
      Rcout << "names:" << VectorToCsv(*names);
      throw LdtException(ErrorType::kLogic, "R-ldt",
                         std::string("invalid number of elements/names."));
    }
    res.names() = wrap(*names);
  }
  return res;
}

NumericMatrix cbind_matrix(NumericMatrix a, NumericMatrix b) {
  int n1 = a.ncol();
  int n2 = b.ncol();
  NumericMatrix result = no_init_matrix(a.nrow(), n1 + n2);
  StringVector nms(n1 + n2);
  StringVector anms = colnames(a);
  StringVector bnms = colnames(b);
  for (int i = 0; i < n1 + n2; i++) {
    if (i < n1) {
      result(_, i) = a(_, i);
      nms[i] = anms[i];
    } else {
      result(_, i) = b(_, i - n1);
      nms[i] = bnms[i - n1];
    }
  }

  colnames(result) = nms;
  return result;
}

NumericMatrix insert_intercept(NumericMatrix a) {
  StringVector nms;
  NumericMatrix result;

  int n1 = a.ncol();
  result = no_init_matrix(a.nrow(), n1 + 1);
  nms = StringVector(n1 + 1);
  StringVector anms = colnames(a);
  for (int i = 0; i < 1 + n1; i++) {
    if (i < 1) {
      result(_, i) = rep(1, a.nrow());
      nms[i] = "Intercept";
    } else {
      result(_, i) = a(_, i - 1);
      if (anms.length() > 0)
        nms[i] = anms[i - 1];
    }
  }

  colnames(result) = nms;
  return result;
}

NumericMatrix cbind_vectormatrix(NumericVector a, NumericMatrix b,
                                 std::string vectorName) {

  int n2 = b.ncol();
  NumericMatrix result = no_init_matrix(a.length(), 1 + n2);
  StringVector nms(1 + n2);
  StringVector bnms = colnames(b);
  for (int i = 0; i < 1 + n2; i++) {
    if (i < 1) {
      result(_, i) = a;
      nms[i] = vectorName;
    } else {
      result(_, i) = b(_, i - 1);
      nms[i] = bnms[i - 1];
    }
  }

  colnames(result) = nms;
  return result;
}

void ReportProgress(bool printMsg, int reportInterval, ModelSet &model,
                    bool &estimating, SearchOptions &options, int allCount) {
  auto start = std::chrono::system_clock::now();
  if (printMsg)
    Rprintf("Calculations Started ...\n");
  int c = 0;
  if (printMsg)
    Rprintf("Expected Number of Models = %i\n", allCount);
  double all = (double)allCount;
  int i = 0;
  while (estimating) {

    std::this_thread::sleep_for(std::chrono::seconds(1));
    // boost::this_thread::sleep_for(boost::chrono::milliseconds(1000));

    try {
      Rcpp::checkUserInterrupt();
    } catch (...) {
      options.RequestCancel = true;
      throw;
    }

    i++;
    if (reportInterval == 0 || i <= reportInterval)
      continue;
    i = 0;

    auto now = std::chrono::system_clock::now();
    double prePecentage = -1;

    c = model.GetNumberOfEstimatedModels();
    double percentage = std::round((c / all) * 10000) / 100;
    if (percentage != prePecentage) {
      std::chrono::duration<double> elapsed_mins = (now - start) / 60.0;
      double remains_mins = (all - c) * elapsed_mins.count() / c;
      if (printMsg)
        Rprintf("    Searched=%i, All=%i  (%.2f%%, %.1f minutes remains)\n", c,
                allCount, percentage > 100 || percentage < 0 ? NAN : percentage,
                remains_mins < 0 ? NAN : remains_mins);
      prePecentage = percentage;
    }
  }

  if (options.RequestCancel)
    throw LdtException(ErrorType::kLogic, "R-ldt", "calculations is canceled");
  else if (printMsg)
    Rprintf("Calculations Ended.\n");
}

static void add_CoefInfo(List &L, std::vector<EstimationKeep *> &list,
                         const char *extra1Label, int exoIndexesPlus,
                         std::vector<std::string> *extra1Names,
                         const char *label, bool addCoefs,
                         std::vector<std::string> &nms, int startIndex = 0) {

  int j = -1;
  for (auto &b : list) {
    j++;
    nms.push_back(std::string(label) + std::to_string(j + 1));

    IntegerVector deps;
    IntegerVector exos;
    if (b->Dependents.Data) {
      deps = as_ivector(b->Dependents);
      deps = deps + 1;
    }
    if (b->Exogenouses.Data) {
      exos = as_ivector(b->Exogenouses);
      exos = exos + exoIndexesPlus;
    }

    List L_j = List::create(
        _["metric"] = b->Metric, _["weight"] = b->Weight,
        _["depIndices"] = b->Dependents.Data ? (SEXP)deps : R_NilValue,
        _["exoIndices"] = b->Exogenouses.Data ? (SEXP)exos : R_NilValue,
        _["mean"] = addCoefs ? wrap(b->Mean) : R_NilValue,
        _["var"] = addCoefs ? wrap(b->Variance) : R_NilValue,
        _[extra1Label] =
            b->Extra.Data ? (SEXP)as_ivector(b->Extra) : R_NilValue);

    if (L_j[6] != R_NilValue)
      ((IntegerVector)L_j[6]).names() = wrap(*extra1Names);

    L[j + startIndex] = L_j;
  }
  L.names() = wrap(nms);
}

static void add_Lengthi(List L, int eIndex, int tIndex, ModelSet &model,
                        std::vector<SearcherSummary *> &list,
                        SearchItems &items, const char *extra1Label,
                        int length1, std::vector<std::string> &length1Names,
                        std::vector<std::string> *extra1Names,
                        int exoIndexesPlus,
                        const char *length1_itemlabel = "item") {

  if (items.KeepBestCount > 0) {
    List L_0 = List(length1);
    std::vector<std::string> L_0_names;
    for (auto i = 0; i < length1; i++) {

      auto L_0_i = List(items.KeepBestCount + 1);
      std::vector<std::string> L_0_i_names;
      L_0_i_names.push_back(std::string("name"));

      L_0_i[0] = wrap(length1Names.at(i));

      auto bests = std::vector<EstimationKeep *>();
      model.CombineBests(eIndex, tIndex, i, list, bests);

      if (bests.size() != 0) {
        add_CoefInfo(L_0_i, bests, extra1Label, exoIndexesPlus, extra1Names,
                     "best", true, L_0_i_names, 1);
        L_0[i] = L_0_i;
      } else {
        L_0[i] = R_NilValue;
      }
      L_0_names.push_back(std::string(length1_itemlabel) +
                          std::to_string(i + 1));
    }

    L_0.names() = L_0_names;
    L[0] = L_0;
  } else
    L[0] = R_NilValue;

  // All ? do we really need all

  if (items.CdfsAt.size() > 0) {

    List L_1 = List(items.CdfsAt.size());
    std::vector<std::string> L_1_names;

    for (int k = 0; k < (int)items.CdfsAt.size(); k++) {
      auto cdf = RunningMoments<1, true, true, Tv>();
      auto mat_d = std::unique_ptr<double[]>(new double[length1 * 3]);
      auto mat = ldt::Matrix<double>(mat_d.get(), length1, 3);
      auto colnames = std::vector<std::string>({"Mean", "Count", "SumWeights"});
      for (auto i = 0; i < length1; i++) {
        model.CombineCdfAt(eIndex, tIndex, i, k, list, cdf);
        mat.Set0(i, 0, cdf.GetMean());
        mat.Set0(i, 1, (double)cdf.Count);
        mat.Set0(i, 2, cdf.SumWeights);
      }
      L_1[k] = as_matrix(mat, &length1Names, &colnames);
      L_1_names.push_back(std::string("cdf") + std::to_string(k + 1));
    }
    L_1.names() = L_1_names;
    L[1] = L_1;
  } else
    L[1] = R_NilValue;

  if (items.ExtremeBoundsMultiplier > 0) {
    double min = 0, max = 0;
    auto mat_d = std::unique_ptr<double[]>(new double[length1 * 2]);
    auto mat = ldt::Matrix<double>(mat_d.get(), length1, 2);
    auto colnames = std::vector<std::string>({"Lower", "Upper"});
    for (auto i = 0; i < length1; i++) {
      model.CombineExtremeBounds(eIndex, tIndex, i, list, min, max);
      mat.Set0(i, 0, min);
      mat.Set0(i, 1, max);
    }

    L[2] = as_matrix(mat, &length1Names, &colnames);
  } else
    L[2] = R_NilValue;

  if (items.KeepMixture) {
    auto mixture = RunningMoments<4, true, true, Tv>();
    auto mat_d = std::unique_ptr<double[]>(new double[length1 * 6]);
    auto mat = ldt::Matrix<double>(mat_d.get(), length1, 6);
    auto colnames = std::vector<std::string>(
        {"Mean", "Variance", "Skewness", "Kurtosis", "Count", "SumWeights"});
    for (auto i = 0; i < length1; i++) {
      model.CombineMixture(eIndex, tIndex, i, list, mixture);
      mat.Set0(i, 0, mixture.GetMean());
      mat.Set0(i, 1, (double)mixture.GetVariance());
      mat.Set0(i, 2, (double)mixture.GetSkewness());
      mat.Set0(i, 3, (double)mixture.GetKurtosis());
      mat.Set0(i, 4, (double)mixture.Count);
      mat.Set0(i, 5, (double)mixture.SumWeights);
    }
    L[3] = as_matrix(mat, &length1Names, &colnames);
  } else
    L[3] = R_NilValue;
}

List GetModelSetResults(ModelSet &model, SearchItems &items,
                        std::vector<std::string> &metricNames, int length1,
                        const char *extra1Label,
                        std::vector<std::string> *extra1Names,
                        int exoIndexesPlus,
                        std::vector<std::string> &length1Names,
                        std::vector<std::string> &inclusionNames,
                        const char *length1Label, const char *length1_itemlabel,
                        bool printMsg) {
  // output structure:

  std::vector<std::string> namesL;
  namesL.push_back("counts");
  for (auto eIndex = 0; eIndex < items.LengthEvals; eIndex++)
    namesL.push_back(metricNames.at(eIndex));
  namesL.push_back("info");

  List L = List(1 + items.LengthEvals + 1);
  L.names() = namesL;

  // general information:
  auto result = SearcherModelingInfo();
  auto list0 = std::vector<SearcherSummary *>();
  auto list1 = std::vector<SearcherSummary *>();
  auto list2 = std::vector<SearcherSummary *>();
  model.CombineInfo(result, list0, list1, list2);
  int fcount = 0;
  List failDetails = List(result.FailsCount.size());
  int i = -1;
  for (const auto &[k, v] : result.FailsCount) {
    i++;
    fcount += v;
    failDetails[i] = List::create(_["message"] = wrap(k), _["count"] = wrap(v));
  }
  L[0] = List::create(_["expectedCount"] = wrap(result.ExpectedCount),
                      _["searchedCount"] = wrap(result.SearchedCount),
                      _["failedCount"] = wrap(fcount),
                      _["failedDetails"] = wrap(failDetails));
  if (fcount > 0 & printMsg)
    Rprintf("** Search process ended successfully. However, there are some "
            "failed estimations. See 'result$counts' for more details.");

  for (auto eIndex = 0; eIndex < items.LengthEvals; eIndex++) {

    List L_i = List(items.LengthTargets);
    std::vector<std::string> L_i_names;
    L[1 + eIndex] = L_i;

    for (auto tIndex = 0; tIndex < items.LengthTargets; tIndex++) {
      L_i_names.push_back(std::string("target") + std::to_string(tIndex + 1));

      List L_i_t = List(3); // when item2 implemented -> 4
      L_i_t.names() = std::vector<std::string>({"name", "model", length1Label});
      L_i[tIndex] = L_i_t;

      L_i_t[0] = wrap(inclusionNames.at(
          tIndex)); // expecting target name in inclusion names

      if (items.KeepModelEvaluations) {

        List L_i_t_m = List(3);
        L_i_t_m.names() =
            std::vector<std::string>({"bests", "all", "inclusion"});
        L_i_t[1] = L_i_t_m;

        if (items.KeepBestCount > 0) {
          auto bests = std::vector<EstimationKeep *>();
          model.CombineBests(eIndex, tIndex, 0, list0, bests);
          List L_i_t_m_0 = List(bests.size());
          std::vector<std::string> nms;
          add_CoefInfo(L_i_t_m_0, bests, extra1Label, exoIndexesPlus,
                       extra1Names, "best", false, nms, 0);
          L_i_t_m[0] = L_i_t_m_0;
        } else
          L_i_t_m[0] = R_NilValue;

        if (items.KeepAll) {
          auto all = std::vector<EstimationKeep *>();
          model.CombineAll(eIndex, tIndex, 0, list0, all);
          List L_i_t_m_1 = List(all.size());
          std::vector<std::string> nms;
          add_CoefInfo(L_i_t_m_1, all, extra1Label, exoIndexesPlus, extra1Names,
                       "model", false, nms, 0);
          L_i_t_m[1] = L_i_t_m_1;
        } else
          L_i_t_m[1] = R_NilValue;

        if (items.KeepInclusionWeights) {
          auto covars =
              items.LengthDependents + items.LengthExogenouses;
          auto incweights = RunningMoments<1, true, false, Tv>();
          auto mat_d = std::unique_ptr<double[]>(new double[covars * 2]);
          auto mat = ldt::Matrix<double>(mat_d.get(), covars, 2);
          auto colnames = std::vector<std::string>({"Mean", "Count"});
          for (auto i = 0; i < covars; i++) {
            model.CombineInclusionWeights(eIndex, tIndex, i, list0, incweights);
            mat.Set0(i, 0, incweights.GetMean());
            mat.Set0(i, 1, (double)incweights.Count);
          }
          L_i_t_m[2] = as_matrix(mat, &inclusionNames, &colnames);

        } else
          L_i_t_m[2] = R_NilValue;
      } else
        L_i_t[1] = R_NilValue; // model is null

      if (items.Length1 > 0) {
        List L_i_t_1 = List(4);
        L_i_t_1.names() = std::vector<std::string>(
            {"bests", "cdfs", "extremeBounds", "mixture"});
        L_i_t[2] = L_i_t_1;

        add_Lengthi(L_i_t_1, eIndex, tIndex, model, list1, items,
                    extra1Label, length1, length1Names, extra1Names,
                    exoIndexesPlus, length1_itemlabel);
        L_i_t[2] = L_i_t_1;
      } else
        L_i_t[2] = R_NilValue;

      if (items.Length2 > 0) {
        throw LdtException(ErrorType::kLogic, "R-ldt",
                           "not implemented: length2>0");
      }
    }
    L_i.names() = L_i_names;
  }

  return L;
}

// Update Options:

void UpdateRocOptions(bool printMsg, List &rocOptionsR, RocOptions &options,
                      const char *startMsg) {

  if (printMsg)
    Rprintf("%s:\n", startMsg);

  options.NormalizePoints = true;
  options.LowerThreshold = as<double>(rocOptionsR["lowerThreshold"]);
  options.UpperThreshold = as<double>(rocOptionsR["upperThreshold"]);
  options.Epsilon = as<double>(rocOptionsR["epsilon"]);

  if (rocOptionsR["costs"] != R_NilValue) {
    auto costs0 = as<NumericVector>(rocOptionsR["costs"]);
    auto costMatrix0 = as<NumericMatrix>(rocOptionsR["costMatrix"]);
    options.Costs.SetData(&costs0[0], costs0.length(), 1);
    options.CostMatrix.SetData(&costMatrix0[0], 2, 2);
  }

  if (printMsg) {
    if ((std::isnan(options.LowerThreshold) || options.LowerThreshold == 0) &&
        (std::isnan(options.UpperThreshold) || options.UpperThreshold == 1))
      Rprintf("    - Not Partial\n");
    else
      Rprintf("    - Partial (%f, %f):\n", options.LowerThreshold,
              options.UpperThreshold);
    Rprintf("    - Epsilon = %f\n", options.Epsilon);
    if (options.Costs.Data) {
      Rprintf("    - Varing Cost\n");
    }
  }
}

void UpdatePcaOptions(List pcaOptionsR,
                      PcaAnalysisOptions &options) {

  options.IgnoreFirstCount = as<int>(pcaOptionsR["ignoreFirst"]);
  options.ExactCount = as<int>(pcaOptionsR["exactCount"]);
  options.CutoffRate = as<double>(pcaOptionsR["cutoffRate"]);
  options.CutoffCountMax = as<int>(pcaOptionsR["max"]);
  if (options.IsEnabled())
    options.CheckValidity();

}

void UpdateLbfgsOptions(bool printMsg, List &lbfgsOptions,
                        LimitedMemoryBfgsbOptions &options) {
  if (printMsg)
    Rprintf("L-BFGS options:\n");
  options.Factor = as<double>(lbfgsOptions["factor"]);
  options.IterationMax = as<int>(lbfgsOptions["maxIterations"]);
  options.ProjectedGradientTol =
      as<double>(lbfgsOptions["projectedGradientTol"]);
  options.mMaxCorrections = as<int>(lbfgsOptions["maxCorrections"]);
  ;

  if (printMsg) {
    Rprintf("    - Maximum Number of Iterations=%i\n", options.IterationMax);
    Rprintf("    - Factor=%f\n", options.Factor);
    Rprintf("    - Projected Gradient Tolerance=%f\n",
            options.ProjectedGradientTol);
    Rprintf("    - Maximum Corrections=%i\n", options.mMaxCorrections);
  }
}

void UpdateNewtonOptions(bool printMsg, List &newtonR, Newton &newton) {

  if (printMsg)
    Rprintf("Newton Optimization Parameters:\n");

  newton.IterationMax = as<int>(newtonR["maxIterations"]);
  newton.TolFunction = as<double>(newtonR["functionTol"]);
  newton.TolGradient = as<double>(newtonR["gradientTol"]);
  newton.UseLineSearch = as<bool>(newtonR["useLineSearch"]);

  if (printMsg) {
    Rprintf("    - Iterations (Maximum)=%i\n", newton.IterationMax);
    Rprintf("    - Function Tolerance=%f\n", newton.TolFunction);
    Rprintf("    - Gradient Tolerance=%f\n", newton.TolGradient);
    Rprintf("    - Use Line Search=%s\n",
            newton.UseLineSearch ? "TRUE" : "FALSE");
  }
}


