// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp17)]]

#pragma once

#include "helpers.h"
#include "matrix.h"
#include "pca.h"
#include "searchers.h"
#include <Rcpp.h>
#include <chrono>
#include <cstring>
#include <exception>
#include <future>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

// #include <Rinternals.h>

using namespace Rcpp;
using namespace ldt;

void UpdateSearchData(List &dataR, ldt::SearchData &data);

void UpdateSearchCombinations(List combinationsR,
                              ldt::SearchCombinations &combinations);

void UpdateSearchOptions(List &optionsR, ldt::SearchOptions &options);

void UpdateSearchItems(List &itemsR, ldt::SearchItems &items, Ti length1,
                       Ti length2, const char *length1Informtion,
                       const char *length2Informtion, bool type1NeedsModelEstim,
                       bool type2NeedsModelEstim);

void UpdateModelCheckItems(List &checksR, ldt::SearchModelChecks &checks,
                           const ldt::SearchMetricOptions &metrics,
                           const ldt::SearchItems &items);

void UpdatemetricOptions(List &metricsR, ldt::SearchMetricOptions &metrics,
                         std::vector<std::string> &metricNames,
                         bool isTimeSeries, bool isDc, int numTargets);

void UpdateOptions(List &itemsR, List &metricsR, List &modelChecksR,
                   SearchMetricOptions &res_metric, SearchItems &res_items,
                   SearchModelChecks &res_checks,
                   std::vector<std::string> &metricsNames, int length1,
                   int exoCount, int numTargets, int numDependents,
                   bool isTimeSeries, bool type1NeedsModelEstim,
                   const char *length1Informtion, bool isDc);

void UpdatePcaOptions(List optionsR, ldt::PcaAnalysisOptions &options);

void ReportProgress(int reportInterval, ldt::ModelSet &model, bool &estimating,
                    ldt::SearchOptions &options, int allCount);

NumericMatrix as_matrix(
    const ldt::Matrix<double> &mat,
    const std::vector<std::string> &rowNames = std::vector<std::string>(),
    const std::vector<std::string> &colNames = std::vector<std::string>());

NumericVector
as_vector(const ldt::Matrix<double> &vec,
          const std::vector<std::string> &names = std::vector<std::string>());

IntegerMatrix as_imatrix(
    const ldt::Matrix<int> &mat,
    const std::vector<std::string> &rowNames = std::vector<std::string>(),
    const std::vector<std::string> &colNames = std::vector<std::string>());

IntegerVector
as_ivector(const ldt::Matrix<int> &vec,
           const std::vector<std::string> &names = std::vector<std::string>());

void UpdateRocOptions(bool printMsg, List &rocOptionsR,
                      ldt::RocOptions &options, const char *startMsg);

void UpdateLbfgsOptions(bool printMsg, List &lbfgsOptions,
                        ldt::LimitedMemoryBfgsbOptions &options);

void UpdateNewtonOptions(bool printMsg, List &newtonR, ldt::Newton &newton);

List GetModelSetResults(const ModelSet &model, const SearchItems &items,
                        const std::vector<std::string> &metricNames,
                        const std::vector<std::string> &targetNames,
                        const std::string &extra1Label,
                        const std::vector<std::string> &extra1Names,
                        const std::vector<std::string> &length1Names,
                        const std::vector<std::string> &inclusionNames,
                        const std::string length1Label, const bool &printMsg);
