// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp17)]]

#pragma once

#include "frequency.h"
#include "helpers.h"
#include "variable.h"
#include <Rcpp.h>

#include <boost/date_time/posix_time/posix_time.hpp>

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

std::unique_ptr<ldt::FrequencyWeekBased> GetFreqFromSEXP_week(List f);

std::unique_ptr<ldt::Frequency>
GetFreqFromSEXP(SEXP value, std::vector<std::string> &listItems,
                std::vector<boost::gregorian::date> &listItemsDate);

std::unique_ptr<ldt::Variable<double>> GetVariableFromSEXP(List w);

