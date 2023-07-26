#include "r_tsdata.h"

using namespace Rcpp;
using namespace ldt;

// [[Rcpp::export(.ConvertTo_Daily)]]
List ConvertTo_Daily(SEXP w, SEXP aggregateFun) {
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  Variable<double> v;
  UpdateVariableFromSEXP(w, v, listItems, listItemsDate);

  Variable<double> nv;

  if (aggregateFun == R_NilValue) {
    v.ConvertTo_Daily(nv, nullptr);
  } else if (is<Function>(aggregateFun)) {

    auto F = as<Function>(aggregateFun);
    std::function<double(const std::vector<double> &)> func =
        [&F](const std::vector<double> &x) { return as<double>(F(wrap(x))); };
    v.ConvertTo_Daily(nv, &func);

  } else if (TYPEOF(aggregateFun) == STRSXP) {

    auto descType = FromString_DescriptiveType(as<const char *>(aggregateFun));
    std::function<double(const std::vector<double> &)> func =
        [&descType](const std::vector<double> &x) {
          double v;
          Array<double>::GetDescriptive<true>(&x[0], x.size(), descType,
                                              v); // TODO: true as an option?!
          return v;
        };
    v.ConvertTo_Daily(nv, &func);

  } else {
    throw LdtException(ErrorType::kLogic, "R-freq-convert",
                       "'aggregateFun' should be a character or a function");
  }

  return GetVariableForR(nv);
}

// [[Rcpp::export(.ConvertTo_MultiDaily)]]
List ConvertTo_MultiDaily(SEXP w, int k, SEXP aggregateFun, bool fromEnd) {
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  Variable<double> v;
  UpdateVariableFromSEXP(w, v, listItems, listItemsDate);

  Variable<double> nv;
  if (aggregateFun == R_NilValue) {
    v.ConvertTo_MultiDaily(nv, k, nullptr, fromEnd);
  } else if (is<Function>(aggregateFun)) {

    auto F = as<Function>(aggregateFun);
    std::function<double(const std::vector<double> &)> func =
        [&F](const std::vector<double> &x) { return as<double>(F(wrap(x))); };
    v.ConvertTo_MultiDaily(nv, k, &func, fromEnd);

  } else if (TYPEOF(aggregateFun) == STRSXP) {

    auto descType = FromString_DescriptiveType(as<const char *>(aggregateFun));
    std::function<double(const std::vector<double> &)> func =
        [&descType](const std::vector<double> &x) {
          double v;
          Array<double>::GetDescriptive<true>(&x[0], x.size(), descType,
                                              v); // TODO: true as an option?!
          return v;
        };
    v.ConvertTo_MultiDaily(nv, k, &func, fromEnd);

  } else {
    throw LdtException(
        ErrorType::kLogic, "R-freq-convert",
        "invalid 'aggregateFun'. It should be a character or a function");
  }

  return GetVariableForR(nv);
}

// [[Rcpp::export(.ConvertTo_Weekly)]]
List ConvertTo_Weekly(SEXP w, const char *weekStart, SEXP aggregateFun) {
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  Variable<double> v;
  UpdateVariableFromSEXP(w, v, listItems, listItemsDate);

  auto weekstartType = FromString_DayOfWeek(weekStart);

  Variable<double> nv;
  if (aggregateFun == R_NilValue) {
    v.ConvertTo_Weekly(nv, weekstartType, nullptr);
  } else if (is<Function>(aggregateFun)) {

    auto F = as<Function>(aggregateFun);
    std::function<double(const std::vector<double> &)> func =
        [&F](const std::vector<double> &x) { return as<double>(F(wrap(x))); };
    v.ConvertTo_Weekly(nv, weekstartType, &func);

  } else if (TYPEOF(aggregateFun) == STRSXP) {

    auto descType = FromString_DescriptiveType(as<const char *>(aggregateFun));
    std::function<double(const std::vector<double> &)> func =
        [&descType](const std::vector<double> &x) {
          double v;
          Array<double>::GetDescriptive<true>(&x[0], x.size(), descType,
                                              v); // TODO: true as an option?!
          return v;
        };
    v.ConvertTo_Weekly(nv, weekstartType, &func);

  } else {
    throw LdtException(
        ErrorType::kLogic, "R-freq-convert",
        "invalid 'aggregateFun'. It should be a character or a function");
  }
  return GetVariableForR(nv);
}

// [[Rcpp::export(.ConvertTo_XxYear)]]
List ConvertTo_XxYear(SEXP w, SEXP aggregateFun, int x) {
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  Variable<double> v;
  UpdateVariableFromSEXP(w, v, listItems, listItemsDate);

  Variable<double> nv;
  if (aggregateFun == R_NilValue) {

    if (x == 24)
      v.ConvertTo_XxYear_month_based<24>(nv, nullptr);
    else if (x == 12)
      v.ConvertTo_XxYear_month_based<12>(nv, nullptr);
    else if (x == 6)
      v.ConvertTo_XxYear_month_based<6>(nv, nullptr);
    else if (x == 4)
      v.ConvertTo_XxYear_month_based<4>(nv, nullptr);
    else if (x == 3)
      v.ConvertTo_XxYear_month_based<3>(nv, nullptr);
    else if (x == 2)
      v.ConvertTo_XxYear_month_based<2>(nv, nullptr);
    else if (x == 1)
      v.ConvertTo_XxYear_month_based<1>(nv, nullptr);
    else
      v.ConvertTo_XxYear(nv, x, nullptr);

  } else if (is<Function>(aggregateFun)) {

    auto F = as<Function>(aggregateFun);
    std::function<double(const std::vector<double> &)> func =
        [&F](const std::vector<double> &x) { return as<double>(F(wrap(x))); };

    if (x == 24)
      v.ConvertTo_XxYear_month_based<24>(nv, &func);
    else if (x == 12)
      v.ConvertTo_XxYear_month_based<12>(nv, &func);
    else if (x == 6)
      v.ConvertTo_XxYear_month_based<6>(nv, &func);
    else if (x == 4)
      v.ConvertTo_XxYear_month_based<4>(nv, &func);
    else if (x == 3)
      v.ConvertTo_XxYear_month_based<3>(nv, &func);
    else if (x == 2)
      v.ConvertTo_XxYear_month_based<2>(nv, &func);
    else if (x == 1)
      v.ConvertTo_XxYear_month_based<1>(nv, &func);
    else
      v.ConvertTo_XxYear(nv, x, &func);

  } else if (TYPEOF(aggregateFun) == STRSXP) {

    auto descType = FromString_DescriptiveType(as<const char *>(aggregateFun));
    std::function<double(const std::vector<double> &)> func =
        [&descType](const std::vector<double> &x) {
          double v;
          Array<double>::GetDescriptive<true>(&x[0], x.size(), descType,
                                              v); // TODO: true as an option?!
          return v;
        };

    if (x == 24)
      v.ConvertTo_XxYear_month_based<24>(nv, &func);
    else if (x == 12)
      v.ConvertTo_XxYear_month_based<12>(nv, &func);
    else if (x == 6)
      v.ConvertTo_XxYear_month_based<6>(nv, &func);
    else if (x == 4)
      v.ConvertTo_XxYear_month_based<4>(nv, &func);
    else if (x == 3)
      v.ConvertTo_XxYear_month_based<3>(nv, &func);
    else if (x == 2)
      v.ConvertTo_XxYear_month_based<2>(nv, &func);
    else if (x == 1)
      v.ConvertTo_XxYear_month_based<1>(nv, &func);
    else
      v.ConvertTo_XxYear(nv, x, &func);

  } else {
    throw LdtException(
        ErrorType::kLogic, "R-freq-convert",
        "invalid 'aggregateFun'. It should be a character or a function");
  }
  return GetVariableForR(nv);
}
