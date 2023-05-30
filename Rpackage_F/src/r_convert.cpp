#include "r_tsdata.h"

using namespace Rcpp;
using namespace ldt;

// [[Rcpp::export(.c_DateList_to_Daily)]]
List C_DateList_to_Date(SEXP w) {
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  Variable<double> v;
  UpdateVariableFromSEXP(w, v, listItems, listItemsDate);

  if (listItemsDate.size() == 0)
    throw std::logic_error("Invalid type of variable. Expected frequeny: list of dates.");

  Variable<double> nv;
  v.ConvertTo_Daily(nv);

  return GetVariableForR(nv);
}


// [[Rcpp::export(.c_Daily_to_MultiDaily)]]
List C_Daily_to_MultiDaily(SEXP w, int k, Function aggregateFun, bool fromEnd) {
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  Variable<double> v;
  UpdateVariableFromSEXP(w, v, listItems, listItemsDate);

 if (v.StartFrequency.get()->mClass != FrequencyClass::kDaily)
   throw std::logic_error("Invalid type of frequency. Expected type is daily.");

//TODO: you can make it more efficient for common functions such as mean, etc.
 std::function<double(const std::vector<double>&)> func =
   [&aggregateFun](const std::vector<double>& x) { return as<double>(aggregateFun(wrap(x))); };

   Variable nv;
   v.ConvertTo_MultiDaily(nv, k, func, fromEnd);

  return GetVariableForR(nv);
}
