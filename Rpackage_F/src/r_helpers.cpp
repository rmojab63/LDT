#include "r_tsdata.h"

using namespace Rcpp;
using namespace ldt;

// [[Rcpp::export(.Get_Descriptive)]]
NumericVector Get_Descriptive(NumericVector w, const char* type, bool skipNAN){

  double res = NAN;
  auto descType = FromString_DescriptiveType(type);

  if (skipNAN)
    Array<double>::GetDescriptive<true>(&w[0],w.length(),descType,res);
  else
    Array<double>::GetDescriptive<true>(&w[0],w.length(),descType,res);

  NumericVector V = {res};
  V.names() = CharacterVector({ToString(descType)});

  return V;
}
