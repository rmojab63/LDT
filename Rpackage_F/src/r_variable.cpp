#include "r_tsdata.h"

using namespace Rcpp;
using namespace ldt;

// [[Rcpp::export(.Variable)]]
List Variable(SEXP data, SEXP name, SEXP startFrequency, SEXP fields) {
  List V =
      List::create(_["data"] = data, _["name"] = name,
                   _["startFrequency"] = startFrequency, _["fields"] = fields);
  V.attr("class") = std::vector<std::string>({"ldtv", "list"});
  return V;
}

List GetVariableForR(ldt::Variable<double> &v) {
  auto sf = v.StartFrequency.get()->ToString();
  auto sfc = v.StartFrequency.get()->ToClassString();
  // TODO: fields

  List V = List::create(_["data"] = wrap(v.Data), _["name"] = wrap(v.Name),
                        _["startFrequency"] = Parse_F(sf, sfc),
                        _["fields"] = R_NilValue);

  V.attr("class") = std::vector<std::string>({"ldtv", "list"});
  return V;
}

void UpdateVariableFromSEXP(
    Rcpp::List w, ldt::Variable<double> &variable,
    std::vector<std::string> &listItems,
    std::vector<boost::gregorian::date> &listItemsDate) {

  if (w["name"] != R_NilValue)
    variable.Name = as<std::string>(w["name"]);

  try {
    variable.StartFrequency =
        GetFreqFromSEXP(w["startFrequency"], listItems, listItemsDate);
  } catch (...) {
    throw std::logic_error("Invalid 'startFrequency'.");
  }

  try {
    variable.Data = as<std::vector<double>>(w["data"]);
  } catch (...) {
    throw std::logic_error(
        "Invalid 'data'. It should be a vector of numerics.");
  }

  try {
    if (is<List>(w["fields"])) {
      List F = w["fields"];
      for (int j = 0; j < F.length(); j++) {
        CharacterVector Fj = as<CharacterVector>(F[j]);
        if (Fj.length() < 2)
          throw std::logic_error("Expected a 'key' and a 'value'.");
        variable.Fields.insert(
            {as<std::string>(Fj[0]), as<std::string>(Fj[1])});
      }
    }
  } catch (...) {
    throw std::logic_error(
        "Invalid fields. It should be a 'List' of 'CharacterVector's. The "
        "vector must have 'Key-Value' pairs.");
  }
}

// [[Rcpp::export(.VariableToString)]]
std::string VariableToString(List w) {
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  auto v = ldt::Variable<double>();
  UpdateVariableFromSEXP(w, v, listItems, listItemsDate);
  return v.ToString();
}

// [[Rcpp::export(.BindVariables)]]
List BindVariables(SEXP varList, bool interpolate, bool adjustLeadLags,
                   int numExo, int horizon) {
  auto vars = as<List>(varList);
  int n = vars.length();
  if (n == 0)
    throw std::logic_error("Empty list of variables.");

  int numEndo = n - numExo;
  if (numExo < 0 || numEndo < 0)
    throw std::logic_error(
        "Invalid number of exogenous/exogenous variables (check 'numExo').");
  if (numExo > 0 && horizon < 0)
    throw std::logic_error("Invalid out-of-sample size (check 'horizon').");

  auto info = IntegerMatrix(n, 5);
  colnames(info) =
      wrap(std::vector<const char *>({"Start_Index", "End_Index", "Has_Missing",
                                      "Interpolated_Count", "Lags_Leads"}));

  auto listItems = std::vector<std::vector<std::string>>(n);
  auto listItemsDate = std::vector<std::vector<boost::gregorian::date>>(n);

  auto list = std::vector<ldt::Variable<double>>(n);
  auto list0 = std::vector<ldt::Variable<double> *>(n);
  try {
    for (int i = 0; i < n; i++) {
      if (is<List>(vars[i]) == false)
        throw std::logic_error("Invalid variable type.");
      UpdateVariableFromSEXP(as<List>(vars[i]), list.at(i), listItems.at(i),
                             listItemsDate.at(i));
      list.at(i).Trim();
      list0.at(i) = &list.at(i);

      int count = 0;
      if (interpolate) {
        list.at(i).Interpolate(count);
      }
      info(i, 3) = count;
    }
  } catch (...) {
    throw std::logic_error(
        "Creating a variable failed. Make sure all items of the list are "
        "'ldtv' objects.");
  }
  auto vs = Variables<double>(list0);

  std::vector<std::tuple<int, int>> ranges;

  for (int i = 0; i < n; i++) {
    bool hasMissing = false;
    auto range = vs.GetRange(i, hasMissing);
    info(i, 0) = std::get<0>(range) + 1;
    info(i, 1) = std::get<1>(range) + 1;
    info(i, 2) = hasMissing ? 1 : 0;
    ranges.push_back(range);

    if (adjustLeadLags) {
      // we should update data of the variable
      auto a = &vs.Data[i * vs.NumObs];
      list.at(i).Data = std::vector<double>(a, a + vs.NumObs);
      list.at(i).StartFrequency = vs.StartFrequency.get()->Clone();
      list0.at(i) = &list.at(i);
      list.at(i).Trim();
    }
  }

  // lags, leads
  info(0, 5) = 0;
  int lastIndex = std::get<1>(ranges.at(0));
  int lastIndexHor = lastIndex + horizon;
  int len;
  for (int i = 1; i < n; i++) {
    if (i < numEndo)
      len = lastIndex - std::get<1>(ranges.at(i));
    else
      len = lastIndexHor - std::get<1>(ranges.at(i));
    len *= -1;
    info(i, 4) = len;

    if (adjustLeadLags && len != 0) {
      auto v = list0.at(i);

      v->StartFrequency.get()->Next(-len);
      v->Name = v->Name + std::string("(") +
                (len > 0 ? std::string("+") : std::string("")) +
                std::to_string(len) + std::string(")");
    }
  }

  if (adjustLeadLags) {
    vs = Variables<double>(list0);
    // update
    for (auto j = 0; j < (int)list0.size(); j++) {
      bool hasMissing = false;
      auto range = vs.GetRange(j, hasMissing);
      info(j, 0) = std::get<0>(range) + 1;
      info(j, 1) = std::get<1>(range) + 1;
      info(j, 4) = 0; // it is adjusted now
    }
  }

  auto res = NumericMatrix(vs.NumObs, vs.Names.size(), vs.Data.begin());
  colnames(res) = wrap(vs.Names);
  auto f0 = vs.StartFrequency.get()->Clone();
  auto f = f0.get();
  std::vector<std::string> rns;
  for (int i = 0; i < vs.NumObs; i++) {
    rns.push_back(f->ToString());
    f->Next(1);
  }
  rownames(res) = wrap(rns);

  auto L =
      List::create(_["data"] = res, _["info"] = info,
                   _["startFrequency"] = vs.StartFrequency.get()->ToString(),
                   _["startClass"] = vs.StartFrequency.get()->ToClassString());
  return L;
}
