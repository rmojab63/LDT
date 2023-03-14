#include "r_ldt.h"

using namespace Rcpp;
using namespace ldt;

// clang-format off

//' Creates a Variable
//'
//' @param data Data of the variable
//' @param name Name of the variable
//' @param startFrequency Frequency of the first data-point. It is an \code{ldtf} object. See \code{F_?} functions.
//' @param fields Named list of any other fields
//'
//' @return An object of class \code{ldtv}.
//' @export
//' @examples
//' v1 = ldt::Variable(c(1,2,3,2,3,4,5),"V1",F_Monthly(2022,12),
//'      list(c("key1","value1"), c("key2", "value2")))
// [[Rcpp::export]]
List Variable(std::vector<double> data, std::string name, SEXP startFrequency,
              List fields)
// clang-format on
{
  List V =
      List::create(_["data"] = data, _["name"] = name,
                   _["startFrequency"] = startFrequency, _["fields"] = fields);
  V.attr("class") = std::vector<std::string>({"ldtv", "list"});
  return V;
}

void UpdateVariableFromSEXP(Rcpp::List w, ldt::Variable<double>& variable,
                         std::vector<std::string>& listItems,
                         std::vector<boost::gregorian::date>& listItemsDate) {

  variable.Name = as<std::string>(w["name"]);

  try {
    variable.StartFrequency =
        std::move(GetFreqFromSEXP(w["startFrequency"], listItems, listItemsDate));
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
    List F = w["fields"];
    for (int j = 0; j < F.length(); j++) {
      CharacterVector Fj = as<CharacterVector>(F[j]);
      if (Fj.length() < 2)
        throw std::logic_error("Expected a 'key' and a 'value'.");
      variable.Fields.insert({as<std::string>(Fj[0]), as<std::string>(Fj[1])});
    }
  } catch (...) {
    throw std::logic_error(
        "Invalid fields. It should be a 'List' of 'CharacterVector's. The "
        "vector must have 'Key-Value' pairs.");
  }
}

// clang-format off

//' Converts a Variable to String
//'
//' @param w The variable
//'
//' @return String representation of the variable in compact form
//' @export
// [[Rcpp::export]]
std::string VariableToString(List w)
// clang-format on
{
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  ldt::Variable v;
  UpdateVariableFromSEXP(w,v,listItems,listItemsDate);
  return v.ToString();
}

// clang-format off

//' Binds a List of Variables
//'
//' @param varList A list of variables ((i.e., \code{ldtv} objects)) with similar frequency class
//' @param interpolate If \code{TRUE} missing observations are interpolated.
//' @param adjustLeadLags If \code{TRUE} leads and lags are adjusted with respect to the first variable.
//' @param numEndo (integer) If \code{adjustLeadLags} is \code{TRUE}, this must be
//' the number of endogenous variables. The rest is exogenous.
//' @param horizon (integer) If \code{adjustLeadLags} is \code{TRUE} and there is exogenous variables,
//' this determines the required length of out-of-sample data. It creates lag of exogenous variables
//' or omits \code{NaN}s to make data available.
//'
//' @return (list) results
//' \item{hasMissing}{(logical) If \code{TRUE}, missing observation exists in the original data.}
//' \item{hasMissingIndices}{(numeric matrix) Indices of variables with missing data (first column) and the number of missing data points (second column).}
//' \item{hasLags}{(numeric matrix) Indices of variables with lags (first column) and the value of lag (second column).}
//' \item{hasLeads}{Similar to \code{hasLags} but for leads.}
//' \item{data}{(numeric matrix) Final data after the requested fixes. It is a matrix with variables in the columns and frequencies as the row names.}
//' \item{ranges}{(numeric matrix) Start and end of the final data.}
//' \item{countNansSet}{(integer) Number of generated \code{NAN}s for adjusting out-of-sample data.}
//'
//' @export
//' @examples
//' v1 = ldt::Variable(c(1,2,3,2,3,4,5),"V1",F_Monthly(2022,12), list())
//' v2 = ldt::Variable(c(10,20,30,20,30,40,50),"V2",F_Monthly(2022,8), list())
//' L = ldt::BindVariables(list(v1,v2))
// [[Rcpp::export]]
List BindVariables(SEXP varList, bool interpolate = false,
                            bool adjustLeadLags = false, int numExo = 0,
                            int horizon = 0)
// clang-format on
{
  auto vars = as<List>(varList);
  int n = vars.length();
  if (n == 0)
    throw std::logic_error("Empty list of variables.");

  int numEndo = n - numExo;
  if (numExo < 0 || numEndo < 0)
    throw std::logic_error("Invalid number of exogenous/exogenous variables (check 'numExo').");
  if (numExo >0 && horizon < 0)
    throw std::logic_error("Invalid out-of-sample size (check 'horizon').");

  auto info = IntegerMatrix(n, 5);
  colnames(info) = wrap(std::vector<const char*>({"Start_Index","End_Index", "Has_Missing", "Interpolated_Count", "Lags_Leads"}));

  auto listItems = std::vector<std::vector<std::string>>(n);
  auto listItemsDate = std::vector<std::vector<boost::gregorian::date>>(n);

  auto list =  std::vector<ldt::Variable<double>>(n);
  auto list0 =  std::vector<ldt::Variable<double>*>(n);
  try {
    for (int i = 0; i < n; i++) {
      if (is<List>(vars[i]) == false)
        throw std::logic_error("Invalid variable type.");
      UpdateVariableFromSEXP(as<List>(vars[i]),
                             list.at(i),listItems.at(i),listItemsDate.at(i));
      list.at(i).Trim();
      list0.at(i) = &list.at(i);

      int count = 0;
      if (interpolate){
        auto vec = ldt::Matrix(&list.at(i).Data[0],list.at(i).Data.size(),1);
        vec.InterpolateColumn(count, 0);
      }
      info(i,3) = count;
    }
  } catch (...) {
    throw std::logic_error(
        "Creating a variable failed. Make sure all items of the list are "
        "'ldtv' objects.");
  }
  auto vs = Variables<double>(list0);
  auto vs_data = ldt::Matrix<double>(&vs.Data[0],vs.NumObs, n);

  std::vector<IndexRange> ranges;

  for (int i = 0; i < n; i++) {
    bool hasMissing = false;
    IndexRange range = vs_data.GetRangeColumn(hasMissing, i);
    info(i,0) = range.StartIndex + 1;
    info(i,1) = range.EndIndex + 1;
    info(i,2) = hasMissing ? 1 : 0;
    ranges.push_back(range);

    if (adjustLeadLags){
      // we should update data of the variable
      auto a = &vs_data.Data[i*vs_data.RowsCount];
      list.at(i).Data = std::vector<double>(a,a + vs_data.RowsCount);
      list.at(i).StartFrequency = vs.StartFrequency.get()->Clone();
      list0.at(i) = &list.at(i);
      list.at(i).Trim();
    }
  }

  // lags, leads
  info(0,5) = 0;
  int lastIndex = ranges.at(0).EndIndex;
  int lastIndexHor = lastIndex + horizon;
  int len;
  for (int i = 1; i < n; i++) {
    if (i < numEndo)
      len = lastIndex - ranges.at(i).EndIndex;
    else
      len = lastIndexHor - ranges.at(i).EndIndex;
    len*=-1;
    info(i,4) = len;

    if (adjustLeadLags && len != 0){
      auto v = list0.at(i);
      v->StartFrequency.get()->Next(-len);
      v->Name = v->Name + std::string("(") +
        (len > 0 ? std::string("+") : std::string("")) +
        std::to_string(len) + std::string(")");
    }

  }

  if (adjustLeadLags){
    vs = Variables<double>(list0);
    vs_data = ldt::Matrix<double>(&vs.Data[0],vs.NumObs, n);
    // update
    for (auto j=0;j<vs_data.ColsCount;j++){
      bool hasMissing = false;
      IndexRange range = vs_data.GetRangeColumn(hasMissing, j);
      info(j,0) = range.StartIndex + 1;
      info(j,1) = range.EndIndex + 1;
      info(j,4) = 0; // it is adjusted now
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

  auto L = List::create(
    _["data"] = res,
    _["info"] = info
  );
  return L;
}
