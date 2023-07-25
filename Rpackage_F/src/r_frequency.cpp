#include "r_tsdata.h"

using namespace Rcpp;
using namespace ldt;

// [[Rcpp::export(.F_cross_section)]]
SEXP F_CrossSection(int position) {
  List F = List::create(_["class"] = (int)FrequencyClass::kCrossSection,
                        _["position"] = position);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_yearly)]]
SEXP F_Yearly(int year) {
  List F =
      List::create(_["class"] = (int)FrequencyClass::kYearly, _["year"] = year);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_quarterly)]]
SEXP F_Quarterly(int year, int quarter) {
  List F = List::create(_["class"] = (int)FrequencyClass::kQuarterly,
                        _["year"] = year, _["quarter"] = quarter);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_monthly)]]
SEXP F_Monthly(int year, int month) {
  List F = List::create(_["class"] = (int)FrequencyClass::kMonthly,
                        _["year"] = year, _["month"] = month);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_multi_yearly)]]
SEXP F_MultiYearly(int year, int z) {
  List F = List::create(_["class"] = (int)FrequencyClass::kMultiYear,
                        _["year"] = year, _["z"] = z);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_x_times_a_year)]]
SEXP F_XTimesAYear(int year, int x, int position) {
  List F = List::create(_["class"] = (int)FrequencyClass::kXTimesAYear,
                        _["year"] = year, _["x"] = x, _["position"] = position);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_x_times_z_years)]]
SEXP F_XTimesZYears(int year, int x, int z, int position) {
  List F = List::create(_["class"] = (int)FrequencyClass::kXTimesZYears,
                        _["year"] = year, _["x"] = x, _["z"] = z,
                        _["position"] = position);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_weekly)]]
SEXP F_Weekly(int year, int month, int day)

{
  List F = List::create(_["class"] = (int)FrequencyClass::kWeekly,
                        _["year"] = year, _["month"] = month, _["day"] = day);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_multi_weekly)]]
SEXP F_MultiWeekly(int year, int month, int day, int k)

{
  List F = List::create(_["class"] = (int)FrequencyClass::kMultiWeekly,
                        _["year"] = year, _["month"] = month, _["day"] = day,
                        _["k"] = k);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_daily)]]
SEXP F_Daily(int year, int month, int day)

{
  List F = List::create(_["class"] = (int)FrequencyClass::kDaily,
                        _["year"] = year, _["month"] = month, _["day"] = day);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_multi_daily)]]
SEXP F_MultiDaily(int year, int month, int day, int k)

{
  List F = List::create(_["class"] = (int)FrequencyClass::kMultiDaily,
                        _["year"] = year, _["month"] = month, _["day"] = day,
                        _["k"] = k);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_daily_in_week)]]
SEXP F_DailyInWeek(int year, int month, int day, std::string weekStart,
                   std::string weekEnd, bool forward) {
  auto range = DayOfWeekRange::Parse(weekStart + std::string("-") + weekEnd);
  auto freq = FrequencyWeekBased(boost::gregorian::date(year, month, day),
                                 false, &range, forward);
  List F = List::create(
      _["class"] = (int)FrequencyClass::kDailyInWeek,
      _["year"] = (int)freq.mDay.year(), _["month"] = (int)freq.mDay.month(),
      _["day"] = (int)freq.mDay.day(), _["weekStart"] = (int)range.mStart,
      _["weekEnd"] = (int)range.mEnd);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_list_string)]]
SEXP F_ListString(std::vector<std::string> items, std::string value) {
  List F = List::create(_["class"] = (int)FrequencyClass::kListString,
                        _["items"] = items, _["value"] = value);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_list_date)]]
SEXP F_ListDate(std::vector<std::string> items, std::string value) {
  // auto d = boost::gregorian::date_from_iso_string(value);
  List F = List::create(_["class"] = (int)FrequencyClass::kListDate,
                        _["items"] = items, _["value"] = value);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_hourly)]]
SEXP F_Hourly(SEXP day, int hour)

{
  List F = List::create(_["class"] = (int)FrequencyClass::kHourly,
                        _["day"] = day, _["hour"] = hour);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_minutely)]]
SEXP F_Minutely(SEXP day, int minute) {
  List F = List::create(_["class"] = (int)FrequencyClass::kMinutely,
                        _["day"] = day, _["minute"] = minute);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_secondly)]]
SEXP F_Secondly(SEXP day, int second) {
  List F = List::create(_["class"] = (int)FrequencyClass::kSecondly,
                        _["day"] = day, _["second"] = second);
  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

// [[Rcpp::export(.F_x_times_a_day)]]
SEXP F_XTimesADay(SEXP day, int x, int position)

{
  List F = List::create(_["class"] = (int)FrequencyClass::kXTimesADay,
                        _["day"] = day, _["x"] = x, _["position"] = position);

  F.attr("class") = std::vector<std::string>({"ldtf", "list"});
  return F;
}

void getCh(CharacterVector f, std::vector<std::string> &result) {
  result.resize(f.size());
  for (int i = 0; i < f.size(); i++)
    result.at(i) = std::string(f[i]);
}

std::unique_ptr<FrequencyWeekBased> GetFreqFromSEXP_week(List f) {
  auto fclass = (FrequencyClass)as<int>(f["class"]);

  switch (fclass) {

  case FrequencyClass::kWeekly:
    return FrequencyWeekBased::Weekly(boost::gregorian::date(
        as<int>(f["year"]), as<int>(f["month"]), as<int>(f["day"])));
  case FrequencyClass::kMultiWeekly:
    return FrequencyWeekBased::MultiWeekly(
        boost::gregorian::date(as<int>(f["year"]), as<int>(f["month"]),
                               as<int>(f["day"])),
        as<int>(f["k"]));
  case FrequencyClass::kDaily:
    return FrequencyWeekBased::Daily(boost::gregorian::date(
        as<int>(f["year"]), as<int>(f["month"]), as<int>(f["day"])));
  case FrequencyClass::kMultiDaily:
    return FrequencyWeekBased::MultiDaily(
        boost::gregorian::date(as<int>(f["year"]), as<int>(f["month"]),
                               as<int>(f["day"])),
        as<int>(f["k"]));
  case FrequencyClass::kDailyInWeek: {
    auto range = DayOfWeekRange((DayOfWeek)as<int>(f["weekStart"]),
                                (DayOfWeek)as<int>(f["weekEnd"]));
    return FrequencyWeekBased::DailyInWeek(
        boost::gregorian::date(as<int>(f["year"]), as<int>(f["month"]),
                               as<int>(f["day"])),
        range, true);
  }
  default:
    throw LdtException(ErrorType::kLogic, "R-frequency",
                   "use this class for week-based frequency");
  }
}

std::unique_ptr<Frequency>
GetFreqFromSEXP(SEXP value, std::vector<std::string> &listItems,
                std::vector<boost::gregorian::date> &listItemsDate) {
  auto f = as<List>(value);
  auto fclass = (FrequencyClass)as<int>(f["class"]);

  switch (fclass) {
  case FrequencyClass::kCrossSection:
    return std::unique_ptr<FrequencyCrossSection>(
        new FrequencyCrossSection(as<int>(f["position"])));

  case FrequencyClass::kYearly:
    return FrequencyYearBased::Yearly(as<int>(f["year"]));
  case FrequencyClass::kQuarterly:
    return FrequencyYearBased::Quarterly(as<int>(f["year"]),
                                         as<int>(f["quarter"]));
  case FrequencyClass::kMonthly:
    return FrequencyYearBased::Monthly(as<int>(f["year"]), as<int>(f["month"]));
  case FrequencyClass::kMultiYear:
    return FrequencyYearBased::MultiYearly(as<int>(f["year"]), as<int>(f["z"]));
  case FrequencyClass::kXTimesAYear:
    return FrequencyYearBased::XTimesAYear(as<int>(f["year"]), as<int>(f["x"]),
                                           as<int>(f["position"]));
  case FrequencyClass::kXTimesZYears:
    return FrequencyYearBased::XTimesZYear(as<int>(f["year"]), as<int>(f["x"]),
                                           as<int>(f["position"]),
                                           as<int>(f["z"]));

  case FrequencyClass::kWeekly:
  case FrequencyClass::kMultiWeekly:
  case FrequencyClass::kDaily:
  case FrequencyClass::kMultiDaily:
  case FrequencyClass::kDailyInWeek:
    return GetFreqFromSEXP_week(f);

  case FrequencyClass::kListString: {
    getCh(f["items"], listItems);
    return std::unique_ptr<FrequencyList<std::string>>(
        new FrequencyList<std::string>(as<std::string>(f["value"]),
                                       &listItems));
  }
  case FrequencyClass::kListDate: {
    getCh(f["items"], listItems);
    for (auto const &d : listItems)
      listItemsDate.push_back(boost::gregorian::date_from_iso_string(d));
    return std::unique_ptr<FrequencyList<boost::gregorian::date>>(
        new FrequencyList<boost::gregorian::date>(
            boost::gregorian::date_from_iso_string(as<std::string>(f["value"])),
            &listItemsDate));
  }

  case FrequencyClass::kHourly: {
    auto day = GetFreqFromSEXP_week(as<List>(f["day"]));
    return FrequencyDayBased::Hourly(*day.get(), as<int>(f["hour"]));
  }
  case FrequencyClass::kMinutely: {
    auto day = GetFreqFromSEXP_week(f["day"]);
    return FrequencyDayBased::Minutely(*day.get(), as<int>(f["minute"]));
  }
  case FrequencyClass::kSecondly: {
    auto day = GetFreqFromSEXP_week(f["day"]);
    return FrequencyDayBased::Secondly(*day.get(), as<int>(f["second"]));
  }
  case FrequencyClass::kXTimesADay: {
    auto day = GetFreqFromSEXP_week(f["day"]);
    return FrequencyDayBased::XTimesADay(*day.get(), as<int>(f["x"]),
                                         as<int>(f["position"]));
  }

  default:
    throw LdtException(ErrorType::kLogic, "R-frequency",
                   "not implemeted for this type of frequency");
  }
}

// [[Rcpp::export(.ToString_F)]]
std::string ToString_F(SEXP value)

{
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  auto F = GetFreqFromSEXP(value, listItems, listItemsDate);
  return F.get()->ToString();
}

// [[Rcpp::export(.ToClassString_F)]]
std::string ToClassString_F(SEXP value) {
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  auto F = GetFreqFromSEXP(value, listItems, listItemsDate);
  return F.get()->ToClassString(true);
}

// [[Rcpp::export(.ToString_F0)]]
List ToString_F0(SEXP value) {
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  auto F = GetFreqFromSEXP(value, listItems, listItemsDate);
  List I = List::create(_["value"] = F.get()->ToString(),
                        _["class"] = F.get()->ToClassString(),
                        _["classType"] = ToString(F.get()->mClass));
  return I;
}

SEXP To_SEXP_week(FrequencyClass fClass, Frequency *F) {

  auto f = dynamic_cast<FrequencyWeekBased const &>(*F);
  switch (fClass) {

  case FrequencyClass::kWeekly: {
    return F_Weekly(f.mDay.year(), f.mDay.month(), f.mDay.day());
  }
  case FrequencyClass::kMultiWeekly: {
    return F_MultiWeekly(f.mDay.year(), f.mDay.month(), f.mDay.day(), f.mMulti);
  }
  case FrequencyClass::kDaily: {
    return F_Daily(f.mDay.year(), f.mDay.month(), f.mDay.day());
  }
  case FrequencyClass::kMultiDaily: {
    return F_MultiDaily(f.mDay.year(), f.mDay.month(), f.mDay.day(), f.mMulti);
  }
  case FrequencyClass::kDailyInWeek: {
    return F_DailyInWeek(f.mDay.year(), f.mDay.month(), f.mDay.day(),
                         ToString(f.mRange.mStart), ToString(f.mRange.mEnd),
                         true);
  }

  default:
    throw LdtException(ErrorType::kLogic, "R-frequency",
                   "invalid frequency class. week-based frequency is expected");
  }
}

SEXP To_SEXP(Frequency *F, std::vector<std::string> &listItems,
             std::vector<boost::gregorian::date> &listItemsDate) {
  FrequencyClass fClass = F->mClass;

  switch (fClass) {
  case FrequencyClass::kCrossSection: {
    auto f = dynamic_cast<FrequencyCrossSection const &>(*F);
    return F_CrossSection(f.mPosition);
  }

  case FrequencyClass::kYearly: {
    auto f = dynamic_cast<FrequencyYearBased const &>(*F);
    return F_Yearly(f.mYear);
  }
  case FrequencyClass::kQuarterly: {
    auto f = dynamic_cast<FrequencyYearBased const &>(*F);
    return F_Quarterly(f.mYear, f.mPosition);
  }
  case FrequencyClass::kMonthly: {
    auto f = dynamic_cast<FrequencyYearBased const &>(*F);
    return F_Monthly(f.mYear, f.mPosition);
  }
  case FrequencyClass::kMultiYear: {
    auto f = dynamic_cast<FrequencyYearBased const &>(*F);
    return F_MultiYearly(f.mYear, f.mYearMulti);
  }
  case FrequencyClass::kXTimesAYear: {
    auto f = dynamic_cast<FrequencyYearBased const &>(*F);
    return F_XTimesAYear(f.mYear, f.mPartitionCount, f.mPosition);
  }
  case FrequencyClass::kXTimesZYears: {
    auto f = dynamic_cast<FrequencyYearBased const &>(*F);
    return F_XTimesZYears(f.mYear, f.mPartitionCount, f.mYearMulti,
                          f.mPosition);
  }

  case FrequencyClass::kWeekly:
  case FrequencyClass::kMultiWeekly:
  case FrequencyClass::kDaily:
  case FrequencyClass::kMultiDaily:
  case FrequencyClass::kDailyInWeek:
    return To_SEXP_week(fClass, F);

  case FrequencyClass::kHourly: {
    auto f = dynamic_cast<FrequencyDayBased const &>(*F);
    auto fw = To_SEXP_week(f.mDay.mClass, &f.mDay);
    return F_Hourly(fw, f.mPosition);
  }
  case FrequencyClass::kMinutely: {
    auto f = dynamic_cast<FrequencyDayBased const &>(*F);
    auto fw = To_SEXP_week(f.mDay.mClass, &f.mDay);
    return F_Minutely(fw, f.mPosition);
  }
  case FrequencyClass::kSecondly: {
    auto f = dynamic_cast<FrequencyDayBased const &>(*F);
    auto fw = To_SEXP_week(f.mDay.mClass, &f.mDay);
    return F_Secondly(fw, f.mPosition);
  }
  case FrequencyClass::kXTimesADay: {
    auto f = dynamic_cast<FrequencyDayBased const &>(*F);
    auto fw = To_SEXP_week(f.mDay.mClass, &f.mDay);
    return F_XTimesADay(fw, f.mPartitionCount, f.mPosition);
  }

  case FrequencyClass::kListString: {
    auto f = dynamic_cast<FrequencyList<std::string> const &>(*F);
    auto dd = F_ListString(listItems, f.mValue);
    return dd;
  }

  case FrequencyClass::kListDate: {
    if (listItems.size() == 0) {
      for (auto &d : listItemsDate)
        listItems.push_back(boost::gregorian::to_iso_string(d));
    }
    auto f = dynamic_cast<FrequencyList<boost::gregorian::date> const &>(*F);
    auto dd = F_ListDate(listItems, boost::gregorian::to_iso_string(f.mValue));
    return dd;
  }

  default:
    throw LdtException(ErrorType::kLogic, "R-frequency",
                   "not implemeted for this type of frequency");
  }
}

// [[Rcpp::export(.Parse_F)]]
SEXP Parse_F(std::string str, std::string classStr) {

  FrequencyClass fClass;
  auto F = Frequency::Parse(str, classStr, fClass);

  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  if (F.get()->mClass == FrequencyClass::kListString) {

    auto f0 = FrequencyList<std::string>("", nullptr);
    FrequencyList<std::string>::Parse0(str, classStr, fClass, f0, &listItems);
    return To_SEXP(&f0, listItems, listItemsDate);

  } else if (F.get()->mClass == FrequencyClass::kListDate) {

    auto f0 = FrequencyList<boost::gregorian::date>(boost::gregorian::date(),
                                                    nullptr);
    FrequencyList<boost::gregorian::date>::Parse0(str, classStr, fClass, f0,
                                                  &listItemsDate);
    return To_SEXP(&f0, listItems, listItemsDate);

  } else {
    return To_SEXP(F.get(), listItems, listItemsDate);
  }
}

// [[Rcpp::export(.Sequence_F0)]]
std::vector<std::string> Sequence_F0(SEXP start, int length, int by)

{
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  auto F = GetFreqFromSEXP(start, listItems, listItemsDate);
  auto f = F.get();

  std::vector<std::string> list;
  for (Ti i = 0; i < length; i++) {
    list.push_back(f->ToString());
    f->Next(by);
  }

  return list;
}

// [[Rcpp::export(.Sequence_F)]]
std::vector<std::string> Sequence_F(SEXP from, SEXP to, int by)

{
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  auto F = GetFreqFromSEXP(from, listItems, listItemsDate);
  auto f = F.get();
  auto T = GetFreqFromSEXP(to, listItems, listItemsDate);
  auto t = T.get();

  auto fromNewer = f->IsNewerThan(*t);
  if ((by > 0 && fromNewer) || (by < 0 && fromNewer == false))
    by = -by;

  std::vector<std::string> list;
  while (true) {
    if ((fromNewer == false && f->IsNewerThan(*t)) ||
        (fromNewer && t->IsNewerThan(*f)))
      break;
    list.push_back(f->ToString());
    f->Next(by);
  }

  return list;
}

// [[Rcpp::export(.F_GetClass)]]
int F_GetClass(std::string name) {
  return (int)FromString_FrequencyClass(name.c_str());
}

// [[Rcpp::export(.F_Next)]]
SEXP F_Next(SEXP freq, int steps) {
  std::vector<std::string> listItems;
  std::vector<boost::gregorian::date> listItemsDate;
  auto F = GetFreqFromSEXP(freq, listItems, listItemsDate);
  auto f = F.get();

  f->Next(steps);

  return To_SEXP(f, listItems, listItemsDate);
}

// [[Rcpp::export(.F_Minus)]]
int F_Minus(SEXP freq1, SEXP freq2) {
  std::vector<std::string> listItems1;
  std::vector<boost::gregorian::date> listItemsDate1;
  auto F1 = GetFreqFromSEXP(freq1, listItems1, listItemsDate1);
  auto f1 = F1.get();

  std::vector<std::string> listItems2;
  std::vector<boost::gregorian::date> listItemsDate2;
  auto F2 = GetFreqFromSEXP(freq2, listItems2, listItemsDate2);
  auto f2 = F2.get();
  return f1->Minus(*f2);
}
