
#include "variable.h"

using namespace ldt;

template <typename Tw> Variable<Tw>::Variable() {}

template <typename Tw>
std::unique_ptr<Frequency> Variable<Tw>::GetEndFrequency() const {
  auto end = StartFrequency.get()->Clone();
  auto temp = end.get();
  temp->Next(Data.size());
  return end;
}

template <typename Tw>
bool Variable<Tw>::IsEqualTo(Variable<Tw> &other, const Tv &epsilon) const {
  if (Name != other.Name)
    return false;
  if (Data.size() != other.Data.size())
    return false;
  for (auto i = 0; i < (Ti)Data.size(); i++)
    if (std::abs(Data.at(i) - other.Data.at(i)) > epsilon)
      return false;
  if (StartFrequency.get()->IsEqualTo(*other.StartFrequency.get()) == false)
    return false;
  if (Fields != other.Fields)
    return false;
  return true;
}

template <typename Tw>
void Variable<Tw>::Examples(std::vector<std::unique_ptr<Variable<Tw>>> &values,
                            std::vector<std::unique_ptr<Frequency>> &f_values,
                            std::vector<std::string> &listItemsString,
                            std::vector<boost::gregorian::date> &listItemsDate,
                            std::mt19937 &eng) {

  // if constexpr (std::is_floating_point<Tw>()) {
  std::normal_distribution<Tw> D(0, 1);
  /*} else if constexpr (true) {
    std::uniform_int_distribution<Ti> D((Tw)0, (Tw)1);
  }*/

  Frequency::Examples(f_values, listItemsString, listItemsDate);
  FrequencyClass fc = FrequencyClass::kCrossSection;
  Ti i = 0;
  for (auto &a : f_values) {
    i++;
    Ti length = 20;

    auto v = new Variable<Tw>();
    values.push_back(std::unique_ptr<Variable<Tw>>(v));

    v->Name = std::string("V") + std::to_string(i);

    auto b = a.get();
    if (fc == FrequencyClass::kListString) {
      auto d = dynamic_cast<FrequencyList<std::string> &>(*b);
      d.pItems = &listItemsString;
      length = listItemsString.size();
    } else if (fc == FrequencyClass::kListDate) {
      auto d = dynamic_cast<FrequencyList<boost::gregorian::date> &>(*b);
      d.pItems = &listItemsDate;
      length = listItemsDate.size();
    }

    v->Data = std::vector<Tw>(length);
    for (Ti j = 0; j < length; j++)
      v->Data.at(j) = D(eng);

    v->StartFrequency = std::move(a);

    v->Fields.insert(
        {std::string("Key1"), std::string("V") + std::to_string(i)});
    if (i % 2 == 0)
      v->Fields.insert(
          {std::string("Key2"), std::string("W") + std::to_string(i)});
  }
}

template <typename Tw> std::string Variable<Tw>::ToString() const {

  std::ostringstream ss;
  ss << Name;
  ss << '\t';

  auto sf = StartFrequency.get();
  ss << (sf == nullptr ? std::string("NA") : sf->ToClassString(true));
  ss << '\t';

  ss << (sf == nullptr ? std::string("NA") : sf->ToString());
  ss << '\t';

  ss << std::fixed << std::setprecision(16);
  Ti c = Data.size();
  Ti i = 0;
  for (auto const &d : Data) {
    i++;
    ss << d;
    if (i < c)
      ss << ";";
  }
  ss << '\t';

  c = Fields.size();
  i = 0;
  for (auto const &x : Fields) {
    i++;
    ss << x.first;
    ss << ';'; // assuming that 'key' does not contain it
    ss << x.second;
    if (i < c)
      ss << '\t';
  }

  return ss.str();
}

template <typename Tw>
void Variable<Tw>::Parse(const std::string &str, Variable<Tw> &result,
                         std::vector<std::string> &listItemsString,
                         std::vector<boost::gregorian::date> &listItemsDate) {
  try {

    auto parts = std::vector<std::string>();
    Split(str, "\t", parts);
    if (parts.size() < 5)
      throw std::logic_error("At least 4 tab-separated items is expected.");

    result.Name = parts.at(0);

    // frequency
    FrequencyClass fClass;
    auto frequency = Frequency::Parse(parts.at(2), parts.at(1), fClass);

    if (fClass == FrequencyClass::kListString)
      frequency = FrequencyList<std::string>::ParseList(
          parts.at(2), parts.at(1), fClass, listItemsString);
    else if (fClass == FrequencyClass::kListDate)
      frequency = FrequencyList<boost::gregorian::date>::ParseList(
          parts.at(2), parts.at(1), fClass, listItemsDate);
    result.StartFrequency = std::move(frequency);

    // data
    result.Data.clear();
    std::vector<std::string> parts_d;
    Split(parts.at(3), std::string(";"), parts_d);
    result.Data.reserve(parts_d.size());
    for (auto const &d : parts_d) {
      if constexpr (std::is_same<Tw, double>()) {
        result.Data.push_back(std::stod(d));
      } else if constexpr (std::is_same<Tw, float>()) {
        result.Data.push_back(std::stof(d));
      } else if constexpr (std::is_same<Tw, Ti>()) {
        result.Data.push_back(std::stoi(d));
      } else if constexpr (std::is_same<Tw, long long>()) {
        result.Data.push_back(std::stoll(d));
      } else if constexpr (true) {
        throw std::logic_error(
            "Conversion of the variable's data-type is not implemented.");
      }
    }

    // Fields
    result.Fields.clear();
    for (Ti i = 4; i < (Ti)parts.size(); i++) {
      auto j = parts.at(i).find(std::string(";"));
      if (j < 0)
        throw std::logic_error("Invalid field: Key-Value separator is missing");
      result.Fields.insert(
          {parts.at(i).substr(0, j), parts.at(i).substr(j + 1)});
    }

  } catch (...) {
    Rethrow("Invalid format in parsing 'Variable'.");
  }
}

template <typename Tw> IndexRange Variable<Tw>::GetRange(bool &hasMissing) {

  auto range = Array<Tw>::GetRange(&Data[0], (Ti)Data.size(), hasMissing);
  return range;
}

template <typename Tw> IndexRange Variable<Tw>::Trim() {
  bool hasMissing = false;
  auto range = GetRange(hasMissing);

  if (range.IsNotValid())
    return range;

  auto count = range.EndIndex - range.StartIndex + 1;
  if (count != (Ti)Data.size()) {
    Data = {Data.begin() + range.StartIndex, Data.begin() + range.EndIndex + 1};
    StartFrequency.get()->Next(range.StartIndex);
  }
  return range;
}

template <typename Tw> IndexRange Variable<Tw>::Interpolate(Ti &count) {

  if constexpr (std::numeric_limits<Tw>::has_quiet_NaN) {

    auto range = Array<Tw>::Interpolate(&Data[0], (Ti)Data.size(), count);
    return range;

  } else if constexpr (true)
    throw std::logic_error("invalid operation"); // there is no NAN
}

template <typename Tw>
void Variable<Tw>::ConvertTo_Daily(
    Variable &result,
    const std::function<Tv(const std::vector<Tv> &)> *aggregateFunc) const {

  auto mclass = StartFrequency.get()->mClass;

  if (mclass == FrequencyClass::kListDate) {
    auto startList =
        dynamic_cast<FrequencyList<boost::gregorian::date> const &>(
            *StartFrequency.get());
    result.Data.clear();

    auto dates = startList.pItems;
    auto minmax = std::minmax_element(dates->begin(), dates->end());
    auto min_date = *minmax.first;
    auto duration = (*minmax.second) - min_date;

    for (Ti i = 0; i <= duration.days(); ++i) {
      auto date_to_find = min_date + boost::gregorian::days(i);
      auto it = std::find(dates->begin(), dates->end(), date_to_find);
      if (it != dates->end()) {
        Ti index = std::distance(dates->begin(), it);
        result.Data.push_back(Data.at(index));
      } else { // not found
        result.Data.push_back(NAN);
      }
    }
    result.Name = Name;
    result.StartFrequency = std::move(FrequencyWeekBased::Daily(min_date));
    result.Fields.insert(std::pair("conversion", "from date-list"));

  } else if (mclass == FrequencyClass::kDailyInWeek) {
    auto startF =
        dynamic_cast<FrequencyWeekBased const &>(*StartFrequency.get());
    result.Data.clear();

    auto range = startF.mRange;
    auto startDay = startF.mDay;

    auto weekSize = range.GetLength();
    auto weekEndCount = 7 - weekSize;
    auto numWeeks = std::ceil((Tv)Data.size() / (Tv)weekSize);

    Ti end = 0, pend = 0;
    for (Ti i = 0; i < numWeeks; ++i) {
      end = std::min(end + weekSize, static_cast<Ti>(Data.size()));
      for (Ti p = pend; p < end; p++)
        result.Data.push_back(Data.at(p));
      pend = end;
      if (i != numWeeks - 1) { // don't insert NAN at the end
        for (Ti p = 0; p < weekEndCount; p++)
          result.Data.push_back(NAN);
      }
    }
    result.Name = Name;
    result.StartFrequency = std::move(FrequencyWeekBased::Daily(startDay));
    result.Fields.insert(std::pair("conversion", "from date-list"));

  } else
    throw std::logic_error("Direct conversion from current type of frequency "
                           "to 'Daily' frequency is not "
                           "supported (or not implemented).");
}

template <typename Tw>
void Variable<Tw>::ConvertTo_MultiDaily(
    Variable &result, Ti k,
    const std::function<Tv(const std::vector<Tv> &)> *aggregateFunc,
    bool fromEnd) const {

  if (StartFrequency.get()->mClass == FrequencyClass::kDaily) {

    if (!aggregateFunc)
      throw std::logic_error("Aggregate function is missing.");
    auto aggF = *aggregateFunc;

    std::vector<std::vector<Tv>> partitions;
    Array<Tw>::PartitionEqual(Data, partitions, k, fromEnd);
    std::vector<Tv> newdata;
    for (Ti i = 0; i < (Ti)partitions.size(); i++) {
      newdata.push_back(aggF(partitions.at(i)));
    }

    result.Data = newdata;
    result.Name = Name;

    auto start =
        dynamic_cast<FrequencyWeekBased const &>(*StartFrequency.get());
    result.StartFrequency =
        std::move(FrequencyWeekBased::MultiDaily(start.mDay, k));

  } else {
    throw std::logic_error(
        "Direct conversion from current type of frequency to "
        "'Multi-Day' frequency is not "
        "supported (or not implemented).");
  }
}

template <typename Tw>
void Variable<Tw>::ConvertTo_Weekly(
    Variable &result, DayOfWeek firstDayOfWeek,
    const std::function<Tv(const std::vector<Tv> &)> *aggregateFunc) const {

  if (StartFrequency.get()->mClass == FrequencyClass::kDaily) {

    if (!aggregateFunc)
      throw std::logic_error("Aggregate function is missing.");
    auto aggF = *aggregateFunc;

    auto start =
        dynamic_cast<FrequencyWeekBased const &>(*StartFrequency.get());
    auto wd = start.mDay.day_of_week();
    std::string wd_ss = wd.as_short_string();
    boost::algorithm::to_lower(wd_ss);
    auto pre_end = FromString_DayOfWeek(wd_ss.c_str());
    auto diff = DayOfWeekRange::Distance(firstDayOfWeek, pre_end, true);

    std::vector<std::vector<Tv>> weeks;
    std::vector<Tv> newdata;

    if (diff == 0) {

      Array<Tw>::PartitionEqual(Data, weeks, 7, false);
      for (Ti i = 0; i < (Ti)weeks.size(); i++)
        newdata.push_back(aggF(weeks.at(i)));

    } else { // insert NAN and make it the first day of the week

      auto copydata = Data;
      for (Ti i = 0; i < diff; i++)
        copydata.insert(copydata.begin(), NAN);
      Array<Tw>::PartitionEqual(copydata, weeks, 7, false);
      for (Ti i = 0; i < (Ti)weeks.size(); i++)
        newdata.push_back(aggF(weeks.at(i)));
    }
    result.Data = newdata;
    result.Name = Name;
    boost::gregorian::date_duration dd(diff);
    result.StartFrequency =
        std::move(FrequencyWeekBased::Weekly(start.mDay - dd));

  } else {
    throw std::logic_error(
        "Direct conversion from current type of frequency to "
        "'Weekly' frequency is not "
        "supported (or not implemented).");
  }
}

static Ti get_part(boost::gregorian::date date, Ti x) {
  auto day_of_year = date.day_of_year();
  int partition;
  if (boost::gregorian::gregorian_calendar::is_leap_year(date.year()))
    partition = ((day_of_year - 1) / 366.0) * x + 1;
  else
    partition = ((day_of_year - 1) / 365.0) * x + 1;
  return partition;
}

template <typename Tw>
void Variable<Tw>::ConvertTo_XxYear(
    Variable &result, Ti x,
    const std::function<Tv(const std::vector<Tv> &)> *aggregateFunc) const {

  /*if (mclass == FrequencyClass::kListDate) { // It can be more efficient
    } else*/
  if (StartFrequency.get()->mClass == FrequencyClass::kDaily) {
    auto start =
        dynamic_cast<FrequencyWeekBased const &>(*StartFrequency.get());

    if (!aggregateFunc)
      throw std::logic_error("Aggregate function is missing.");
    auto aggF = *aggregateFunc;

    result.Data.clear();
    std::vector<Tv> partdata;

    auto current_part = get_part(start.mDay, x);
    auto part = current_part;

    for (Ti i = 0; i < (Ti)Data.size(); i++) {

      auto day = start.mDay + boost::gregorian::date_duration(i);
      part = get_part(day, x);

      if (part != current_part) {
        result.Data.push_back(aggF(partdata));
        partdata.clear();
      }
      partdata.push_back(Data.at(i));
      current_part = part;
    }
    if (partdata.size() > 0)
      result.Data.push_back(aggF(partdata));

    result.Name = Name;
    result.StartFrequency = std::move(FrequencyYearBased::XTimesAYear(
        start.mDay.year(), x, get_part(start.mDay, x)));

  } else {
    throw std::logic_error(
        "Direct conversion from current type of frequency to "
        "'x times a year' frequency is not "
        "supported (or not implemented).");
  }
}

template class ldt::Variable<Tv>;
// template class ldt::Variable < Ti>;
