
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
      } else if constexpr (std::is_same<Tw, int>()) {
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

template <typename Tw>
std::tuple<Ti, Ti> Variable<Tw>::GetRange(bool &hasMissing) {

  if constexpr (std::numeric_limits<Tw>::has_quiet_NaN) {
    Ti start, end;
    hasMissing = false;
    for (start = 0; start < (Ti)Data.size(); start++)
      if (std::isnan(Data.at(start)) == false)
        break;

    for (end = Data.size(); end > 0; end--)
      if (std::isnan(Data.at(end - 1)) == false) {
        end--;
        break;
      }

    if (start > end)
      return std::tuple(start, end);

    for (Ti i = start; i <= end; i++)
      if (std::isnan(Data.at(i))) {
        hasMissing = true;
        break;
      }

    return std::tuple(start, end);

  } else if constexpr (true) {
    throw std::logic_error("invalid operation"); // there is no NAN
  }
}

template <typename Tw> std::tuple<Ti, Ti> Variable<Tw>::Trim() {
  bool hasMissing = false;
  auto range = GetRange(hasMissing);
  auto start = std::get<0>(range);
  auto end = std::get<1>(range);
  if (start > end)
    return range;

  auto count = end - start + 1;
  if (count != (Ti)Data.size()) {
    Data = {Data.begin() + start, Data.begin() + end + 1};
    StartFrequency.get()->Next(start);
  }
  return range;
}

template <typename Tw> std::tuple<Ti, Ti> Variable<Tw>::Interpolate(Ti &count) {

  if constexpr (std::numeric_limits<Tw>::has_quiet_NaN) {

    bool hasMissing = false;
    auto range = GetRange(hasMissing);
    count = 0;
    if (hasMissing) {
      bool inMissing = false;
      Tw first = NAN, last = NAN;
      Ti length = 1;
      auto start = std::get<0>(range);
      auto end = std::get<1>(range);

      Tw *col = &Data[0];
      for (Ti i = start; i <= end; i++) {
        auto isNaN = std::isnan(col[i]);
        if (isNaN)
          length++;
        if (isNaN == false && inMissing) {
          last = col[i];
          // calculate and set
          Tw step = (last - first) / length;
          for (int p = 1; p < length; p++) {
            col[i - p] = col[i] - p * step;
            count++;
          }
          length = 1;
          inMissing = false;
        }
        if (isNaN && inMissing == false) {
          first = col[i - 1];
          inMissing = true;
        }
      }
    }
    return range;
  } else if constexpr (true)
    throw std::logic_error("invalid operation"); // there is no NAN
}

template <typename Tw> void Variable<Tw>::ConvertTo_Daily(Variable &result) {

  if (StartFrequency.get()->mClass == FrequencyClass::kListDate) {
    auto startList =
        dynamic_cast<FrequencyList<boost::gregorian::date> const &>(
            *StartFrequency.get());
    result.Data.clear();

    auto dates = startList.pItems;
    auto minmax = std::minmax_element(dates->begin(), dates->end());
    auto min_date = *minmax.first;
    auto duration = (*minmax.second) - min_date;

    for (int i = 0; i <= duration.days(); ++i) {
      auto date_to_find = min_date + boost::gregorian::days(i);
      auto it = std::find(dates->begin(), dates->end(), date_to_find);
      if (it != dates->end()) {
        int index = std::distance(dates->begin(), it);
        result.Data.push_back(Data.at(index));
      } else { // not found
        result.Data.push_back(NAN);
      }
    }
    result.Name = Name;
    result.StartFrequency = std::move(FrequencyWeekBased::Daily(min_date));
    result.Fields.insert(std::pair("conversion", "from date-list"));

  } else
    throw std::logic_error(
        "Converting from current type of frequency to 'Daily' frequency is not "
        "supported (or not implemented).");
}

template <typename Tw>
void Variable<Tw>::PartitionData(std::vector<Tv> &data,
                                 std::vector<std::vector<Tv>> &result, Ti size,
                                 bool fromEnd) {
  result.clear();
  if (fromEnd) {
    for (Ti i = (Ti)data.size(); i >= 0; i -= size) {
      int start = std::max(i - size, 0);
      result.insert(result.begin(),
                    std::vector<Tv>(data.begin() + start, data.begin() + i));
    }
  } else {
    for (Ti i = 0; i < (Ti)data.size(); i += size) {
      int end = std::min(i + size, static_cast<int>(data.size()));
      result.emplace_back(data.begin() + i, data.begin() + end);
    }
  }
}

template <typename Tw>
void Variable<Tw>::ConvertTo_MultiDaily(
    Variable &result, int k,
    const std::function<double(const std::vector<double> &)> &aggregateFunc,
    bool fromEnd) {

  if (StartFrequency.get()->mClass == FrequencyClass::kDaily) {

    std::vector<std::vector<double>> partitions;
    PartitionData(Data, partitions, k, fromEnd);
    std::vector<double> newdata;
    for (int i = 0; i < (int)partitions.size(); i++) {
      newdata.push_back(aggregateFunc(partitions.at(i)));
    }

    result.Data = newdata;
    result.Name = Name;

    auto start =
        dynamic_cast<FrequencyWeekBased const &>(*StartFrequency.get());
    result.StartFrequency =
        std::move(FrequencyWeekBased::MultiDaily(start.mDay, k));

  } else
    throw std::logic_error("Converting from current type of frequency to "
                           "'Multi-Day' frequency is not "
                           "supported (or not implemented).");
}

template class ldt::Variable<Tv>;
// template class ldt::Variable < Ti>;
