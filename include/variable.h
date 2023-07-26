#pragma once

#include "array.h"
#include "frequency.h"
#include "helpers.h"
#include "ldt_base.h"
#include <map>
#include <random>
#include <string>
#include <vector>

namespace ldt {

/// @brief A class for representing vector of data with name, frequency and
/// field
/// @tparam Tw Type of data
template <class Tw = Tv> class LDT_EXPORT Variable {

public:
  /// @brief Data of the variable.
  /// array
  std::vector<Tw> Data;

  /// @brief A pointer to the frequency of the first data-point.
  std::unique_ptr<Frequency> StartFrequency;

  /// @brief Name of the variable
  std::string Name = std::string("V");

  /// @brief Fields of the variable. Keys must not contain ';' or parsing fails
  std::map<std::string, std::string> Fields;

  /// @brief Initializes a new instance of this class
  Variable();

  /// @brief Gets the end frequency
  /// @return End frequency
  std::unique_ptr<Frequency> GetEndFrequency() const;

  bool IsEqualTo(Variable<Tw> &other, const Tv &epsilon = 1e-14) const;

  /// @brief Converts this variable to a string (compact form)
  /// @return The string
  std::string ToString() const;

  /// @brief Parses a string (compact form) to a variable
  /// @param str The string
  /// @param variable On exit, the updated variable
  /// @param listItemsString On exit, if frequency is of type string list, it is
  /// its items.
  /// @param listItemsDate n exit, if frequency is of type string list, it is
  /// its items.
  static void Parse(const std::string &str, Variable<Tw> &result,
                    std::vector<std::string> &listItemsString,
                    std::vector<boost::gregorian::date> &listItemsDate);

  /// @brief Generates a list of variables with all available frequencies
  /// (generally to be used in testing)
  /// @param values On exit, it contains the examples
  /// @param f_values On exit, it contains the frequency of the examples
  /// @param listItemsString On exit, it contains the list items for the string
  /// list frequency
  /// @param listItemsDate On exit, it contains the list items for the date list
  /// frequency
  /// @param end rnd
  static void Examples(std::vector<std::unique_ptr<Variable<Tw>>> &values,
                       std::vector<std::unique_ptr<Frequency>> &f_values,
                       std::vector<std::string> &listItemsString,
                       std::vector<boost::gregorian::date> &listItemsDate,
                       std::mt19937 &eng);

  /// @brief Removes leading and trailing  NAN and adjusts the \ref
  /// StartFrequency
  /// @return Range of original data. It might be invalid (start>end).
  IndexRange Trim();

  IndexRange GetRange(bool &hasMissing);

  IndexRange Interpolate(Ti &count);

  void ConvertTo_Daily(
      Variable &result,
      const std::function<Tv(const std::vector<Tv> &)> *aggregateFunc) const;

  void ConvertTo_MultiDaily(
      Variable &result, Ti k,
      const std::function<Tv(const std::vector<Tv> &)> *aggregateFunc,
      bool fromEnd = true) const;

  void ConvertTo_Weekly(
      Variable &result, DayOfWeek firstDayOfWeek,
      const std::function<Tv(const std::vector<Tv> &)> *aggregateFunc) const;

  void ConvertTo_XxYear(
      Variable &result, Ti x,
      const std::function<Tv(const std::vector<Tv> &)> *aggregateFunc) const;

  template <Ti x>
  void ConvertTo_XxYear_month_based(
      Variable &result,
      const std::function<Tv(const std::vector<Tv> &)> *aggregateFunc) const {

    /*if (mclass == FrequencyClass::kListDate) { // It can be more efficient
      } else*/
    if (StartFrequency.get()->mClass == FrequencyClass::kDaily) {
      auto start =
          dynamic_cast<FrequencyWeekBased const &>(*StartFrequency.get());

      if (!aggregateFunc)
        throw LdtException(ErrorType::kLogic, "variable.h",
                           "aggregate function is missing");
      auto aggF = *aggregateFunc;

      result.Data.clear();
      std::vector<Tv> partdata;

      auto current_part = get_part_month_based<x>(start.mDay);
      auto part = current_part;

      for (Ti i = 0; i < (Ti)Data.size(); i++) {

        auto day = start.mDay + boost::gregorian::date_duration(i);
        part = get_part_month_based<x>(day);

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
      if constexpr (x == 1) { // it should be different because 'get_part...'
                              // returns the year.
        result.StartFrequency =
            std::move(FrequencyYearBased::Yearly(start.mDay.year()));
      } else {
        result.StartFrequency = std::move(FrequencyYearBased::XTimesAYear(
            start.mDay.year(), x, get_part_month_based<x>(start.mDay)));
      }
    } else {
      throw LdtException(ErrorType::kLogic, "variable.h",
                         "direct conversion from current type of frequency to "
                         "'x times a year' frequency is not "
                         "supported (or not implemented)");
    }
  }

private:
  template <Ti x> static Ti get_part_month_based(boost::gregorian::date date) {

    if constexpr (x == 1) {
      Ti year = date.year();
      return year;
    } else if constexpr (true) {
      // I divide year based on the month

      Ti month = date.month();
      if (x == 12) {
        return month;
      } else if constexpr (x == 6) {
        return (month - 1) / 2 + 1;
      } else if constexpr (x == 4) {
        return month > 9 ? 4 : month > 6 ? 3 : month > 3 ? 2 : 1;
      } else if constexpr (x == 3) {
        return month > 8 ? 3 : month > 4 ? 2 : 1;
      } else if constexpr (x == 2) {
        return month > 6 ? 2 : 1;
      } else if constexpr (x == 24) {
        auto day = date.day();
        if constexpr (x == 24) {
          return 2 * (month - 1) + (day > 15 ? 2 : 1);
        }
      } else if constexpr (true) {
        throw LdtException(ErrorType::kLogic, "variable.h",
                           "it is not implemented for this type of "
                           "partition. Use the non-templated method");
      }
    }
  }
};

/// @brief A class for grouping a list of variables with similar frequencies
/// @tparam Tw Type of data
template <class Tw = Tv> class LDT_EXPORT Variables {

public:
  /// @brief Data of the variable. You can create a
  /// column major matrix from it array
  std::vector<Tw> Data;

  /// @brief A pointer to the frequency of the first
  /// data-point.
  std::unique_ptr<Frequency> StartFrequency;

  /// @brief Number of observations
  Ti NumObs = 0;

  /// @brief Names of the variables
  std::vector<std::string> Names;

  /// @brief Initializes a new instance of this class
  Variables();

  /// @brief Initializes a new instance of this class
  /// @param vars List of variables
  Variables(const std::vector<Variable<Tw> *> vars);

  /// @brief Gets the range of data in a variable
  /// @param j the variable index
  /// @param hasMissing on exit, it determines if the
  /// variable contains missing observations
  /// @return a tuple in which the first value is the
  /// start and the second value is the end (it might
  /// be invalid: start>end). Its length will be
  /// 'end-start+1'
  std::tuple<Ti, Ti> GetRange(Ti j, bool &hasMissing);
};

extern template class ldt::Variable<Tv>;
// extern template class ldt::Variable < Ti>;

extern template class ldt::Variables<Tv>;
// extern template class ldt::Variable < Ti>;

} // namespace ldt
