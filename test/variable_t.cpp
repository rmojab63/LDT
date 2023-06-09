/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "matrix.h"
#include "variable.h"
#include <gtest/gtest.h>
#include <random>

using namespace ldt;

TEST(Variable_t, import_export) {

  std::vector<std::unique_ptr<Variable<Tv>>> values;
  std::vector<std::unique_ptr<Frequency>> f_values;
  std::vector<std::string> listItemsString;
  std::vector<boost::gregorian::date> listItemsDate;

  std::mt19937 eng(340);
  Variable<Tv>::Examples(values, f_values, listItemsString, listItemsDate, eng);
  for (auto &v : values) {
    auto b = v.get();
    auto str = b->ToString();

    std::vector<std::string> listItemsString0;
    std::vector<boost::gregorian::date> listItemsDate0;
    Variable<Tv> w;
    b->Parse(str, w, listItemsString0, listItemsDate0);

    EXPECT_EQ(true, w.IsEqualTo(*b));
  }
}

TEST(Variable_t, variables_init) {
  Variable<Tv> v1;
  v1.Name = std::string("V1");
  v1.StartFrequency = std::unique_ptr<Frequency>(new FrequencyCrossSection(2));
  v1.Data.insert(v1.Data.end(), {1, 2, 3});

  Variable<Tv> v2;
  v2.Name = std::string("V2");
  v2.StartFrequency = std::unique_ptr<Frequency>(new FrequencyCrossSection(1));
  v2.Data.insert(v2.Data.end(), {4, 5, 6, 7, 8});

  auto vs = Variables<Tv>(std::vector<Variable<Tv> *>({&v1, &v2}));
  auto d = Matrix<Tv>(&vs.Data, vs.NumObs, (Ti)vs.Names.size());
  auto ex = Matrix<Tv>(new Tv[10]{NAN, 1, 2, 3, NAN, 4, 5, 6, 7, 8}, 5, 2);
  EXPECT_EQ(true, d.Equals(ex, 1e-14, true, true));
}

TEST(Variable_t, trim) {
  Variable<Tv> v1;
  v1.Name = std::string("V1");
  v1.StartFrequency = std::unique_ptr<Frequency>(new FrequencyCrossSection(1));
  v1.Data.insert(v1.Data.end(), {NAN, NAN, 1, 2, 3, NAN});

  Variable<Tv> v2;
  v2.Name = std::string("V2");
  v2.StartFrequency = std::unique_ptr<Frequency>(new FrequencyCrossSection(3));
  v2.Data.insert(v2.Data.end(), {1, 2, 3});

  v1.Trim();
  EXPECT_EQ(v1.Data, v2.Data);
  EXPECT_EQ(true, v1.StartFrequency.get()->IsEqualTo(*v2.StartFrequency.get()));
}

TEST(Variable_t, range) {

  Variable<Tv> v1;
  v1.Name = std::string("V1");
  v1.StartFrequency = std::unique_ptr<Frequency>(new FrequencyCrossSection(1));
  v1.Data.insert(v1.Data.end(), {NAN, 1, 2, 3, 4, 5, 6, NAN, NAN});

  Variable<Tv> v2;
  v2.Name = std::string("V2");
  v2.StartFrequency = std::unique_ptr<Frequency>(new FrequencyCrossSection(1));
  v2.Data.insert(v2.Data.end(), {NAN, NAN, 2, 3, NAN, 5, NAN, NAN, NAN});

  auto vs = Variables<Tv>(std::vector<Variable<Tv> *>({&v1, &v2}));

  bool hasmissing;
  auto range = vs.GetRange(1, hasmissing);
  ASSERT_EQ(std::get<0>(range), (Ti)2);
  ASSERT_EQ(std::get<1>(range), (Ti)5);
  ASSERT_EQ(hasmissing, true);

  range = vs.GetRange(0, hasmissing);
  ASSERT_EQ(std::get<0>(range), (Ti)1);
  ASSERT_EQ(std::get<1>(range), (Ti)6);
  ASSERT_EQ(hasmissing, false);
}

TEST(Variable_t, con_daylist_daily) {

  std::vector<boost::gregorian::date> dates;
  dates.push_back(boost::gregorian::from_string("2022-06-10"));
  dates.push_back(boost::gregorian::from_string("2022-06-08"));
  dates.push_back(boost::gregorian::from_string("2022-06-12"));
  Variable<Tv> v;
  v.Name = std::string("V1");
  v.StartFrequency = std::unique_ptr<Frequency>(
      new FrequencyList<boost::gregorian::date>(dates.at(0), &dates));
  v.Data.insert(v.Data.end(), {2, 1, 3});

  Variable<Tv> w;
  v.ConvertTo_Daily(w, nullptr);

  ASSERT_EQ(5, w.Data.size());
  ASSERT_EQ(1, w.Data.at(0));
  ASSERT_EQ(2, w.Data.at(2));
  ASSERT_EQ(3, w.Data.at(4));
  ASSERT_TRUE(std::isnan(w.Data.at(1)));
  ASSERT_TRUE(std::isnan(w.Data.at(3)));
}

TEST(Variable_t, con_dayinweek_daily) {

  auto range = DayOfWeekRange(DayOfWeek::kSat, DayOfWeek::kWed);
  Variable<Tv> v;
  v.Name = std::string("V1");
  v.StartFrequency = std::move(FrequencyWeekBased::DailyInWeek(
      boost::gregorian::from_string("2023-06-07"), range, true));
  v.Data.insert(v.Data.end(), {1, 2, 3, 4, 5, 6, 7, 8, 9, 10});

  Variable<Tv> w;
  v.ConvertTo_Daily(w, nullptr);
  ASSERT_EQ(12, w.Data.size());
  ASSERT_EQ(1, w.Data.at(0));
  ASSERT_EQ(5, w.Data.at(4));
  ASSERT_TRUE(std::isnan(w.Data.at(5)));
  ASSERT_TRUE(std::isnan(w.Data.at(6)));
  ASSERT_EQ(6, w.Data.at(7));
  ASSERT_EQ(10, w.Data.at(11));

  // smaller week
  range = DayOfWeekRange(DayOfWeek::kSat, DayOfWeek::kTue);
  v.StartFrequency = std::move(FrequencyWeekBased::DailyInWeek(
      boost::gregorian::from_string("2023-06-07"), range, true));
  v.ConvertTo_Daily(w, nullptr);
  ASSERT_EQ(16, w.Data.size());
  ASSERT_EQ(1, w.Data.at(0));
  ASSERT_EQ(4, w.Data.at(3));
  ASSERT_TRUE(std::isnan(w.Data.at(4)));
  ASSERT_TRUE(std::isnan(w.Data.at(6)));
  ASSERT_EQ(5, w.Data.at(7));
  ASSERT_EQ(10, w.Data.at(15));
  ASSERT_TRUE(std::isnan(w.Data.at(11)));
  ASSERT_TRUE(std::isnan(w.Data.at(13)));
}

TEST(Variable_t, con_daily_weekly) {

  Variable<Tv> v;
  v.Name = std::string("V1");
  v.StartFrequency =
      FrequencyWeekBased::Daily(boost::gregorian::from_string("2023-05-29"));
  v.Data.insert(v.Data.end(), {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});

  std::function<double(const std::vector<double> &)> func =
      [](const std::vector<double> &v) {
        Tv m;
        Array<Tv>::Mean<true>(&v[0], v.size(), m);
        return m;
      };
  Variable<Tv> w;

  v.ConvertTo_Weekly(w, DayOfWeek::kSat, &func);
  ASSERT_EQ((1 + 2 + 3 + 4 + 5) / 5.0, w.Data.at(0));
  ASSERT_EQ((6 + 7.0 + 8 + 9 + 10 + 11 + 12) / 7.0, w.Data.at(1));
  ASSERT_EQ((13.0 + 14) / 2.0, w.Data.at(2));

  v.ConvertTo_Weekly(w, DayOfWeek::kMon, &func);
  ASSERT_EQ(4, w.Data.at(0));
  ASSERT_EQ(11, w.Data.at(1));
}

TEST(Variable_t, con_daily_monthly) {

  Variable<Tv> v;
  v.Name = std::string("V1");
  v.StartFrequency =
      FrequencyWeekBased::Daily(boost::gregorian::from_string("2023-06-01"));
  v.Data = std::vector<double>(70);
  std::iota(v.Data.begin(), v.Data.begin() + 70, 1);

  std::function<double(const std::vector<double> &)> func =
      [](const std::vector<double> &v) {
        Tv m;
        Array<Tv>::Mean<true>(&v[0], v.size(), m);
        return m;
      };

  Variable<Tv> w;

  v.ConvertTo_XxYear_month_based<12>(w, &func);
  ASSERT_EQ(15.5, w.Data.at(0));
  ASSERT_EQ(46, w.Data.at(1));
  ASSERT_EQ(66, w.Data.at(2));
}

TEST(Variable_t, con_daily_5xyear) {

  Variable<Tv> v;
  v.Name = std::string("V1");
  v.StartFrequency =
      FrequencyWeekBased::Daily(boost::gregorian::from_string("2023-01-01"));
  v.Data = std::vector<double>(5 * 365);
  std::iota(v.Data.begin(), v.Data.begin() + 5 * 365, 1);

  std::function<double(const std::vector<double> &)> func =
      [](const std::vector<double> &v) {
        Tv m;
        Array<Tv>::Mean<true>(&v[0], v.size(), m);
        return m;
      };

  Variable<Tv> w;

  v.ConvertTo_XxYear(w, 5, &func);
  ASSERT_EQ(25, w.Data.size());
}
