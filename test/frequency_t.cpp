/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "frequency.h"
#include <gtest/gtest.h>
#include <random>

using namespace ldt;

TEST(Frequency_t, cross_section) {
  // Compare
  auto y0 = FrequencyCrossSection(9);
  auto y1 = FrequencyCrossSection(9);
  auto y2 = FrequencyCrossSection(10);
  EXPECT_EQ(true, y0.IsEqualTo(y1));
  EXPECT_NE(true, y0.IsEqualTo(y2));
  EXPECT_EQ(true, y0.IsEqualOrOlderThan(y1));
  EXPECT_EQ(true, y0.IsOlderThan(y2));
  EXPECT_EQ(true, y2.IsNewerThan(y1));

  // Next
  y0.Next(4);
  y2.Next(3);
  EXPECT_EQ(true, y0.IsEqualTo(y2));

  // Minus
  EXPECT_EQ(0, y0.Minus(y2));
  EXPECT_EQ(-4, y1.Minus(y0));

  // Parse
  FrequencyClass fc;
  auto str = std::string("40");
  auto fStr = std::string("cs");
  auto y3 = Frequency::Parse(str, fStr, fc);
  auto y4 = FrequencyCrossSection(40);
  EXPECT_EQ(true, y3.get()->IsEqualTo(y4));
  EXPECT_EQ(fc, FrequencyClass::kCrossSection);

  // ToString
  auto y4str = y4.ToString();
  auto y4cstr = y4.ToClassString();
  auto y5 = Frequency::Parse(y4str, y4cstr, fc);
  EXPECT_EQ(true, y5.get()->IsEqualTo(y4));
}

TEST(Frequency_t, year_based_yearly) {
  // Compare
  auto y0 = *FrequencyYearBased::Yearly(9).get();
  auto y1 = *FrequencyYearBased::Yearly(9).get();
  auto y2 = *FrequencyYearBased::Yearly(10).get();
  EXPECT_EQ(true, y0.IsEqualTo(y1));
  EXPECT_NE(true, y0.IsEqualTo(y2));
  EXPECT_EQ(true, y0.IsEqualOrOlderThan(y1));
  EXPECT_EQ(true, y0.IsOlderThan(y2));
  EXPECT_EQ(true, y2.IsNewerThan(y1));

  // Next
  y0.Next(4);
  y2.Next(3);
  EXPECT_EQ(true, y0.IsEqualTo(y2));

  // Minus
  EXPECT_EQ(0, y0.Minus(y2));
  EXPECT_EQ(-4, y1.Minus(y0));

  // Parse
  FrequencyClass fc;
  auto str = std::string("40");
  auto fStr = std::string("y");
  auto y3 = Frequency::Parse(str, fStr, fc);
  auto y4 = *FrequencyYearBased::Yearly(40).get();
  EXPECT_EQ(true, y3.get()->IsEqualTo(y4));
  EXPECT_EQ(fc, FrequencyClass::kYearly);

  // ToString
  auto y4str = y4.ToString();
  auto y4cstr = y4.ToClassString();
  auto y5 = Frequency::Parse(y4str, y4cstr, fc);
  EXPECT_EQ(true, y5.get()->IsEqualTo(y4));
}

TEST(Frequency_t, year_based_xTimesZyears) {
  // we will test 3 times each 2 years

  // Compare
  auto q0 = *FrequencyYearBased::XTimesZYear(2000, 3, 1, 2).get();
  auto q1 = *FrequencyYearBased::XTimesZYear(2000, 3, 1, 2).get();
  auto q2 = *FrequencyYearBased::XTimesZYear(2000, 3, 2, 2).get();
  auto q3 = *FrequencyYearBased::XTimesZYear(2002, 3, 1, 2).get();
  EXPECT_EQ(true, q0.IsEqualTo(q1));
  EXPECT_EQ(false, q0.IsEqualTo(q2));
  EXPECT_EQ(true, q0.IsEqualOrOlderThan(q1));
  EXPECT_EQ(true, q0.IsOlderThan(q2));
  EXPECT_EQ(true, q3.IsNewerThan(q2));

  // Next
  q3 = *FrequencyYearBased::XTimesZYear(2000, 3, 1, 2).get();
  q3.Next(2);
  EXPECT_EQ(true, q3.IsEqualTo(
                      *FrequencyYearBased::XTimesZYear(2000, 3, 3, 2).get()));

  q3 = *FrequencyYearBased::XTimesZYear(2000, 3, 1, 2).get();
  q3.Next(3);
  EXPECT_EQ(true, q3.IsEqualTo(
                      *FrequencyYearBased::XTimesZYear(2002, 3, 1, 2).get()));

  q3 = *FrequencyYearBased::XTimesZYear(2000, 3, 1, 2).get();
  q3.Next(-1);
  EXPECT_EQ(true, q3.IsEqualTo(
                      *FrequencyYearBased::XTimesZYear(1998, 3, 3, 2).get()));

  q3 = *FrequencyYearBased::XTimesZYear(2000, 3, 1, 2).get();
  q3.Next(-4);
  EXPECT_EQ(true, q3.IsEqualTo(
                      *FrequencyYearBased::XTimesZYear(1996, 3, 3, 2).get()));

  // Minus
  q2 = *FrequencyYearBased::XTimesZYear(2000, 3, 1, 2).get();
  q3 = *FrequencyYearBased::XTimesZYear(2000, 3, 2, 2).get();
  EXPECT_EQ(1, q3.Minus(q2));
  EXPECT_EQ(-1, q2.Minus(q3));
  for (Ti i = -6; i < 6; i++) {
    q2 = *FrequencyYearBased::XTimesZYear(2000, 3, 1, 2).get();
    q3 = *FrequencyYearBased::XTimesZYear(2000, 3, 1, 2).get();
    q3.Next(i);
    EXPECT_EQ(i, q3.Minus(q2));
  }

  // Parse
  FrequencyClass fc;
  auto str = std::string("40:1");
  auto fStr = std::string("x3z2");
  auto y3 = Frequency::Parse(str, fStr, fc);
  auto y4 = *FrequencyYearBased::XTimesZYear(40, 3, 1, 2).get();
  EXPECT_EQ(true, y3.get()->IsEqualTo(y4));
  EXPECT_EQ(fc, FrequencyClass::kXTimesZYears);

  // ToString
  auto y4str = y4.ToString();
  auto y4cstr = y4.ToClassString();
  auto y5 = Frequency::Parse(y4str, y4cstr, fc);
  EXPECT_EQ(true, y5.get()->IsEqualTo(y4));
}

TEST(Frequency_t, year_based_quarterly) {
  // Compare
  auto q0 = *FrequencyYearBased::Quarterly(2000, 1).get();
  auto q1 = *FrequencyYearBased::Quarterly(2000, 1).get();
  auto q2 = *FrequencyYearBased::Quarterly(2000, 2).get();
  auto q3 = *FrequencyYearBased::Quarterly(2001, 1).get();
  EXPECT_EQ(true, q0.IsEqualTo(q1));
  EXPECT_EQ(false, q0.IsEqualTo(q2));
  EXPECT_EQ(true, q0.IsEqualOrOlderThan(q1));
  EXPECT_EQ(true, q0.IsOlderThan(q2));
  EXPECT_EQ(true, q3.IsNewerThan(q2));

  // Next
  q3 = *FrequencyYearBased::Quarterly(2000, 2).get();
  q3.Next(1);
  EXPECT_EQ(true, q3.IsEqualTo(*FrequencyYearBased::Quarterly(2000, 3).get()));

  q3 = *FrequencyYearBased::Quarterly(2000, 2).get();
  q3.Next(11);
  EXPECT_EQ(true, q3.IsEqualTo(*FrequencyYearBased::Quarterly(2003, 1).get()));

  q3 = *FrequencyYearBased::Quarterly(2000, 2).get();
  q3.Next(16);
  EXPECT_EQ(true, q3.IsEqualTo(*FrequencyYearBased::Quarterly(2004, 2).get()));

  q3 = *FrequencyYearBased::Quarterly(2000, 2).get();
  q3.Next(-1);
  EXPECT_EQ(true, q3.IsEqualTo(*FrequencyYearBased::Quarterly(2000, 1).get()));

  q3 = *FrequencyYearBased::Quarterly(2000, 2).get();
  q3.Next(-11);
  EXPECT_EQ(true, q3.IsEqualTo(*FrequencyYearBased::Quarterly(1997, 3).get()));

  q3 = *FrequencyYearBased::Quarterly(2000, 2).get();
  q3.Next(-16);
  EXPECT_EQ(true, q3.IsEqualTo(*FrequencyYearBased::Quarterly(1996, 2).get()));

  // Minus
  q0 = *FrequencyYearBased::Quarterly(2000, 2).get();
  q1 = *FrequencyYearBased::Quarterly(2000, 2).get();
  q2 = *FrequencyYearBased::Quarterly(2000, 3).get();
  q3 = *FrequencyYearBased::Quarterly(2001, 1).get();
  auto q4 = *FrequencyYearBased::Quarterly(2001, 3).get();
  EXPECT_EQ(0, q0.Minus(q1));
  EXPECT_EQ(-1, q0.Minus(q2));
  EXPECT_EQ(3, q3.Minus(q0));
  EXPECT_EQ(5, q4.Minus(q0));

  // Parse
  FrequencyClass fc;
  auto str = std::string("40q1");
  auto fStr = std::string("q");
  auto y3 = Frequency::Parse(str, fStr, fc);
  auto y4 = *FrequencyYearBased::Quarterly(40, 1).get();
  EXPECT_EQ(true, y3.get()->IsEqualTo(y4));
  EXPECT_EQ(fc, FrequencyClass::kQuarterly);

  // ToString
  auto y4str = y4.ToString();
  auto y4cstr = y4.ToClassString();
  auto y5 = Frequency::Parse(y4str, y4cstr, fc);
  EXPECT_EQ(true, y5.get()->IsEqualTo(y4));
}

boost::gregorian::date getDate(Ti year, Ti month, Ti day) {
  return boost::gregorian::date(year, month, day);
}

TEST(Frequency_t, week_based_weekly) {

  // Compare
  auto q0 = *FrequencyWeekBased::Weekly(getDate(2015, 4, 2)).get();
  auto q1 = *FrequencyWeekBased::Weekly(getDate(2015, 4, 2)).get();
  auto q2 = *FrequencyWeekBased::Weekly(getDate(2015, 4, 9)).get();
  auto q3 = *FrequencyWeekBased::Weekly(getDate(2015, 3, 2)).get();
  EXPECT_EQ(true, q0.IsEqualTo(q1));
  EXPECT_EQ(false, q0.IsEqualTo(q2));
  EXPECT_EQ(true, q0.IsEqualOrOlderThan(q1));
  EXPECT_EQ(true, q0.IsOlderThan(q2));
  EXPECT_EQ(false, q3.IsNewerThan(q2));

  // Next
  q3 = *FrequencyWeekBased::Weekly(getDate(2015, 4, 2)).get();
  q3.Next(1);
  EXPECT_EQ(true, q3.IsEqualTo(
                      *FrequencyWeekBased::Weekly(getDate(2015, 4, 9)).get()));

  q3 = *FrequencyWeekBased::Weekly(getDate(2015, 4, 2)).get();
  q3.Next(3);
  EXPECT_EQ(true, q3.IsEqualTo(
                      *FrequencyWeekBased::Weekly(getDate(2015, 4, 23)).get()));

  q3 = *FrequencyWeekBased::Weekly(getDate(2015, 4, 12)).get();
  q3.Next(-1);
  EXPECT_EQ(true, q3.IsEqualTo(
                      *FrequencyWeekBased::Weekly(getDate(2015, 4, 5)).get()));

  q3 = *FrequencyWeekBased::Weekly(getDate(2015, 4, 22)).get();
  q3.Next(-3);
  EXPECT_EQ(true, q3.IsEqualTo(
                      *FrequencyWeekBased::Weekly(getDate(2015, 4, 1)).get()));

  // Minus
  q0 = *FrequencyWeekBased::Weekly(getDate(2015, 4, 22)).get();
  q1 = *FrequencyWeekBased::Weekly(getDate(2015, 4, 15)).get();
  q2 = *FrequencyWeekBased::Weekly(getDate(2015, 4, 8)).get();
  q3 = *FrequencyWeekBased::Weekly(getDate(2015, 4, 1)).get();
  EXPECT_EQ(0, q0.Minus(q0));
  EXPECT_EQ(1, q0.Minus(q1));
  EXPECT_EQ(-1, q1.Minus(q0));
  EXPECT_EQ(-3, q3.Minus(q0));

  // Parse
  FrequencyClass fc;
  auto str = std::string("20150401");
  auto fStr = std::string("w");
  auto y3 = Frequency::Parse(str, fStr, fc);
  auto y4 = *FrequencyWeekBased::Weekly(getDate(2015, 4, 1)).get();
  EXPECT_EQ(true, y3.get()->IsEqualTo(y4));
  EXPECT_EQ(fc, FrequencyClass::kWeekly);

  // ToString
  auto y4str = y4.ToString();
  auto y4cstr = y4.ToClassString();
  auto y5 = Frequency::Parse(y4str, y4cstr, fc);
  EXPECT_EQ(true, y5.get()->IsEqualTo(y4));
}

TEST(Frequency_t, week_based_dailyInWeek) {

  DayOfWeekRange range;

  for (Ti i = 0; i < 2; i++) {
    if (i == 0)
      range = DayOfWeekRange(DayOfWeek::kMon, DayOfWeek::kWed); // 1-2-3
    else
      range = DayOfWeekRange(DayOfWeek::kSat, DayOfWeek::kMon); // 6-0-1

    // Compare
    auto q0 = *FrequencyWeekBased::DailyInWeek(getDate(2015, 4, 2), range, true)
                   .get();
    auto q1 = *FrequencyWeekBased::DailyInWeek(getDate(2015, 4, 2), range, true)
                   .get();
    auto q2 = *FrequencyWeekBased::DailyInWeek(getDate(2015, 4, 9), range, true)
                   .get();
    auto q3 = *FrequencyWeekBased::DailyInWeek(getDate(2015, 3, 2), range, true)
                   .get();
    EXPECT_EQ(true, q0.IsEqualTo(q1));
    EXPECT_EQ(false, q0.IsEqualTo(q2));
    EXPECT_EQ(true, q0.IsEqualOrOlderThan(q1));
    EXPECT_EQ(true, q0.IsOlderThan(q2));
    EXPECT_EQ(false, q3.IsNewerThan(q2));

    // Initialize
    q0 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 17), range,
                                          true)
              .get(); // Thursday and forward
    q1 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 17), range,
                                          false)
              .get(); // Thursday and backward
    q2 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 16), range,
                                          true)
              .get(); // Wednesday
    q3 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 19), range,
                                          true)
              .get(); // Saturday
    auto q4 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 21), range,
                                               true)
                   .get(); // Monday
    auto q5 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 14), range,
                                               true)
                   .get(); // Monday this week

    if (i == 0) {
      EXPECT_EQ(true, q0.IsEqualTo(q4));
      EXPECT_EQ(true, q1.IsEqualTo(q2));
    } else {
      EXPECT_EQ(true, q0.IsEqualTo(q3));
      EXPECT_EQ(true, q1.IsEqualTo(q5));
    }

    // Next
    q5 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 14), range, true)
              .get();
    q5.Next(1);
    if (i == 0)
      EXPECT_EQ(true, q5.IsEqualTo(*FrequencyWeekBased::DailyInWeek(
                                        getDate(2022, 11, 15), range, true)
                                        .get()));
    else
      EXPECT_EQ(true, q5.IsEqualTo(*FrequencyWeekBased::DailyInWeek(
                                        getDate(2022, 11, 19), range, true)
                                        .get()));

    q5 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 14), range, true)
              .get();
    q5.Next(4);
    if (i == 0)
      EXPECT_EQ(true, q5.IsEqualTo(*FrequencyWeekBased::DailyInWeek(
                                        getDate(2022, 11, 22), range, true)
                                        .get()));
    else
      EXPECT_EQ(true, q5.IsEqualTo(*FrequencyWeekBased::DailyInWeek(
                                        getDate(2022, 11, 26), range, true)
                                        .get()));

    q5 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 14), range, true)
              .get();
    q5.Next(-1);
    if (i == 0)
      EXPECT_EQ(true, q5.IsEqualTo(*FrequencyWeekBased::DailyInWeek(
                                        getDate(2022, 11, 9), range, true)
                                        .get()));
    else
      EXPECT_EQ(true, q5.IsEqualTo(*FrequencyWeekBased::DailyInWeek(
                                        getDate(2022, 11, 13), range, true)
                                        .get()));

    q5 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 14), range, true)
              .get();
    q5.Next(-5);
    if (i == 0)
      EXPECT_EQ(true, q5.IsEqualTo(*FrequencyWeekBased::DailyInWeek(
                                        getDate(2022, 11, 1), range, true)
                                        .get()));
    else
      EXPECT_EQ(true, q5.IsEqualTo(*FrequencyWeekBased::DailyInWeek(
                                        getDate(2022, 11, 5), range, true)
                                        .get()));

    // Minus
    q5 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 14), range,
                                          true)
              .get(); // Monday
    EXPECT_EQ(0, q4.Minus(q4));
    EXPECT_EQ(3, q4.Minus(q5));
    EXPECT_EQ(-3, q5.Minus(q4));
    q1 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 13), range,
                                          false)
              .get(); // 0:Wednesday, 1:Sunday
    EXPECT_EQ(1, q5.Minus(q1));
    EXPECT_EQ(-1, q1.Minus(q5));
    q1 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 12), range, false)
              .get();
    if (i == 0) {
      EXPECT_EQ(1, q5.Minus(q1));
      EXPECT_EQ(-1, q1.Minus(q5));
    } else {
      EXPECT_EQ(2, q5.Minus(q1));
      EXPECT_EQ(-2, q1.Minus(q5));
    }
    q1 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 7), range, false)
              .get();
    EXPECT_EQ(3, q5.Minus(q1));
    EXPECT_EQ(-3, q1.Minus(q5));
    q1 = *FrequencyWeekBased::DailyInWeek(getDate(2022, 10, 31), range, false)
              .get();
    EXPECT_EQ(6, q5.Minus(q1));
    EXPECT_EQ(-6, q1.Minus(q5));

    // Parse
    FrequencyClass fc;
    auto str = std::string("20221114");
    auto fStr = std::string(i == 0 ? "i:mon-wed" : "i:sat-mon");
    auto y3 = Frequency::Parse(str, fStr, fc);
    auto y4 =
        *FrequencyWeekBased::DailyInWeek(getDate(2022, 11, 14), range, true)
             .get();
    EXPECT_EQ(true, y3.get()->IsEqualTo(y4));
    EXPECT_EQ(fc, FrequencyClass::kDailyInWeek);

    // ToString
    auto y4str = y4.ToString();
    auto y4cstr = y4.ToClassString();
    auto y5 = Frequency::Parse(y4str, y4cstr, fc);
    EXPECT_EQ(true, y5.get()->IsEqualTo(y4));
  }
}

TEST(Frequency_t, list_string) {
  // Compare
  auto items = std::vector<std::string>(
      {std::string("A"), std::string("B"), std::string("C"), std::string("D")});
  auto q0 = FrequencyList<std::string>(std::string("A"), &items);
  auto q1 = FrequencyList<std::string>(std::string("A"), &items);
  auto q2 = FrequencyList<std::string>(std::string("B"), &items);
  auto q3 = FrequencyList<std::string>(std::string("C"), &items);
  auto q4 = FrequencyList<std::string>(std::string("D"), &items);
  EXPECT_EQ(true, q0.IsEqualTo(q1));
  EXPECT_EQ(false, q0.IsEqualTo(q2));
  EXPECT_EQ(true, q0.IsEqualOrOlderThan(q1));
  EXPECT_EQ(true, q0.IsOlderThan(q2));
  EXPECT_EQ(true, q3.IsNewerThan(q2));

  // Next
  q3.Next(1);
  EXPECT_EQ(true, q3.IsEqualTo(q4));
  q2.Next(2);
  EXPECT_EQ(true, q3.IsEqualTo(q4));

  //     out
  //          forward
  q2 = FrequencyList<std::string>(std::string("B"), &items);
  q2.Next(2);
  EXPECT_EQ(q2.mValue, std::string("D"));
  q2.Next(1);
  EXPECT_EQ(q2.ToString(), std::string("out_item:1"));
  q2 = FrequencyList<std::string>(std::string("B"), &items);
  q2.Next(3);
  EXPECT_EQ(q2.ToString(), std::string("out_item:1"));
  q2 = FrequencyList<std::string>(std::string("B"), &items);
  q2.Next(4);
  EXPECT_EQ(q2.ToString(), std::string("out_item:2"));
  q2.Next(1);
  EXPECT_EQ(q2.ToString(), std::string("out_item:3"));
  q2.Next(10);
  EXPECT_EQ(q2.ToString(), std::string("out_item:13"));

  //          backward
  q2 = FrequencyList<std::string>(std::string("B"), &items);
  q2.Next(-1);
  EXPECT_EQ(q2.ToString(), std::string("A"));
  q2.Next(-1);
  EXPECT_EQ(q2.ToString(), std::string("out_item:-1"));
  q2 = FrequencyList<std::string>(std::string("B"), &items);
  q2.Next(-3);
  EXPECT_EQ(q2.ToString(), std::string("out_item:-2"));
  q2 = FrequencyList<std::string>(std::string("B"), &items);
  q2.Next(-4);
  EXPECT_EQ(q2.ToString(), std::string("out_item:-3"));
  q2.Next(-1);
  EXPECT_EQ(q2.ToString(), std::string("out_item:-4"));
  q2.Next(-4);
  EXPECT_EQ(q2.ToString(), std::string("out_item:-8"));
  q2.Next(4);
  EXPECT_EQ(q2.ToString(), std::string("out_item:-4"));
  q2.Next(4);
  EXPECT_EQ(q2.ToString(), std::string("A"));
  q2.Next(-1);
  EXPECT_EQ(q2.ToString(), std::string("out_item:-1"));
  q2.Next(4);
  EXPECT_EQ(q2.ToString(), std::string("D"));

  q2 = FrequencyList<std::string>(std::string("B"), &items);
  q3 = FrequencyList<std::string>(std::string("C"), &items);
  q3.Next(-1);
  EXPECT_EQ(true, q3.IsEqualTo(q2));
  q1 = FrequencyList<std::string>(std::string("A"), &items);
  q4 = FrequencyList<std::string>(std::string("D"), &items);
  q4.Next(-3);
  EXPECT_EQ(true, q4.IsEqualTo(q1));

  // Minus
  q1 = FrequencyList<std::string>(std::string("A"), &items);
  q2 = FrequencyList<std::string>(std::string("B"), &items);
  q3 = FrequencyList<std::string>(std::string("C"), &items);
  q4 = FrequencyList<std::string>(std::string("D"), &items);
  EXPECT_EQ(0, q1.Minus(q1));
  EXPECT_EQ(-1, q1.Minus(q2));
  EXPECT_EQ(1, q2.Minus(q1));
  EXPECT_EQ(3, q4.Minus(q1));
  EXPECT_EQ(-3, q1.Minus(q4));

  //     out
  q2.Next(10);
  EXPECT_EQ(11, q2.Minus(q1));
  q3.Next(-10);
  EXPECT_EQ(-8, q3.Minus(q1));
  q4.Next(2);
  q1.Next(-3);
  EXPECT_EQ(8, q4.Minus(q1));

  // Parse
  FrequencyClass fc;
  auto str = std::string("A");
  auto fStr = std::string("Ls:A;B;C;D");
  auto items0 = std::vector<std::string>();
  auto y3 = FrequencyList<std::string>::ParseList(str, fStr, fc, items0);
  auto y4 = FrequencyList<std::string>(std::string("A"), &items);
  EXPECT_EQ(true, y3.get()->IsEqualTo(y4));
  EXPECT_EQ(fc, FrequencyClass::kListString);

  // ToString
  auto y4str = y4.ToString();
  auto y4cstr = y4.ToClassString(true);
  items0.clear();
  auto y5 = FrequencyList<std::string>::ParseList(y4str, y4cstr, fc, items0);
  EXPECT_EQ(true, y5.get()->IsEqualTo(y4));
}

TEST(Frequency_t, day_based_xTimes_inWeek) {

  DayOfWeekRange range = DayOfWeekRange(DayOfWeek::kMon, DayOfWeek::kWed);
  auto day1 =
      *FrequencyWeekBased::DailyInWeek(getDate(2015, 4, 2), range, true).get();
  auto day2 =
      *FrequencyWeekBased::DailyInWeek(getDate(2015, 5, 1), range, true).get();

  // Compare
  auto q0 = *FrequencyDayBased::XTimesADay(day1, 12, 5).get();
  auto q1 = *FrequencyDayBased::XTimesADay(day1, 12, 5).get();
  auto q2 = *FrequencyDayBased::XTimesADay(day1, 12, 12).get();
  auto q3 = *FrequencyDayBased::XTimesADay(day1, 12, 1).get();
  auto q4 = *FrequencyDayBased::XTimesADay(day2, 12, 5).get();

  EXPECT_EQ(true, q0.IsEqualTo(q1));
  EXPECT_EQ(false, q0.IsEqualTo(q2));
  EXPECT_EQ(true, q0.IsEqualOrOlderThan(q1));
  EXPECT_EQ(true, q0.IsOlderThan(q2));
  EXPECT_EQ(false, q3.IsNewerThan(q2));
  EXPECT_EQ(true, q4.IsNewerThan(q0));

  // Next
  q0 = *FrequencyDayBased::XTimesADay(day1, 12, 5).get();
  q0.Next(1);
  EXPECT_EQ(true,
            q0.IsEqualTo(*FrequencyDayBased::XTimesADay(day1, 12, 6).get()));

  q0 = *FrequencyDayBased::XTimesADay(day1, 12, 5).get();
  q0.Next(8);
  day1.Next(1);
  EXPECT_EQ(true,
            q0.IsEqualTo(*FrequencyDayBased::XTimesADay(day1, 12, 1).get()));

  q0 = *FrequencyDayBased::XTimesADay(day1, 12, 5).get();
  q0.Next(-1);
  EXPECT_EQ(true,
            q0.IsEqualTo(*FrequencyDayBased::XTimesADay(day1, 12, 4).get()));

  q0 = *FrequencyDayBased::XTimesADay(day1, 12, 5).get();
  q0.Next(-5);
  day1.Next(-1);
  EXPECT_EQ(true,
            q0.IsEqualTo(*FrequencyDayBased::XTimesADay(day1, 12, 12).get()));

  // Minus
  q0 = *FrequencyDayBased::XTimesADay(day1, 12, 5).get();
  q1 = *FrequencyDayBased::XTimesADay(day1, 12, 6).get();
  EXPECT_EQ(0, q0.Minus(q0));
  EXPECT_EQ(1, q1.Minus(q0));
  EXPECT_EQ(-1, q0.Minus(q1));

  q0 = *FrequencyDayBased::XTimesADay(day1, 12, 5).get();
  day1.Next(1);
  q1 = *FrequencyDayBased::XTimesADay(day1, 12, 6).get();
  EXPECT_EQ(13, q1.Minus(q0));
  EXPECT_EQ(-13, q0.Minus(q1));

  q0 = *FrequencyDayBased::XTimesADay(day1, 12, 6).get();
  day1.Next(10);
  q1 = *FrequencyDayBased::XTimesADay(day1, 12, 7).get();
  EXPECT_EQ(10 * 12 + 1, q1.Minus(q0));
  EXPECT_EQ(-10 * 12 - 1, q0.Minus(q1));

  // Parse
  FrequencyClass fc;
  auto str = std::string("20150402:5/12");
  auto fStr = std::string("da12|i:mon-wed");
  auto y3 = Frequency::Parse(str, fStr, fc);
  day1 =
      *FrequencyWeekBased::DailyInWeek(getDate(2015, 4, 2), range, true).get();
  auto y4 = *FrequencyDayBased::XTimesADay(day1, 12, 5).get();
  /*
  EXPECT_EQ(true, y3.get()->IsEqualTo(y4));
  EXPECT_EQ(fc, FrequencyClass::kXTimesADay);*/

  // ToString
  auto y4str = y4.ToString();
  auto y4cstr = y4.ToClassString();
  auto y5 = Frequency::Parse(y4str, y4cstr, fc);
  EXPECT_EQ(true, y5.get()->IsEqualTo(y4));
}

TEST(Frequency_t, parse_general) {
  auto values = std::vector<std::unique_ptr<Frequency>>();
  auto listS = std::vector<std::string>();
  auto listD = std::vector<boost::gregorian::date>();
  Frequency::Examples(values, listS, listD);
  FrequencyClass fc;
  for (auto &a : values) {
    auto b = a.get();
    auto str = b->ToString();
    auto classStr = b->ToClassString();
    auto c = Frequency::Parse(str, classStr, fc);
    if (fc == FrequencyClass::kListString) {
      auto d = dynamic_cast<FrequencyList<std::string> &>(*c.get());
      d.pItems = &listS;
      EXPECT_EQ(true, d.IsEqualTo(*b));
    } else if (fc == FrequencyClass::kListDate) {
      auto d = dynamic_cast<FrequencyList<boost::gregorian::date> &>(*c.get());
      d.pItems = &listD;
      EXPECT_EQ(true, d.IsEqualTo(*b));
    } else
      EXPECT_EQ(true, c.get()->IsEqualTo(*b));
  }
}

TEST(Frequency_t, next_general) {
  auto values = std::vector<std::unique_ptr<Frequency>>();
  auto listS = std::vector<std::string>();
  auto listD = std::vector<boost::gregorian::date>();
  Frequency::Examples(values, listS, listD);
  FrequencyClass fc;
  for (auto &a : values) {
    auto b = a.get();
    if (fc == FrequencyClass::kListString) {
      auto d = dynamic_cast<FrequencyList<std::string> &>(*b);
      d.pItems = &listS;
    } else if (fc == FrequencyClass::kListDate) {
      auto d = dynamic_cast<FrequencyList<boost::gregorian::date> &>(*b);
      d.pItems = &listD;
    }

    for (Ti i = -4; i < 4; i++) {
      auto s1 = b->ToString();
      b->Next(i);
      b->Next(-i);
      auto s2 = b->ToString();
      EXPECT_EQ(true, AreEqual_i(s1.c_str(), s2.c_str()));
    }
  }
}
