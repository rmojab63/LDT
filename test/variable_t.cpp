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

TEST(Variable_t, convert_datelist) {

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
  v.ConvertTo_Daily(w);

  ASSERT_EQ(5, w.Data.size());
  ASSERT_EQ(1, w.Data.at(0));
  ASSERT_EQ(2, w.Data.at(2));
  ASSERT_EQ(3, w.Data.at(4));
  ASSERT_TRUE(std::isnan(w.Data.at(1)));
  ASSERT_TRUE(std::isnan(w.Data.at(3)));
}

TEST(Variable_t, partition) {

  auto v = std::vector<Tv>({1, 2, 3, 4, 5, 6, 7});
  std::vector<std::vector<Tv>> partitions;
  Variable<Tv>::PartitionData(v, partitions, 2, false);
  ASSERT_EQ(4, partitions.size());
  ASSERT_EQ(std::vector<Tv>({1, 2}), partitions.at(0));
  ASSERT_EQ(std::vector<Tv>({7}), partitions.at(3));

  // backward
  Variable<Tv>::PartitionData(v, partitions, 2, true);
  ASSERT_EQ(4, partitions.size());
  ASSERT_EQ(std::vector<Tv>({1}), partitions.at(0));
  ASSERT_EQ(std::vector<Tv>({6, 7}), partitions.at(3));
}