/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "matrix.h"
#include "searchers.h"
#include <gtest/gtest.h>

using namespace ldt;

TEST(Searcher_T, indexes1) {

  auto searchItems = SearchItems();
  searchItems.KeepAll = true;
  searchItems.LengthEvals = 1;
  searchItems.LengthTargets = 6;
  searchItems.Length1 = 1;

  auto gi = std::vector<std::vector<Ti>>({
      std::vector<Ti>({1, 2}),
      std::vector<Ti>({3}),
      std::vector<Ti>({4, 5, 6}),
  });
  auto searchOptions = SearchOptions();
  auto measures = SearchMeasureOptions();
  auto checks = SearchModelChecks();

  auto s1 =
      SearcherTest(searchOptions, searchItems, measures, checks, 3, gi, 0);
  auto sers = std::vector<Searcher *>({&s1});
  auto ms = ModelSet(sers, gi, searchOptions, searchItems, measures, checks);

  ms.Start(nullptr, nullptr);

  auto list0 = std::vector<SearcherSummary *>();
  auto list1 = std::vector<SearcherSummary *>();
  auto list2 = std::vector<SearcherSummary *>();
  auto d = SearcherModelingInfo();
  ms.CombineInfo(d, list0, list1, list2);
  auto v = std::vector<std::string>();
  for (auto &a : list0.at(0)->All)
    v.push_back(a->Dependents.ToString0('\t', ','));

  ASSERT_EQ((Ti)6, d.SearchedCount);
  ASSERT_EQ(d.ExpectedCount, d.SearchedCount);
  ASSERT_EQ(d.FailsCount.size(), (Ti)0);
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "1,3,4") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "2,3,4") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "2,3,5") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "1,3,5") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "1,3,6") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "2,3,4") != v.end());
}

TEST(Searcher_T, indexes2) {

  auto searchItems = SearchItems();
  searchItems.KeepAll = true;
  searchItems.LengthEvals = 1;
  searchItems.LengthTargets = 7;
  searchItems.Length1 = 1;

  auto gi = std::vector<std::vector<Ti>>(
      {std::vector<Ti>({1, 2}), std::vector<Ti>({3}),
       std::vector<Ti>({4, 5, 6}), std::vector<Ti>({7, 8})});
  auto searchOptions = SearchOptions();
  auto measures = SearchMeasureOptions();
  auto checks = SearchModelChecks();

  auto s1 =
      SearcherTest(searchOptions, searchItems, measures, checks, 3, gi, 0);
  auto sers = std::vector<Searcher *>({&s1});
  auto ms = ModelSet(sers, gi, searchOptions, searchItems, measures, checks);

  ms.Start(nullptr, nullptr);

  auto list0 = std::vector<SearcherSummary *>();
  auto list1 = std::vector<SearcherSummary *>();
  auto list2 = std::vector<SearcherSummary *>();
  auto d = SearcherModelingInfo();
  ms.CombineInfo(d, list0, list1, list2);

  auto v = std::vector<std::string>();
  for (auto &a : list0.at(0)->All)
    v.push_back(a->Dependents.ToString0('\t', ','));

  ASSERT_EQ((Ti)28, d.SearchedCount);
  ASSERT_EQ(d.ExpectedCount, d.SearchedCount);
  ASSERT_EQ(d.FailsCount.size(), (Ti)0);

  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "3,4,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "2,4,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "1,6,8") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "2,3,8") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "2,3,6") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "2,5,7") != v.end());
}

TEST(Searcher_T, indexes3_skip_target) {

  auto searchItems = SearchItems();
  searchItems.KeepAll = true;
  searchItems.LengthEvals = 1;
  searchItems.LengthTargets = 1;
  searchItems.Length1 = 1;

  auto gi = std::vector<std::vector<Ti>>(
      {std::vector<Ti>({0, 1, 2}), std::vector<Ti>({3}),
       std::vector<Ti>({4, 5, 6}), std::vector<Ti>({7, 8})});
  auto searchOptions = SearchOptions();
  auto measures = SearchMeasureOptions();
  auto checks = SearchModelChecks();

  auto s1 = SearcherTest(searchOptions, searchItems, measures, checks, 3, gi, 0,
                         searchItems.LengthTargets);
  auto sers = std::vector<Searcher *>({&s1});
  auto ms = ModelSet(sers, gi, searchOptions, searchItems, measures, checks);

  ms.Start(nullptr, nullptr);

  auto list0 = std::vector<SearcherSummary *>();
  auto list1 = std::vector<SearcherSummary *>();
  auto list2 = std::vector<SearcherSummary *>();
  auto d = SearcherModelingInfo();
  ms.CombineInfo(d, list0, list1, list2);

  auto v = std::vector<std::string>();
  for (auto &a : list0.at(0)->All)
    v.push_back(a->Dependents.ToString0('\t', ','));

  ASSERT_EQ((Ti)11, d.SearchedCount);
  ASSERT_EQ(d.ExpectedCount, d.SearchedCount);
  ASSERT_EQ(d.FailsCount.size(), (Ti)0);

  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,3,4") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,3,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,3,5") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,3,8") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,3,6") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,4,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,4,8") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,5,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,5,8") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,6,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,6,8") != v.end());
}

TEST(Searcher_T, indexes4_skip_target) {

  auto searchItems = SearchItems();
  searchItems.KeepAll = true;
  searchItems.LengthEvals = 1;
  searchItems.LengthTargets = 2;
  searchItems.Length1 = 1;

  auto gi = std::vector<std::vector<Ti>>(
      {std::vector<Ti>({0}), std::vector<Ti>({1, 2, 3}), std::vector<Ti>({4}),
       std::vector<Ti>({7, 8})});
  auto searchOptions = SearchOptions();
  auto measures = SearchMeasureOptions();
  auto checks = SearchModelChecks();

  auto s1 = SearcherTest(searchOptions, searchItems, measures, checks, 3, gi, 0,
                         searchItems.LengthTargets);
  auto sers = std::vector<Searcher *>({&s1});
  auto ms = ModelSet(sers, gi, searchOptions, searchItems, measures, checks);

  ms.Start(nullptr, nullptr);

  auto list0 = std::vector<SearcherSummary *>();
  auto list1 = std::vector<SearcherSummary *>();
  auto list2 = std::vector<SearcherSummary *>();
  auto d = SearcherModelingInfo();
  ms.CombineInfo(d, list0, list1, list2);

  auto v = std::vector<std::string>();
  for (auto &a : list0.at(0)->All)
    v.push_back(a->Dependents.ToString0('\t', ','));
  for (auto &a : list0.at(1)->All)
    v.push_back(a->Dependents.ToString0('\t', ','));

  ASSERT_EQ((Ti)13, d.SearchedCount);
  ASSERT_EQ(d.ExpectedCount, d.SearchedCount);
  ASSERT_EQ(d.FailsCount.size(), (Ti)0);

  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,1,4") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,2,4") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,3,4") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,1,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,1,8") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,2,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,2,8") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,3,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,3,8") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,4,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "0,4,8") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "1,4,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "1,4,8") != v.end());
}

TEST(Searcher_T, indexes_fixed1) {

  auto searchItems = SearchItems();
  searchItems.KeepAll = true;
  searchItems.LengthEvals = 1;
  searchItems.LengthTargets = 6;
  searchItems.Length1 = 1;

  auto gi = std::vector<std::vector<Ti>>(
      {std::vector<Ti>({1, 2}), std::vector<Ti>({3}),
       std::vector<Ti>({4, 5, 6}), std::vector<Ti>({7, 8})});

  auto searchOptions = SearchOptions();
  auto measures = SearchMeasureOptions();
  auto checks = SearchModelChecks();

  auto s1 =
      SearcherTest(searchOptions, searchItems, measures, checks, 3, gi, 1);
  auto sers = std::vector<Searcher *>({&s1});
  auto ms = ModelSet(sers, gi, searchOptions, searchItems, measures, checks);

  ms.Start(nullptr, nullptr);

  auto list0 = std::vector<SearcherSummary *>();
  auto list1 = std::vector<SearcherSummary *>();
  auto list2 = std::vector<SearcherSummary *>();
  auto d = SearcherModelingInfo();
  ms.CombineInfo(d, list0, list1, list2);

  auto v = std::vector<std::string>();
  for (auto &a : list0.at(0)->All)
    v.push_back(a->Dependents.ToString0('\t', ','));

  ASSERT_EQ((Ti)22, d.SearchedCount);
  ASSERT_EQ(d.ExpectedCount, d.SearchedCount);
  ASSERT_EQ(d.FailsCount.size(), (Ti)0);

  ASSERT_EQ(false,
            std::find(v.begin(), v.end(), (std::string) "3,4,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "2,4,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "1,6,8") != v.end());
}

TEST(Searcher_T, indexes_fixed2) {

  auto searchItems = SearchItems();
  searchItems.KeepAll = true;
  searchItems.LengthEvals = 1;
  searchItems.LengthTargets = 6;
  searchItems.Length1 = 1;

  auto gi = std::vector<std::vector<Ti>>(
      {std::vector<Ti>({1, 2}), std::vector<Ti>({3}),
       std::vector<Ti>({4, 5, 6}), std::vector<Ti>({7, 8})});

  auto searchOptions = SearchOptions();
  auto measures = SearchMeasureOptions();
  auto checks = SearchModelChecks();

  auto s1 =
      SearcherTest(searchOptions, searchItems, measures, checks, 3, gi, 2);
  auto sers = std::vector<Searcher *>({&s1});
  auto ms = ModelSet(sers, gi, searchOptions, searchItems, measures, checks);

  ms.Start(nullptr, nullptr);

  auto list0 = std::vector<SearcherSummary *>();
  auto list1 = std::vector<SearcherSummary *>();
  auto list2 = std::vector<SearcherSummary *>();
  auto d = SearcherModelingInfo();
  ms.CombineInfo(d, list0, list1, list2);

  auto v = std::vector<std::string>();
  for (auto &a : list0.at(0)->All)
    v.push_back(a->Dependents.ToString0('\t', ','));

  ASSERT_EQ((Ti)10, d.SearchedCount);
  ASSERT_EQ(d.ExpectedCount, d.SearchedCount);
  ASSERT_EQ(d.FailsCount.size(), (Ti)0);

  ASSERT_EQ(false,
            std::find(v.begin(), v.end(), (std::string) "1,4,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "2,3,7") != v.end());
  ASSERT_EQ(true,
            std::find(v.begin(), v.end(), (std::string) "1,3,8") != v.end());
}
