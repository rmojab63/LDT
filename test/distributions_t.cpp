/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "distributions.h"
#include "matrix.h"
#include "statistics.h"
#include <gtest/gtest.h>

using namespace ldt;

TEST(Distributions_T, empirical103_quantile) {

  double probs[103];
  probs[0] = 0.0001;
  probs[1] = 0.001;
  int j = 2;
  for (double d = 0.01; d < 0.999; d += 0.01)
    probs[j++] = d;
  probs[101] = 0.999;
  probs[102] = 0.9999;

  auto normal_cdfs = std::vector<double>(103);
  auto dis = Distribution<DistributionType::kNormal>(0.0, 1.0);
  for (int i = 0; i < 103; i++)
    normal_cdfs[i] = dis.GetQuantile(probs[i]);

  auto emp103 = DistributionEmpirical103(&normal_cdfs[0]);

  ASSERT_NEAR(dis.GetQuantile(0.0001), emp103.GetQuantileApprox(0.0001), 1e-14);
  ASSERT_NEAR(dis.GetQuantile(0.0001 - 0.0000000001),
              emp103.GetQuantileApprox(0.0001 - 0.0000000001), 1e-6);
  ASSERT_NEAR(dis.GetQuantile(0.0001 + 0.0000000001),
              emp103.GetQuantileApprox(0.0001 + 0.0000000001), 1e-6);

  ASSERT_NEAR(dis.GetQuantile(0.9999), emp103.GetQuantileApprox(0.9999), 1e-14);
  ASSERT_NEAR(dis.GetQuantile(0.9999 + 0.0000000001),
              emp103.GetQuantileApprox(0.9999 + 0.0000000001), 1e-6);
  ASSERT_NEAR(dis.GetQuantile(0.9999 - 0.0000000001),
              emp103.GetQuantileApprox(0.9999 - 0.0000000001), 1e-6);

  ASSERT_NEAR(dis.GetQuantile(0.001), emp103.GetQuantileApprox(0.001), 1e-14);
  ASSERT_NEAR(dis.GetQuantile(0.25), emp103.GetQuantileApprox(0.25), 1e-14);
  ASSERT_NEAR(dis.GetQuantile(0.5), emp103.GetQuantileApprox(0.5), 1e-14);
  ASSERT_NEAR(dis.GetQuantile(0.75), emp103.GetQuantileApprox(0.75), 1e-14);
  ASSERT_NEAR(dis.GetQuantile(0.999), emp103.GetQuantileApprox(0.999), 1e-14);

  ASSERT_NEAR(dis.GetQuantile(0.2500001), emp103.GetQuantileApprox(0.2500001),
              1e-6);
  ASSERT_NEAR(dis.GetQuantile(0.89999999), emp103.GetQuantileApprox(0.89999999),
              1e-6);
}

TEST(Distributions_T, empirical103_cdf) {

  double probs[103];
  probs[0] = 0.0001;
  probs[1] = 0.001;
  int j = 2;
  for (double d = 0.01; d < 0.999; d += 0.01)
    probs[j++] = d;
  probs[101] = 0.999;
  probs[102] = 0.9999;

  auto normal_qs = std::vector<double>(103);
  auto dis = Distribution<DistributionType::kNormal>(0.0, 1.0);
  for (int i = 0; i < 103; i++)
    normal_qs[i] = dis.GetQuantile(probs[i]);

  auto emp103 = DistributionEmpirical103(&normal_qs[0]);

  ASSERT_NEAR(dis.GetCdf(normal_qs[0]), emp103.GetCDFApprox(normal_qs[0]),
              1e-14);
  ASSERT_NEAR(dis.GetCdf(normal_qs[0] - 0.0000000001),
              emp103.GetCDFApprox(normal_qs[0] - 0.0000000001), 1e-6);
  ASSERT_NEAR(dis.GetCdf(normal_qs[0] + 0.0000000001),
              emp103.GetCDFApprox(normal_qs[0] + 0.0000000001), 1e-6);

  ASSERT_NEAR(dis.GetCdf(normal_qs[102]), emp103.GetCDFApprox(normal_qs[102]),
              1e-14);
  ASSERT_NEAR(dis.GetCdf(normal_qs[102] + 0.0000000001),
              emp103.GetCDFApprox(normal_qs[102] + 0.0000000001), 1e-6);
  ASSERT_NEAR(dis.GetCdf(normal_qs[102] - 0.0000000001),
              emp103.GetCDFApprox(normal_qs[102] - 0.0000000001), 1e-6);

  ASSERT_NEAR(dis.GetCdf(normal_qs[1]), emp103.GetCDFApprox(normal_qs[1]),
              1e-14);
  ASSERT_NEAR(dis.GetCdf(normal_qs[1] - 0.01),
              emp103.GetCDFApprox(normal_qs[1] - 0.01), 1e-4);
  ASSERT_NEAR(dis.GetCdf(normal_qs[1] + 0.01),
              emp103.GetCDFApprox(normal_qs[1] + 0.01), 1e-4);

  ASSERT_NEAR(dis.GetCdf(normal_qs[25]), emp103.GetCDFApprox(normal_qs[25]),
              1e-14);
  ASSERT_NEAR(dis.GetCdf(normal_qs[52]), emp103.GetCDFApprox(normal_qs[52]),
              1e-14);
  ASSERT_NEAR(dis.GetCdf(normal_qs[75]), emp103.GetCDFApprox(normal_qs[75]),
              1e-14);
  ASSERT_NEAR(dis.GetCdf(normal_qs[101]), emp103.GetCDFApprox(normal_qs[101]),
              1e-14);

  ASSERT_NEAR(dis.GetCdf(normal_qs[101] + 0.001),
              emp103.GetCDFApprox(normal_qs[101] + 0.001), 1e-5);
  ASSERT_NEAR(dis.GetCdf(normal_qs[101] - 0.001),
              emp103.GetCDFApprox(normal_qs[101] - 0.001), 1e-5);

  ASSERT_NEAR(dis.GetCdf(2.434), emp103.GetCDFApprox(2.434), 1e-2);
  ASSERT_NEAR(dis.GetCdf(0.89999999), emp103.GetCDFApprox(0.89999999), 1e-4);
  ASSERT_NEAR(dis.GetCdf(0.12), emp103.GetCDFApprox(0.12), 1e-5);
}

TEST(Distributions_T, empirical103_combine) {

  double probs[103];
  probs[0] = 0.0001;
  probs[1] = 0.001;
  int j = 2;
  for (double d = 0.01; d < 0.999; d += 0.01)
    probs[j++] = d;
  probs[101] = 0.999;
  probs[102] = 0.9999;

  auto normal_qs1 = std::vector<double>(103);
  auto normal_qs2 = std::vector<double>(103);
  auto dis1 = Distribution<DistributionType::kNormal>(0.0, 1.0);
  auto dis2 = Distribution<DistributionType::kNormal>(0.0, 1.0);
  for (int i = 0; i < 103; i++) {
    normal_qs1[i] = dis1.GetQuantile(probs[i]);
    normal_qs2[i] = dis2.GetQuantile(probs[i]);
  }

  auto emp103_1 = DistributionEmpirical103(&normal_qs1[0]);
  auto emp103_2 = DistributionEmpirical103(&normal_qs2[0]);

  auto dists = std::vector<DistributionEmpirical103>({emp103_1, emp103_2});
  auto weights = std::vector<Tv>();
  auto result = std::vector<double>(103);
  DistributionEmpirical103::Combine(dists, nullptr, &result[0]);
  auto res_dis = DistributionEmpirical103(&result[0]);

  // mixture
  auto a_weights = std::vector<double>({0.5, 0.5});
  auto a_dists = std::vector<DistributionBase *>({&dis1, &dis2});
  auto mix = DistributionMixture(a_weights, a_dists);

  ASSERT_NEAR(mix.GetCdf(-1.0), res_dis.GetCDFApprox(-1.0), 1e-4);
  ASSERT_NEAR(mix.GetCdf(1.0), res_dis.GetCDFApprox(1.0), 1e-4);
  ASSERT_NEAR(mix.GetCdf(-2.0), res_dis.GetCDFApprox(-2.0), 1e-3);
  ASSERT_NEAR(mix.GetCdf(2.0), res_dis.GetCDFApprox(2.0), 1e-3);
}

TEST(Distributions_T, mnormal_sample) {
  Ti N = 100000;
  auto storage = Matrix<Tv>(NAN, new Tv[2 * N], (Ti)N, (Ti)2);
  Tv W[2 * 2 + 4];

  // identity var
  auto mn1 = NormalM(2, new Matrix<Tv>(new Tv[2]{0, 0}, 2), nullptr, N, true,
                     false, 0.0, true, 1.0, true, 0.0);
  mn1.GetSample(storage.Data, W, 99);
  auto varm = Matrix<Tv>(0.0, new Tv[4]{0, 0, 0, 0}, (Ti)2, (Ti)2);
  auto sv = std::vector<Ti>();
  storage.ColumnsVariance(varm, sv, false);
  ASSERT_NEAR(varm.Get(0, 0), 1.0, 1e-2); // variance 1
  ASSERT_NEAR(varm.Get(0, 1), 0.0, 1e-2); // covariance 0

  // general
  auto var = Matrix<Tv>(new Tv[4]{1.04, 4.04, 4.04, 16.04}, 2, 2);
  auto var_copy = Matrix<Tv>(new Tv[4]{1.04, 4.04, 4.04, 16.04}, 2, 2);
  auto mn2 = NormalM(2, new Matrix<Tv>(new Tv[2]{0, 0}, (Ti)2), &var,
                     N); // var gets destroyed
  mn2.GetSample(storage.Data, W, 340);
  auto su = std::vector<Ti>();
  storage.ColumnsVariance(varm, su, false);
  ASSERT_EQ(varm.Equals(var_copy, 1e-1), true);
}

TEST(Distributions_T, mnormal_density) {

  // density
  Tv WW[3 + 4];
  //     zero variance%
  auto mn3 = NormalM(2, new Matrix<Tv>(new Tv[2]{1, 1}, 2), nullptr, 0, true,
                     false, 0, true, 0, true, 0);
  auto x = Matrix<Tv>(new Tv[6]{1, 2, 1, 1, 2, 2}, 2, 3);
  auto storage = Matrix<Tv>(-10.0, new Tv[3]{0, 0, 0}, 3);
  mn3.GetDensity(&x, &storage, WW, false);
  ASSERT_EQ(storage.Get(0), 0.0);
  ASSERT_EQ(true, std::isinf(storage.Get(1)));
  ASSERT_EQ(storage.Get(2), 0.0);
  //     log

  mn3.GetDensity(&x, &storage, WW, true);
  ASSERT_EQ(true, std::isinf(-storage.Get(0)));
  ASSERT_EQ(true, std::isinf(storage.Get(1)));
  ASSERT_EQ(true, std::isinf(-storage.Get(2)));

  //     identity variance
  auto mn4 = NormalM(2, new Matrix<Tv>(new Tv[2]{1, 1}, 2), nullptr, 0, true,
                     false, 0, true, 1, true, 0);
  x = Matrix<Tv>(new Tv[6]{1, 2, 1, 1, 2, 3}, 2, 3);
  storage = Matrix<Tv>(-10.0, new Tv[3]{0, 0, 0}, 3);
  mn4.GetDensity(&x, &storage, WW, false);
  ASSERT_EQ(
      storage.Equals(Matrix<Tv>(new Tv[3]{0.096532, 0.159154, 0.013064}, 3, 1),
                     1e-5),
      true);
  //     log
  x = Matrix<Tv>(new Tv[6]{1, 2, 1, 1, 2, 3}, 2, 3);
  mn4.GetDensity(&x, &storage, WW, true);
  ASSERT_EQ(
      storage.Equals(
          Matrix<Tv>(new Tv[3]{-2.337877, -1.8378770, -4.337877}, 3, 1), 1e-5),
      true);

  //     general variance
  auto P = Matrix<Tv>(new Tv[4]{1, 0.2, 0.2, 4}, 2, 2);
  auto mn5 = NormalM(2, new Matrix<Tv>(new Tv[2]{1, 1}, 2), &P, 0);
  x = Matrix<Tv>(new Tv[6]{1, 2, 1, 1, 2, 3}, 2, 3);
  storage = Matrix<Tv>(-10.0, new Tv[3], 3);
  mn5.GetDensity(&x, &storage, WW, false);
  ASSERT_EQ(
      storage.Equals(
          Matrix<Tv>(new Tv[3]{0.070491, 0.0799783, 0.0322225}, 3, 1), 1e-5),
      true);
  //     log
  x = Matrix<Tv>(new Tv[6]{1, 2, 1, 1, 2, 3}, 2, 3);
  mn5.GetDensity(&x, &storage, WW, true);
  ASSERT_EQ(
      storage.Equals(
          Matrix<Tv>(new Tv[3]{-2.6522617, -2.5259990, -3.43508998}, 3, 1),
          1e-5),
      true);
}

TEST(Distributions_T, normal) {
  auto dis = Distribution<DistributionType::kNormal>(0.0, 1.0);
  ASSERT_EQ(0.0, dis.GetMean());
  ASSERT_EQ(1.0, dis.GetVariance());
  ASSERT_EQ(0.0, dis.GetMedian());
  ASSERT_EQ(0.0, dis.GetMode());
  ASSERT_EQ(0.0, dis.GetSkewness());
  ASSERT_EQ(0.0, dis.GetKurtosis());

  ASSERT_NEAR(0.24197072451914337, dis.GetPdfOrPmf(1.0), 1e-15);
  ASSERT_NEAR(0.0, dis.GetPdfOrPmf(INFINITY), 1e-15);
  ASSERT_NEAR(0.0, dis.GetPdfOrPmf(-INFINITY), 1e-15);
  ASSERT_NEAR(std::log(0.24197072451914337), dis.GetPdfOrPmfLog(1.0), 1e-15);

  ASSERT_NEAR(0.5, dis.GetCdf(0.0), 1e-15);
  ASSERT_NEAR(0.8413447460685429, dis.GetCdf(1.0), 1e-15);
  ASSERT_NEAR(1.0, dis.GetCdf(INFINITY), 1e-15);
  ASSERT_NEAR(0.0, dis.GetCdf(-INFINITY), 1e-15);

  ASSERT_EQ(true, std::isinf(-dis.GetQuantile(0.0)));
  ASSERT_EQ(true, std::isinf(dis.GetQuantile(1.0)));
  ASSERT_NEAR(0.0, dis.GetQuantile(0.5), 1e-15);

  dis = Distribution<DistributionType::kNormal>(1.0, 2.0);
  Ti N = 10000;
  auto sample = new Tv[N];
  dis.GetSample(sample, N, 555);
  auto mat = Matrix<Tv>(sample, N, 1);
  auto desc = Descriptive(&mat);
  auto mv = desc.MeanVariance(false);
  ASSERT_NEAR(1.0, std::get<0>(mv), 1e-2);
  ASSERT_NEAR(4.0, std::get<1>(mv), 0.2);

  // export
  /*
  Tv args[4]{ 1.0,2.0,NAN,NAN };
  dist_random(N, sample, args, 1, 555);
  mat = Matrix<Tv>(sample, N, 1);
  desc = Descriptive(&mat);
  auto mv1 = desc.MeanVariance(false);
  ASSERT_NEAR(std::get<0>(mv), std::get<0>(mv1), 1e-16);
  ASSERT_NEAR(std::get<1>(mv), std::get<1>(mv1), 1e-16);
  */

  // bounds
  N = 10;
  dis = Distribution<DistributionType::kNormal>(INFINITY, 2.0);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kNormal>(-INFINITY, 2.0);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kNormal>(0.0, INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kNormal>(INFINITY, INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);
}

TEST(Distributions_T, uniform) {
  auto dis = Distribution<DistributionType::kUniformCon>(-1.0, 1.0);
  ASSERT_EQ(0.0, dis.GetMean());
  ASSERT_EQ(1.0 / 3.0, dis.GetVariance());
  ASSERT_EQ(0.0, dis.GetMedian());
  ASSERT_EQ(0.0, dis.GetMode());
  ASSERT_EQ(0.0, dis.GetSkewness());
  ASSERT_EQ(-1.2, dis.GetKurtosis());

  ASSERT_NEAR(0.5, dis.GetPdfOrPmf(1.0), 1e-15);
  ASSERT_NEAR(0.0, dis.GetPdfOrPmf(INFINITY), 1e-15);
  ASSERT_NEAR(0.0, dis.GetPdfOrPmf(-INFINITY), 1e-15);
  ASSERT_NEAR(std::log(0.5), dis.GetPdfOrPmfLog(1.0), 1e-15);

  ASSERT_NEAR(0.25, dis.GetCdf(-0.5), 1e-15);
  ASSERT_NEAR(0.5, dis.GetCdf(0.0), 1e-15);
  ASSERT_NEAR(1.0, dis.GetCdf(1.0), 1e-15);
  ASSERT_NEAR(1.0, dis.GetCdf(INFINITY), 1e-15);
  ASSERT_NEAR(0.0, dis.GetCdf(-1.0), 1e-15);

  ASSERT_NEAR(0.0, dis.GetQuantile(0.5), 1e-15);
  ASSERT_NEAR(-0.5, dis.GetQuantile(0.25), 1e-15);
  ASSERT_NEAR(-1, dis.GetQuantile(0.0), 1e-15);
  ASSERT_NEAR(1, dis.GetQuantile(1.0), 1e-15);

  Ti N = 10000;
  auto sample = new Tv[N];
  dis.GetSample(sample, N, 555);
  auto mat = Matrix<Tv>(sample, N, 1);
  auto desc = Descriptive(&mat);
  auto mv = desc.MeanVariance(false);
  ASSERT_NEAR(0.0, std::get<0>(mv), 1e-2);
  ASSERT_NEAR(1.0 / 3.0, std::get<1>(mv), 1e-2);

  // bounds
  N = 10;
  dis = Distribution<DistributionType::kUniformCon>(2.0, INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kUniformCon>(-INFINITY, 2.0);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  // dis.GetSample(sample, N, 555);    throws exception ?!

  dis = Distribution<DistributionType::kUniformCon>(0.0, INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kUniformCon>(-INFINITY, INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  // dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kUniformCon>(INFINITY, INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kUniformCon>(2.0, 2.0);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);
}

TEST(Distributions_T, chi2) {
  auto dis = Distribution<DistributionType::kChi2>(10.0);
  ASSERT_EQ(10.0, dis.GetMean());
  ASSERT_EQ(20.0, dis.GetVariance());
  // ASSERT_NEAR(0.0, dis.GetMedian());
  ASSERT_EQ(8.0, dis.GetMode());
  ASSERT_NEAR(0.894427, dis.GetSkewness(), 1e-5);
  ASSERT_EQ(1.2, dis.GetKurtosis());

  ASSERT_NEAR(0.045111761078870897298, dis.GetPdfOrPmf(4.0), 1e-15);
  ASSERT_NEAR(std::log(0.045111761078870897298), dis.GetPdfOrPmfLog(4.0),
              1e-15);

  ASSERT_NEAR(0.05265, dis.GetCdf(4), 1e-5);
  ASSERT_NEAR(0.0, dis.GetCdf(0.0), 1e-15);
  ASSERT_NEAR(1.0, dis.GetCdf(INFINITY), 1e-15);

  ASSERT_NEAR(0.0, dis.GetQuantile(0.0), 1e-15);
  ASSERT_EQ(true, std::isinf(dis.GetQuantile(1.0)));
  ASSERT_NEAR(3.94030, dis.GetQuantile(0.05), 1e-5);

  Ti N = 100000;
  auto sample = new Tv[N];
  dis.GetSample(sample, N, 340);
  auto mat = Matrix<Tv>(sample, N, 1);
  auto desc = Descriptive(&mat);
  auto mv = desc.MeanVariance(false);
  ASSERT_NEAR(10.0, std::get<0>(mv), 1e-1);
  ASSERT_NEAR(20.0, std::get<1>(mv), 1e-1);

  // export
  Ti M = 10;
  auto result = new Tv[M];
  auto Xs = new Tv[M];
  Tv args[4]{10.0, NAN, NAN, NAN};
  // dist_function(result, M, Xs, 2, args, 0, NAN, 0.05, NAN, 0.95);

  // special params

  // bounds
  N = 10;
  /*
  dis = Distribution<DistributionType::kChi2>( 0.0);
  auto p = dis.GetPdfOrPmf(2);
  auto c = dis.GetCdf(0.0);
  auto q = dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kChi2>( INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);*/
}

TEST(Distributions_T, t) {
  auto dis = Distribution<DistributionType::kT>(10.0);
  ASSERT_EQ(0.0, dis.GetMean());
  ASSERT_EQ(1.25, dis.GetVariance());
  ASSERT_EQ(0.0, dis.GetMedian());
  ASSERT_EQ(0.0, dis.GetMode());
  ASSERT_NEAR(0.0, dis.GetSkewness(), 1e-5);
  ASSERT_EQ(1.0, dis.GetKurtosis());

  ASSERT_NEAR(0.3703984615527454378394, dis.GetPdfOrPmf(0.3), 1e-15);
  ASSERT_NEAR(0.0, dis.GetPdfOrPmf(INFINITY), 1e-15);
  ASSERT_NEAR(0.0, dis.GetPdfOrPmf(-INFINITY), 1e-15);
  ASSERT_NEAR(std::log(0.3703984615527454378394), dis.GetPdfOrPmfLog(0.3),
              1e-15);

  ASSERT_NEAR(0.6148396962171006977095, dis.GetCdf(0.3), 1e-5);
  ASSERT_NEAR(0.5, dis.GetCdf(0.0), 1e-15);
  ASSERT_NEAR(1.0, dis.GetCdf(INFINITY), 1e-15);
  ASSERT_NEAR(0.1704465661510299363374, dis.GetCdf(-1.0), 1e-5);

  ASSERT_EQ(true, std::isinf(dis.GetQuantile(1.0)));
  ASSERT_EQ(true, std::isinf(dis.GetQuantile(0.0)));
  ASSERT_EQ(true, dis.GetQuantile(0.0) < 0);
  ASSERT_EQ(true, dis.GetQuantile(1.0) > 0);
  ASSERT_NEAR(0.0, dis.GetQuantile(0.5), 1e-10);
  ASSERT_NEAR(0.3, dis.GetQuantile(0.6148396962171006977095), 1e-10);
  ASSERT_NEAR(-1.0, dis.GetQuantile(0.1704465661510299363374), 1e-10);

  // more on GetCdf/inv
  dis = Distribution<DistributionType::kT>(1.0);
  ASSERT_NEAR(0.1591549430918953357689, dis.GetPdfOrPmf(1.0), 1e-15);
  ASSERT_NEAR(std::log(0.1591549430918953357689), dis.GetPdfOrPmfLog(1.0),
              1e-15);
  ASSERT_NEAR(0.75, dis.GetCdf(1.0), 1e-15);
  ASSERT_NEAR(1.0, dis.GetQuantile(0.75), 1e-15);

  dis = Distribution<DistributionType::kT>(2.0);
  ASSERT_NEAR(0.7886751345948128822546, dis.GetCdf(1.0), 1e-15);
  ASSERT_NEAR(1.0, dis.GetQuantile(0.7886751345948128822546), 1e-15);

  dis = Distribution<DistributionType::kT>(3.0);
  ASSERT_NEAR(0.8044988905221146790445, dis.GetCdf(1.0), 1e-15);
  ASSERT_NEAR(1.0, dis.GetQuantile(0.8044988905221146790445), 1e-15);

  dis = Distribution<DistributionType::kT>(4.0);
  ASSERT_NEAR(0.8130495168499705574973, dis.GetCdf(1.0), 1e-15);
  ASSERT_NEAR(1.0, dis.GetQuantile(0.8130495168499705574973), 1e-15);

  dis = Distribution<DistributionType::kT>(5.0);
  ASSERT_NEAR(0.8183912661754386871999, dis.GetCdf(1.0), 1e-15);
  ASSERT_NEAR(1.0, dis.GetQuantile(0.8183912661754386871999), 1e-15);

  dis = Distribution<DistributionType::kT>(6.0);
  ASSERT_NEAR(0.8220411581252089129838, dis.GetCdf(1.0), 1e-15);
  ASSERT_NEAR(1.0, dis.GetQuantile(0.8220411581252089129838), 1e-15);

  Ti N = 10000;
  auto sample = new Tv[N];
  dis.GetSample(sample, N, 555);
  auto mat = Matrix<Tv>(sample, N, 1);
  auto desc = Descriptive(&mat);
  auto mv = desc.MeanVariance(false);
  ASSERT_NEAR(dis.GetMean(), std::get<0>(mv), 1e-1);
  ASSERT_NEAR(dis.GetVariance(), std::get<1>(mv), 1e-1);

  // bounds
  N = 10;
  /*
  dis = Distribution<DistributionType::kT>( 0.0);
  auto p = dis.GetPdfOrPmf(2);
  auto c = dis.GetCdf(0.0);
  auto q = dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kT>( INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);
  */
}

TEST(Distributions_T, exponential) {
  auto dis = Distribution<DistributionType::kExponential>(10.0);
  ASSERT_EQ(0.1, dis.GetMean());
  ASSERT_EQ(0.01, dis.GetVariance());
  ASSERT_NEAR(0.0693147, dis.GetMedian(), 1e-5);
  ASSERT_EQ(0.0, dis.GetMode());
  ASSERT_EQ(2.0, dis.GetSkewness());
  ASSERT_EQ(6.0, dis.GetKurtosis());

  ASSERT_NEAR(0.49787068367863946, dis.GetPdfOrPmf(0.3), 1e-15);
  ASSERT_NEAR(std::log(0.49787068367863946), dis.GetPdfOrPmfLog(0.3), 1e-15);
  ASSERT_NEAR(0.95021293163213605, dis.GetCdf(0.3), 1e-5);

  ASSERT_NEAR(0.3, dis.GetQuantile(0.95021293163213605), 1e-10);

  Ti N = 10000;
  auto sample = new Tv[N];
  dis.GetSample(sample, N, 555);
  auto mat = Matrix<Tv>(sample, N, 1);
  auto desc = Descriptive(&mat);
  auto mv = desc.MeanVariance(false);
  ASSERT_NEAR(0.1, std::get<0>(mv), 1e-1);
  ASSERT_NEAR(0.01, std::get<1>(mv), 1e-1);

  // bounds
  N = 10;
  /*
  dis = Distribution<DistributionType::kExponential>( 0.0);
  auto p = dis.GetPdfOrPmf(2);
  auto c = dis.GetCdf(0.0);
  auto q = dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kT>( INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);
  */
}

TEST(Distributions_T, lognormal) {
  auto dis = Distribution<DistributionType::kLogNormal>(0.1, 0.04);
  ASSERT_NEAR(1.10606, dis.GetMean(), 1e-4);
  ASSERT_NEAR(0.00195894, dis.GetVariance(), 1e-4);
  ASSERT_NEAR(1.10517, dis.GetMedian(), 1e-4);
  ASSERT_NEAR(1.1034, dis.GetMode(), 1e-4);
  ASSERT_NEAR(0.120112, dis.GetSkewness(), 1);
  ASSERT_NEAR(0.025659, dis.GetKurtosis(), 1e-4);

  dis = Distribution<DistributionType::kLogNormal>(2, 4);

  ASSERT_NEAR(0.088016331691074867, dis.GetPdfOrPmf(1.0), 1e-15);
  ASSERT_NEAR(std::log(0.088016331691074867), dis.GetPdfOrPmfLog(1.0), 1e-15);
  ASSERT_NEAR(0.30853753872598694, dis.GetCdf(1.0), 1e-5);

  ASSERT_NEAR(1.0, dis.GetQuantile(0.30853753872598694), 1e-10);

  Ti N = 10000;
  auto sample = new Tv[N];
  dis.GetSample(sample, N, 555);
  auto mat = Matrix<Tv>(sample, N, 1);
  auto desc = Descriptive(&mat);
  auto mv = desc.MeanVariance(false);
  // error is too large. test is trivial
  // ASSERT_NEAR(dis.GetMean(), std::get<0>(mv), 1000);
  // ASSERT_NEAR(dis.GetVariance(), std::get<1>(mv), 1000);

  // bounds
  N = 10;
  /*
  dis = Distribution<DistributionType::kLogNormal>( 2.0, 0.0);
  auto p = dis.GetPdfOrPmf(2);
  auto c = dis.GetCdf(0.0);
  auto q = dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kLogNormal>( -2.0, INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);
  */
}

TEST(Distributions_T, f) {
  auto dis = Distribution<DistributionType::kF>(9, 10);
  ASSERT_NEAR(5.0 / 4.0, dis.GetMean(), 1e-4);
  ASSERT_NEAR(425.0 / 432.0, dis.GetVariance(), 1e-4);
  // ASSERT_NEAR(0.992321, dis.GetMedian(), 1e-4);
  ASSERT_NEAR(70.0 / (12.0 * 9.0), dis.GetMode(), 1e-4);
  ASSERT_NEAR(26.0 / std::sqrt(51.0), dis.GetSkewness(), 1e-4);
  ASSERT_NEAR(778.0 / 17.0, dis.GetKurtosis(), 1e-4);

  dis = Distribution<DistributionType::kF>(3, 5);
  ASSERT_NEAR(0.39924120425900500, dis.GetPdfOrPmf(0.90), 1e-15);
  ASSERT_NEAR(std::log(0.39924120425900500), dis.GetPdfOrPmfLog(0.9), 1e-15);

  ASSERT_NEAR(0.5, dis.GetCdf(0.9071462198190186999416), 1e-5);
  ASSERT_NEAR(0.9071462198190186999416, dis.GetQuantile(0.5), 1e-10);

  dis = Distribution<DistributionType::kF>(5, 2);
  ASSERT_NEAR(0.6, dis.GetCdf(1.764421462571285798743), 1e-5);
  ASSERT_NEAR(1.764421462571285798743, dis.GetQuantile(0.6), 1e-10);

  // bounds
  auto N = 10;
  auto sample = new Tv[N];
  /*
  dis = Distribution<DistributionType::kF>( 0.0, 0.0);
  auto p = dis.GetPdfOrPmf(2);
  auto c = dis.GetCdf(0.0);
  auto q = dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kF>( 0.0, 1.0);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kF>( 2.0, 0.0);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kF>( INFINITY, INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kF>( 1.0, INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kF>( INFINITY, 3.0);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);*/
}

TEST(Distributions_T, gamma) {
  auto dis = Distribution<DistributionType::kGamma>(4, 0.5);
  ASSERT_NEAR(2.0, dis.GetMean(), 1e-4);
  ASSERT_NEAR(1.0, dis.GetVariance(), 1e-4);
  // ASSERT_NEAR(1.83603, dis.GetMedian(), 1e-4);
  ASSERT_NEAR(1.5, dis.GetMode(), 1e-4);
  ASSERT_NEAR(1.0, dis.GetSkewness(), 1e-4);
  ASSERT_NEAR(1.5, dis.GetKurtosis(), 1e-4);

  ASSERT_NEAR(0.3907336296263291795993, dis.GetPdfOrPmf(2.0), 1e-15);
  ASSERT_NEAR(std::log(0.3907336296263291795993), dis.GetPdfOrPmfLog(2.0),
              1e-15);
  ASSERT_NEAR(0.566529879633291066382, dis.GetCdf(2.0), 1e-5);
  ASSERT_NEAR(2.0, dis.GetQuantile(0.566529879633291066382), 1e-10);

  dis = Distribution<DistributionType::kGamma>(5, 3);
  ASSERT_NEAR(0.01178192253733708830453, dis.GetCdf(4.0), 1e-5);
  ASSERT_NEAR(4.0, dis.GetQuantile(0.01178192253733708830453), 1e-10);

  Ti N = 1000000;
  auto sample = new Tv[N];
  dis.GetSample(sample, N, 66);
  auto mat = Matrix<Tv>(sample, N, 1);
  auto desc = Descriptive(&mat);
  auto mv = desc.MeanVariance(false);
  ASSERT_NEAR(dis.GetMean(), std::get<0>(mv), 1e-1);
  ASSERT_NEAR(dis.GetVariance(), std::get<1>(mv), 1e-1);

  // bounds
  N = 10;
  /*
  dis = Distribution<DistributionType::kGamma>( 0.0, 1.0);
  auto p = dis.GetPdfOrPmf(2);
  auto c = dis.GetCdf(0.0);
  auto q = dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kGamma>( 1.0, 0.0);
  p = dis.GetPdfOrPmf(2);
  c = dis.GetCdf(0.0);
  q = dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kGamma>( 0.0, 0.0);
  p = dis.GetPdfOrPmf(2);
  c = dis.GetCdf(0.0);
  q = dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kGamma>( 2.0, INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);
  */
}

TEST(Distributions_T, beta) {
  auto dis = Distribution<DistributionType::kBeta>(2, 3);
  ASSERT_NEAR(0.4, dis.GetMean(), 1e-15);
  ASSERT_NEAR(0.04, dis.GetVariance(), 1e-14);
  // ASSERT_NEAR(1.83603, dis.GetMedian(), 1e-4);
  ASSERT_NEAR(1.0 / 3.0, dis.GetMode(), 1e-15);
  ASSERT_NEAR(0.285714, dis.GetSkewness(), 1e-4);
  ASSERT_NEAR(-0.642857, dis.GetKurtosis(), 1e-4);

  ASSERT_NEAR(1.536, dis.GetPdfOrPmf(0.2), 1e-15);
  ASSERT_NEAR(std::log(1.536), dis.GetPdfOrPmfLog(0.2), 1e-15);
  ASSERT_NEAR(0.1808, dis.GetCdf(0.2), 1e-5);
  ASSERT_NEAR(0.2, dis.GetQuantile(0.1808), 1e-10);

  dis = Distribution<DistributionType::kBeta>(0.5, 0.6);
  ASSERT_NEAR(0.3317880199648356919427, dis.GetCdf(0.2), 1e-5);
  ASSERT_NEAR(0.2, dis.GetQuantile(0.3317880199648356919427), 1e-10);

  /*
          Ti N = 10000;
          auto sample = new Tv[N];
          dis.GetSample(sample, N, 555);
          auto mat = Matrix<Tv>(sample, N, 1);
          auto desc = Descriptive(&mat);
          auto mv = desc.MeanVariance(false);
          ASSERT_NEAR(dis.GetMean(), std::get<0>(mv), 1e-1);
          ASSERT_NEAR(dis.GetVariance(), std::get<1>(mv), 1e-1);


          // bounds
          N = 10;
  */
  /*dis = Distribution<DistributionType::kBeta>( 0.0, 1.0);
  auto p = dis.GetPdfOrPmf(2);
  auto c = dis.GetCdf(0.0);
  auto q = dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kBeta>( 1.0, 0.0);
  auto p = dis.GetPdfOrPmf(2);
  auto c = dis.GetCdf(0.0);
  auto q = dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kBeta>( 0.0, 0.0);
  p = dis.GetPdfOrPmf(2);
  c = dis.GetCdf(0.0);
  q = dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);

  dis = Distribution<DistributionType::kBeta>( 2.0, INFINITY);
  dis.GetPdfOrPmf(2);
  dis.GetCdf(0.0);
  dis.GetQuantile(0.33);
  dis.GetSample(sample, N, 555);*/
}

TEST(Distributions_T, duniform) {
  auto dis = Distribution<DistributionType::kUniformDis>(10, 100);
  ASSERT_NEAR(55, dis.GetMean(), 1e-15);
  ASSERT_NEAR(690, dis.GetVariance(), 1e-15);
  ASSERT_NEAR(-4141.0 / 3450.0, dis.GetKurtosis(), 1e-15);
  ASSERT_NEAR(55, dis.GetMedian(), 1e-15);
  ASSERT_NEAR(0.010989010989010990, dis.GetPdfOrPmf(20), 1e-15);
  ASSERT_NEAR(std::log(0.010989010989010990), dis.GetPdfOrPmfLog(20), 1e-15);
  ASSERT_NEAR(0.12087912087912088, dis.GetCdf(20), 1e-5);
  ASSERT_NEAR(20, dis.GetQuantile(0.12087912087912088), 1e-10);

  ASSERT_NEAR(0.98901098901098905, dis.GetCdf(99), 1e-15);
  ASSERT_NEAR(99, dis.GetQuantile(0.98901098901098905), 1e-15);

  Ti N = 100000;
  auto sample = new Tv[N];
  dis.GetSample(sample, N, 55);
  auto mat = Matrix<Tv>(sample, N, 1);
  auto desc = Descriptive(&mat);
  auto mv = desc.MeanVariance(false);
  ASSERT_NEAR(dis.GetMean(), std::get<0>(mv), 1e-1);
  ASSERT_NEAR(dis.GetVariance(), std::get<1>(mv), 1.0);
}

TEST(Distributions_T, poisson) {
  auto dis = Distribution<DistributionType::kPoisson>(10.0);
  ASSERT_NEAR(10, dis.GetMean(), 1e-15);
  ASSERT_NEAR(10, dis.GetVariance(), 1e-15);
  ASSERT_NEAR(0.316228, dis.GetSkewness(), 1e-6);
  ASSERT_NEAR(0.1, dis.GetKurtosis(), 1e-15);
  ASSERT_NEAR(10, dis.GetMedian(), 1e-15);
  ASSERT_NEAR(0.12511003572113372, dis.GetPdfOrPmf(10), 1e-15);
  ASSERT_NEAR(std::log(0.12511003572113372), dis.GetPdfOrPmfLog(10), 1e-15);
  ASSERT_NEAR(0.58303975019298537, dis.GetCdf(10), 1e-5);
  ASSERT_NEAR(10, dis.GetQuantile(0.58303975019298537), 1e-10);

  ASSERT_NEAR(0.010336050675925690, dis.GetCdf(3), 1e-15);
  ASSERT_NEAR(3, dis.GetQuantile(0.010336050675925690), 1e-15);

  Ti N = 10000;
  auto sample = new Tv[N];
  dis.GetSample(sample, N, 555);
  auto mat = Matrix<Tv>(sample, N, 1);
  auto desc = Descriptive(&mat);
  auto mv = desc.MeanVariance(false);
  ASSERT_NEAR(dis.GetMean(), std::get<0>(mv), 1e-1);
  ASSERT_NEAR(dis.GetVariance(), std::get<1>(mv), 1.0);

  dis = Distribution<DistributionType::kPoisson>(3.0);
  ASSERT_NEAR(3, dis.GetMean(), 1e-15);
  ASSERT_NEAR(3, dis.GetVariance(), 1e-15);
  ASSERT_NEAR(0.33333333, dis.GetKurtosis(), 1e-6);
  ASSERT_NEAR(3, dis.GetMedian(), 1e-15);
  ASSERT_NEAR(0.22404180765538775, dis.GetPdfOrPmf(2), 1e-15);
  ASSERT_NEAR(std::log(0.22404180765538775), dis.GetPdfOrPmfLog(2), 1e-15);

  // export
  Ti M = 10;
  auto result = new Tv[M];
  auto Xs = new Tv[M];
  Tv args[4]{3.0, NAN, NAN, NAN};
  auto min = std::numeric_limits<Tv>().max();
  auto max = std::numeric_limits<Tv>().max();
  // dist_function(result, M, Xs, 2, args, 0, min, 0.05, max, 0.95);
  auto dd = Matrix<Tv>(Xs, M, 1);
  auto dds = dd.ToString();
  auto mm = Matrix<Tv>(result, M, 1);
  auto mms = mm.ToString();

  // bounds
  dis = Distribution<DistributionType::kPoisson>(0.0);
  dis.GetPdfOrPmf(2);
  dis.GetPdfOrPmfLog(1);
  dis.GetCdf(2);
  dis.GetQuantile(0.3);
}

TEST(Distributions_T, geometric) {
  auto dis = Distribution<DistributionType::kGeometric>(0.3);
  ASSERT_NEAR(2.33333333333333, dis.GetMean(), 1e-9);
  ASSERT_NEAR(7.77777777777777, dis.GetVariance(), 1e-9);
  ASSERT_NEAR(1.0, dis.GetMedian(), 1e-9);
  ASSERT_NEAR(6.12857, dis.GetKurtosis(), 1e-5);
  ASSERT_NEAR(2.03189, dis.GetSkewness(), 1e-5);
  ASSERT_NEAR(0.0, dis.GetMode(), 1e-9);
  ASSERT_NEAR(0.035294699999999984, dis.GetPdfOrPmf(6), 1e-15);
  ASSERT_NEAR(std::log(0.035294699999999984), dis.GetPdfOrPmfLog(6), 1e-15);
  ASSERT_NEAR(0.91764570000000001, dis.GetCdf(6), 1e-5);
  ASSERT_NEAR(6, dis.GetQuantile(0.91764570000000001), 1e-10);

  ASSERT_NEAR(0.75990000000000002, dis.GetCdf(3), 1e-15);
  ASSERT_NEAR(3, dis.GetQuantile(0.75990000000000002), 1e-15);
  ASSERT_NEAR(4, dis.GetQuantile(0.75990000000000002 + 1e-10), 1e-15);

  Ti N = 10000;
  auto sample = new Tv[N];
  dis.GetSample(sample, N, 555);
  auto mat = Matrix<Tv>(sample, N, 1);
  auto desc = Descriptive(&mat);
  auto mv = desc.MeanVariance(false);
  ASSERT_NEAR(dis.GetMean(), std::get<0>(mv), 1e-1);
  ASSERT_NEAR(dis.GetVariance(), std::get<1>(mv), 1.0);
}

TEST(Distributions_T, binomial) {
  auto dis = Distribution<DistributionType::kBinomial>(0.3, 100);
  ASSERT_NEAR(30, dis.GetMean(), 1e-9);
  ASSERT_NEAR(21, dis.GetVariance(), 1e-9);
  ASSERT_NEAR(30.0, dis.GetMedian(), 1e-9);
  ASSERT_NEAR(-0.01238, dis.GetKurtosis(), 1e-5);
  ASSERT_NEAR(0.0872872, dis.GetSkewness(), 1e-5);
  ASSERT_NEAR(0.086783864753427600, dis.GetPdfOrPmf(30), 1e-15);
  ASSERT_NEAR(std::log(0.086783864753427600), dis.GetPdfOrPmfLog(30), 1e-10);
  ASSERT_NEAR(0.54912360076878874, dis.GetCdf(30), 1e-5);
  ASSERT_NEAR(30, dis.GetQuantile(0.54912360076878874), 1e-10);

  ASSERT_NEAR(0.98750159283356176, dis.GetCdf(40), 1e-15);
  ASSERT_NEAR(40, dis.GetQuantile(0.98750159283356176), 1e-15);
  ASSERT_NEAR(41, dis.GetQuantile(0.98750159283356176 + 1e-10), 1e-15);

  Ti N = 10000;
  auto sample = new Tv[N];
  dis.GetSample(sample, N, 555);
  auto mat = Matrix<Tv>(sample, N, 1);
  auto desc = Descriptive(&mat);
  auto mv = desc.MeanVariance(false);
  ASSERT_NEAR(dis.GetMean(), std::get<0>(mv), 1e-1);
  ASSERT_NEAR(dis.GetVariance(), std::get<1>(mv), 1.0);
}

TEST(Distributions_T, bernoulli) {
  auto dis = Distribution<DistributionType::kBernoulli>(0.3);
  ASSERT_NEAR(0.30, dis.GetMean(), 1e-9);
  ASSERT_NEAR(0.21, dis.GetVariance(), 1e-9);
  ASSERT_NEAR(0.0, dis.GetMedian(), 1e-9);
  ASSERT_NEAR(-1.2381, dis.GetKurtosis(), 1e-5);
  ASSERT_NEAR(0.872872, dis.GetSkewness(), 1e-5);
  ASSERT_NEAR(0.7, dis.GetPdfOrPmf(0), 1e-15);
  ASSERT_NEAR(std::log(0.7), dis.GetPdfOrPmfLog(0), 1e-10);
  ASSERT_NEAR(0.3, dis.GetPdfOrPmf(1), 1e-15);
  ASSERT_NEAR(std::log(0.3), dis.GetPdfOrPmfLog(1), 1e-10);

  ASSERT_NEAR(0.0, dis.GetCdf(-1), 1e-5);
  ASSERT_NEAR(1.0, dis.GetCdf(1), 1e-5);

  Ti N = 10000;
  auto sample = new Tv[N];
  dis.GetSample(sample, N, 555);
  auto mat = Matrix<Tv>(sample, N, 1);
  auto desc = Descriptive(&mat);
  auto mv = desc.MeanVariance(false);
  ASSERT_NEAR(dis.GetMean(), std::get<0>(mv), 1e-1);
  ASSERT_NEAR(dis.GetVariance(), std::get<1>(mv), 1.0);
};

TEST(Distributions_T, gld_FKML) {

  // none zero
  auto dis = Distribution<DistributionType::kGldFkml>(0.2, 0.3, 0.4, 0.5);
  ASSERT_NEAR(0.04126984126984035, dis.GetMean(), 1e-12);
  ASSERT_NEAR(10.717108899261879, dis.GetVariance(), 1e-12);
  ASSERT_NEAR(-0.11152089036689010, dis.GetSkewness(), 1e-12);
  ASSERT_NEAR(-0.83612928785873786, dis.GetKurtosis(), 1e-12);
  ASSERT_NEAR(-0.8544374353170032, dis.GetQuantile(0.4), 1e-12);
  ASSERT_NEAR(0.09921105528066751, dis.GetDensityQuantile(0.4), 1e-12);

  // both zero
  dis = Distribution<DistributionType::kGldFkml>(0.2, 0.3, 0.0, 0.0);
  ASSERT_NEAR(0.2, dis.GetMean(), 1e-12);
  ASSERT_NEAR(36.554090374405035, dis.GetVariance(), 1e-12);
  ASSERT_NEAR(0.0, dis.GetSkewness(), 1e-12);
  ASSERT_NEAR(1.2, dis.GetKurtosis(), 1e-12);
  ASSERT_NEAR(-1.151550360360548, dis.GetQuantile(0.4), 1e-12);
  ASSERT_NEAR(0.072, dis.GetDensityQuantile(0.4), 1e-12);

  // L3 zero
  dis = Distribution<DistributionType::kGldFkml>(0.2, 0.3, 0.0, 0.5);
  ASSERT_NEAR(-0.9111111111111112, dis.GetMean(), 1e-12);
  ASSERT_NEAR(21.8875744853365610, dis.GetVariance(), 1e-12);
  ASSERT_NEAR(-1.0427032034996768, dis.GetSkewness(), 1e-12);
  ASSERT_NEAR(1.7297501612443469, dis.GetKurtosis(), 1e-12);
  ASSERT_NEAR(-1.351613567857073, dis.GetQuantile(0.4), 1e-12);
  ASSERT_NEAR(0.07913490881002, dis.GetDensityQuantile(0.4), 1e-12);

  // L4 zero
  dis = Distribution<DistributionType::kGldFkml>(0.2, 0.3, 0.4, 0.0);
  ASSERT_NEAR(1.1523809523809525, dis.GetMean(), 1e-12);
  ASSERT_NEAR(23.3921394431432503, dis.GetVariance(), 1e-12);
  ASSERT_NEAR(0.9223339139764826, dis.GetSkewness(), 1e-12);
  ASSERT_NEAR(1.4603763611867260, dis.GetKurtosis(), 1e-12);
  ASSERT_NEAR(-0.6543742278204780, dis.GetQuantile(0.4), 1e-12);
  ASSERT_NEAR(0.08824752484682569, dis.GetDensityQuantile(0.4), 1e-12);

  // inf
  dis = Distribution<DistributionType::kGldFkml>(0.2, 0.3, INFINITY, 0.5);
  ASSERT_NEAR(1.702688871723444, dis.GetQuantile(0.4), 1e-12);
  ASSERT_NEAR(0.232379000772445, dis.GetDensityQuantile(0.4), 1e-12);

  dis = Distribution<DistributionType::kGldFkml>(0.2, 0.3, INFINITY, 0.0);
  ASSERT_NEAR(1.902752079219969, dis.GetQuantile(0.4), 1e-12);
  ASSERT_NEAR(0.18, dis.GetDensityQuantile(0.4), 1e-12);

  dis = Distribution<DistributionType::kGldFkml>(0.2, 0.3, 0.4, INFINITY);
  ASSERT_NEAR(-2.357126307040447, dis.GetQuantile(0.4), 1e-12);
  ASSERT_NEAR(0.1731239887088656, dis.GetDensityQuantile(0.4), 1e-12);

  dis = Distribution<DistributionType::kGldFkml>(0.2, 0.3, 0.0, INFINITY);
  ASSERT_NEAR(-2.854302439580517, dis.GetQuantile(0.4), 1e-12);
  ASSERT_NEAR(0.12, dis.GetDensityQuantile(0.4), 1e-12);
}

TEST(Distributions_T, gld_FKML_Ms) {

  // check Ms
  Tv L3 = 0.4;
  Tv L4 = 0.5;
  Tv m1, m2, m3, m4;
  DistributionGld::GetMs(L3, L4, m1, m2, m3, m4);
  Tv b1 = DistributionGld::GetMk(1, L3, L4);
  Tv b2 = DistributionGld::GetMk(2, L3, L4);
  Tv b3 = DistributionGld::GetMk(3, L3, L4);
  Tv b4 = DistributionGld::GetMk(4, L3, L4);
  ASSERT_NEAR(b1, m1, 1e-16);
  ASSERT_NEAR(b2, m2, 1e-16);
  ASSERT_NEAR(b3, m3, 1e-16);
  ASSERT_NEAR(b4, m4, 1e-16);

  // both zero
  L3 = 0.0;
  L4 = 0.0;
  DistributionGld::GetMs(L3, L4, m1, m2, m3, m4);
  b1 = DistributionGld::GetMk(1, L3, L4);
  b2 = DistributionGld::GetMk(2, L3, L4);
  b3 = DistributionGld::GetMk(3, L3, L4);
  b4 = DistributionGld::GetMk(4, L3, L4);
  ASSERT_NEAR(b1, m1, 1e-16);
  ASSERT_NEAR(b2, m2, 1e-16);
  ASSERT_NEAR(b3, m3, 1e-16);
  ASSERT_NEAR(b4, m4, 1e-16);

  // L3 zero
  L3 = 0.0;
  L4 = 0.5;
  DistributionGld::GetMs(L3, L4, m1, m2, m3, m4);
  b1 = DistributionGld::GetMk(1, L3, L4);
  b2 = DistributionGld::GetMk(2, L3, L4);
  b3 = DistributionGld::GetMk(3, L3, L4);
  b4 = DistributionGld::GetMk(4, L3, L4);
  ASSERT_NEAR(b1, m1, 1e-16);
  ASSERT_NEAR(b2, m2, 1e-16);
  ASSERT_NEAR(b3, m3, 1e-16);
  ASSERT_NEAR(b4, m4, 1e-16);

  // L4 zero
  L3 = 0.4;
  L4 = 0.0;
  DistributionGld::GetMs(L3, L4, m1, m2, m3, m4);
  b1 = DistributionGld::GetMk(1, L3, L4);
  b2 = DistributionGld::GetMk(2, L3, L4);
  b3 = DistributionGld::GetMk(3, L3, L4);
  b4 = DistributionGld::GetMk(4, L3, L4);
  ASSERT_NEAR(b1, m1, 1e-16);
  ASSERT_NEAR(b2, m2, 1e-16);
  ASSERT_NEAR(b3, m3, 1e-16);
  ASSERT_NEAR(b4, m4, 1e-16);
}

TEST(Distributions_T, gld_FKML_Ms_T_from_moments) {

  Ti type;
  Tv mean, variance, skewness, kurtosis;

  // no constrain
  type = 0;
  mean = 0;
  variance = 1;
  skewness = 1.0;
  kurtosis = 100.0;
  auto optim = NelderMead(2);
  auto Ls = DistributionGld::GetFromMoments(mean, variance, skewness, kurtosis,
                                            type, optim);
  auto dis = Distribution<DistributionType::kGldFkml>(
      std::get<0>(Ls), std::get<1>(Ls), std::get<2>(Ls), std::get<3>(Ls));
  ASSERT_NEAR(mean, dis.GetMean(), 1e-5);
  ASSERT_NEAR(variance, dis.GetVariance(), 1e-5);
  ASSERT_NEAR(skewness, dis.GetSkewness(), 1e-3);
  ASSERT_NEAR(kurtosis, dis.GetKurtosis(), 1e-3);
}

TEST(Distributions_T, mixture_sample) {
  auto weight = std::vector<Tv>({0.8, 1.1, 4.0, 2.0, 0.3, 0.5});
  auto dists = std::vector<DistributionBase *>(
      {new Distribution<DistributionType::kNormal>(0.1, 3.0),
       new Distribution<DistributionType::kT>(6.0),
       new Distribution<DistributionType::kGamma>(2, 4),
       new Distribution<DistributionType::kUniformDis>(4, 6),
       new Distribution<DistributionType::kPoisson>(4.0),
       new Distribution<DistributionType::kBinomial>(0.33, 70)});

  auto N = 10000;
  auto dis = DistributionMixture(weight, dists);
  Tv storage[N];
  dis.GetSample(storage, N, 99);
  Matrix<Tv> mat = Matrix<Tv>(storage, N, 1);
  Descriptive desc = Descriptive(&mat);
  auto mean = desc.MeanVariance(true);
  ASSERT_NEAR(std::get<0>(mean), 6.2973532287391443,
              1e-16); // I use it (hard-coded) in c#

  // export
  //	Ti n = 10;
  //	Tv output[10];
  //	Tv args[8]{ 0,1,0,0,2,3,0,0 };
  //	Ti ids[2]{ 1,1 };
  //	Tv weigts[2]{ 1,1 };
  // Ti iscont[2]{ 1,1 };
  // dist_random_mix(n, output, args, ids, 2, weigts, 44);
}

TEST(Distributions_T, mixture_dis_sample) {
  auto disweight = std::vector<Tv>({0.7, 0.3, 0.5});
  auto disdists = std::vector<DistributionBase *>(
      {new Distribution<DistributionType::kUniformDis>(4, 6),
       new Distribution<DistributionType::kPoisson>(4.0),
       new Distribution<DistributionType::kBinomial>(0.33, 70)});
  auto dis = DistributionMixture(disweight, disdists);

  auto N = 10;
  Tv storage[N];

  dis.GetSample(storage, N, 166);

  ASSERT_EQ(storage[0], (Tv)6); // I use it (hard-coded) in c#
  ASSERT_EQ(storage[1], (Tv)24);
  ASSERT_EQ(storage[2], (Tv)18);

  // export
  //	Ti n = 10;
  //	Tv output[10];
  //	Tv args[8]{ 0,1,0,0,2,3,0,0 };
  //	Ti ids[2]{ 1,1 };
  //	Tv weigts[2]{ 1,1 };
  //	Ti iscont[2]{ 1,1 };
  // dist_random_mix(n, output, args, ids, 2, weigts, 44);
}

TEST(Distributions_T, mixture_function) {
  auto weight = std::vector<Tv>({0.8, 1.1, 4.0, 2.0, 0.3, 0.5});
  auto dists = std::vector<DistributionBase *>(
      {new Distribution<DistributionType::kNormal>(0.1, 3.0),
       new Distribution<DistributionType::kT>(6.0),
       new Distribution<DistributionType::kGamma>(2, 4),
       new Distribution<DistributionType::kUniformDis>(4, 6),
       new Distribution<DistributionType::kPoisson>(4.0),
       new Distribution<DistributionType::kBinomial>(0.33, 70)});
  DistributionMixture dis = DistributionMixture(weight, dists);

  auto GetCdf = dis.GetCdf(3);
  auto max = dis.GetMaximum();
  auto min = dis.GetMinimum();

  // export
  //	Ti n = 10;
  //	Tv output[10];
  //	Tv X[10];
  //	Tv args[8]{ 0,1,0,0,2,3,0,0 };
  //	Ti ids[2]{ 1,1 };
  //	Tv weigts[2]{ 1,1 };
  // Ti iscont[2]{ 1,1 };
  // dist_function_mix(output, n, X, args, ids, 2, weigts, 0, NAN, NAN);
  // dist_function_mix(output, n, X, args, ids, 2, weigts, 1, NAN, NAN);
  // dist_function_mix(output, n, X, args, ids, 2, weigts, 2, NAN, NAN);
}

TEST(Distributions_T, mixture_function_pdf) {
  auto weight = std::vector<Tv>({0.8, 1.1, 4.0});
  auto dists = std::vector<DistributionBase *>(
      {new Distribution<DistributionType::kNormal>(0.1, 3.0),
       new Distribution<DistributionType::kT>(6.0),
       new Distribution<DistributionType::kGamma>(2, 4)});
  DistributionMixture dis = DistributionMixture(weight, dists);

  auto pdf = dis.GetPdfOrPmf(10);
  auto pdf1 = dis.GetPdfOrPmf(5);
  auto pdfln = dis.GetPdfOrPmfLog(10);
}

TEST(Distributions_T, mixture_function_pmf_T) {
  auto weight = std::vector<Tv>({2.0, 0.3, 0.5});
  auto dists = std::vector<DistributionBase *>(
      {new Distribution<DistributionType::kUniformDis>(4, 6),
       new Distribution<DistributionType::kPoisson>(4.0),
       new Distribution<DistributionType::kBinomial>(0.33, 70)});
  DistributionMixture dis = DistributionMixture(weight, dists);

  auto pdf = dis.GetPdfOrPmf(10);
  auto pdf1 = dis.GetPdfOrPmf(5);
  auto pdfln = dis.GetPdfOrPmfLog(10);
}

TEST(Distributions_T, mixture_moments) {

  // just one normal
  auto conweight = std::vector<Tv>({1});
  auto condists = std::vector<DistributionBase *>(
      {new Distribution<DistributionType::kNormal>(0.0, 1.0)});
  DistributionMixture dis = DistributionMixture(conweight, condists);

  Tv mean = 0.0;
  Tv variance = 0.0;
  Tv kurtosis = 0.0;
  Tv skewness = 0.0;
  dis.GetMoments(mean, variance, skewness, kurtosis);
  ASSERT_NEAR(0.0, mean, 1e-16);
  ASSERT_NEAR(1.0, variance, 1e-16);
  ASSERT_NEAR(0.0, skewness, 1e-16);
  ASSERT_NEAR(0.0, kurtosis, 1e-16);

  // two similar normals
  conweight = std::vector<Tv>({1, 4});
  condists = std::vector<DistributionBase *>(
      {new Distribution<DistributionType::kNormal>(2.0, 2.0),
       new Distribution<DistributionType::kNormal>(2.0, 2.0)});
  dis = DistributionMixture(conweight, condists);

  dis.GetMoments(mean, variance, skewness, kurtosis);
  ASSERT_NEAR(2.0, mean, 1e-16);
  ASSERT_NEAR(4.0, variance, 1e-16);
  ASSERT_NEAR(0.0, skewness, 1e-16);
  ASSERT_NEAR(0.0, kurtosis, 1e-16);

  // three different normals
  conweight = std::vector<Tv>({1, 4, 6});
  condists = std::vector<DistributionBase *>({
      new Distribution<DistributionType::kNormal>(2.0, 3.0),
      new Distribution<DistributionType::kNormal>(3.0, 4.0),
      new Distribution<DistributionType::kNormal>(-5, 9.0),
  });
  dis = DistributionMixture(conweight, condists);
  dis.GetMoments(mean, variance, skewness, kurtosis);
  // check by the normal version
  Tv mean_n = 0.0;
  Tv variance_n = 0.0;
  Tv kurtosis_n = 0.0;
  Tv skewness_n = 0.0;
  dis.GetMomentsNormal(mean_n, variance_n, skewness_n, kurtosis_n);
  ASSERT_NEAR(mean_n, mean, 1e-16);
  ASSERT_NEAR(variance_n, variance, 1e-16);
  ASSERT_NEAR(skewness_n, skewness, 1e-16);
  ASSERT_NEAR(kurtosis_n, kurtosis, 1e-16);

  // normal, kChi2
  conweight = std::vector<Tv>({1, 4});
  condists = std::vector<DistributionBase *>(
      {new Distribution<DistributionType::kNormal>(1.0, 2.0),
       new Distribution<DistributionType::kChi2>(4)});
  dis = DistributionMixture(conweight, condists);

  dis.GetMoments(mean, variance, skewness, kurtosis);

  // test it by sampling

  Ti N = 1000000;
  Tv *storage = new Tv[N];
  dis.GetSample(storage, N, 434);
  Matrix<Tv> mat = Matrix<Tv>(storage, N, 1);
  Descriptive desc = Descriptive(&mat);
  auto res = desc.MeanVarianceKurtosisSkewness(false);
  ASSERT_NEAR(std::get<0>(res), mean, 1e-2);
  ASSERT_NEAR(std::get<1>(res), variance, 1e-1);
  ASSERT_NEAR(std::get<2>(res), skewness, 1e-1);
  ASSERT_NEAR(std::get<3>(res), kurtosis, 1e-1);
}

TEST(Distributions_T, mixture_dis_moments) {

  // just one duniform
  auto disweight = std::vector<Tv>({1});
  auto disdists = std::vector<DistributionBase *>(
      {new Distribution<DistributionType::kUniformDis>(0, 3)});
  auto dis = DistributionMixture(disweight, disdists);
  auto dis0 = Distribution<DistributionType::kUniformDis>(0, 3);
  Tv mean = 0.0;
  Tv variance = 0.0;
  Tv kurtosis = 0.0;
  Tv skewness = 0.0;
  dis.GetMoments(mean, variance, skewness, kurtosis);
  ASSERT_NEAR(dis0.GetMean(), mean, 1e-16);
  ASSERT_NEAR(dis0.GetVariance(), variance, 1e-16);
  ASSERT_NEAR(dis0.GetSkewness(), skewness, 1e-14);
  ASSERT_NEAR(dis0.GetKurtosis(), kurtosis, 1e-14);

  // 3 similar duniforms
  disweight = std::vector<Tv>({1, 4, 6});
  disdists = std::vector<DistributionBase *>(
      {new Distribution<DistributionType::kUniformDis>(0, 3),
       new Distribution<DistributionType::kUniformDis>(0, 3),
       new Distribution<DistributionType::kUniformDis>(0, 3)});
  dis = DistributionMixture(disweight, disdists);
  dis.GetMoments(mean, variance, skewness, kurtosis);
  ASSERT_NEAR(dis0.GetMean(), mean, 1e-16);
  ASSERT_NEAR(dis0.GetVariance(), variance, 1e-16);
  ASSERT_NEAR(dis0.GetSkewness(), skewness, 1e-14);
  ASSERT_NEAR(dis0.GetKurtosis(), kurtosis, 1e-14);

  // multiple dists
  disweight = std::vector<Tv>({1, 4, 3, 5});
  disdists = std::vector<DistributionBase *>({
      new Distribution<DistributionType::kUniformDis>(0, 3),
      new Distribution<DistributionType::kBinomial>(0.6, 30),
      new Distribution<DistributionType::kGeometric>(0.3),
      new Distribution<DistributionType::kPoisson>(5.0),
  });
  dis = DistributionMixture(disweight, disdists);
  dis.GetMoments(mean, variance, skewness, kurtosis);

  // test it by sampling

  Ti N = 100000;
  Tv *storage = new Tv[N];
  dis.GetSample(storage, N, 434);
  Matrix<Tv> mat = Matrix<Tv>(storage, N, 1);
  Descriptive desc = Descriptive(&mat);
  auto res = desc.MeanVarianceKurtosisSkewness(false);
  ASSERT_NEAR(std::get<0>(res), mean, 3e-1);
  ASSERT_NEAR(std::get<1>(res), variance, 3e-1);
  ASSERT_NEAR(std::get<2>(res), skewness, 3e-1);
  ASSERT_NEAR(std::get<3>(res), kurtosis, 3e-1);
}
