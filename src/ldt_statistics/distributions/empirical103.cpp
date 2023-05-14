/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "distributions.h"
#include "functionx.h"

using namespace ldt;

DistributionEmpirical103::DistributionEmpirical103() {}

DistributionEmpirical103::DistributionEmpirical103(std::vector<Tv> &cdfs) {
  CDFs = &cdfs;
}

Tv DistributionEmpirical103::GetCDFApprox(Tv x) {
  if (std::isinf(x))
    return x > 0 ? 1 : 0;

  auto cdfs = *CDFs;
  Tv shift;
  Tv y;

  // deal with tails
  if (x <= cdfs[0]) { // calculate last change and use it
    shift = (x - cdfs[0]) * (0.001 - 0.0001) / (cdfs[1] - cdfs[0]);
    y = 0.0001;
  } else if (x >= cdfs[102]) {
    shift = (x - cdfs[102]) * (0.9999 - 0.999) / (cdfs[102] - cdfs[101]);
    y = 0.9999;
  }

  else { // find the position and calculate the shift
    Ti i = 1;
    Tv d, h;
    for (i; i < 103; i++) {
      d = cdfs[i];
      if (x <= d)
        break;
    }

    if (i == 1) { // cdf[0] <x< cdf[1]
      h = 0.001 - 0.0001;
      y = 0.001;
    } else if (i == 102) { // cdf[101] <x< cdf[102]
      h = 0.9999 - 0.999;
      y = 0.9999;
    } else if (i == 2) { // cdf[1] <x< cdf[2]
      h = 0.01 - 0.001;
      y = 0.01;
    } else if (i == 101) { // cdf[100] <x< cdf[101]
      h = 0.999 - 0.99;
      y = 0.999;
    } else { //    e.g. cdf[24] <x< cdf[25]
      h = 0.01;
      y = (i - 1) / 100.0;
    }
    shift = (d - x) * h / (cdfs.at(i) - cdfs.at(i - 1));
  }
  Tv res = y - shift;
  return res;
}

Tv DistributionEmpirical103::GetQuantileApprox(Tv p) {
  if (p <= 0)
    return -INFINITY;
  if (p >= 1)
    return INFINITY;

  auto cdfs = *CDFs;
  Tv shift;
  Tv y;

  // deal with tails
  if (p <= 0.0001) { // calculate last change and use it
    shift = (p - 0.0001) * (cdfs.at(1) - cdfs.at(0)) / (0.001 - 0.0001);
    y = cdfs.at(0);
  } else if (p >= 0.9999) {
    shift = (p - 0.9999) * (cdfs.at(102) - cdfs.at(101)) / (0.9999 - 0.999);
    y = cdfs.at(102);
  }

  else { // find the position and calculate the shift
    Tv d, h;
    Ti i;

    if (p < 0.001) { // 0.0001<p<0.001
      i = 0;
      d = 0.0001;
      h = 0.001 - 0.0001;
    } else if (p > 0.999) { // 0.999<p<0.9999
      i = 101;
      d = 0.999;
      h = 0.9999 - 0.999;
    } else if (p < 0.01) { // 0.001<p<0.01
      i = 1;
      d = 0.001;
      h = 0.01 - 0.001;
    } else if (p > 0.99) { // 0.99<p<0.999
      i = 100;
      d = 0.99;
      h = 0.999 - 0.99;
    } else {                         // e.g. p=0.244   =>    0.24<p<0.25
      d = std::floor(p * 100) / 100; // floor
      i = (Ti)(d * 100) + 1;
      h = 0.01;
    }
    shift = (p - d) * (cdfs.at(i + 1) - cdfs.at(i)) / h;
    y = cdfs.at(i);
  }
  Tv res = y + shift;
  return res;
}

void DistributionEmpirical103::Combine(
    const std::vector<DistributionEmpirical103 *> &dists,
    const std::vector<Tv> &weights, std::vector<Tv> &result) {

  if (weights.size() > 0 && dists.size() != weights.size())
    throw std::logic_error("Inconsistent number of weights and distributions.");
  Tv weight;
  if (weights.size() == 0)
    weight = 1.0 / dists.size();

  result.resize(103);

  if (dists.size() == 0)
    throw std::logic_error("Empty list of distributions!");
  else if (dists.size() == 1) {
    for (Ti i = 0; i < 103; i++)
      result.at(i) = dists.at(0)->CDFs->at(i);
    return;
  }

  // we find an x which minimizes the distance between prob and weighted
  // average of the CDFs
  Tv prob = NAN;
  std::function<Tv(const Tv &)> objective = [&weight, &prob, &dists,
                                             &weights](const Tv &x) -> Tv {
    Tv wa = 0;
    if (weights.size() == 0)
      for (Ti j = 0; j < dists.size(); j++)
        wa += weight * dists.at(j)->GetCDFApprox(x);
    else
      for (Ti j = 0; j < dists.size(); j++)
        wa += weights.at(j) * dists.at(j)->GetCDFApprox(x);
    return std::abs(wa - prob);
  };

  // x cannot be less than start
  Tv start;

  for (Ti i = 0; i < 103; i++) {
    prob = i == 0     ? 0.0001
           : i == 1   ? 0.001
           : i == 101 ? 0.999
           : i == 102 ? 0.9999
                      : (i - 1) / 100.0;

    if (i == 0) { // let start to be the minimum of all distributions
      Tv mi = INFINITY;
      for (Ti j = 0; j < dists.size(); j++)
        if (mi > dists.at(j)->CDFs->at(0))
          mi = dists.at(j)->CDFs->at(0);
      start = mi;
    }
    auto optim_x =
        NelderMead::Minimize1(objective, start, 0.1, 500, 1e-8, start, NAN);
    start = std::get<0>(optim_x);

    /*for (int t = 0; t < 1000000; t++) {
      auto diff = objective(mstart);
      if (diff < 0.00001)
        break;
      mstart.Data[0] += 0.00001;
    }*/

    result.at(i) = start;
  }
}