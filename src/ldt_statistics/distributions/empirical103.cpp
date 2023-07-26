/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "distributions.h"
#include "functionx.h"

using namespace ldt;

DistributionEmpirical103::DistributionEmpirical103() {}

DistributionEmpirical103::DistributionEmpirical103(Tv *quantiles) {
  Quantiles = quantiles;
}

Tv DistributionEmpirical103::GetCDFApprox(Tv x) const {
  if (std::isinf(x))
    return x > 0 ? 1 : 0;

  Tv shift;
  Tv y;

  // deal with tails
  if (x <= Quantiles[0]) { // calculate last change and use it
    shift =
        (x - Quantiles[0]) * (0.001 - 0.0001) / (Quantiles[1] - Quantiles[0]);
    y = 0.0001;
  } else if (x >= Quantiles[102]) {
    shift = (x - Quantiles[102]) * (0.9999 - 0.999) /
            (Quantiles[102] - Quantiles[101]);
    y = 0.9999;
  }

  else { // find the position and calculate the shift
    Ti i = 1;
    Tv d, h;
    for (; i < 103; i++) {
      d = Quantiles[i];
      if (x <= d)
        break;
    }

    if (i == 1) { // Quantiles[0] <x< Quantiles[1]
      h = 0.001 - 0.0001;
      y = 0.001;
    } else if (i == 102) { // Quantiles[101] <x< Quantiles[102]
      h = 0.9999 - 0.999;
      y = 0.9999;
    } else if (i == 2) { // Quantiles[1] <x< Quantiles[2]
      h = 0.01 - 0.001;
      y = 0.01;
    } else if (i == 101) { // Quantiles[100] <x< Quantiles[101]
      h = 0.999 - 0.99;
      y = 0.999;
    } else { //    e.g. Quantiles[24] <x< Quantiles[25]
      h = 0.01;
      y = (i - 1) / 100.0;
    }
    shift = (d - x) * h / (Quantiles[i] - Quantiles[i - 1]);
  }
  Tv res = y - shift;
  return res;
}

Tv DistributionEmpirical103::GetQuantileApprox(Tv p) const {
  if (p <= 0)
    return -INFINITY;
  if (p >= 1)
    return INFINITY;

  Tv shift;
  Tv y;

  // deal with tails
  if (p <= 0.0001) { // calculate last change and use it
    shift = (p - 0.0001) * (Quantiles[1] - Quantiles[0]) / (0.001 - 0.0001);
    y = Quantiles[0];
  } else if (p >= 0.9999) {
    shift = (p - 0.9999) * (Quantiles[102] - Quantiles[101]) / (0.9999 - 0.999);
    y = Quantiles[102];
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
    shift = (p - d) * (Quantiles[i + 1] - Quantiles[i]) / h;
    y = Quantiles[i];
  }
  Tv res = y + shift;
  return res;
}

void DistributionEmpirical103::Combine(
    const std::vector<DistributionEmpirical103> &dists, const Tv *weights,
    Tv *result) {

  Tv weight;
  if (!weights)
    weight = 1.0 / dists.size();

  if (dists.size() == 0)
    throw LdtException(ErrorType::kLogic, "empirical103",
                       "empty list of distributions!");
  else if (dists.size() == 1) {
    for (Ti i = 0; i < 103; i++)
      result[i] = dists.at(0).Quantiles[i];
    return;
  }

  // we find an x which minimizes the distance between prob and weighted
  // average of the quantiles
  Tv prob = NAN;
  std::function<Tv(const Tv &)> objective = [&weight, &prob, &dists,
                                             &weights](const Tv &x) -> Tv {
    Tv wa = 0;
    if (!weights)
      for (Ti j = 0; j < (Ti)dists.size(); j++)
        wa += weight * dists.at(j).GetCDFApprox(x);
    else
      for (Ti j = 0; j < (Ti)dists.size(); j++)
        wa += weights[j] * dists.at(j).GetCDFApprox(x);
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
      for (Ti j = 0; j < (Ti)dists.size(); j++)
        if (mi > dists.at(j).Quantiles[0])
          mi = dists.at(j).Quantiles[0];
      start = mi;
    }
    auto optim_x =
        NelderMead::Minimize1(objective, start, 0.1, 500, 1e-8, start, NAN);
    start = std::get<0>(optim_x);

    result[i] = start;
  }
}