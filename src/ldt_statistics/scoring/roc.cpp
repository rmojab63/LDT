/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "functionx.h"
#include "scoring.h"

using namespace ldt;

template <bool hasWeight, bool hasCost> ROC<hasWeight, hasCost>::ROC(){};

template <bool hasWeight, bool hasCost> ROC<hasWeight, hasCost>::ROC(Ti n){};

template <bool hasWeight, bool hasCost>
void ROC<hasWeight, hasCost>::Calculate(Matrix<Tv> &y, Matrix<Tv> &scores,
                                        Matrix<Tv> *weights,
                                        bool normalizePoints, Tv lowerThreshold,
                                        Tv upperThreshold, Tv epsilon,
                                        bool pessimistic, Matrix<Tv> *costs,
                                        Matrix<Tv> *frequencyCost) {
  bool isPartial = std::isnan(lowerThreshold) == false &&
                   std::isnan(upperThreshold) == false;
  if (isPartial) {
    normalizePoints =
        true; // we should normalize the points before using thresholds
    if (upperThreshold < lowerThreshold || upperThreshold > 1 ||
        lowerThreshold < 0)
      throw std::logic_error("Invalid bounds in partial AUC.");
  }

  if constexpr (hasCost) {
    if (!frequencyCost || frequencyCost->RowsCount != 2 ||
        frequencyCost->ColsCount != 2)
      throw std::logic_error("Missing or invalid cost matrix.");
  }

  Ti n = y.length();
  if (n == 0)
    throw std::logic_error("zero number of observations in calculating ROC.");

  // first column must be the scores (i.e., probability of negative
  // observations) score=1 => It is definitely negative
  auto sortedIndexes = std::vector<Ti>();
  SortIndexes(scores.Data, n, sortedIndexes, true);

  bool isNeg = false;
  Ti ind;
  Tv yi, th, thPre = scores.Data[sortedIndexes[0]], w = 1, sumFP = 0, horiz = 0,
             sumTP = 0, verti = 0, area = 0;
  Points.clear();
  Points.push_back(std::make_tuple<Tv, Tv>(0, 0)); // start from origin
  for (Ti i = 0; i < n; i++) {
    ind = sortedIndexes[i];
    th = scores.Data[ind];
    yi = y.Data[ind];
    w = 1;
    if constexpr (hasWeight) {
      w = weights->Data[ind];
    }

    if (std::abs(th - thPre) > epsilon) { // push new point(s)
      sumTP += verti; // overall vertical move (for the rectangle)
      sumFP += horiz;
      if (pessimistic)
        Points.push_back(std::make_tuple(sumFP, 0)); // no vertical move
      Points.push_back(std::make_tuple(sumFP, sumTP));
      thPre = th;
      horiz = 0;
      verti = 0;
    }

    //    At the current threshold: this observation and all
    //    previous observations are predicted to be positive

    isNeg = yi == 0;

    if constexpr (hasCost) {
      Tv xi = costs ? costs->Data[ind] : 1;
      Tv b_tp = frequencyCost->Data[0] * xi - frequencyCost->Data[2];
      if (b_tp < 0)
        throw std::logic_error("Invalid cost matrix: benefit of TP is "
                               "negative. Check the first row.");
      Tv c_fp = -(frequencyCost->Data[1] * xi - frequencyCost->Data[3]);
      if (c_fp < 0)
        throw std::logic_error("Invalid cost matrix: cost of FP is negative. "
                               "Check the second row.");
      if (isNeg)
        horiz += w * c_fp;
      else
        verti += w * b_tp;

    } else if constexpr (true) {

      if (isNeg) // A false-positive -> a horizontal move in ROC
        horiz += w;
      else // a true-positive -> a vertical move
        verti += w;
    }

    // for pessimistic: see
    // https://link.springer.com/article/10.1007/s00357-019-09345-1
  }

  sumTP += verti;
  sumFP += horiz;
  Points.push_back(std::make_tuple(sumFP, sumTP));

  if (normalizePoints) {
    for (auto &p : Points) {
      std::get<0>(p) /= sumFP;
      std::get<1>(p) /= sumTP;
    }
    if (isPartial) { // TODO: create a Bounded AUC class and move this logic
                     // to that class
      // throw std::logic_error("Partial AUC is not implemented.");
      std::vector<std::tuple<Tv, Tv>> newPoints;
      Tv x, y, x0 = 0, y0 = 0, slope;
      for (auto &p : Points) {
        x = std::get<0>(p);
        y = std::get<1>(p);
        slope = (y - y0) / (x - x0);

        if (x >= lowerThreshold && x0 <= upperThreshold) {
          if (x > lowerThreshold &&
              x0 < lowerThreshold) // don't miss the first point
            newPoints.push_back(std::make_tuple(
                lowerThreshold, y0 + (lowerThreshold - x0) * slope));

          if (x >= lowerThreshold && x <= upperThreshold)
            newPoints.push_back(std::make_tuple(x, y));

          if (x > upperThreshold &&
              x0 < upperThreshold) // don't miss the last point
            newPoints.push_back(std::make_tuple(
                upperThreshold, y - (x - upperThreshold) * slope));
        }

        x0 = x;
        y0 = y;
      }
      AucPoints<false> a(newPoints, 0);
      Result = a.Result / (upperThreshold - lowerThreshold);
    } else {
      AucPoints<false> a(Points, 0);
      Result = a.Result;
    }
  } else {
    AucPoints<false> a(Points, 0);
    Result = a.Result / (sumFP * sumTP); // normalize
  }
}

template class ldt::ROC<true, true>;
template class ldt::ROC<true, false>;
template class ldt::ROC<false, true>;
template class ldt::ROC<false, false>;
