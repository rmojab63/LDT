/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "functionx.h"
#include "scoring.h"

using namespace ldt;

template <bool hasWeight, bool isBinary> ROC<hasWeight, isBinary>::ROC(){};

template <bool hasWeight, bool isBinary> ROC<hasWeight, isBinary>::ROC(Ti n){};

template <bool hasWeight, bool isBinary>
void ROC<hasWeight, isBinary>::Calculate(Matrix<Tv> &y, Matrix<Tv> &scores,
                                         Matrix<Tv> *weights,
                                         bool normalizePoints,
                                         Tv lowerThreshold, Tv upperThreshold) {
  bool isPartial = std::isnan(lowerThreshold) == false &&
                   std::isnan(upperThreshold) == false;
  if (isPartial)
    normalizePoints =
        true; // we should normalize the points before using thresholds

  Ti n = y.length();
  if (n == 0)
    throw std::logic_error("zero number of observations in calculating ROC.");

  if constexpr (isBinary) {
    // first column must be the scores (i.e., probability of negative
    // observations) score=1 => It is definitely negative
    auto sortedIndexes = std::vector<Ti>();
    SortIndexes(scores.Data, n, sortedIndexes, true);

    bool isNeg = false;
    Ti ind;
    Tv yi, th, thPre = scores.Data[sortedIndexes[0]], w = 1, sumFP = 0,
               horiz = 0, sumTP = 0, verti = 0, area = 0;
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

      if (th != thPre) { // push a new point
        sumTP += verti;  // overall vertical move (for the rectangle)
        sumFP += horiz;
        Points.push_back(std::make_tuple(sumFP, sumTP));
        thPre = th;
        horiz = 0;
        verti = 0;
      }

      //    At the current threshold: this observation and all
      //    previous observations are predicted to be positive
      isNeg = yi == 0;

      if (isNeg) { // A false-positive -> a horizontal move in ROC
        horiz += w;
      } else { // a true-positive -> a vertical move
        verti += w;
      }
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

          if (x >= lowerThreshold &&
              x0 < lowerThreshold) // we should add the first point
            newPoints.push_back(std::make_tuple(
                lowerThreshold, y0 + (lowerThreshold - x0) * slope));

          if (x > upperThreshold &&
              x0 <= upperThreshold) // we should add the last point
            newPoints.push_back(std::make_tuple(
                upperThreshold, y - (x - upperThreshold) * slope));

          else if (x >= lowerThreshold && x0 >= lowerThreshold) {
            newPoints.push_back(std::make_tuple(x, y));
          }
          x0 = x;
          y0 = y;
        }
        AucPoints<false> a(newPoints, 0);
        Result = a.Result;
      } else {
        AucPoints<false> a(Points, 0);
        Result = a.Result;
      }
    } else {
      AucPoints<false> a(Points, 0);
      Result = a.Result / (sumFP * sumTP); // normalize
    }

  } else {
    throw std::logic_error("Not implemented (ROC).");
  }
}

template class ldt::ROC<true, true>;
template class ldt::ROC<true, false>;
template class ldt::ROC<false, true>;
template class ldt::ROC<false, false>;
