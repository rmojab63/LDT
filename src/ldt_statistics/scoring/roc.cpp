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
                                         bool normalizePoints) {
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
      AucPoints<false> a(Points, 0);
      Result = a.Result; // normalize
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
