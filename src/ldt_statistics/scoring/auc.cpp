/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "scoring.h"

using namespace ldt;

template <bool hasWeight, bool isBinary> AUC<hasWeight, isBinary>::AUC(){};

template <bool hasWeight, bool isBinary> AUC<hasWeight, isBinary>::AUC(Ti n){};

template <bool hasWeight, bool isBinary>
void AUC<hasWeight, isBinary>::Calculate(Matrix<Tv> &y, Matrix<Tv> &scores,
                                         Matrix<Tv> *weights,
                                         Matrix<Tv> *multi_class_weights_) {
  Ti n = y.length();
  if (n == 0)
    throw std::logic_error("zero number of observations in calculating AUC.");

  if constexpr (isBinary) {
    // first column must be the scores (i.e., probability of negative
    // observations) score=1 => It is definitely negative
    auto sortedIndexes = std::vector<Ti>();
    SortIndexes(scores.Data, n, sortedIndexes, true);

    bool isNeg = false;
    Ti ind;
    Tv yi, w = 1, sumFP = 0, sumTP = 0, sumTpPre = 0, area = 0;
    for (Ti i = 0; i < n; i++) {
      ind = sortedIndexes[i];
      yi = y.Data[ind];
      w = 1;
      if constexpr (hasWeight) {
        w = weights->Data[ind];
      }
      isNeg = yi == 0;

      //    Note that at the current threshold: this observation and all
      //    previous observations are predicted to be positive

      if (isNeg) { // A false-positive (a horizontal move in ROC) Let's increase
                   // the area: add a rectangle and subtract a triangle
        sumFP += w;
        area += w * (sumTpPre + (sumTP - sumTpPre) / 2);
        sumTpPre = sumTP;
      } else { // a true-positive -> a vertical move. wait
        sumTP += w;
      }
    }
    if (isNeg == false) // add the last horizontal move
      area += w * (sumTpPre + (sumTP - sumTpPre) / 2);
    Result = area / (sumFP * sumTP); // normalize
  } else {
    throw std::logic_error("Not implemented (AUC).");
  }
}

template class ldt::AUC<true, true>;
template class ldt::AUC<true, false>;
template class ldt::AUC<false, true>;
template class ldt::AUC<false, false>;
