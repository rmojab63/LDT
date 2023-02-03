/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "discrete_choice.h"

using namespace ldt;

DiscreteChoiceExtended::DiscreteChoiceExtended(
    DiscreteChoiceModelType modelType, DiscreteChoiceDistType distType, Ti rows,
    Ti cols, bool hasWeight, bool checkNan, Ti numChoices, bool doDetails,
    Ti numForecast, PcaAnalysisOptions *pcaOptions,
    std::vector<Matrix<Tv>> *costMatrices, bool weightedEval) {

  mModelType = modelType;
  mDoDetails = doDetails;
  mNumChoices = numChoices;
  mHasWeight = hasWeight;
  mCheckNan = checkNan;
  mWeightedEval = hasWeight && weightedEval;
  pCostMatrices = costMatrices;

  Ti numExo = cols - (hasWeight ? 2 : 1);

  StorageSize = 0;
  WorkSize = 0;

  Data = Dataset<Tv>(rows, cols, checkNan, false /*for this, you should find a way to select column in forecast too*/);
  StorageSize += Data.StorageSize;

  if (pcaOptions && pcaOptions->IsEnabled()) {
    pcaOptions->CheckValidity();
    pPcaOptions = pcaOptions;
    Pca = PcaAnalysis(rows, numExo, numForecast, true, true, true, true);
    StorageSize += Pca.StorageSize;
    WorkSize = std::max(WorkSize, Pca.WorkSize);

    numExo = std::min(
        numExo, pcaOptions->CutoffCountMax +
                    pcaOptions->IgnoreFirstCount); // for estimating the model
  }

  Model = DiscreteChoiceBase::GetFromType(modelType, distType, rows, numExo,
                                          numChoices, doDetails);
  StorageSize += Model->StorageSize;
  WorkSize = std::max(WorkSize, Model->WorkSize);

  if (numForecast > 0) {
    StorageSize +=
        mNumChoices * numForecast +
        numForecast * cols; // for copying forecast matrix (TODO: this is work)
    WorkSize = std::max(WorkSize, rows + numChoices - 2);
  }

  if (doDetails) {                     // Auc and Cost
    StorageSize += mNumChoices * rows; // for projections
    WorkSize = std::max(WorkSize, rows + mNumChoices - 2);
    if (costMatrices) {
      if (hasWeight && weightedEval) {
        auto cost = CostMatrix<true>(costMatrices->size());
        StorageSize += cost.StorageSize;
      } else {
        auto cost = CostMatrix<false>(costMatrices->size());
        StorageSize += cost.StorageSize;
      }
    }
  }
}

DiscreteChoiceExtended::~DiscreteChoiceExtended() { delete Model; }

void DiscreteChoiceExtended::Calculate(const Matrix<Tv> &data, Tv *storage,
                                       Tv *work, bool olsInitial,
                                       const Matrix<Tv> *xForecast) {
  // you can add colIndexes as an argument, but as I wrote in the constructor,
  // you should handle forecast indexation too

  // if (PredProbs && !xForecast)

  Ti numForecast = xForecast ? xForecast->RowsCount : 0;
  Ti numObs = data.RowsCount;
  Ti numExo = data.ColsCount - (mHasWeight ? 2 : 1);

  auto temp = DiscreteChoiceExtended(
      Model->mModelType, Model->mDistType, data.RowsCount, data.ColsCount,
      mHasWeight, mCheckNan, mNumChoices, mDoDetails, numForecast, pPcaOptions,
      pCostMatrices, mWeightedEval);
  if (temp.WorkSize > WorkSize || temp.StorageSize > StorageSize)
    throw std::logic_error("Inconsistent arguments.");

  Ti p = 0;

  // Prepare
  Data.Calculate(data, nullptr, storage);
  p += Data.StorageSize;
  auto useMat = &Data.Result;

  // create matrixes but note that we might change X by PCA

  numExo = useMat->ColsCount - (mHasWeight ? 2 : 1);
  Ti count = useMat->RowsCount;
  Tv *d = useMat->Data;
  Y.SetData(d, count, 1);
  if (mHasWeight) {
    W.SetData(&d[count], count, 1);
    X.SetData(&d[2 * count], count, numExo);
  } else {
    X.SetData(&d[count], count, numExo);
  }

  auto useForecast = Matrix<Tv>();
  if (pPcaOptions) { // update useMat with PCs
    if (xForecast) {
      if (xForecast->ColsCount != X.ColsCount)
        throw std::logic_error(
            "Data and forecast data has different number of columns.");
      useForecast.SetData(&storage[p], numForecast, X.ColsCount);
      p += useForecast.length();
      useForecast.CopyFrom00(*xForecast);
    }

    pPcaOptions->CalculateForModel(Pca, X, work, &storage[p],
                                   xForecast ? &useForecast : nullptr);
    p += Pca.StorageSize;
  }

  // Estimate
  Model->Calculate(Y, X, mHasWeight ? &W : nullptr, &storage[p], work,
                   mNumChoices, olsInitial);
  p += Model->StorageSize;

  // predict probabilities
  if (xForecast) {
    PredProbs.SetData(&storage[p], xForecast->RowsCount, mNumChoices);
    p += xForecast->RowsCount * mNumChoices;
    if (pPcaOptions)
      Model->GetProbabilities(useForecast, PredProbs, work);
    else
      Model->GetProbabilities(*xForecast, PredProbs, work);
  }

  if (mDoDetails) {
    // Projections
    Projections = Matrix<Tv>(&storage[p], numObs, mNumChoices);
    p += numObs * mNumChoices;
    Model->GetProbabilities(X, Projections, work);

    if (mHasWeight && mWeightedEval) {
      if (mModelType == DiscreteChoiceModelType::kBinary) {
        auto auc = ROC<true, true>(numObs);
        auc.Calculate(Y, Projections, &W, false);
        Auc = auc.Result;
      } else {
        auto auc = ROC<true, false>(numObs);
        auc.Calculate(Y, Projections, &W, false);
        Auc = auc.Result;
      }
    } else {
      if (mModelType == DiscreteChoiceModelType::kBinary) {
        auto auc = ROC<false, true>(numObs);
        auc.Calculate(Y, Projections, nullptr, false);
        Auc = auc.Result;
      } else {
        auto auc = ROC<false, false>(numObs);
        auc.Calculate(Y, Projections, nullptr, false);
        Auc = auc.Result;
      }
    }

    if (pCostMatrices) {
      if (mHasWeight && mWeightedEval) {
        auto cost = CostMatrix<true>(pCostMatrices->size());
        cost.Calculate(*pCostMatrices, Y, Projections, &W, &storage[p]);
        p += cost.StorageSize;
        CostRatioAvg = cost.AverageRatio;
      } else {
        auto cost = CostMatrix<false>(pCostMatrices->size());
        cost.Calculate(*pCostMatrices, Y, Projections, nullptr, &storage[p]);
        p += cost.StorageSize;
        CostRatioAvg = cost.AverageRatio;
      }
    }
  }
}