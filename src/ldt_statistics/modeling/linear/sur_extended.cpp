/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "sur.h"

using namespace ldt;

SurExtended::SurExtended(Ti N, Ti m, Ti k, bool isRestricted, bool checkNan,
                         bool doDetails, Ti numProjections, Ti maxSigSearchIter,
                         bool doForcVariance, PcaAnalysisOptions *pcaOptionsY,
                         PcaAnalysisOptions *pcaOptionsX) {
  mCheckNan = checkNan;
  mHasPcaY = pcaOptionsY && pcaOptionsY->IsEnabled();
  mHasPcaX = pcaOptionsX && pcaOptionsX->IsEnabled();

  StorageSize = 0;
  WorkSize = 0;

  if (checkNan || mHasPcaY || mHasPcaX) {
    Data = Dataset<Tv>(N, m + k, checkNan, false);
    StorageSize += Data.StorageSize;
  }

  if (mHasPcaY) {
    pcaOptionsY->CheckValidity();
    pPcaOptionsY = pcaOptionsY;
    PcaY = PcaAnalysis(N, m, 0, true, true, true, true);
    StorageSize += PcaY.StorageSize;
    WorkSize = std::max(WorkSize, PcaY.WorkSize);
    // update number of equations
    m = std::min(m,
                 pcaOptionsY->CutoffCountMax +
                     pcaOptionsY->IgnoreFirstCount); // for estimating the model
  }

  if (mHasPcaX) {
    pcaOptionsX->CheckValidity();
    pPcaOptionsX = pcaOptionsX;
    PcaX = PcaAnalysis(N, k, numProjections, true, true, true, true);
    StorageSize += PcaX.StorageSize;
    WorkSize = std::max(WorkSize, PcaX.WorkSize);
    // update number of exogenous data
    k = std::min(k,
                 pcaOptionsX->CutoffCountMax +
                     pcaOptionsX->IgnoreFirstCount); // for estimating the model

    if (numProjections > 0)
      StorageSize += numProjections * k; // copy X
  }

  Model = Sur(N, m, k, isRestricted, doDetails, maxSigSearchIter);
  StorageSize += Model.StorageSize;
  WorkSize = std::max(WorkSize, Model.WorkSize);

  if (numProjections > 0) {
    Projections =
        SurProjection(numProjections, m, k, isRestricted, doForcVariance);
    StorageSize += Projections.StorageSize;
    WorkSize = std::max(WorkSize, Projections.WorkSize);
  }
}

void SurExtended::Calculate(const Matrix<Tv> &data, Ti m, Tv *storage, Tv *work,
                            Matrix<Tv> *R, Tv sigSearchMaxProb,
                            const Matrix<Tv> *newX,
                            const SearchModelChecks *checks) {
  Ti N = data.RowsCount;
  Ti k = data.ColsCount - m;
  if (k < 0)
    throw LdtException(ErrorType::kLogic, "sur-extended",
                       "invalid number of equations in SUR extended");

  Ti numProjections = newX ? newX->RowsCount : 0;

  auto temp =
      SurExtended(N, m, k, Model.mIsRestricted, mCheckNan, Model.mDoDetails,
                  numProjections, Model.mSigSearchMaxIter,
                  Projections.mDoVariance, pPcaOptionsY, pPcaOptionsX);
  if (temp.WorkSize > WorkSize || temp.StorageSize > StorageSize)
    throw LdtException(ErrorType::kLogic, "sur-extended",
                       "inconsistent arguments (in SurExtended)");

  Ti p = 0;

  // Prepare
  Ti count = Data.Result.RowsCount;
  Tv *d = nullptr;
  if (mCheckNan || mHasPcaY || mHasPcaX) {
    Data.Calculate(data, nullptr, storage);
    p += Data.StorageSize;
    d = Data.Result.Data;
    count = Data.Result.RowsCount;
  } else {
    d = data.Data;
    count = data.RowsCount;
  }

  if (checks) {
    if (checks->MinObsCount > 0 && checks->MinObsCount > count)
      throw LdtException(ErrorType::kLogic, "sur-extended",
                         "model check: minimum no. obs");
    if (checks->MinDof > 0 && checks->MinDof > count - k)
      throw LdtException(ErrorType::kLogic, "sur-extended",
                         "model check: minimum dof");
  }

  // create matrixes but note that we might update them by PCA
  Y.SetData(d, count,
            m); // since variables are in columns and matrix is column-major
  X.SetData(&d[m * count], count, k);

  if (mHasPcaY) { // change data for Y
    pPcaOptionsY->CalculateForModel(PcaY, Y, work, &storage[p], nullptr,
                                    true); // throw if endogenous is constant
    p += PcaY.StorageSize;
  }

  auto useNewX = Matrix<Tv>();
  if (mHasPcaX) { // change data for X
    if (numProjections > 0) {
      useNewX.SetData(&storage[p], numProjections, k);
      p += numProjections * k;
      newX->CopyTo00(useNewX);
    }
    pPcaOptionsX->CalculateForModel(PcaX, X, work, &storage[p],
                                    numProjections > 0 ? &useNewX : nullptr,
                                    false); // ignore constant variables ?!
    p += PcaX.StorageSize;
  } else if (numProjections > 0) {
    if (newX->ColsCount != k)
      throw LdtException(ErrorType::kLogic, "sur-extended",
                         "invalid number of variables in new exogenous data");
    useNewX.SetData(newX->Data, numProjections, k);
  }

  if (Model.mSigSearchMaxIter > 0) {
    if (!R)
      throw LdtException(
          ErrorType::kLogic, "sur-extended",
          "Restriction matrix cannot be null when significance search is "
          "enabled");
    auto km0 = Y.ColsCount * X.ColsCount;
    R->Restructure0(km0, km0);
  }
  // else If PCA changes the dimension and R dimension is not predicted
  // correctly, I expect Model to throw an exception.

  // Estimate
  Model.Calculate(Y, X, &storage[p], work, R, sigSearchMaxProb);
  p += Model.StorageSize;

  if (checks) {
    if (checks->mCheckCN_all &&
        Model.condition_number > checks->MaxConditionNumber)
      throw LdtException(ErrorType::kLogic, "sur-extended",
                         "model check: maximum cn");
    if (checks->MaxAic < Model.Aic)
      throw LdtException(ErrorType::kLogic, "sur-extended",
                         "model check: maximum aic");
    if (checks->MaxSic < Model.Sic)
      throw LdtException(ErrorType::kLogic, "sur-extended",
                         "model check: maximum sic");
    if (checks->MinR2 > Model.r2)
      throw LdtException(ErrorType::kLogic, "sur-extended",
                         "model check: maximum R2");
  }

  // predict
  if (numProjections > 0)
    Projections.Calculate(Model, useNewX, &storage[p], work);
}
