/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "varma.h"

using namespace ldt;

/*


VarmaSimulationDetail::VarmaSimulationDetail() {
        Actuals = new Matrix<Tv>();
        Forecasts = new Matrix<Tv>();
        ForecastsStd = new Matrix<Tv>();
        IsValid = false;
        SampleEnds = -1;
}


VarmaSimulationDetail::~VarmaSimulationDetail() {
        delete Actuals;
        delete Forecasts;
        delete ForecastsStd;
}
*/

// #pragma region Simulation

VarmaSimulation::VarmaSimulation(const VarmaSizes &sizes, Ti count,
                                 const std::vector<Ti> &horizons,
                                 const std::vector<ScoringType> &metrics,
                                 LimitedMemoryBfgsbOptions *optimOptions,
                                 bool isExtended,
                                 VarmaRestrictionType restriction,
                                 bool checkNan, PcaAnalysisOptions *pcaOptionsY,
                                 PcaAnalysisOptions *pcaOptionsX) {
  pSizes = &sizes;
  pHorizons = &horizons;
  pMetrics = &metrics;
  mCount = count;
  auto count0 = mCount;
  if (count == 0 || count0 >= sizes.T)
    throw std::logic_error(
        std::string(
            "Invalid number of simulations. It is zero or larger than the "
            "number of observations.") +
        std::to_string(count0) + std::string("...") + std::to_string(sizes.T));

  IsExtended = isExtended;
  auto horizonMax = horizons.at(horizons.size() - 1);

  mDoForecastVar = false;
  for (auto &s : metrics)
    if (Scoring::RequiresVariance(s)) {
      mDoForecastVar = true;
      break;
    }

  StorageSize = (Ti)(metrics.size() * horizons.size() *
                     sizes.EqsCount); // for each metric and horizon and
                                      // variable, there is a value
  Results = std::vector<Matrix<Tv>>(metrics.size());
  StorageSize += (Ti)metrics.size() * sizes.EqsCount;
  WorkSize = 0;
  // calculate work size
  if (IsExtended) {
    EModel = VarmaExtended(sizes, restriction, checkNan, true, true, horizonMax,
                           pcaOptionsY, pcaOptionsX, optimOptions);
    WorkSize += EModel.WorkSize + EModel.StorageSize;
  } else {
    Model = Varma(sizes, true, true, true,
                  optimOptions); // we save space for details and variance
                                 // of coefficients
    WorkSize +=
        Model.Result.StorageSize; // TODO: it can be more efficient. we
                                  // need a partition of these observations
    Forecast = VarmaForecast(Model.Sizes, horizonMax, mDoForecastVar,
                             false); // TODO: coefficients variance
    WorkSize += Forecast.StorageSize +
                std::max(Model.Result.WorkSize, Forecast.WorkSize);
  }

  //    others
  WorkSize += 5 * (sizes.EqsCount * (Ti)horizons.size()) + sizes.EqsCount;
}

void GetScore(const ScoringType &type, Matrix<Tv> &result,
              const Matrix<Tv> &act, const Matrix<Tv> &forc,
              const Matrix<Tv> &err, const Matrix<Tv> &std,
              const Matrix<Tv> &last_m) {
  Ti i, j;
  switch (type) {
  case ScoringType::kDirection:
    Tv d, f, a, l;
    for (i = 0; i < act.RowsCount; i++) {
      for (j = 0; j < act.ColsCount; j++) {
        d = 0;
        f = forc.Get0(i, j);
        a = act.Get0(i, j);
        l = last_m.Data[i];

        if (std::isnan(f))
          d = NAN;
        else if (a > l) {
          if (f > l)
            d = 1;
        } else if (a < l) {
          if (f < l)
            d = 1;
        } else if (a == l) {
          if (f == l)
            d = 1;
        }
        result.Set0(i, j, d);
      }
    }
    break;

  case ScoringType::kSign: {
    std::function<Tv(Tv, Tv)> g = [](Tv a, Tv f) -> Tv {
      if (std::isnan(f))
        return NAN;
      else if (a == 0 || f == 0)
        return 0.5; // award of 0 is 0.5
      else if ((a > 0 && f < 0) || (a < 0 && f > 0))
        return 0.0;
      return 1.0;
    };
    act.Apply0(forc, g, result);
  } break;

  case ScoringType::kMae: {
    std::function<Tv(Tv)> q = [](Tv x) -> Tv { return std::abs(x); };
    err.Apply0(q, result);
  } break;
  case ScoringType::kMape: {
    for (i = 0; i < act.length(); i++) {
      if (act.Data[i] == 0)
        result.Data[i] = NAN; //??
      else
        result.Data[i] = std::abs(err.Data[i]) / std::abs(act.Data[i]);
    }
  } break;

  case ScoringType::kRmse: {
    std::function<Tv(Tv)> k = [](Tv x) -> Tv { return x * x; };
    err.Apply0(k,
               result); // note that a sqrt is required for aggregating
  } break;
  case ScoringType::kRmspe: {
    for (i = 0; i < act.length(); i++) {
      if (act.Data[i] == 0)
        result.Data[i] = NAN; //??
      else
        result.Data[i] =
            err.Data[i] * err.Data[i] / (act.Data[i] * act.Data[i]);
    }
  } break;

  case ScoringType::kCrps: {
    std::function<Tv(Tv, Tv)> v = [](Tv e, Tv s) -> Tv {
      return Scoring::GetScoreCrpsNormal(e, 0, s);
    };
    err.Apply0(std, v, result);
  } break;

  default:
    throw std::logic_error("not implemented");
  }
}

void VarmaSimulation::Calculate(
    Tv *storage, Tv *work, Matrix<Tv> &data, bool &cancel, Matrix<Tv> *exo,
    const Matrix<Tv> *R, const Matrix<Tv> *r, bool usePreviousEstim,
    Tv maxCondNum, Tv stdMultipler, bool coefUncer, Ti maxInvalidSim,
    const std::function<void(Tv &)> *transformForMetrics) {

  if (cancel)
    return;

  auto horizons = *pHorizons;
  auto metrics = *pMetrics;
  auto sizes = *pSizes;
  Ti hh = (Ti)horizons.size();
  Ti hMin = horizons.at(0);
  Ti hMax = horizons.at(hh - 1);
  auto T = data.ColsCount;
  Ti mm = (Ti)metrics.size();
  Ti yy = sizes.EqsCount; // number of endogenous data
  bool isRestricted = R && R->length() > 0;

  // check size
  auto tempm = VarmaSimulation(sizes, mCount, horizons, metrics,
                               &Model.Result.Optim.Options);
  if (WorkSize < tempm.WorkSize || StorageSize < tempm.StorageSize)
    throw std::logic_error("inconsistent arguments in VARMA simulation");
  tempm = VarmaSimulation();

  auto count = mCount;

  if (T - count <= 0)
    throw std::logic_error(
        "invalid number of simulations. It is larger than available data");
  if (count == 0)
    throw std::logic_error(
        "invalid number of simulations. It is zero of negative");

  if (cancel)
    return;

  std::function<void(Tv &)> tfm;
  if (transformForMetrics)
    tfm = *transformForMetrics;

  // set storage
  Ti pos = 0;
  for (int i = 0; i < mm; i++) {
    Results.at(i).SetData(0, &storage[pos], yy, hh);
    pos += yy * hh; // horizons are in columns
  }
  ResultAggs.SetData(0, &storage[pos], mm, yy);
  pos += mm * yy;

  // set work
  pos = 0;
  double *estimStorage = nullptr;
  double *forecastStorage = nullptr;
  double *estimForecastWork = nullptr;

  estimStorage = &work[pos];
  pos += Model.Result.StorageSize;
  forecastStorage = &work[pos];
  pos += Forecast.StorageSize;
  estimForecastWork = &work[pos];
  pos += std::max(Forecast.WorkSize, Model.Result.WorkSize);

  //      others
  Ti co = yy * hh;
  auto act = Matrix<Tv>(&work[pos], yy, hh);
  pos += co;
  auto forc = Matrix<Tv>(&work[pos], yy, hh);
  pos += co;
  auto err = Matrix<Tv>(&work[pos], yy, hh);
  pos += co;
  auto std = Matrix<Tv>(&work[pos], yy, hh);
  pos += co;
  auto temp = Matrix<Tv>(&work[pos], yy, hh);
  pos += co;
  auto last = Matrix<Tv>(&work[pos], yy, 1);
  pos += yy;

  ValidCounts = 0;
  Ti k, invalidCounts = 0;
  Ti forIndex;
  Ti colc;
  bool useCurrentEstims;
  Ti actIndex;
  Ti effectiveH;
  bool success = false;
  Ti row = 0;
  auto counters = std::vector<Ti>(horizons.size());
  Ti counter = 0;
  for (Ti se = count; se > 0; se--) {
    if (cancel)
      return;
    useCurrentEstims = usePreviousEstim && success;
    invalidCounts++;
    actIndex = T - se - 1;
    if (actIndex < 0)
      break;

    effectiveH = std::min(se, hMax);
    if (effectiveH < hMin)
      break;

    colc = 0;
    for (auto &h : horizons) {
      if (h <= effectiveH)
        colc++;
    }
    if (colc <= 0)
      break; // no more forecasts

    // estimate and forecast
    // update sizes by current sampleEnd. se is different which means Xt,
    // resid, etc are different
    auto forecast0 = &Forecast;
    try {
      Model.EstimateMl(data, exo, estimForecastWork, estimStorage,
                       isRestricted ? R : nullptr, isRestricted ? r : nullptr,
                       se, useCurrentEstims, stdMultipler, maxCondNum);
      if (cancel)
        return;
      Forecast.Calculate(Model, exo, &data, forecastStorage, estimForecastWork,
                         effectiveH);
      if (cancel)
        return;

    } catch (...) {
      row++;
      continue; // will be invalid
    }
    invalidCounts--;
    if (invalidCounts > maxInvalidSim)
      throw std::logic_error("Model check failed: Minimum Valid Simulations");

    auto forecast = *forecast0;

    forIndex = forecast.StartIndex - 1;
    data.GetColumn0(actIndex, last);

    k = 0;
    for (auto &h : horizons) {
      if (h <= effectiveH) {
        if (cancel)
          return;
        counter++;
        counters[k]++;
        act.SetColumnFromColumn0(k, data, actIndex + h);

        forc.SetColumnFromColumn0(k, forecast.Forecast, forIndex + h);
        if (forecast.Variance.Data)
          std.SetColumnFromColumn0(k, forecast.Variance, forIndex + h);
      } else {
        act.SetColumn0(k, NAN);
        forc.SetColumn0(k, NAN);
        std.SetColumn0(k, NAN);
      }
      k++;
    }

    // Convert to STD
    if (mDoForecastVar) { // if you are going to check NAN, you should consider
                          // the effective horizon
      for (int i = 0; i < co; i++)
        std.Data[i] = std::sqrt(std.Data[i]);
    }

    if (transformForMetrics) {

      for (int i = 0; i < yy; i++)
        tfm(last.Data[i]);

      for (int i = 0; i < co; i++) {
        tfm(act.Data[i]);
        tfm(forc.Data[i]);
      }

      if (mDoForecastVar) { // transform STDs for CRPS too
        for (int i = 0; i < co; i++) {
          tfm(std.Data[i]);
        }
      }
    }

    act.Subtract0(forc, err);

    if (cancel)
      return;

    /*if (keepDetails) {
        details->sampleEnd = se;
        details->isValid = true;
        act.CopyTo00(details->actual);
        forc.CopyTo00(details->forecast);
        std.CopyTo00(details->stdforecast);
    }*/
    success = true;
    ValidCounts++;

    Ti c = 0;
    for (auto &eval : metrics) {
      if (cancel)
        return;
      GetScore(eval, temp, act, forc, err, std, last);

      // summation for calculating mean
      auto me = Results.at(c);
      k = 0;
      for (auto &h : horizons) {
        if (h > effectiveH)
          continue;
        for (Ti i = 0; i < temp.RowsCount; i++) {
          me.Set0(i, k, me.Get0(i, k) + temp.Get0(i, k));
          ResultAggs.Set0(c, i, ResultAggs.Get0(c, i) + temp.Get0(i, k));
        }
        k++;
      }

      // if (keepDetails) {
      //     temp.CopyTo(details->metrics.at(c));
      // }
      c++;
    }

    row++;
  }

  if (counter == 0 || invalidCounts > maxInvalidSim)
    throw std::logic_error("Model check failed: Minimum Valid Simulations");

  if (cancel)
    return;

  // average:
  for (Ti j = 0; j < hh; j++) {
    auto o = (Tv)1 / (Tv)counters[j];
    for (Ti i = 0; i < mm; i++) {
      auto m = &Results.at(i);
      if (metrics.at(i) == ScoringType::kRmse ||
          metrics.at(i) == ScoringType::kRmspe) {
        for (k = 0; k < m->RowsCount; k++)
          m->Set0(k, j, std::sqrt(m->Get0(k, j) * o));
      } else {
        for (k = 0; k < m->RowsCount; k++)
          m->Set0(k, j, m->Get0(k, j) * o);
      }

      // convert to percentage
      if (metrics.at(i) == ScoringType::kMape ||
          metrics.at(i) == ScoringType::kRmspe) {
        for (k = 0; k < m->RowsCount; k++)
          m->Set0(k, j, m->Get0(k, j) * 100);
      }
    }
  }

  if (cancel)
    return;

  // average (aggregate over horizons)
  for (Ti i = 0; i < mm; i++) {
    if (metrics.at(i) == ScoringType::kRmse ||
        metrics.at(i) == ScoringType::kRmspe) {
      for (k = 0; k < yy; k++)
        ResultAggs.Set0(i, k, std::sqrt(ResultAggs.Get0(i, k) / (Tv)counter));
    } else {
      for (k = 0; k < yy; k++)
        ResultAggs.Set0(i, k, ResultAggs.Get0(i, k) / (Tv)counter);
    }
  }
}

void VarmaSimulation::CalculateE(
    Tv *storage, Tv *work, Matrix<Tv> &data, Tv maxCondNum, Tv stdMultipler,
    bool coefUncer, bool usePreviousEstim,
    const std::function<void(Tv &)> *transformForMetrics) {
  auto horizons = *pHorizons;
  auto metrics = *pMetrics;
  auto sizes = *pSizes;
  Ti hh = (Ti)horizons.size();
  Ti hMin = horizons.at(0);
  Ti hMax = horizons.at(hh - 1);

  std::function<void(Tv &)> tfm;
  if (transformForMetrics)
    tfm = *transformForMetrics;

  // we should get t:
  auto D = DatasetTs<false>(data.RowsCount, data.ColsCount, true, false);
  D.Data(data);
  D.Update(nullptr, work);
  if (D.HasMissingData)
    throw std::logic_error("Missing data is found in VARMA data.");
  if (D.Start > D.End)
    throw std::logic_error("Data is not valid.");

  auto T = D.End - D.Start + 1;
  Ti mm = (Ti)metrics.size();
  Ti yy = sizes.EqsCount; // number of endogenous data

  // check size
  auto tempm = VarmaSimulation(sizes, mCount, horizons, metrics,
                               &EModel.Model.Result.Optim.Options, true,
                               EModel.mRestriction, EModel.mCheckNan,
                               EModel.pPcaOptionsY, EModel.pPcaOptionsX);
  if (WorkSize < tempm.WorkSize || StorageSize < tempm.StorageSize)
    throw std::logic_error("inconsistent arguments in VARMA simulation");
  tempm = VarmaSimulation();

  auto count = mCount;

  if (T - count <= 0)
    throw std::logic_error(
        "invalid number of simulations. It is larger than available data");
  if (count == 0)
    throw std::logic_error(
        "invalid number of simulations. It is zero of negative");

  // set storage
  Ti pos = 0;
  for (int i = 0; i < mm; i++) {
    Results.at(i).SetData(0, &storage[pos], yy, hh);
    pos += yy * hh; // horizons are in columns
  }
  ResultAggs.SetData(0, &storage[pos], mm, yy);
  pos += mm * yy;

  // set work
  pos = 0;
  double *estimStorage = nullptr;
  double *estimForecastWork = nullptr;

  estimStorage = &work[pos];
  pos += EModel.StorageSize;
  estimForecastWork = &work[pos];
  pos += EModel.WorkSize;

  //      others
  Ti co = yy * hh;
  auto act = Matrix<Tv>(&work[pos], yy, hh);
  pos += co;
  auto forc = Matrix<Tv>(&work[pos], yy, hh);
  pos += co;
  auto err = Matrix<Tv>(&work[pos], yy, hh);
  pos += co;
  auto std = Matrix<Tv>(&work[pos], yy, hh);
  pos += co;
  auto temp = Matrix<Tv>(&work[pos], yy, hh);
  pos += co;
  auto last = Matrix<Tv>(&work[pos], yy, 1);
  pos += yy;
  auto lastIndex = D.End;

  bool useCurrentEstims;
  Ti k;
  Ti forIndex;
  Ti colc;
  Ti actIndex;
  Ti effectiveH;
  bool success = false;
  Ti row = 0;
  auto counters = std::vector<Ti>(horizons.size());
  Ti counter = 0;
  for (Ti se = count; se > 0; se--) {
    useCurrentEstims = usePreviousEstim && success;
    actIndex = lastIndex - se;
    if (actIndex < 0)
      break;

    effectiveH = std::min(se, hMax);
    if (effectiveH < hMin)
      break;

    colc = 0;
    for (auto &h : horizons) {
      if (h <= effectiveH)
        colc++;
    }
    if (colc <= 0)
      break; // no more forecasts

    // estimate and forecast
    // update sizes by current sampleEnd. se is different which means Xt,
    // resid, etc are different
    auto forecast0 = &Forecast;
    try {
      EModel.Calculate(data, estimStorage, estimForecastWork, useCurrentEstims,
                       effectiveH, se, maxCondNum, stdMultipler);
      forecast0 = &EModel.Forecasts;
    } catch (std::exception &ex) {
      AddError(ex.what());
      row++;
      continue;
    } catch (std::string &ex) {
      AddError(ex.c_str());
      row++;
      continue;
    } catch (const char *ex) {
      AddError(ex);
      row++;
      continue;
    } catch (...) {
      AddError("unknown error!");
      row++;
      continue;
    }
    auto forecast = *forecast0;

    forIndex = forecast.StartIndex - 1;
    for (Ti b = 0; b < sizes.EqsCount; b++) // copy just the endogenous part
      last.Data[b] = data.Get0(actIndex, b);

    k = 0;
    for (auto &h : horizons) // set columns
    {
      if (h <= effectiveH) {
        counter++;
        counters[k]++;
        for (Ti b = 0; b < sizes.EqsCount; b++) // copy just the endogenous part
          act.Set0(b, k, data.Get0(actIndex + h, b));

        forc.SetColumnFromColumn0(k, forecast.Forecast, forIndex + h);
        if (mDoForecastVar)
          std.SetColumnFromColumn0(k, forecast.Variance, forIndex + h);
      } else {
        act.SetColumn0(k, NAN);
        forc.SetColumn0(k, NAN);
        std.SetColumn0(k, NAN);
      }
      k++;
    }

    // Convert to STD
    if (mDoForecastVar) {
      for (int i = 0; i < co; i++)
        std.Data[i] = std::sqrt(std.Data[i]);
    }

    if (transformForMetrics) {

      for (int i = 0; i < yy; i++)
        tfm(last.Data[i]);

      for (int i = 0; i < co; i++) {
        tfm(act.Data[i]);
        tfm(forc.Data[i]);
      }

      if (mDoForecastVar) { // transform STDs for CRPS too
        for (int i = 0; i < co; i++) {
          tfm(std.Data[i]);
        }
      }
    }

    act.Subtract0(forc, err);

    /*if (keepDetails) {
        details->sampleEnd = se;
        details->isValid = true;
        act.CopyTo00(details->actual);
        forc.CopyTo00(details->forecast);
        std.CopyTo00(details->stdforecast);
    }*/
    success = true;
    ValidCounts++;

    Ti c = 0;
    for (auto &eval : metrics) {
      GetScore(eval, temp, act, forc, err, std, last);

      // summation for calculating mean
      auto me = Results.at(c);
      k = 0;
      for (auto &h : horizons) {
        if (h > effectiveH)
          continue;
        if (KeepDetails) {
          for (int vi = 0; vi < yy; vi++)
            Details.push_back(std::make_tuple(
                se, c, h, vi, last.Data[vi], act.Get0(vi, k), forc.Get0(vi, k),
                err.Get0(vi, k), std.Get0(vi, k)));
        }
        for (Ti i = 0; i < temp.RowsCount; i++) {
          me.Set0(i, k, me.Get0(i, k) + temp.Get0(i, k));
          ResultAggs.Set0(c, i, ResultAggs.Get0(c, i) + temp.Get0(i, k));
        }
        k++;
      }

      // if (keepDetails) {
      //     temp.CopyTo(details->metrics.at(c));
      // }
      c++;
    }

    row++;
  }

  // average:
  for (Ti j = 0; j < hh; j++) {
    auto o = (Tv)1 / (Tv)counters[j];
    for (Ti i = 0; i < mm; i++) {
      auto m = &Results.at(i);
      if (metrics.at(i) == ScoringType::kRmse ||
          metrics.at(i) == ScoringType::kRmspe) {
        for (k = 0; k < m->RowsCount; k++)
          m->Set0(k, j, std::sqrt(m->Get0(k, j) * o));
      } else {
        for (k = 0; k < m->RowsCount; k++)
          m->Set0(k, j, m->Get0(k, j) * o);
      }
    }
  }

  // average (aggregate ove horizons)
  for (Ti i = 0; i < mm; i++) {
    if (metrics.at(i) == ScoringType::kRmse ||
        metrics.at(i) == ScoringType::kRmspe) {
      for (k = 0; k < yy; k++)
        ResultAggs.Set0(i, k, std::sqrt(ResultAggs.Get0(i, k) / (Tv)counter));
    } else {
      for (k = 0; k < yy; k++)
        ResultAggs.Set0(i, k, ResultAggs.Get0(i, k) / (Tv)counter);
    }
  }
}

void VarmaSimulation::AddError(std::string state) {
  if (state.empty())
    return;
  if (Errors.find(state) != Errors.end())
    Errors.at(state)++;
  else {
    Errors.insert(std::pair<std::string, Ti>(state, 1));
  }
}

// #pragma endregion
