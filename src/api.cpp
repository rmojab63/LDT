/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "api.h"

using namespace ldt;

// #pragma region Exception

void LDT_GetLastError(char *buffer) { std::strcpy(buffer, LastErrorMsg); }

void LDT_GetLastError_TEST() {
  API_BEGIN()
  throw std::logic_error("LDT ERROR TEST");
  API_END()
}

// #pragma endregion

// #pragma region Distribution

double LDT_GetDistributionProperty(int distributionType, int propertyType,
                                   double param1, double param2, double param3,
                                   double param4) {
  API_BEGIN()

  auto distType = static_cast<ldt::DistributionType>(distributionType);
  auto dist = ldt::DistributionBase::GetDistributionFromType(
      distType, param1, param2, param3, param4);
  auto propType = static_cast<ldt::DistributionProperty>(propertyType);
  return dist.get()->GetProperty(propType);

  API_END()
  return NAN;
}

void LDT_GetDistributionCDFs(int distributionType, double *probs,
                             int probsLength, double *result, double param1,
                             double param2, double param3, double param4) {

  API_BEGIN()

  auto distType = static_cast<ldt::DistributionType>(distributionType);
  auto dist = ldt::DistributionBase::GetDistributionFromType(
      distType, param1, param2, param3, param4);

  for (int i = 0; i < probsLength; i++)
    result[i] = dist.get()->GetQuantile(probs[i]);

  API_END()
}

// #pragma endregion

// #pragma region Frequency

ldt::Frequency *LDT_FrequencyCreate(const char *str, const char *classStr,
                                    int &freqClass) {
  API_BEGIN()

  auto fClass = ldt::FrequencyClass::kCrossSection;
  auto str0 = std::string(str);
  auto str1 = std::string(classStr);
  auto f = ldt::Frequency::Parse(str0, str1, fClass);
  freqClass = (int)fClass;
  if (fClass == ldt::FrequencyClass::kListString) {
    auto fList = new ldt::FrequencyList<std::string>("", nullptr);
    auto items = new std::vector<std::string>();
    try {
      ldt::FrequencyList<std::string>::Parse0(str0, str1, fClass, *fList,
                                              items);
    } catch (...) {
      delete fList;
      delete items;
      throw;
    }
    fList->pItems = items;
    return fList;
  } else if (fClass == ldt::FrequencyClass::kListDate) {
    auto fList = new ldt::FrequencyList<boost::gregorian::date>(
        boost::gregorian::date(), nullptr);
    auto items = new std::vector<boost::gregorian::date>();
    try {
      ldt::FrequencyList<boost::gregorian::date>::Parse0(str0, str1, fClass,
                                                         *fList, items);
    } catch (...) {
      delete fList;
      delete items;
      throw;
    }
    fList->pItems = items;
    return fList;
  } else {
    return f.release();
  }

  API_END()
  return nullptr;
}

ldt::Frequency *LDT_FrequencyCopy(ldt::Frequency *f) {
  API_BEGIN()

  auto newf = f->Clone();
  return newf.release();

  API_END()
  return nullptr;
}

void LDT_FrequencyNext(ldt::Frequency *f, int steps) {
  API_BEGIN()

  f->Next(steps);

  API_END()
}

void LDT_FrequencyToString(ldt::Frequency *f, char *asString,
                           char *asClassString, int &fClass, int &length) {
  API_BEGIN()

  auto a1 = f->ToString();
  auto a2 = f->ToClassString(true);
  auto reqLength = std::max(a1.length(), a2.length()) + 1;
  if (length < (int)reqLength) { // return to get new buffer
    length = reqLength;
    return;
  }

  fClass = (int)f->mClass;
  std::strcpy(asString, a1.c_str());
  std::strcpy(asClassString, a2.c_str());

  API_END()
}

void LDT_FrequencyToString0(ldt::Frequency *f, char *asString, int &length) {
  API_BEGIN()

  auto a1 = f->ToString();
  auto reqLength = a1.length() + 1;
  if (length < (int)reqLength) { // return to get new buffer
    length = reqLength;
    return;
  }

  std::strcpy(asString, a1.c_str());

  API_END()
}

int LDT_FrequencyMinus(ldt::Frequency *f, ldt::Frequency *g) {
  API_BEGIN()

  return f->Minus(*g);

  API_END()
  return -1;
}

void LDT_DisposeFrequency(ldt::Frequency *f) {
  API_BEGIN()

  if (f != NULL) {
    if (f->mClass == ldt::FrequencyClass::kListString) {
      auto fList = dynamic_cast<ldt::FrequencyList<std::string> const &>(*f);
      if (fList.pItems)
        delete fList.pItems;
    } else if (f->mClass == ldt::FrequencyClass::kListDate) {
      auto fList =
          dynamic_cast<ldt::FrequencyList<boost::gregorian::date> const &>(*f);
      if (fList.pItems)
        delete fList.pItems;
    }
    delete f;
    f = nullptr;
  }

  API_END()
}

// #pragma endregion

// #pragma region Variable

ldt::Frequency *LDT_VariableCombine(double *data1, int dataLength1,
                                    ldt::Frequency *startFrequency1,
                                    double *data2, int dataLength2,
                                    ldt::Frequency *startFrequency2,
                                    double *result, int &resultLength,
                                    double &maxDiff_perc) {
  API_BEGIN()
  maxDiff_perc = 0;

  auto v1 = ldt::Variable<double>();
  v1.Data.assign(data1, data1 + dataLength1);
  v1.StartFrequency = std::move(startFrequency1->Clone());
  v1.Trim();

  auto v2 = ldt::Variable<double>();
  v2.Data.assign(data2, data2 + dataLength2);
  v2.StartFrequency = std::move(startFrequency2->Clone());
  v2.Trim();

  auto vars = std::vector<Variable<double> *>({&v1, &v2});

  auto vs = ldt::Variables(vars);

  if (vs.NumObs > resultLength) {
    resultLength = vs.NumObs;
    return nullptr;
  }
  resultLength = vs.NumObs;
  auto mat = ldt::Matrix<double>(&vs.Data.at(0), vs.NumObs, 2);
  for (int i = 0; i < mat.RowsCount; i++) {
    auto d1 = mat.Get0(i, 0);
    auto d2 = mat.Get0(i, 1);

    if (std::isnan(d1) == false && std::isnan(d2) == false) {
      double diff = std::fabs((d1 - d2) / d1) *
                    100; // how the second one changed relative to the first one
      if (diff > maxDiff_perc)
        maxDiff_perc = diff;
    }
    result[i] =
        std::isnan(d2) ? d1 : d2; // priority is with the second variable. we
                                  // select it unless it is NAN
  }

  auto f = vs.StartFrequency.get();
  vs.StartFrequency.release();

  return f;

  API_END()

  return nullptr;
}

// #pragma endregion

// #pragma region Scoring

double LDT_GetScoreCrpsNormal(double y, double mean, double variance) {
  API_BEGIN()

  return Scoring::GetScoreCrpsNormal(y, mean, std::sqrt(variance));

  API_END()
  return NAN;
}

double LDT_GetScoreCrpsLogNormal(double y, double meanLog, double varianceLog) {
  API_BEGIN()

  return Scoring::GetScoreCrpsLogNormal(y, meanLog, std::sqrt(varianceLog));

  API_END()
  return NAN;
}

// #pragma endregion
