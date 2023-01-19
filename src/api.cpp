/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "api.h"

using namespace ldt;

//#pragma region Exception

void LDT_GetLastError(char *buffer) { std::strcpy(buffer, LastErrorMsg); }

void LDT_GetLastError_TEST() {
  API_BEGIN()
  throw std::logic_error("LDT ERROR TEST");
  API_END()
}

//#pragma endregion

//#pragma region Distribution

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

//#pragma endregion

//#pragma region Frequency

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

//#pragma endregion

//#pragma region Variable

ldt::Variable<double> *LDT_VariableCreate(const char *name, const double *data,
                                          int dataLength,
                                          ldt::Frequency *startFrequency) {
  API_BEGIN()

  auto v = new ldt::Variable<double>();
  try {
    v->Name = std::string(name);
    v->Data = std::vector<double>(data, data + dataLength);
    v->StartFrequency = startFrequency->Clone(); // copy, new owner
  } catch (...) {
    delete v;
    throw;
  }
  return v;

  API_END()
  return nullptr;
}

void LDT_VariableAddNewField(ldt::Variable<double> *v, const char *key,
                             const char *value) {
  API_BEGIN()

  v->Fields.insert({std::string(key), std::string(value)});

  API_END()
}

void LDT_VariableGetField(ldt::Variable<double> *v, int index, char *key,
                          char *value) {
  API_BEGIN()

  int i = -1;
  for (auto const &x : v->Fields) {
    i++;
    if (i == index) {
      std::strcpy(key, x.first.c_str());
      std::strcpy(value, x.second.c_str());
      return;
    }
  }

  API_END()
}

ldt::Frequency *LDT_VariableGetData(ldt::Variable<double> *v, char *name,
                                    double *data, int dataLength) {

  API_BEGIN()

  std::strcpy(name, v->Name.c_str());

  for (int i = 0; i < dataLength; i++)
    data[i] = v->Data.at(i);

  return v->StartFrequency.get();

  API_END()

  return nullptr;
}

void LDT_VariableToString(ldt::Variable<double> *v, char *result, int &length) {
  API_BEGIN()

  auto a1 = v->ToString();
  auto reqLength = a1.length() + 1;
  if (length < (int)reqLength) { // return to get new buffer
    length = reqLength;
    return;
  }

  std::strcpy(result, a1.c_str());

  API_END()
}

ldt::Variable<double> *LDT_VariableParse(const char *str, int &dataLength,
                                         int &nameLength, int &fieldsCount,
                                         int &maxFieldKeyLength,
                                         int &maxFieldValueLength) {
  API_BEGIN()

  auto v = new ldt::Variable<double>();
  auto listItemsString = new std::vector<std::string>();
  auto listItemsDate = new std::vector<boost::gregorian::date>();
  try {
    ldt::Variable<double>::Parse(str, *v, *listItemsString, *listItemsDate);
    if (v->StartFrequency.get()->mClass != ldt::FrequencyClass::kListString)
      delete listItemsString;
    if (v->StartFrequency.get()->mClass == ldt::FrequencyClass::kListDate)
      delete listItemsDate;

    dataLength = v->Data.size();
    nameLength = v->Name.length();
    fieldsCount = v->Fields.size();

    return v;
  } catch (...) {
    delete v;
    delete listItemsString;
    delete listItemsDate;
  }

  API_END()

  return nullptr;
}

void LDT_VariableDispose(ldt::Variable<double> *v) {
  API_BEGIN()

  if (v != NULL) {
    auto f = v->StartFrequency.get();
    if (f) {
      if (f->mClass != ldt::FrequencyClass::kListString) {
        auto fList = dynamic_cast<ldt::FrequencyList<std::string> const &>(*f);
        if (fList.pItems)
          delete fList.pItems;
      } else if (f->mClass != ldt::FrequencyClass::kListDate) {
        auto fList =
            dynamic_cast<ldt::FrequencyList<boost::gregorian::date> const &>(
                *f);
        if (fList.pItems)
          delete fList.pItems;
      }
    }

    delete v;
    v = nullptr;
  }

  API_END()
}

//#pragma endregion

//#pragma region Scoring

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

//#pragma endregion
