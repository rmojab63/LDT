#pragma once

#include "distributions.h"
#include "frequency.h"
#include "ldt_base.h"
#include "scoring.h"
#include "variable.h"
#include "varma.h"

#define LDT_C_EXPORT extern "C" LDT_EXPORT
// #pragma region Exception

static char LastErrorMsg[2048] = ""; // thread_local?

LDT_C_EXPORT void LDT_GetLastError(char *buffer);

LDT_C_EXPORT void LDT_GetLastError_TEST();

#define API_BEGIN()                                                            \
  try {                                                                        \
    snprintf(LastErrorMsg, 2048, "%s", "");
#define API_END()                                                              \
  }                                                                            \
  catch (std::exception & ex) {                                                \
    snprintf(LastErrorMsg, 2048, "%s", ex.what());                             \
  }                                                                            \
  catch (std::string & ex) {                                                   \
    snprintf(LastErrorMsg, 2048, "%s", ex.c_str());                            \
  }                                                                            \
  catch (const char *ex) {                                                     \
    snprintf(LastErrorMsg, 2048, "%s", ex);                                    \
  }                                                                            \
  catch (...) {                                                                \
    snprintf(LastErrorMsg, 2048, "%s", "Unknown error occurred");              \
  }

// #pragma endregion

// #pragma region Distribution

LDT_C_EXPORT double LDT_GetDistributionProperty(int distributionType,
                                                int propertyType, double param1,
                                                double param2, double param3,
                                                double param4);

// #pragma endregion

// #pragma region Frequency

LDT_C_EXPORT ldt::Frequency *
LDT_FrequencyCreate(const char *str, const char *classStr, int &freqClass);

LDT_C_EXPORT ldt::Frequency *LDT_FrequencyCopy(ldt::Frequency *f);

LDT_C_EXPORT void LDT_FrequencyNext(ldt::Frequency *f, int steps);

LDT_C_EXPORT void LDT_FrequencyToString(ldt::Frequency *f, char *asString,
                                        char *asClassString, int &fClass,
                                        int &length);

LDT_C_EXPORT void LDT_FrequencyToString0(ldt::Frequency *f, char *asString,
                                         int &length);

LDT_C_EXPORT int LDT_FrequencyMinus(ldt::Frequency *f, ldt::Frequency *g);

LDT_C_EXPORT void LDT_DisposeFrequency(ldt::Frequency *f);

// #pragma endregion

// #pragma region Variable

LDT_C_EXPORT ldt::Variable<double> *
LDT_VariableCreate(const char *name, const double *data, int dataLength,
                   ldt::Frequency *startFrequency);

LDT_C_EXPORT ldt::Frequency *LDT_VariableGetFrequency(ldt::Variable<double> *v);

LDT_C_EXPORT void LDT_VariableAddNewField(ldt::Variable<double> *v,
                                          const char *key, const char *value);

LDT_C_EXPORT void LDT_VariableGetField(ldt::Variable<double> *v, int index,
                                       char *key, char *value);

LDT_C_EXPORT ldt::Frequency *LDT_VariableGetData(ldt::Variable<double> *v,
                                                 char *name, double *data,
                                                 int dataLength);

LDT_C_EXPORT void LDT_VariableToString(ldt::Variable<double> *v, char *result,
                                       int &length);

LDT_C_EXPORT ldt::Variable<double> *
LDT_VariableParse(const char *str, int &dataLength, int &nameLength,
                  int &fieldsCount, int &maxFieldKeyLength,
                  int &maxFieldValueLength);

LDT_C_EXPORT void LDT_VariableDispose(ldt::Variable<double> *v);

// #pragma endregion

LDT_C_EXPORT double LDT_GetScoreCrpsNormal(double y, double mean,
                                           double variance);

LDT_C_EXPORT double LDT_GetScoreCrpsLogNormal(double y, double meanLog,
                                              double varianceLog);

// #pragma endregion

// #pragma region Forecast

LDT_C_EXPORT void LDT_Forecast(double *data, int dataRows, int dataCols,
                               int exoStart, int horizon, bool &cancel) {
  throw std::logic_error("Not implemented");

  /*
    auto source = Matrix<Tv>(data, dataRows, dataCols);
    source.Transpose();

    auto options = SearchOptions();
    options.Parallel = true;

    // measures
    auto measures = SearchMeasureOptions();
    for (auto i = 1; i <= horizon; i++)
      measures.Horizons.push_back(i);
    measures.MeasuresIn.push_back(ldt::GoodnessOfFitType::kAic);
    measures.MeasuresOut = std::vector<ldt::ScoringType>(
        {ScoringType::kDirection, ScoringType::kScaledRmse,
    ScoringType::kCrps}); measures.SimFixSize = 0; measures.TrainRatio = 0.7;
    */
}

// #pragma endregion