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

LDT_C_EXPORT void LDT_GetDistributionQuantiles(int distributionType,
                                               double *probs, int probsLength,
                                               double *result, double param1,
                                               double param2, double param3,
                                               double param4);

LDT_C_EXPORT std::vector<double *> *
LDT_CombineDistribution103(std::vector<double *> *dists, double *dist,
                           double *weights, double *result);

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

LDT_C_EXPORT ldt::Frequency *
LDT_VariableCombine(double *data1, int dataLength1,
                    ldt::Frequency *startFrequency1, double *data2,
                    int dataLength2, ldt::Frequency *startFrequency2,
                    double *result, int &resultLength, double &maxDiff_perc);

// #pragma endregion

LDT_C_EXPORT double LDT_GetScoreCrpsNormal(double y, double mean,
                                           double variance);

LDT_C_EXPORT double LDT_GetScoreCrpsLogNormal(double y, double meanLog,
                                              double varianceLog);

// #pragma endregion

// #pragma region Forecast

LDT_C_EXPORT void LDT_Forecast(double *data, int dataRows, int dataCols,
                               int exoStart, int horizon, bool &cancel) {
  throw ldt::LdtException(ldt::ErrorType::kLogic, "api", "not implemented");

  /*
    auto source = Matrix<Tv>(data, dataRows, dataCols);
    source.Transpose();

    auto options = SearchOptions();
    options.Parallel = true;

    // metrics
    auto metrics = SearchMetricOptions();
    for (auto i = 1; i <= horizon; i++)
      metrics.Horizons.push_back(i);
    metrics.MetricsIn.push_back(ldt::GoodnessOfFitType::kAic);
    metrics.MetricsOut = std::vector<ldt::ScoringType>(
        {ScoringType::kDirection, ScoringType::kRmspe,
    ScoringType::kCrps}); metrics.SimFixSize = 0; metrics.TrainRatio = 0.7;
    */
}

// #pragma endregion