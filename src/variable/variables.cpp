
#include "variable.h"
#include <random>

using namespace ldt;

template <typename Tw>
Variables<Tw>::Variables(const std::vector<Variable<Tw> *> vars) {

  if (vars.size() == 0)
    throw LdtException(ErrorType::kLogic, "variables",
                       "no variable is available");

  StartFrequency = vars.at(0)->StartFrequency.get()->Clone();
  auto maxEnds = vars.at(0)->GetEndFrequency();

  for (auto const &v : vars) {
    try {
      if (StartFrequency.get()->IsNewerThan(*v->StartFrequency.get()))
        StartFrequency = std::move(v->StartFrequency.get()->Clone());
      auto temp = v->GetEndFrequency();
      if (maxEnds->IsOlderThan(*temp.get()))
        maxEnds = std::move(temp);
    } catch (...) {

      try {
        std::rethrow_exception(std::current_exception());
      } catch (const std::exception &e) {
        throw LdtException(ErrorType::kLogic, "variables",
                           "mixed frequency is not supported in 'Variables'",
                           &e);
      }
    }
    Names.push_back(v->Name);
  }
  NumObs = maxEnds.get()->Minus(*StartFrequency.get());
  if (NumObs == 0)
    throw LdtException(ErrorType::kLogic, "variables",
                       "no observation is available");

  Data.resize(NumObs * vars.size());
  Ti i = 0;
  std::unique_ptr<Frequency> temp2;
  for (auto const &v : vars) {
    auto start = v->StartFrequency.get()->Minus(*StartFrequency.get());
    temp2 = v->GetEndFrequency();
    auto e = temp2.get();
    auto end = maxEnds->Minus(*e);
    for (Ti j = 0; j < start; j++)
      Data.at(i++) = NAN;
    for (Ti j = 0; j < (Ti)v->Data.size(); j++)
      Data.at(i++) = v->Data.at(j);
    for (Ti j = 0; j < end; j++)
      Data.at(i++) = NAN;
  }
}

template <typename Tw>
std::tuple<Ti, Ti> Variables<Tw>::GetRange(Ti j, bool &hasMissing) {

  if constexpr (std::numeric_limits<Tw>::has_quiet_NaN) {
    hasMissing = false;
    Ti start;
    Ti end;
    Tw *col = &Data[NumObs * j];
    for (start = 0; start < NumObs; start++)
      if (std::isnan(col[start]) == false)
        break;

    for (end = NumObs; end > 0; end--)
      if (std::isnan(col[end - 1]) == false) {
        end--;
        break;
      }

    auto range = std::tuple(start, end);

    if (start > end)
      return range;

    for (Ti i = start; i <= end; i++)
      if (std::isnan(col[i])) {
        hasMissing = true;
        break;
      }
    return range;
  } else if constexpr (true) {
    throw LdtException(ErrorType::kLogic, "variables",
                       "invalid operation"); // there is no NAN
  }
}

template class ldt::Variables<Tv>;
// template class ldt::Variable < Ti>;