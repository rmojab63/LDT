
#include "variable.h"
#include <random>

using namespace ldt;

template <typename Tw>
Variables<Tw>::Variables(const std::vector<Variable<Tw> *> vars) {

  if (vars.size() == 0)
    throw std::logic_error("Variables: No variable is available.");

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
      Rethrow("Mixed frequency is not supported in 'Variables'.");
    }
    Names.push_back(v->Name);
  }
  NumObs = maxEnds.get()->Minus(*StartFrequency.get());
  if (NumObs == 0)
    throw std::logic_error("Variables: No observation is available.");

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

template class ldt::Variables<Tv>;
// template class ldt::Variable < Ti>;