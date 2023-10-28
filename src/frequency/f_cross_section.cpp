/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "frequency.h"

using namespace ldt;

FrequencyCrossSection::FrequencyCrossSection(Ti position) {
  this->mClass = FrequencyClass::kCrossSection;
  mPosition = position;
}

std::unique_ptr<Frequency> FrequencyCrossSection::Clone() const {
  return std::make_unique<FrequencyCrossSection>(*this);
}

void FrequencyCrossSection::Next(Ti steps) { mPosition += steps; }

Ti FrequencyCrossSection::CompareTo(Frequency const &other) {
  CheckClassEquality(*this, other);
  auto second = dynamic_cast<FrequencyCrossSection const &>(other);
  if (mPosition > second.mPosition)
    return 1;
  else if (mPosition < second.mPosition)
    return -1;
  return 0;
}

Ti FrequencyCrossSection::Minus(Frequency const &other) {
  CheckClassEquality(*this, other);
  auto second = dynamic_cast<FrequencyCrossSection const &>(other);
  return mPosition - second.mPosition;
}

void FrequencyCrossSection::Parse0(const std::string &str,
                                   FrequencyCrossSection &result) {
  try {
    result.mPosition = std::stoi(str, nullptr, 10);
  } catch (...) {

    try {
      std::rethrow_exception(std::current_exception());
    } catch (const std::exception &e) {
      throw LdtException(
          ErrorType::kLogic, "freq-cs",
          "Parsing cross-section frequency failed. Invalid format.", &e);
    }
  }
}

std::string FrequencyCrossSection::ToString() const {
  return std::to_string(mPosition);
}

std::string FrequencyCrossSection::ToClassString(bool details) const {
  return "cs";
}