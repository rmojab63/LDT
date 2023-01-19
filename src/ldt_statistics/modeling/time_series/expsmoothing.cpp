/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "scoring.h"
#include "smoothing.h"

using namespace ldt;

//#pragma region export

//#pragma endregion

template <ExpSmoothingType type> int expsmoothing<type>::check_getexpclass() {
  if (_additive_season != Null && _seasonCount == 0)
    throw std::logic_error("zero number of seasons in a seasonal model");

  if constexpr (type == ExpSmoothingType::aaa) {
    return 1;
  } else if constexpr (type == ExpSmoothingType::aam) {
    return 5;
  } else if constexpr (type == ExpSmoothingType::aan) {
    return 1;
  } else if constexpr (type == ExpSmoothingType::aa_da) {
    return 1;
  } else if constexpr (type == ExpSmoothingType::aa_dm) {
    return 5;
  } else if constexpr (type == ExpSmoothingType::aa_dn) {
    return 1;
  } else if constexpr (type == ExpSmoothingType::ama) {
    return 5;
  } else if constexpr (type == ExpSmoothingType::amm) {
    return 5;
  } else if constexpr (type == ExpSmoothingType::amn) {
    return 5;
  } else if constexpr (type == ExpSmoothingType::am_da) {
    return 5;
  } else if constexpr (type == ExpSmoothingType::am_dm) {
    return 5;
  } else if constexpr (type == ExpSmoothingType::am_dn) {
    return 5;
  } else if constexpr (type == ExpSmoothingType::ana) {
    return 1;
  } else if constexpr (type == ExpSmoothingType::anm) {
    return 5;
  } else if constexpr (type == ExpSmoothingType::ann) {
    return 1;
  } else if constexpr (type == ExpSmoothingType::maa) {
    return 2;
  } else if constexpr (type == ExpSmoothingType::mam) {
    return 3;
  } else if constexpr (type == ExpSmoothingType::man) {
    return 2;
  } else if constexpr (type == ExpSmoothingType::ma_da) {
    return 2;
  } else if constexpr (type == ExpSmoothingType::ma_dm) {
    return 3;
  } else if constexpr (type == ExpSmoothingType::ma_dn) {
    return 2;
  } else if constexpr (type == ExpSmoothingType::mma) {
    return 5;
  } else if constexpr (type == ExpSmoothingType::mmm) {
    return 4;
  } else if constexpr (type == ExpSmoothingType::mmn) {
    return 4;
  } else if constexpr (type == ExpSmoothingType::mm_da) {
    return 5;
  } else if constexpr (type == ExpSmoothingType::mm_dm) {
    return 4;
  } else if constexpr (type == ExpSmoothingType::mm_dn) {
    return 4;
  } else if constexpr (type == ExpSmoothingType::mna) {
    return 2;
  } else if constexpr (type == ExpSmoothingType::mnm) {
    return 3;
  } else if constexpr (type == ExpSmoothingType::mnn) {
    return 2;
  } else if constexpr (true) {
    throw std::logic_error("not implemented");
  }
}

template <ExpSmoothingType type>
expsmoothing<type>::expsmoothing(unsigned short seasonCount) {
  if (seasonCount < 2)
    seasonCount = 0;
  _type = type;

  if constexpr (type == ExpSmoothingType::aaa) {
    _additive_season = True;
    _additive_trend = True;
  } else if constexpr (type == ExpSmoothingType::aam) {
    _additive_season = False;
    _additive_trend = True;
  } else if constexpr (type == ExpSmoothingType::aan) {
    _additive_trend = True;
  } else if constexpr (type == ExpSmoothingType::aa_da) {
    _additive_season = True;
    _additive_trend = True;
    _damped = true;
  } else if constexpr (type == ExpSmoothingType::aa_dm) {
    _additive_season = False;
    _additive_trend = True;
    _damped = true;
  } else if constexpr (type == ExpSmoothingType::aa_dn) {
    _additive_trend = True;
    _damped = true;
  } else if constexpr (type == ExpSmoothingType::ama) {
    _additive_season = True;
    _additive_trend = False;
  } else if constexpr (type == ExpSmoothingType::amm) {
    _additive_season = False;
    _additive_trend = False;
  } else if constexpr (type == ExpSmoothingType::amn) {
    _additive_trend = False;
  } else if constexpr (type == ExpSmoothingType::am_da) {
    _additive_season = True;
    _additive_trend = False;
    _damped = true;
  } else if constexpr (type == ExpSmoothingType::am_dm) {
    _additive_season = False;
    _additive_trend = False;
    _damped = true;
  } else if constexpr (type == ExpSmoothingType::am_dn) {
    _additive_trend = False;
    _damped = true;
  } else if constexpr (type == ExpSmoothingType::ana) {
    _additive_season = True;
  } else if constexpr (type == ExpSmoothingType::anm) {
    _additive_season = False;
  } else if constexpr (type == ExpSmoothingType::ann) {
    //_additive_trend = False;
    //_additive_season = False;
    _additive_error = true;
  } else if constexpr (type == ExpSmoothingType::maa) {
    _additive_season = True;
    _additive_trend = True;
    _additive_error = false;
  } else if constexpr (type == ExpSmoothingType::mam) {
    _additive_season = False;
    _additive_trend = True;
    _additive_error = false;
  } else if constexpr (type == ExpSmoothingType::man) {
    _additive_trend = True;
    _additive_error = false;
  } else if constexpr (type == ExpSmoothingType::ma_da) {
    _additive_season = True;
    _additive_trend = True;
    _additive_error = false;
    _damped = true;
  } else if constexpr (type == ExpSmoothingType::ma_dm) {
    _additive_season = False;
    _additive_trend = True;
    _additive_error = false;
    _damped = true;
  } else if constexpr (type == ExpSmoothingType::ma_dn) {
    _additive_trend = True;
    _additive_error = false;
    _damped = true;
  } else if constexpr (type == ExpSmoothingType::mma) {
    _additive_season = True;
    _additive_trend = False;
    _additive_error = false;
  } else if constexpr (type == ExpSmoothingType::mmm) {
    _additive_season = False;
    _additive_trend = False;
    _additive_error = false;
  } else if constexpr (type == ExpSmoothingType::mmn) {
    _additive_trend = False;
    _additive_error = false;
  } else if constexpr (type == ExpSmoothingType::mm_da) {
    _additive_season = True;
    _additive_trend = False;
    _additive_error = false;
    _damped = true;
  } else if constexpr (type == ExpSmoothingType::mm_dm) {
    _additive_season = False;
    _additive_trend = False;
    _additive_error = false;
    _damped = true;
  } else if constexpr (type == ExpSmoothingType::mm_dn) {
    _additive_trend = False;
    _additive_error = false;
    _damped = true;
  } else if constexpr (type == ExpSmoothingType::mna) {
    _additive_season = True;
    _additive_error = false;
  } else if constexpr (type == ExpSmoothingType::mnm) {
    _additive_season = False;
    _additive_error = false;
  } else if constexpr (type == ExpSmoothingType::mnn) {
    _additive_error = false;
  } else if constexpr (true) {
    throw std::logic_error("not implemented");
  }

  _seasonCount = _additive_season != Null ? seasonCount : 0;
  _m = (unsigned short)1 +
       (_additive_trend != Null ? (unsigned short)1 : (unsigned short)0) +
       (_seasonCount == 0 ? 0 : (_seasonCount - 1));
  _k = (unsigned short)1 +
       (_additive_trend != Null
            ? (_damped ? (unsigned short)2 : (unsigned short)1)
            : (unsigned short)0) +
       (_additive_season != Null ? (unsigned short)1 : (unsigned short)0);

  _class = check_getexpclass();
}

expsmoothing_base *expsmoothing_base::get_model(ExpSmoothingType type,
                                                unsigned short seasonCount) {

  switch (type) {
  case ldt::ExpSmoothingType::ann:
    return new expsmoothing<ExpSmoothingType::ann>(seasonCount);
  case ldt::ExpSmoothingType::aan:
    return new expsmoothing<ExpSmoothingType::aan>(seasonCount);
  case ldt::ExpSmoothingType::amn:
    return new expsmoothing<ExpSmoothingType::amn>(seasonCount);
  case ldt::ExpSmoothingType::ana:
    return new expsmoothing<ExpSmoothingType::ana>(seasonCount);
  case ldt::ExpSmoothingType::aaa:
    return new expsmoothing<ExpSmoothingType::aaa>(seasonCount);
  case ldt::ExpSmoothingType::ama:
    return new expsmoothing<ExpSmoothingType::ama>(seasonCount);
  case ldt::ExpSmoothingType::anm:
    return new expsmoothing<ExpSmoothingType::anm>(seasonCount);
  case ldt::ExpSmoothingType::aam:
    return new expsmoothing<ExpSmoothingType::aam>(seasonCount);
  case ldt::ExpSmoothingType::amm:
    return new expsmoothing<ExpSmoothingType::amm>(seasonCount);
  case ldt::ExpSmoothingType::aa_dn:
    return new expsmoothing<ExpSmoothingType::aa_dn>(seasonCount);
  case ldt::ExpSmoothingType::am_dn:
    return new expsmoothing<ExpSmoothingType::am_dn>(seasonCount);
  case ldt::ExpSmoothingType::aa_da:
    return new expsmoothing<ExpSmoothingType::aa_da>(seasonCount);
  case ldt::ExpSmoothingType::am_da:
    return new expsmoothing<ExpSmoothingType::am_da>(seasonCount);
  case ldt::ExpSmoothingType::aa_dm:
    return new expsmoothing<ExpSmoothingType::aa_dm>(seasonCount);
  case ldt::ExpSmoothingType::am_dm:
    return new expsmoothing<ExpSmoothingType::am_dm>(seasonCount);
  case ldt::ExpSmoothingType::mnn:
    return new expsmoothing<ExpSmoothingType::mnn>(seasonCount);
  case ldt::ExpSmoothingType::man:
    return new expsmoothing<ExpSmoothingType::man>(seasonCount);
  case ldt::ExpSmoothingType::mmn:
    return new expsmoothing<ExpSmoothingType::mmn>(seasonCount);
  case ldt::ExpSmoothingType::mna:
    return new expsmoothing<ExpSmoothingType::mna>(seasonCount);
  case ldt::ExpSmoothingType::maa:
    return new expsmoothing<ExpSmoothingType::maa>(seasonCount);
  case ldt::ExpSmoothingType::mma:
    return new expsmoothing<ExpSmoothingType::mma>(seasonCount);
  case ldt::ExpSmoothingType::mnm:
    return new expsmoothing<ExpSmoothingType::mnm>(seasonCount);
  case ldt::ExpSmoothingType::mam:
    return new expsmoothing<ExpSmoothingType::mam>(seasonCount);
  case ldt::ExpSmoothingType::mmm:
    return new expsmoothing<ExpSmoothingType::mmm>(seasonCount);
  case ldt::ExpSmoothingType::ma_dn:
    return new expsmoothing<ExpSmoothingType::ma_dn>(seasonCount);
  case ldt::ExpSmoothingType::mm_dn:
    return new expsmoothing<ExpSmoothingType::mm_dn>(seasonCount);
  case ldt::ExpSmoothingType::ma_da:
    return new expsmoothing<ExpSmoothingType::ma_da>(seasonCount);
  case ldt::ExpSmoothingType::mm_da:
    return new expsmoothing<ExpSmoothingType::mm_da>(seasonCount);
  case ldt::ExpSmoothingType::ma_dm:
    return new expsmoothing<ExpSmoothingType::ma_dm>(seasonCount);
  case ldt::ExpSmoothingType::mm_dm:
    return new expsmoothing<ExpSmoothingType::mm_dm>(seasonCount);
  default:
    throw std::logic_error("not implemented");
  }
}

expsmoothing_baseTv *expsmoothing_base::get_model(bool additiveError,
                                                  bool3 additiveTrend,
                                                  bool3 additiveSeason,
                                                  unsigned short seasonCount,
                                                  bool damped) {
  if (seasonCount < 2)
    seasonCount = 0;

  ExpSmoothingType type = ExpSmoothingType::aaa;

  if (additiveError) {
    if (additiveTrend == Null) {
      if (additiveSeason == Null) {
        type = ExpSmoothingType::ann;
      } else if (additiveSeason == True) {
        type = ExpSmoothingType::ana;
      } else {
        type = ExpSmoothingType::anm;
      }
    } else if (additiveTrend == True) {
      if (additiveSeason == Null) {
        type = damped ? ExpSmoothingType::aa_dn : ExpSmoothingType::aan;
      } else if (additiveSeason == True) {
        type = damped ? ExpSmoothingType::aa_da : ExpSmoothingType::aaa;
      } else {
        type = damped ? ExpSmoothingType::aa_dm : ExpSmoothingType::aam;
      }
    } else {
      if (additiveSeason == Null) {
        type = damped ? ExpSmoothingType::am_dn : ExpSmoothingType::amn;
      } else if (additiveSeason == True) {
        type = damped ? ExpSmoothingType::am_da : ExpSmoothingType::ama;
      } else {
        type = damped ? ExpSmoothingType::am_dm : ExpSmoothingType::amm;
      }
    }
  } else {
    if (additiveTrend == Null) {
      if (additiveSeason == Null) {
        type = ExpSmoothingType::mnn;
      } else if (additiveSeason == True) {
        type = ExpSmoothingType::mna;
      } else {
        type = ExpSmoothingType::mnm;
      }
    } else if (additiveTrend == True) {
      if (additiveSeason == Null) {
        type = damped ? ExpSmoothingType::ma_dn : ExpSmoothingType::man;
      } else if (additiveSeason == True) {
        type = damped ? ExpSmoothingType::ma_da : ExpSmoothingType::maa;
      } else {
        type = damped ? ExpSmoothingType::ma_dm : ExpSmoothingType::mam;
      }
    } else {
      if (additiveSeason == Null) {
        type = damped ? ExpSmoothingType::mm_dn : ExpSmoothingType::mmn;
      } else if (additiveSeason == True) {
        type = damped ? ExpSmoothingType::mm_da : ExpSmoothingType::mma;
      } else {
        type = damped ? ExpSmoothingType::mm_dm : ExpSmoothingType::mmm;
      }
    }
  }

  return get_model(type, seasonCount);
}

void setLastSeason(bool additive, Tv *seasons, unsigned short &seasonsCount) {
  if (seasonsCount == 2) {
    seasons[1] = additive ? (-seasons[0]) : ((Tv)1 / seasons[0]);
  } else {
    Tv sum = 0;
    for (Ti i = 0; i < seasonsCount - 1; i++)
      sum += seasons[i];
    seasons[seasonsCount - 1] = additive ? (-sum) : (Tv)1 / sum;
  }
}

//#pragma region recursives A

// table 2.2, hyndman, et ct (2008)

void rec_ANN(bool sim, Tv &y, Tv alpha, Tv &level, Tv &e, Tv &mu) {
  mu = level;
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  level += alpha * e;
}

void rec_AAN(bool sim, Tv &y, Tv alpha, Tv beta, Tv &level, Tv &trend, Tv &e,
             Tv &mu) {
  mu = level + trend;
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  level += trend + alpha * e;
  trend += beta * e;
}

void rec_AA_dN(bool sim, Tv &y, Tv alpha, Tv beta, Tv phi, Tv &level, Tv &trend,
               Tv &e, Tv &mu) {
  Tv p = phi * trend;
  mu = level + phi * trend;
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  level = mu + alpha * e;
  trend = p + beta * e;
}

void rec_AMN(bool sim, Tv &y, Tv alpha, Tv beta, Tv &level, Tv &trend, Tv &e,
             Tv &mu) {
  mu = level * trend;
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  Tv temp = level;
  level = mu + alpha * e;
  trend += beta * e / temp;
}

void rec_AM_dN(bool sim, Tv &y, Tv alpha, Tv beta, Tv phi, Tv &level, Tv &trend,
               Tv &e, Tv &mu) {
  Tv p = std::pow(trend, phi);
  mu = level * p;
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  Tv temp = level;
  level = mu + alpha * e;
  trend = p + beta * e / temp;
}

// Next Column

void rec_ANA(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
             Tv gamma, Tv &level, Tv *seasons, Tv &e, Tv &mu) {
  Ti rem = t % seasonCount;
  mu = level + seasons[rem];
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  level += alpha * e;
  seasons[rem] += gamma * e;
}

void rec_AAA(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
             Tv beta, Tv gamma, Tv &level, Tv &trend, Tv *seasons, Tv &e,
             Tv &mu) {
  Ti rem = t % seasonCount;
  mu = level + trend + seasons[rem];
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  level += trend + alpha * e;
  trend += beta * e;
  seasons[rem] += gamma * e;
}

void rec_AA_dA(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
               Tv beta, Tv phi, Tv gamma, Tv &level, Tv &trend, Tv *seasons,
               Tv &e, Tv &mu) {
  Ti rem = t % seasonCount;
  Tv p = phi * trend;
  mu = level + p + seasons[rem];
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  level += p + alpha * e;
  trend = p + beta * e;
  seasons[rem] += gamma * e;
}

void rec_AMA(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
             Tv beta, Tv gamma, Tv &level, Tv &trend, Tv *seasons, Tv &e,
             Tv &mu) {
  Ti rem = t % seasonCount;
  mu = level * trend + seasons[rem];
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  Tv temp = level;
  level = level * trend + alpha * e;
  trend += beta * e / temp;
  seasons[rem] += gamma * e;
}

void rec_AM_dA(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
               Tv beta, Tv phi, Tv gamma, Tv &level, Tv &trend, Tv *seasons,
               Tv &e, Tv &mu) {
  Ti rem = t % seasonCount;
  Tv p = std::pow(trend, phi);
  mu = level * p + seasons[rem];
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  Tv temp = level;
  level = level * p + alpha * e;
  trend = p + beta * e / temp;
  seasons[rem] += gamma * e;
}

// THIRD COLUMN

void rec_ANM(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
             Tv gamma, Tv &level, Tv *seasons, Tv &e, Tv &mu) {
  Ti rem = t % seasonCount;
  mu = level * seasons[rem];
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  Tv temp = level;
  level += alpha * e / seasons[rem];
  seasons[rem] += gamma * e / temp;
}

void rec_AAM(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
             Tv beta, Tv gamma, Tv &level, Tv &trend, Tv *seasons, Tv &e,
             Tv &mu) {
  Ti rem = t % seasonCount;
  mu = (level + trend) * seasons[rem];
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  Tv temp = level + trend;
  level += trend + alpha * e / seasons[rem];
  trend += beta * e / seasons[rem];
  seasons[rem] += gamma * e / temp;
}

void rec_AA_dM(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
               Tv beta, Tv phi, Tv gamma, Tv &level, Tv &trend, Tv *seasons,
               Tv &e, Tv &mu) {
  Ti rem = t % seasonCount;
  Tv p = phi * trend;
  Tv temp = level + p;
  mu = temp * seasons[rem];
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  level = temp + alpha * e / seasons[rem];
  trend = p + beta * e / seasons[rem];
  seasons[rem] += gamma * e / temp;
}

void rec_AMM(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
             Tv beta, Tv gamma, Tv &level, Tv &trend, Tv *seasons, Tv &e,
             Tv &mu) {
  Ti rem = t % seasonCount;
  Tv temp = level * trend;
  mu = temp * seasons[rem];
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update
  Tv temp0 = level;
  level = temp + alpha * e / seasons[rem];
  trend += beta * e / (seasons[rem] * temp0);
  seasons[rem] += gamma * e / temp;
}

void rec_AM_dM(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
               Tv beta, Tv phi, Tv gamma, Tv &level, Tv &trend, Tv *seasons,
               Tv &e, Tv &mu) {
  Ti rem = t % seasonCount;
  Tv p = std::pow(trend, phi);
  Tv temp = level * p;
  mu = temp * seasons[rem];
  if (sim)
    y = mu + e;
  else
    e = y - mu;
  // update

  Tv temp0 = level;
  level = temp + alpha * e / seasons[rem];
  trend = p + beta * e / (seasons[rem] * temp0);
  seasons[rem] += gamma * e / temp;
}

//#pragma endregion

//#pragma region recursives M

// table 2.2, hyndman, et ct (2008)
// see also table 2.3
// I replaced e with mu * e

void calc_err(bool sim, Tv &y, Tv &e, Tv &mu, Tv &e0) {
  if (sim) {
    e0 = mu * e;
    y = mu + e0;
  } else {
    e0 = y - mu;
    e = e0 / mu;

    // in a multiplicative error model, mu might be very large, i.e. infinity.
    // In these cases, e (and therefore sigma2) will be NAN.
  }
}

void rec_MNN(bool sim, Tv &y, Tv alpha, Tv &level, Tv &e, Tv &mu) {
  mu = level;

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  level += alpha * e0;
}

void rec_MAN(bool sim, Tv &y, Tv alpha, Tv beta, Tv &level, Tv &trend, Tv &e,
             Tv &mu) {
  mu = level + trend;

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  level += trend + alpha * e0;
  trend += beta * e0;
}

void rec_MA_dN(bool sim, Tv &y, Tv alpha, Tv beta, Tv phi, Tv &level, Tv &trend,
               Tv &e, Tv &mu) {
  Tv p = phi * trend;
  mu = level + phi * trend;

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  level = mu + alpha * e0;
  trend = p + beta * e0;
}

void rec_MMN(bool sim, Tv &y, Tv alpha, Tv beta, Tv &level, Tv &trend, Tv &e,
             Tv &mu) {
  mu = level * trend;

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  Tv temp = level;
  level = mu + alpha * e0;
  trend += beta * e0 / temp;
}

void rec_MM_dN(bool sim, Tv &y, Tv alpha, Tv beta, Tv phi, Tv &level, Tv &trend,
               Tv &e, Tv &mu) {
  Tv p = std::pow(trend, phi);
  mu = level * p;

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  Tv temp = level;
  level = mu + alpha * e0;
  trend = p + beta * e0 / temp;
}

// Next Column

void rec_MNA(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
             Tv gamma, Tv &level, Tv *seasons, Tv &e, Tv &mu) {
  Ti rem = t % seasonCount;
  mu = level + seasons[rem];

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  level += alpha * e0;
  seasons[rem] += gamma * e0;
}

void rec_MAA(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
             Tv beta, Tv gamma, Tv &level, Tv &trend, Tv *seasons, Tv &e,
             Tv &mu) {
  Ti rem = t % seasonCount;
  mu = level + trend + seasons[rem];

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  level += trend + alpha * e0;
  trend += beta * e0;
  seasons[rem] += gamma * e0;
}

void rec_MA_dA(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
               Tv beta, Tv phi, Tv gamma, Tv &level, Tv &trend, Tv *seasons,
               Tv &e, Tv &mu) {
  Ti rem = t % seasonCount;
  Tv p = phi * trend;
  mu = level + p + seasons[rem];

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  level += p + alpha * e0;
  trend = p + beta * e0;
  seasons[rem] += gamma * e0;
}

void rec_MMA(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
             Tv beta, Tv gamma, Tv &level, Tv &trend, Tv *seasons, Tv &e,
             Tv &mu) {
  Ti rem = t % seasonCount;
  mu = level * trend + seasons[rem];

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  Tv temp = level;
  level = level * trend + alpha * e0;
  trend += beta * e0 / temp;
  seasons[rem] += gamma * e0;
}

void rec_MM_dA(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
               Tv beta, Tv phi, Tv gamma, Tv &level, Tv &trend, Tv *seasons,
               Tv &e, Tv &mu) {
  Ti rem = t % seasonCount;
  Tv p = std::pow(trend, phi);
  mu = level * p + seasons[rem];

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  Tv temp = level;
  level = level * p + alpha * e0;
  trend = p + beta * e0 / temp;
  seasons[rem] += gamma * e0;
}

// THIRD COLUMN

void rec_MNM(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
             Tv gamma, Tv &level, Tv *seasons, Tv &e, Tv &mu) {
  Ti rem = t % seasonCount;
  mu = level * seasons[rem];

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  Tv temp = level;
  level += alpha * e0 / seasons[rem];
  seasons[rem] += gamma * e0 / temp;
}

void rec_MAM(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
             Tv beta, Tv gamma, Tv &level, Tv &trend, Tv *seasons, Tv &e,
             Tv &mu) {
  Ti rem = t % seasonCount;
  mu = (level + trend) * seasons[rem];

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  Tv temp = level + trend;
  level += trend + alpha * e0 / seasons[rem];
  trend += beta * e0 / seasons[rem];
  seasons[rem] += gamma * e0 / temp;
}

void rec_MA_dM(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
               Tv beta, Tv phi, Tv gamma, Tv &level, Tv &trend, Tv *seasons,
               Tv &e, Tv &mu) {
  Ti rem = t % seasonCount;
  Tv p = phi * trend;
  Tv temp = level + p;
  mu = temp * seasons[rem];

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  level = temp + alpha * e0 / seasons[rem];
  trend = p + beta * e0 / seasons[rem];
  seasons[rem] += gamma * e0 / temp;
}

void rec_MMM(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
             Tv beta, Tv gamma, Tv &level, Tv &trend, Tv *seasons, Tv &e,
             Tv &mu) {
  Ti rem = t % seasonCount;
  Tv temp = level * trend;
  mu = temp * seasons[rem];

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update
  Tv temp0 = level;
  level = temp + alpha * e0 / seasons[rem];
  trend += beta * e0 / (seasons[rem] * temp0);
  seasons[rem] += gamma * e0 / temp;
}

void rec_MM_dM(bool sim, Ti t, unsigned short seasonCount, Tv &y, Tv alpha,
               Tv beta, Tv phi, Tv gamma, Tv &level, Tv &trend, Tv *seasons,
               Tv &e, Tv &mu) {
  Ti rem = t % seasonCount;
  Tv p = std::pow(trend, phi);
  Tv temp = level * p;
  mu = temp * seasons[rem];

  Tv e0;
  calc_err(sim, y, e, mu, e0);
  // update

  Tv temp0 = level;
  level = temp + alpha * e0 / seasons[rem];
  trend = p + beta * e0 / (seasons[rem] * temp0);
  seasons[rem] += gamma * e0 / temp;
}

//#pragma endregion

void initializestates_size(Ti countUse, unsigned short seasonCount,
                           Ti &workSize) {

  if (seasonCount == 0) {
    workSize = 0;
  } else {
    workSize = seasonCount + 2 * countUse;
  }
}

void initializestates(Matrix<Tv> *data, ExpSmoothingResultTv *storage, Tv *WORK,
                      bool3 additiveTrend, bool3 additiveSeason,
                      unsigned short seasonCount, Ti countUse) {
  if (countUse >= data->length())
    throw std::logic_error("initialization size is larger than the array size");

  // parameters
  if (std::isnan(storage->alpha))
    storage->alpha = 0.2;
  if (std::isnan(storage->beta))
    storage->beta = 0.1;
  if (std::isnan(storage->gamma))
    storage->gamma = 0.05;
  if (std::isnan(storage->phi))
    storage->phi = 0.99;

  // initial states
  Ti T = data->length();
  data->Restructure0(countUse, 1);
  auto stat = Descriptive(data);

  if (additiveSeason == Null) {
    if ((std::isnan(storage->init_Level)) ||
        (additiveTrend != Null && std::isnan(storage->init_Trend))) {
      Tv trend[2];
      stat.RegressionTrend(trend);

      if (additiveTrend == Null && std::isnan(storage->init_Level))
        storage->init_Level = trend[0];
      else if (additiveTrend == True) {
        if (std::isnan(storage->init_Level))
          storage->init_Level = trend[0];
        if (std::isnan(storage->init_Trend))
          storage->init_Trend = trend[1];
      } else {
        if (std::isnan(storage->init_Level))
          storage->init_Level = trend[0];
        if (std::isnan(storage->init_Trend))
          storage->init_Trend = 1 + trend[1] / trend[0];
      }
    }
  } else {
    Ti q = 0;
    auto s_means = Matrix<Tv>(&WORK[q], seasonCount, 1);
    q += seasonCount;
    auto s_seasonal = Matrix<Tv>(&WORK[q], countUse, 1);
    q += countUse;
    auto s_trend = Matrix<Tv>(&WORK[q], countUse, 1);
    q += countUse;
    stat.SeasonalDecompositionMa(&s_trend, &s_seasonal, {}, &s_means, {},
                                 seasonCount, additiveSeason == False, false,
                                 false);
    s_trend.RemoveNanVector_in(true);
    auto vstat = Descriptive(&s_trend);
    Tv trend[2];
    vstat.RegressionTrend(trend);

    if ((std::isnan(storage->init_Level)) ||
        (additiveTrend != Null && std::isnan(storage->init_Trend))) {
      if (additiveTrend == Null && std::isnan(storage->init_Level))
        storage->init_Level = trend[0];
      else if (additiveTrend == True) {
        if (std::isnan(storage->init_Level))
          storage->init_Level = trend[0];
        if (std::isnan(storage->init_Trend))
          storage->init_Trend = trend[1];
      } else {
        if (std::isnan(storage->init_Level))
          storage->init_Level = trend[0];
        if (std::isnan(storage->init_Trend))
          storage->init_Trend = 1 + trend[1] / trend[0];
      }
    }
    // set seasons
    Ti j;
    for (j = 0; j < seasonCount - 1; j++)
      if (std::isnan(storage->init_seasons->Data[j]))
        storage->init_seasons->Data[j] = s_means.Data[j];
    setLastSeason(additiveSeason == True, storage->init_seasons->Data,
                  seasonCount); // restrict last season
  }

  data->Restructure0(T, 1);
}

void adjustInitialLength(Ti &initialUpdateLength, unsigned short seasonCount,
                         Ti T) {
  if (initialUpdateLength <= 0) {
    if (seasonCount == 0)
      initialUpdateLength = std::min((Ti)10, T);
    else
      initialUpdateLength = std::min((Ti)(6 * seasonCount), T);
  }
  if (initialUpdateLength >= T)
    initialUpdateLength = T - 1;
}

template <ExpSmoothingType type>
void expsmoothing<type>::estimate_size(Matrix<Tv> *data,
                                       LimitedMemoryBFGSBTv *optim,
                                       ldt::DerivativeTv *derv, Ti &workSize,
                                       Ti &storageSize, bool keepAll,
                                       Ti initialUpdateLength) {

  Ti T = data->length();
  Ti x0Length = _m + _k;
  workSize =
      4 * x0Length + _seasonCount; // x0, gradient, lower and upper bounds, temp
                                   // for seasons in optim function

  adjustInitialLength(initialUpdateLength, _seasonCount, T);

  Ti wS = 0;
  initializestates_size(initialUpdateLength, _seasonCount, wS);
  workSize += wS;

  workSize += derv->first_nWORK(x0Length);
  workSize += optim->worksizes(x0Length);

  // storage
  storageSize = 2 * _seasonCount;
  if (keepAll)
    storageSize += (3 + _seasonCount) * T;
}

void ExpSmoothingResult::setStorage(Tv *S, unsigned short seasoncount,
                                    Ti length, bool keepAll) {
  if (seasoncount < 2)
    seasoncount = 0;
  Ti q = 0;
  init_seasons = new Matrix<Tv>(NAN, &S[q], seasoncount, 1);
  q += seasoncount;
  last_seasons = new Matrix<Tv>(NAN, &S[q], seasoncount, 1);
  q += seasoncount;

  if (keepAll) {
    resid = new Matrix<Tv>(&S[q], length, 1);
    q += length;
    levels = new Matrix<Tv>(&S[q], length, 1);
    q += length;
    trends = new Matrix<Tv>(&S[q], length, 1);
    q += length;
    seasons = new std::vector<Matrix<Tv> *>();
    for (int i = 0; i < seasoncount; i++) {
      seasons->push_back(new Matrix<Tv>(&S[q], length, 1));
      q += length;
    }
  }
}

ExpSmoothingResult::~ExpSmoothingResult() {
  if (init_seasons)
    delete init_seasons;
  if (last_seasons)
    delete last_seasons;
  if (resid) {
    delete resid;
    delete levels;
    delete trends;
    for (std::vector<int>::size_type i = 0; i < seasonsCount; i++)
      delete seasons->at(i);
    delete seasons;
  }
}

void keepStates(Ti t, ExpSmoothingResultTv *res, Tv level, Tv trend,
                Tv *seasons, Ti seasonsCount, Tv e) {
  res->levels->Data[t] = level;
  res->trends->Data[t] = trend;
  for (auto i = 0; i < seasonsCount; i++)
    res->seasons->at(i)->Data[t] = seasons[i];
  res->resid->Data[t] = e;
}

template <ExpSmoothingType type>
void expsmoothing<type>::estimate(Matrix<Tv> *data, ExpSmoothingResultTv *res,
                                  Tv *WORK, LimitedMemoryBFGSBTv *optim,
                                  ldt::DerivativeTv *derv,
                                  ExpSmoothingResultTv *initials,
                                  Ti &initialUpdateLength, bool keepAll) {

  if (data->length() <= 0)
    throw std::logic_error("invalid length");

  if (_additive_error == false)
    for (Ti i = 0; i < data->length(); i++)
      if (data->Data[i] <= 0)
        throw std::logic_error(
            "Negative data for a model with multiplicative error");

  if (_seasonCount != 0) {
    if (!initials->init_seasons)
      throw std::logic_error("seasons are not initialized");
    else if (initials->init_seasons->length() != _seasonCount)
      throw std::logic_error("Inconsistent number of seasons");
  }

  auto T = data->length();
  Ti i, j, r;
  res->seasonsCount = _seasonCount;

  // Update NAN initials
  adjustInitialLength(initialUpdateLength, _seasonCount, T);

  try {
    Ti wS = 0;
    initializestates_size(initialUpdateLength, _seasonCount, wS);
    auto W = WORK;
    initializestates(data, initials, W, _additive_trend, _additive_season,
                     _seasonCount, initialUpdateLength);
  } catch (...) {
    throw std::logic_error("Initialization failed");
  }

  Ti q = 0;

  // set initial vector
  // states: level, growth (if has trend), seasons (if is seasonal)
  // params: alpha, beta, gamma, phi
  auto x0Length = _m + _k;

  auto x0 = &WORK[q];
  q += x0Length;

  r = 0;
  if (_additive_error)
    x0[r++] = initials->init_Level;
  if (_additive_trend)
    x0[r++] = initials->init_Trend;
  if (_seasonCount != 0)
    for (i = 0; i < _seasonCount - 1; i++) // -1, the last season is restricted
      x0[r++] = initials->init_seasons->Data[i];

  // Initial Pars
  x0[r++] = initials->alpha;

  if (_additive_trend != Null) {
    x0[r++] = initials->beta;
    if (_damped)
      x0[r++] = initials->phi; // phi is before gamma
  }
  if (_additive_season != Null)
    x0[r++] = initials->gamma;

  // optim

  // restrictions
  auto lower = Matrix<Tv>(&WORK[q], x0Length, 1);
  q += x0Length;
  auto upper = Matrix<Tv>(&WORK[q], x0Length, 1);
  q += x0Length;
  for (i = 0; i < _m; i++) {
    lower.Data[i] = -INFINITY;
    upper.Data[i] = INFINITY;
  }
  for (i = _m; i < x0Length; i++) {
    lower.Data[i] = -0.0001;
    upper.Data[i] = 0.9999;
  }

  optim->setbounds(&lower, &upper);
  bool keep = false;

  Tv alpha = NAN;
  Tv beta = NAN;
  Tv gamma = NAN;
  Tv phi = NAN;
  Tv sigma2 = NAN;

  Tv level = NAN;
  Tv trend = NAN;
  auto seasons = &WORK[q];
  q += _seasonCount;

  Tv e = 0;
  Tv e2 = 0;
  Tv mu = 1;
  Tv lnabsr = 0;

  std::function<(const Matrix<Tv> *)> f =
      [&](const Matrix<Tv> *x) -> Tv {
    auto xd = x->Data;
    level = xd[0];
    trend = 0;
    i = 1;
    if (_additive_trend != Null)
      trend = xd[i++];
    if (_additive_season != Null) {
      for (j = 0; j < _seasonCount - 1; j++)
        seasons[j] = xd[i++];
      setLastSeason(_additive_season == True, seasons, _seasonCount);
    }

    alpha = xd[_m];
    i = _m + 1;
    beta = NAN;
    gamma = NAN;
    phi = NAN;

    if (_additive_trend != Null) {
      beta = xd[i++];
      if (_damped)
        phi = xd[i++];
    }
    if (_additive_season != Null)
      gamma = xd[i++];

    e = 0;
    e2 = 0;
    mu = 1;
    lnabsr = 0;

    for (Ti t = 0; t < T; t++) {

      if constexpr (type == ExpSmoothingType::aaa) {
        rec_AAA(false, t, _seasonCount, data->Data[t], alpha, beta, gamma,
                level, trend, seasons, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::aam) {
        rec_AAM(false, t, _seasonCount, data->Data[t], alpha, beta, gamma,
                level, trend, seasons, e, mu);
        e2 += e * e;
      } else if constexpr (type == ExpSmoothingType::aan) {
        rec_AAN(false, data->Data[t], alpha, beta, level, trend, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::aa_da) {
        rec_AA_dA(false, t, _seasonCount, data->Data[t], alpha, beta, phi,
                  gamma, level, trend, seasons, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::aa_dm) {
        rec_AA_dM(false, t, _seasonCount, data->Data[t], alpha, beta, phi,
                  gamma, level, trend, seasons, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::aa_dn) {
        rec_AA_dN(false, data->Data[t], alpha, beta, phi, level, trend, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::ama) {
        rec_AMA(false, t, _seasonCount, data->Data[t], alpha, beta, gamma,
                level, trend, seasons, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::amm) {
        rec_AMM(false, t, _seasonCount, data->Data[t], alpha, beta, gamma,
                level, trend, seasons, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::amn) {
        rec_AMN(false, data->Data[t], alpha, beta, level, trend, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::am_da) {
        rec_AM_dA(false, t, _seasonCount, data->Data[t], alpha, beta, phi,
                  gamma, level, trend, seasons, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::am_dm) {
        rec_AM_dM(false, t, _seasonCount, data->Data[t], alpha, beta, phi,
                  gamma, level, trend, seasons, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::am_dn) {
        rec_AM_dN(false, data->Data[t], alpha, beta, phi, level, trend, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::ana) {
        rec_ANA(false, t, _seasonCount, data->Data[t], alpha, gamma, level,
                seasons, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::anm) {
        rec_ANM(false, t, _seasonCount, data->Data[t], alpha, gamma, level,
                seasons, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::ann) {
        rec_ANN(false, data->Data[t], alpha, level, e, mu);
        e2 += e * e;
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::maa) {
        rec_MAA(false, t, _seasonCount, data->Data[t], alpha, beta, gamma,
                level, trend, seasons, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::mam) {
        rec_MAM(false, t, _seasonCount, data->Data[t], alpha, beta, gamma,
                level, trend, seasons, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::man) {
        rec_MAN(false, data->Data[t], alpha, beta, level, trend, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::ma_da) {
        rec_MA_dA(false, t, _seasonCount, data->Data[t], alpha, beta, phi,
                  gamma, level, trend, seasons, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::ma_dm) {
        rec_MA_dM(false, t, _seasonCount, data->Data[t], alpha, beta, phi,
                  gamma, level, trend, seasons, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::ma_dn) {
        rec_MA_dN(false, data->Data[t], alpha, beta, phi, level, trend, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::mma) {
        rec_MMA(false, t, _seasonCount, data->Data[t], alpha, beta, gamma,
                level, trend, seasons, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::mmm) {
        rec_MMM(false, t, _seasonCount, data->Data[t], alpha, beta, gamma,
                level, trend, seasons, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::mmn) {
        rec_MMN(false, data->Data[t], alpha, beta, level, trend, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::mm_da) {
        rec_MM_dA(false, t, _seasonCount, data->Data[t], alpha, beta, phi,
                  gamma, level, trend, seasons, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::mm_dm) {
        rec_MM_dM(false, t, _seasonCount, data->Data[t], alpha, beta, phi,
                  gamma, level, trend, seasons, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::mm_dn) {
        rec_MM_dN(false, data->Data[t], alpha, beta, phi, level, trend, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::mna) {
        rec_MNA(false, t, _seasonCount, data->Data[t], alpha, gamma, level,
                seasons, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::mnm) {
        rec_MNM(false, t, _seasonCount, data->Data[t], alpha, gamma, level,
                seasons, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (type == ExpSmoothingType::mnn) {
        rec_MNN(false, data->Data[t], alpha, level, e, mu);
        e2 += e * e;
        lnabsr += std::log(std::abs(mu));
        if (keep)
          keepStates(t, res, level, trend, seasons, _seasonCount, e);
      } else if constexpr (true) {
        throw std::logic_error("not implemented");
      }
    }

    sigma2 = e2 / (T - 1); // unbiased
    if (_additive_error == false)
      return T * std::log(e2) + 2 * lnabsr;
    else
      return T * std::log(e2);
  }

                                  // DerivativeTv

                                  derv->setfunction(&f);
  auto dervW = &WORK[q];
  q += derv->first_nWORK(x0Length);
  std::function<void(const Matrix<Tv> *, Matrix<Tv> *)> g =
      [derv, dervW](const Matrix<Tv> *x,
                    Matrix<Tv> *storage) -> void {
    derv->first(x, storage, dervW);
  }

                                            optim->setfunctions(&f, &g);
  auto optimW = &WORK[q];
  q += optim->worksizes(x0Length);
  auto gradient = Matrix<Tv>(&WORK[q], x0Length, 1);
  q += x0Length;
  auto mat_x0 = Matrix<Tv>(x0, x0Length, 1);

  optim->minimize0(&mat_x0, &gradient, optimW);

  // last run
  keep = keepAll;
  res->logL_con = f(&mat_x0) / -2;
  res->logL =
      (T * (c_ln_2Pi_plus_one - std::log(T))) / -2 + res->logL_con; // p. 69
  res->aic = -2 * res->logL_con + 2 * (x0Length);
  res->sic = -2 * res->logL_con + std::log(T) * (x0Length);
  auto ystat = Descriptive(data);
  res->r2 = sigma2 * (T - 1) / ystat.SumOfSquares(true);

  res->type = _type;

  res->sigma2 = sigma2;
  res->T = T;

  res->alpha = alpha;
  res->beta = beta;
  res->phi = phi;
  res->gamma = gamma;

  r = 0;
  res->init_Level = mat_x0.Data[r++];
  if (_additive_trend != Null)
    res->init_Trend = mat_x0.Data[r++];

  if (_additive_season != Null) {
    for (i = r; i < _m; i++)
      res->init_seasons->Data[i - r] = mat_x0.Data[i];
    setLastSeason(_additive_season == True, res->init_seasons->Data,
                  _seasonCount);
  }

  res->last_Level = level;
  res->last_Trend = trend;
  for (i = 0; i < _seasonCount; i++)
    res->last_seasons->Data[i] = seasons[i];
}

//#pragma region Simulate

template <ExpSmoothingType type>
void expsmoothing<type>::simulate(ExpSmoothingResultTv *pars, Tv *storage, Ti n,
                                  Ti burnin, std::function<(Ti)> shock) {

  if (shock == NULL) {
    std::random_device rdev{};
    auto eng = std::mt19937(rdev());
    std::normal_distribution<Tv> dist(0, std::sqrt(pars->sigma2));
    shock = [&dist, &eng](Ti) -> Tv { return dist(eng); };
  }

  // TODO: use constexpr

  Tv mu = 0;
  Tv y = 0;
  Tv e;
  for (Ti t = 0; t < n + burnin; t++) {
    e = shock(t);

    switch (pars->type) {
    case ldt::ExpSmoothingType::ann:
      rec_ANN(true, y, pars->alpha, pars->init_Level, e, mu);
      break;
    case ldt::ExpSmoothingType::aan:
      rec_AAN(true, y, pars->alpha, pars->beta, pars->init_Level,
              pars->init_Trend, e, mu);
      break;
    case ldt::ExpSmoothingType::amn:
      rec_AMN(true, y, pars->alpha, pars->beta, pars->init_Level,
              pars->init_Trend, e, mu);
      break;
    case ldt::ExpSmoothingType::ana:
      rec_ANA(true, t, pars->seasonsCount, y, pars->alpha, pars->gamma,
              pars->init_Level, pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::aaa:
      rec_AAA(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
              pars->gamma, pars->init_Level, pars->init_Trend,
              pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::ama:
      rec_AMA(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
              pars->gamma, pars->init_Level, pars->init_Trend,
              pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::anm:
      rec_ANM(true, t, pars->seasonsCount, y, pars->alpha, pars->gamma,
              pars->init_Level, pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::aam:
      rec_AAM(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
              pars->gamma, pars->init_Level, pars->init_Trend,
              pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::amm:
      rec_AMM(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
              pars->gamma, pars->init_Level, pars->init_Trend,
              pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::aa_dn:
      rec_AA_dN(true, y, pars->alpha, pars->beta, pars->phi, pars->init_Level,
                pars->init_Trend, e, mu);
      break;
    case ldt::ExpSmoothingType::am_dn:
      rec_AM_dN(true, y, pars->alpha, pars->beta, pars->phi, pars->init_Level,
                pars->init_Trend, e, mu);
      break;
    case ldt::ExpSmoothingType::aa_da:
      rec_AA_dA(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
                pars->phi, pars->gamma, pars->init_Level, pars->init_Trend,
                pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::am_da:
      rec_AM_dA(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
                pars->phi, pars->gamma, pars->init_Level, pars->init_Trend,
                pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::aa_dm:
      rec_AA_dM(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
                pars->phi, pars->gamma, pars->init_Level, pars->init_Trend,
                pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::am_dm:
      rec_AM_dM(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
                pars->phi, pars->gamma, pars->init_Level, pars->init_Trend,
                pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::mnn:
      rec_MNN(true, y, pars->alpha, pars->init_Level, e, mu);
      break;
    case ldt::ExpSmoothingType::man:
      rec_MAN(true, y, pars->alpha, pars->beta, pars->init_Level,
              pars->init_Trend, e, mu);
      break;
    case ldt::ExpSmoothingType::mmn:
      rec_MMN(true, y, pars->alpha, pars->beta, pars->init_Level,
              pars->init_Trend, e, mu);
      break;
    case ldt::ExpSmoothingType::mna:
      rec_MNA(true, t, pars->seasonsCount, y, pars->alpha, pars->gamma,
              pars->init_Level, pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::maa:
      rec_MAA(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
              pars->gamma, pars->init_Level, pars->init_Trend,
              pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::mma:
      rec_MMA(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
              pars->gamma, pars->init_Level, pars->init_Trend,
              pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::mnm:
      rec_MNM(true, t, pars->seasonsCount, y, pars->alpha, pars->gamma,
              pars->init_Level, pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::mam:
      rec_MAM(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
              pars->gamma, pars->init_Level, pars->init_Trend,
              pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::mmm:
      rec_MMM(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
              pars->gamma, pars->init_Level, pars->init_Trend,
              pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::ma_dn:
      rec_MA_dN(true, y, pars->alpha, pars->beta, pars->phi, pars->init_Level,
                pars->init_Trend, e, mu);
      break;
    case ldt::ExpSmoothingType::mm_dn:
      rec_MM_dN(true, y, pars->alpha, pars->beta, pars->phi, pars->init_Level,
                pars->init_Trend, e, mu);
      break;
    case ldt::ExpSmoothingType::ma_da:
      rec_MA_dA(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
                pars->phi, pars->gamma, pars->init_Level, pars->init_Trend,
                pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::mm_da:
      rec_MM_dA(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
                pars->phi, pars->gamma, pars->init_Level, pars->init_Trend,
                pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::ma_dm:
      rec_MA_dM(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
                pars->phi, pars->gamma, pars->init_Level, pars->init_Trend,
                pars->init_seasons->Data, e, mu);
      break;
    case ldt::ExpSmoothingType::mm_dm:
      rec_MM_dM(true, t, pars->seasonsCount, y, pars->alpha, pars->beta,
                pars->phi, pars->gamma, pars->init_Level, pars->init_Trend,
                pars->init_seasons->Data, e, mu);
      break;
    default:
      break;
    }

    if (t >= burnin)
      storage[t - burnin] = y;
  }
}

//#pragma endregion

//#pragma region Forecast

template <ExpSmoothingType type>
void forecastSim(Ti h, ExpSmoothingResultTv *estim, Matrix<Tv> *forecast,
                 Matrix<Tv> *varianceF, Tv *WORK, Ti count, std::mt19937 &eng) {
  Ti t, counter;
  Tv y, mu, level, e, trend;

  auto W_shocks = WORK;    // h
  auto seasons = &WORK[h]; // seasonsCount

  std::normal_distribution<Tv> dist(0, std::sqrt(estim->sigma2));
  auto mvs = std::vector<RunningWeightedVariance *>();
  for (t = 0; t < h; t++)
    mvs.push_back(new RunningWeightedVariance());

  counter = 0;
  while (counter < count) {
    counter++;

    for (t = 0; t < h; t++)
      W_shocks[t] = dist(eng);

    y = 0;
    mu = 0;
    level = estim->last_Level;
    trend = estim->last_Trend;
    if (estim->seasonsCount != 0)
      for (t = 0; t < estim->seasonsCount; t++)
        seasons[t] = estim->last_seasons->Data[t];

    for (t = 0; t < h; t++) {
      if constexpr (type == ExpSmoothingType::aaa) {
        e = W_shocks[t];
        rec_AAA(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                estim->beta, estim->gamma, level, trend, seasons, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::aam) {
        e = W_shocks[t];
        rec_AAM(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                estim->beta, estim->gamma, level, trend, seasons, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::aan) {
        e = W_shocks[t];
        rec_AAN(true, y, estim->alpha, estim->beta, level, trend, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::aa_da) {
        e = W_shocks[t];
        rec_AA_dA(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                  estim->beta, estim->phi, estim->gamma, level, trend, seasons,
                  e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::aa_dm) {
        e = W_shocks[t];
        rec_AA_dM(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                  estim->beta, estim->phi, estim->gamma, level, trend, seasons,
                  e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::aa_dn) {
        e = W_shocks[t];
        rec_AA_dN(true, y, estim->alpha, estim->beta, estim->phi, level, trend,
                  e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::ama) {
        e = W_shocks[t];
        rec_AMA(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                estim->beta, estim->gamma, level, trend, seasons, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::amm) {
        e = W_shocks[t];
        rec_AMM(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                estim->beta, estim->gamma, level, trend, seasons, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::amn) {
        e = W_shocks[t];
        rec_AMN(true, y, estim->alpha, estim->beta, level, trend, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::am_da) {
        e = W_shocks[t];
        rec_AM_dA(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                  estim->beta, estim->phi, estim->gamma, level, trend, seasons,
                  e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::am_dm) {
        e = W_shocks[t];
        rec_AM_dM(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                  estim->beta, estim->phi, estim->gamma, level, trend, seasons,
                  e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::am_dn) {
        e = W_shocks[t];
        rec_AM_dN(true, y, estim->alpha, estim->beta, estim->phi, level, trend,
                  e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::ana) {
        e = W_shocks[t];
        rec_ANA(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                estim->gamma, level, seasons, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::anm) {
        e = W_shocks[t];
        rec_ANM(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                estim->gamma, level, seasons, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::ann) {
        e = W_shocks[t];
        rec_ANN(true, y, estim->alpha, level, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::maa) {
        e = W_shocks[t];
        rec_MAA(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                estim->beta, estim->gamma, level, trend, seasons, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::mam) {
        e = W_shocks[t];
        rec_MAM(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                estim->beta, estim->gamma, level, trend, seasons, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::man) {
        e = W_shocks[t];
        rec_MAN(true, y, estim->alpha, estim->beta, level, trend, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::ma_da) {
        e = W_shocks[t];
        rec_MA_dA(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                  estim->beta, estim->phi, estim->gamma, level, trend, seasons,
                  e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::ma_dm) {
        e = W_shocks[t];
        rec_MA_dM(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                  estim->beta, estim->phi, estim->gamma, level, trend, seasons,
                  e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::ma_dn) {
        e = W_shocks[t];
        rec_MA_dN(true, y, estim->alpha, estim->beta, estim->phi, level, trend,
                  e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::mma) {
        e = W_shocks[t];
        rec_MMA(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                estim->beta, estim->gamma, level, trend, seasons, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::mmm) {
        e = W_shocks[t];
        rec_MMM(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                estim->beta, estim->gamma, level, trend, seasons, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::mmn) {
        e = W_shocks[t];
        rec_MMN(true, y, estim->alpha, estim->beta, level, trend, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::mm_da) {
        e = W_shocks[t];
        rec_MM_dA(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                  estim->beta, estim->phi, estim->gamma, level, trend, seasons,
                  e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::mm_dm) {
        e = W_shocks[t];
        rec_MM_dM(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                  estim->beta, estim->phi, estim->gamma, level, trend, seasons,
                  e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::mm_dn) {
        e = W_shocks[t];
        rec_MM_dN(true, y, estim->alpha, estim->beta, estim->phi, level, trend,
                  e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::mna) {
        e = W_shocks[t];
        rec_MNA(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                estim->gamma, level, seasons, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::mnm) {
        e = W_shocks[t];
        rec_MNM(true, estim->T + t, estim->seasonsCount, y, estim->alpha,
                estim->gamma, level, seasons, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (type == ExpSmoothingType::mnn) {
        e = W_shocks[t];
        rec_MNN(true, y, estim->alpha, level, e, mu);
        mvs[t]->PushNew(y, 1);
      } else if constexpr (true) {
        throw std::logic_error("not implemented");
      }
    }
  }

  for (t = 0; t < h; t++) {
    forecast->Data[t] = mvs[t]->GetMean();
    varianceF->Data[t] = mvs[t]->GetVariancePopulation();

    delete mvs[t];
  }
}

Tv phj(Ti j, bool _damped, Tv phi) {
  if (_damped == false)
    return static_cast<Tv>(j);
  Tv s = 0;
  for (Ti jj = 1; jj <= j; jj++)
    s += std::pow(phi, jj);
  return s;
}

template <ExpSmoothingType type>
void forecastC12(expsmoothing<type> *exp, ExpSmoothingResultTv *estim, Ti T,
                 Ti h, bool c1, Tv *thetas, Tv &mean, Tv &var,
                 std::function<Tv(Ti, Ti)> meanj,
                 std::function<(Ti, bool)> cj) {
  Ti j;
  Ti rem = 0;
  if (exp->_additive_season != Null)
    rem = (T + h - 1) % exp->_seasonCount;

  mean = meanj(h, rem);
  var = NAN;

  if (c1) {
    // eq. 6.1, p. 81
    if (h > 1) {
      var = 0;
      if (exp->_seasonCount == 0)
        for (j = 1; j < h; j++)
          var += std::pow(cj(j, false), 2);
      else {
        Ti p = 0;
        for (j = 1; j < h; j++) {
          p++;
          if (p == exp->_seasonCount)
            p = 0;
          var += std::pow(cj(j, p == 0), 2);
        }
      }
      var = (1 + var) * estim->sigma2;
    } else
      var = estim->sigma2;
  } else {
    // eq. 6.2, p. 83
    auto meanmean = mean * mean;

    if (exp->_seasonCount == 0) {
      auto meanh =
          meanj(1, 0); // no season, therefore second argument is not used
      thetas[0] = meanh * meanh;
      for (Ti H = 1; H < h; H++) {
        var = 0;
        for (j = 1; j <= H; j++)
          var += std::pow(cj(j, false), 2) * thetas[H - j];
        meanh = meanj(H + 1, 0);
        thetas[H] = meanh * meanh + estim->sigma2 * var;
      }
    } else {
      rem = T % exp->_seasonCount;
      auto meanh = meanj(1, rem);
      thetas[0] = meanh * meanh;
      Ti p;
      for (Ti H = 1; H < h; H++) {
        p = 0;
        var = 0;
        for (j = 1; j <= H; j++) {
          p++;
          if (p == exp->_seasonCount)
            p = 0;
          var += std::pow(cj(j, p == 0), 2) * thetas[H - j];
        }
        rem++;
        if (rem == exp->_seasonCount)
          rem = 0;
        meanh = meanj(H + 1, rem);
        thetas[H] = meanh * meanh + estim->sigma2 * var;
      }
    }

    var = (1 + estim->sigma2) * thetas[h - 1] - meanmean;
  }
}

template <ExpSmoothingType type>
void expsmoothing<type>::forecast(Ti horizon, ExpSmoothingResultTv *estim,
                                  Tv *WORK, bool &simulatedVar,
                                  Matrix<Tv> *forecast, Matrix<Tv> *varianceF,
                                  bool forceSimV, Ti simFixSize,
                                  unsigned int seed) {
  // WORK size is horizon + seasonsCount

  Ti h, q;
  Tv mean, variance;

  // forecast, varianceF horizon x 1
  simulatedVar = forceSimV;
  q = 0;
  if (forceSimV || _class > 2) {
    if (std::isnan(estim->sigma2))
      throw std::logic_error(
          "In forecasting with the exponential model, variance of residual "
          "is not a number. One possibility is that 'mu' becomes very large "
          "(infinity( while updating a multiplicative error model");

    simulatedVar = true;

    std::mt19937 eng;
    if (seed != 0)
      eng = std::mt19937(seed);
    else {
      std::random_device rdev{};
      eng = std::mt19937(rdev());
    }

    forecastSim<type>(horizon, estim, forecast, varianceF, WORK, simFixSize,
                      eng);
    return;
  }

  simulatedVar = false;

  std::function<(Ti, bool)> cj = NULL;
  std::function<(Ti, Ti)> meanj = NULL;

  switch (_type) {
  case ExpSmoothingType::ann:
  case ExpSmoothingType::mnn:
    meanj = [estim](Ti, Ti) -> Tv { return estim->last_Level; };
    cj = [estim](Ti, bool) -> Tv { return estim->alpha; };
    break;
  case ExpSmoothingType::aan:
  case ExpSmoothingType::man:
  case ExpSmoothingType::aa_dn:
  case ExpSmoothingType::ma_dn:
    meanj = [estim, this](Ti j, Ti) -> Tv {
      return estim->last_Level +
             phj(j, _damped, estim->phi) * estim->last_Trend;
    };
    cj = [estim, this](Ti j, bool) -> Tv {
      return estim->alpha + phj(j, _damped, estim->phi) * estim->beta;
    };
    break;
  case ExpSmoothingType::ana:
  case ExpSmoothingType::mna:
    meanj = [estim](Ti, Ti r) -> Tv {
      return estim->last_Level + estim->last_seasons->Data[r];
    };
    cj = [estim](Ti, bool g) -> Tv {
      return estim->alpha + (g ? estim->gamma : 0);
    };
    break;
  case ExpSmoothingType::aaa:
  case ExpSmoothingType::aa_da:
  case ExpSmoothingType::maa:
  case ExpSmoothingType::ma_da:
    meanj = [estim, this](Ti j, Ti r) -> Tv {
      return estim->last_Level +
             phj(j, _damped, estim->phi) * estim->last_Trend +
             estim->last_seasons->Data[r];
    };
    cj = [estim, this](Ti j, bool g) -> Tv {
      return estim->alpha + phj(j, _damped, estim->phi) * estim->beta +
             (g ? estim->gamma : 0);
    };
    break;
  };

  for (h = 1; h <= horizon; h++) {

    forecastC12(this, estim, estim->T, h, _class == 1, WORK, mean, variance,
                meanj, cj);
    forecast->Data[h - 1] = mean;
    varianceF->Data[h - 1] = variance;
  }
}

void expsmoothing_base::interval(Tv &lower, Tv &upper, Tv forecast, Tv sd,
                                 Tv confidence) {
  Tv prob = (1 - confidence) / 2;
  auto norm = Distribution<DistributionType::kNormal>(forecast, sd);
  lower = norm.GetQuantile(prob);
  upper = norm.GetQuantile(1 - prob);
}

//#pragma endregion

//#pragma region Out - of - Sample

Ti expsmoothing_outsample_res::getStorageSize() {

  Ti s = measures->length();
  if (hasDetails)
    for (auto &m : *details)
      s += m->act_for_sd->length() + m->measures->length();
  return s;
}

void expsmoothing_outsample_res::setStorage(Tv *storage) {
  Ti q = 0;
  measures->SetData(0, &storage[q]);
  q += measures->length();
  if (hasDetails) {
    for (auto &m : *details) {
      m->act_for_sd->Data = &storage[q];
      q += m->act_for_sd->length();
      m->measures->SetData(0, &storage[q]);
      q += m->measures->length();
    }
  }
}

expsmoothing_outsample_res::~expsmoothing_outsample_res() {
  delete measures;
  if (hasDetails) {
    for (Ti i = 0; i < (Ti)details->size(); i++) {
      auto m = details->at(i);
      delete m->act_for_sd;
      delete m->measures;
      delete m;
    }
    delete details;
  }
}

template <ExpSmoothingType type>
void expsmoothing<type>::calculatemeasure(
    bool forsize, std::vector<Ti> *measures, std::vector<Ti> *horizons,
    Ti outOfSampleCount, Tv outOfSamplePercentage, bool keepDetails,
    expsmoothing_outsample_res *result, Matrix<Tv> *data,
    ExpSmoothingResultTv *initials, Ti initializeLength, Tv *WORK,
    LimitedMemoryBFGSBTv *optim, DerivativeTv *derv, bool usePreviousEstim,
    bool &simulateForecast, bool forceSimFor, Ti simFixSize,
    unsigned int seed) {
  result->hasDetails = keepDetails;
  auto T = data->length();
  Ti maxSampleEnd = outOfSampleCount;
  if (maxSampleEnd == 0)
    maxSampleEnd = static_cast<Ti>(std::floor(T * outOfSamplePercentage / 100));
  if (T - maxSampleEnd <= 0 || maxSampleEnd == 0)
    throw std::logic_error(
        "invalid number of simulations.It is larger than the available data, "
        "zero or negative");

  Ti hh = horizons->size();
  Ti hMin = horizons->at(0);
  Ti hMax = horizons->at(hh - 1);

  // size for estimation and forecast and keeping the details
  Ti estim_w = 0;
  Ti estim_s = 0;
  estimate_size(data, optim, derv, estim_w, estim_s, false,
                initializeLength); // no need to keep the states
  Ti fore_w = hMax + _seasonCount;
  Ti fore_s = 2 * hMax;

  result->WORKsize = estim_w + estim_s + fore_w + fore_s + 5 * hh;

  Ti mm = measures->size();
  Ti row = 0;
  if (forsize) {
    result->measures = new Matrix<Tv>(
        hh, mm); // row i is for horizons[i] in different measures

    if (keepDetails) {
      result->details =
          new std::vector<expsmoothing_outsample_res_detailsTv *>(maxSampleEnd);
      for (Ti se = maxSampleEnd; se > 0; se--) {
        auto re = new expsmoothing_outsample_res_details();
        re->sampleEnd = se;
        re->isValid = false; // default is false
        re->act_for_sd = new Matrix<Tv>(hh, 3);
        re->measures = new Matrix<Tv>(hh, mm);
        result->details->at(row) = re;
        row++;
      }
    }

    return;
  }

  std::mt19937 eng;
  if (seed != 0)
    eng = std::mt19937(seed);
  else {
    std::random_device rdev{};
    eng = std::mt19937(rdev());
  }
  std::uniform_int_distribution<unsigned int> seeder(10, 1000000);

  auto estim = ExpSmoothingResult();
  estim.seasonsCount = _seasonCount;
  estim.T = T;
  Ti pos = 0;
  auto WorkEstim = &WORK[pos];
  pos += estim_w; // to be used in the loop
  estim.setStorage(&WORK[pos], _seasonCount, T, false);
  pos += estim_s;
  auto forecast = Matrix<Tv>(&WORK[pos], hMax, 1);
  pos += hMax;
  auto variancef = Matrix<Tv>(&WORK[pos], hMax, 1);
  pos += hMax;
  auto WORKfore = &WORK[pos];
  pos += fore_w; // to be used in the loop

  auto act = Matrix<Tv>(&WORK[pos], hh, 1);
  pos += hh;
  auto forc = Matrix<Tv>(&WORK[pos], hh, 1);
  pos += hh;
  auto std = Matrix<Tv>(&WORK[pos], hh, 1);
  pos += hh;
  auto err = Matrix<Tv>(&WORK[pos], hh, 1);
  pos += hh;
  auto temp = Matrix<Tv>(&WORK[pos], hh, 1);
  pos += hh;

  Tv last;
  Ti k, colc, actIndex, effectiveH;
  bool useCurrentEstims;
  expsmoothing_outsample_res_detailsTv *details = nullptr;
  bool success = false;
  row = 0;
  auto counters = std::vector<Ti>(horizons->size());
  for (Ti se = maxSampleEnd; se > 0; se--) {
    useCurrentEstims = usePreviousEstim && success;

    details = nullptr;
    if (keepDetails)
      details = result->details->at(row);

    actIndex = T - se - 1;
    if (actIndex < 0)
      break;

    effectiveH = std::min(se, hMax);
    if (effectiveH < hMin)
      break;

    colc = 0;
    for (auto &h : *horizons) {
      if (h <= effectiveH)
        colc++;
    }
    if (colc <= 0)
      break; // no more forecasts

    // estimate and forecast
    data->Restructure0(T - se, 1);
    if (useCurrentEstims) {
      initials->alpha = estim.alpha;
      initials->beta = estim.beta;
      initials->phi = estim.phi;
      initials->gamma = estim.gamma;
      initials->init_Level = estim.init_Level;
      initials->init_Trend = estim.init_Trend;
      if (_seasonCount > 1)
        initials->init_seasons->CopyFrom(estim.init_seasons);
    }
    // todo: consider adding an option to estimate, in order to skip
    // initialization whenever it is not needed

    try {
      estimate(data, &estim, WorkEstim, optim, derv, initials, initializeLength,
               false);

      unsigned int newseed = seeder(eng); // TODO: use the same engine?!
      expsmoothing<type>::forecast(effectiveH, &estim, WORKfore,
                                   simulateForecast, &forecast, &variancef,
                                   forceSimFor, simFixSize, newseed);
    } catch (...) {
      row++;
      continue;
    }

    last = data->Data[actIndex];

    k = 0;
    for (auto &h : *horizons) {
      if (h <= effectiveH) {
        counters[k]++;
        act.Data[k] = data->Data[actIndex + h];
        forc.Data[k] = forecast.Data[h - 1];
        std.Data[k] = variancef.Data[h - 1];
      } else {
        act.Data[k] = NAN;
        forc.Data[k] = NAN;
        std.Data[k] = NAN;
      }
      k++;
    }
    std.Apply_in([](Tv x) -> Tv { return std::sqrt(x); });
    act.Subtract0(&forc, &err);

    if (keepDetails) {
      details->sampleEnd = se;
      details->isValid = true;
      details->act_for_sd->SetColumn(0, &act);
      details->act_for_sd->SetColumn(1, &forc);
      details->act_for_sd->SetColumn(2, &std);
    }
    success = true;
    result->validCount++;

    Ti c = 0;
    for (auto &eval : *measures) {
      Scoring::calculate_forecast_measure(false, eval, &err, &temp, &act, &forc,
                                          &std, last, {});

      // summation for calculating mean
      k = 0;
      for (auto &h : *horizons) {
        if (h > effectiveH)
          continue;
        result->measures->Set_Plus0(k, c, temp.Data[k]);
        k++;
      }

      if (keepDetails) {
        auto stdd = temp.ToString();
        details->measures->SetColumn(c, &temp);
      }

      c++;
    }
    row++;
  }

  // average:
  /// see also <see cref="EvaluationMeasureHelper.Aggregate(EvaluationMeasureID,
  /// Tv[])"/>
  for (Ti j = 0; j < hh; j++) {
    Tv o = (Tv)1 / counters[j]; // denominator for horizon j
    for (Ti i = 0; i < (Ti)measures->size(); i++) {
      if (measures->at(i) == 0 || measures->at(i) == 9) // RMSE or Scaled RMSE
        result->measures->Set0(j, i,
                               std::sqrt(result->measures->Get0(j, i) * o));
      else
        result->measures->Set0(j, i, result->measures->Get0(j, i) * o);
    }
  }

  data->Restructure0(T, 1);
}

//#pragma endregion

template class ldt::expsmoothing<ExpSmoothingType::aaa>;
template class ldt::expsmoothing<ExpSmoothingType::aam>;
template class ldt::expsmoothing<ExpSmoothingType::aan>;
template class ldt::expsmoothing<ExpSmoothingType::aa_da>;
template class ldt::expsmoothing<ExpSmoothingType::aa_dm>;
template class ldt::expsmoothing<ExpSmoothingType::aa_dn>;
template class ldt::expsmoothing<ExpSmoothingType::ama>;
template class ldt::expsmoothing<ExpSmoothingType::amm>;
template class ldt::expsmoothing<ExpSmoothingType::amn>;
template class ldt::expsmoothing<ExpSmoothingType::am_da>;
template class ldt::expsmoothing<ExpSmoothingType::am_dm>;
template class ldt::expsmoothing<ExpSmoothingType::am_dn>;
template class ldt::expsmoothing<ExpSmoothingType::ana>;
template class ldt::expsmoothing<ExpSmoothingType::anm>;
template class ldt::expsmoothing<ExpSmoothingType::ann>;
template class ldt::expsmoothing<ExpSmoothingType::maa>;
template class ldt::expsmoothing<ExpSmoothingType::mam>;
template class ldt::expsmoothing<ExpSmoothingType::man>;
template class ldt::expsmoothing<ExpSmoothingType::ma_da>;
template class ldt::expsmoothing<ExpSmoothingType::ma_dm>;
template class ldt::expsmoothing<ExpSmoothingType::ma_dn>;
template class ldt::expsmoothing<ExpSmoothingType::mma>;
template class ldt::expsmoothing<ExpSmoothingType::mmm>;
template class ldt::expsmoothing<ExpSmoothingType::mmn>;
template class ldt::expsmoothing<ExpSmoothingType::mm_da>;
template class ldt::expsmoothing<ExpSmoothingType::mm_dm>;
template class ldt::expsmoothing<ExpSmoothingType::mm_dn>;
template class ldt::expsmoothing<ExpSmoothingType::mna>;
template class ldt::expsmoothing<ExpSmoothingType::mnm>;
template class ldt::expsmoothing<ExpSmoothingType::mnn>;
