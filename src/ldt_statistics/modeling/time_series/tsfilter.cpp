/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "smoothing.h"

using namespace ldt;

TimeSeriesFilters::TimeSeriesFilters(const Matrix<Tv> *source, bool byrow,
                                     bool checkNAN) {
  _source = source;
  _byrow = byrow;
  auto indexes = MatrixIndexes(source, byrow, true, false);

  if (indexes.getAnyHasMissing() || indexes.getAnyIsInvalid())
    throw std::logic_error("missing data or invalid row/column is found.");

  auto rq = _byrow ? _source->ColsCount : source->RowsCount;
  if (checkNAN) {
    auto r = indexes.biggestWithoutNaN(nullptr);
    _ranges = new IndexRange(r.StartIndex, r.EndIndex);
  } else
    _ranges = new IndexRange(0, rq);

  if (_ranges->Count() != rq) {
    _delete = true;
    auto c = _ranges->Count();

    if (_byrow) {
      Data =
          new Matrix<Tv>(new Tv[c * source->RowsCount], source->RowsCount, c);
      source->GetSub((Ti)0, _ranges->StartIndex, Data->RowsCount,
                     _ranges->Count(), Data);
    } else {
      Data =
          new Matrix<Tv>(new Tv[c * source->ColsCount], c, source->ColsCount);
      source->GetSub(_ranges->StartIndex, (Ti)0, _ranges->Count(),
                     Data->ColsCount, Data);
    }
  } else
    Data = (Matrix<Tv> *)source;
}

//#pragma region HP

void TimeSeriesFilters::calculate_hp_size(Ti &workSize, Ti &storageRow,
                                          Ti &storageCol) {
  Ti T = _ranges->Count();
  workSize = 2 * T * T;
  if (_byrow) {
    storageRow = Data->RowsCount;
    storageCol = T;
  } else {
    storageRow = T;
    storageCol = Data->ColsCount;
  }
}

void TimeSeriesFilters::calculate_hp(Tv lambda, Matrix<Tv> *storage_filt,
                                     Tv *WORK) {

  Ti T = _ranges->Count();

  Ti q = 0;
  Matrix<Tv> *mat = nullptr;
  //    if (T <= 25)
  mat = new Matrix<Tv>(0.0, &WORK[q], T, T);
  q += T * T;
  //     else
  //        mat = Matrix.Sparse(T, T);  // TODO

  auto a = 6 * lambda + 1;
  auto b = -4 * lambda;
  auto c = lambda;

  auto l2 = -2 * c;
  auto l5 = 5 * c + 1;
  mat->Set0(0, 0, 1 + c);
  mat->Set0(T - 1, T - 1, 1 + c);
  mat->Set0(0, 1, l2);
  mat->Set0(T - 2, T - 1, l2);
  mat->Set0(1, 0, l2);
  mat->Set0(T - 1, T - 2, l2);
  mat->Set0(1, 1, l5);
  mat->Set0(T - 2, T - 2, l5);
  Ti row = 0;
  for (Ti col = 2; col < T - 2; col++) {
    mat->Set0(row, col, c);
    mat->Set0(row + 1, col, b);
    mat->Set0(row + 2, col, a);
    mat->Set0(row + 3, col, b);
    mat->Set0(row + 4, col, c);
    row++;
  }
  mat->Set0(2, 0, c);
  mat->Set0(2, 1, b);
  mat->Set0(3, 1, c);
  mat->Set0(T - 3, T - 1, c);
  mat->Set0(T - 3, T - 2, b);
  mat->Set0(T - 4, T - 2, c);

  auto Wi = new Ti[T + 1];
  mat->Inv00(Wi, &WORK[q]);

  if (_byrow)
    Data->Dot(mat, storage_filt);
  else
    mat->Dot(Data, storage_filt);

  storage_filt->Restructure0(Data->RowsCount, Data->ColsCount);
  delete mat;
}

//#pragma endregion

//#pragma region BK

void TimeSeriesFilters::calculate_bk_size(Ti &storageRow, Ti &storageCol) {
  Ti T = _ranges->Count();
  if (_byrow) {
    storageRow = Data->RowsCount;
    storageCol = T;
  } else {
    storageRow = T;
    storageCol = Data->ColsCount;
  }
}

void TimeSeriesFilters::calculate_bk(Tv minf, Tv maxf, Ti k,
                                     Matrix<Tv> *storage_cycle) {
  Tv res, a, b;
  Ti i, j, p, q;

  throw std::logic_error(
      "not implemented"); // the indexation seems to be wrong. 'p' parameter in
                          // the loop below exceeds its valid values

  Ti T = _ranges->Count();
  q = 2 * k + 1;
  Tv *WORK = new Tv[q];
  auto ws = Matrix<Tv>(WORK, q);

  a = 2 * c_pi / maxf;
  b = 2 * c_pi / minf;

  ws.Data[k] = (b - a) / c_pi;
  for (i = 1; i <= k; i++) {
    res = (std::sin(b * i) - std::sin(a * i)) / (c_pi * i);
    ws.Data[k + i] = res;
    ws.Data[k - i] = res;
  }
  res = ws.Sum();
  ws.Add_in(-res / ws.length());

  if (_byrow) {
    for (j = 0; j < T; j++) {
      for (p = k; p < T - k; p++) {
        res = 0.0;
        q = p - k;
        for (i = 0; i < ws.length(); i++) {
          res += ws.Data[i] * Data->Get0(j, q); // transpose
          q++;
        }
        storage_cycle->Set0(j, p, res); // transpose
      }
    }
  } else {
    for (j = 0; j < T; j++) {
      for (p = k; p < T - k; p++) {
        res = 0.0;
        q = p - k;
        for (i = 0; i < ws.length(); i++) {
          res += ws.Data[i] * Data->Get0(q, j); // transpose
          q++;
        }
        storage_cycle->Set0(p, j, res); // transpose
      }
    }
  }
  delete[] WORK;
}

//#pragma endregion
