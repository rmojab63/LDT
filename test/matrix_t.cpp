/*
 Copyright (C) 2022-2023 Ramin Mojab
 Licensed under the GPL 3.0.
 See accompanying file LICENSE
*/

#include "matrix.h"
#include "matrix_utils.h"
#include <gtest/gtest.h>

using namespace ldt;

class matrix_t {
public:
  template <class Tw> static void initialize0() {
    Tw W1[20];
    Tw W2[20];

    auto vec1 = Matrix<Tw>(W1, 2, 1);
    vec1.Data[1] = 96;
    ASSERT_EQ(vec1.Data[1], (Tw)96);
    ASSERT_EQ(W1[1], (Tw)96);

    auto vec2 = Matrix<Tw>(W2, 2);
    vec2.SetValue(5);
    ASSERT_EQ(vec2.Get(1, 0), (Tw)5);
    ASSERT_EQ(true, vec2.Equals(Matrix<Tw>(new Tw[]{5, 5}, (Ti)2, (Ti)1)));

    Tw data[] = {1, 2, 3, 4};
    auto vec3 = Matrix<Tw>(data, (Ti)4);
    ASSERT_EQ(vec3.Get(1, 0), (Tw)2);
    vec3.Set(3, (Tw)6);
    ASSERT_EQ(vec3.Get(3, 0), (Tw)6);
    ASSERT_EQ(true,
              vec3.Equals(Matrix<Tw>(new Tw[]{1, 2, 3, 6}, (Ti)4, (Ti)1)));

    auto mat1 = Matrix<Tw>(W1, (Ti)2, (Ti)2);
    auto mat2 = Matrix<Tw>(W2, (Ti)2, (Ti)2);
    mat2.SetValue((Tw)7);
    ASSERT_EQ(mat2.Get(1, 0), (Tw)7);
    ASSERT_EQ(true,
              mat2.Equals(Matrix<Tw>(new Tw[]{7, 7, 7, 7}, (Ti)2, (Ti)2)));

    Tw data2[] = {1, 2, 3, 4, 5, 6};
    auto mat3 = Matrix<Tw>(data2, (Ti)2, (Ti)3);
    ASSERT_EQ(mat3.Get(1, 0), (Tw)2);
    ASSERT_EQ(mat3.Get(0, 1), (Tw)3);
    mat3.Set((Ti)1, (Ti)2, (Tw)9);
    ASSERT_EQ(mat3.Get(1, 2), (Tw)9);
    ASSERT_EQ(true, mat3.Equals(
                        Matrix<Tw>(new Tw[]{1, 2, 3, 4, 5, 9}, (Ti)2, (Ti)3)));

    // vector
    std::vector<Tw> V = {1, 2, 3, 4};
    auto mV = Matrix<Tw>(&V, (Ti)2, (Ti)2);
    V[1] = (Tw)100;
    ASSERT_EQ(mV.Data[1], (Tw)100);

    // VMatrix
    auto vmat = VMatrix<Tw>(2, 2);
    vmat.Mat.Set(0, 0, 2);
    ASSERT_EQ(2, vmat.Vec.at(0));

    vmat = VMatrix<Tw>({1, 2, 3, 4}, 4);
    ASSERT_EQ(3, vmat.Mat.Get(2, 0));
  }

  template <class Tw> static void iterators0() {
    Tw data[] = {1, 2, 3, 4, 5, 6};
    auto m = Matrix<Tw>(data, (Ti)2, (Ti)3);

    Tw sum = 0;
    for (auto &d : m)
      sum += d;
    ASSERT_EQ(21, sum);

    // columns
    sum = 0;
    for (auto it = m.ColBegin(1); it != m.ColEnd(1); ++it)
      sum += *it;
    ASSERT_EQ(7, sum);

    // generally you should use nested loops to iterate over multiple rows and
    // columns, however I guess the following might work:
    sum = 0;
    for (auto it = m.ColBegin(1); it != m.ColEnd(2); ++it)
      sum += *it;
    ASSERT_EQ(18, sum);

    // rows
    sum = 0;
    for (auto it = m.RowBegin(1); it != m.RowEnd(1); ++it)
      sum += *it;
    ASSERT_EQ(12, sum);

    sum = 0;
    for (int row = 0; row <= 1; ++row) {
      for (auto it = m.RowBegin(row); it != m.RowEnd(row); ++it)
        sum += *it;
    }
    ASSERT_EQ(21, sum);
  }

  template <class Tw> static void getdata0() {
    Tw data[] = {1, 2, 3, 4, 5, 6};
    auto mat = Matrix<Tw>(data, (Ti)2, (Ti)3);

    data[3] = 100;
    ASSERT_EQ(mat.Get(1, 1), (Tw)100); // its not a copy
  }

  template <class Tw> static void equality0() {

    Tw data0[] = {1, 1, 1, 1};
    auto mat0 = Matrix<Tw>(data0, (Ti)2, (Ti)2);
    Tw data1[] = {1, 2, 3, 4, 5, 6};
    auto mat1 = Matrix<Tw>(data1, (Ti)2, (Ti)3);
    Tw data2[] = {1, 2, 3, 4, 5, 6};
    auto mat2 = Matrix<Tw>(data2, (Ti)2, (Ti)3);
    Tw data3[] = {1, 2, 3, 4, 5, 6};
    auto mat3 = Matrix<Tw>(data3, (Ti)3, (Ti)1);

    ASSERT_EQ(mat0.EqualsValue((Tw)1), true);
    ASSERT_EQ(mat1.EqualsValue((Tw)2), false);
    ASSERT_EQ(mat1.Equals(mat1), true);
    ASSERT_EQ(mat1.Equals(mat2), true);
    ASSERT_EQ(mat1.Equals(mat3, (Tw)1e-16, false), false);

    // NaN
    if constexpr (std::numeric_limits<Tw>::has_quiet_NaN) {
      mat1.Set((Ti)1, (Ti)1, std::numeric_limits<Tw>::quiet_NaN());
      ASSERT_EQ(mat1.Equals(mat2), false);
      ASSERT_EQ(mat1.Equals(mat2, (Tw)1e-16, true),
                false); // one of them is NaN, the other is a number
      mat2.Set((Ti)1, (Ti)1, std::numeric_limits<Tw>::quiet_NaN());
      ASSERT_EQ(mat1.Equals(mat2, (Tw)1e-16, true, true),
                true); // both are NaN, so it must be true
    }
  }

  template <class Tw> static void issomething0() {
    // symmetric
    Tw data[] = {1, 3, 4, 3, 5, 2, 4, 2, 4};
    auto mat = Matrix<Tw>(data, 3, 3);
    ASSERT_EQ(mat.IsSymmetric(), true);

    if constexpr (std::is_same<Tw, double>()) {
      mat.Set(1, 0, mat.Get(1, 0) + (Tw)1e-8);
      ASSERT_EQ(mat.IsSymmetric(), false);
      ASSERT_EQ(mat.IsSymmetric((Tw)1e-6), true);
    }
  }

  template <class Tw> static void applycopy0() {

    Tw data[] = {1, 2, 3, 4, 5, 6};
    auto mat = Matrix<Tw>(data, (Ti)2, (Ti)3);
    std::function<Tw(Tw)> f = [](Tw x) -> Tw { return x * x; };
    mat.Apply_in(f);

    ASSERT_EQ(mat.Get(1, 1), (Tw)16);
    Tw W[6];
    auto mat2 = Matrix<Tw>(W, (Ti)2, (Ti)3);
    std::function<Tw(Tw)> g = [](Tw x) -> Tw { return (Tw)sqrt(x); };
    mat.Apply(g, mat2);
    ASSERT_EQ(mat2.Get((Ti)1, (Ti)2), (Tw)6);

    // copy
    Tw data2[] = {1, 2, 3, 4, 5, 6};
    mat = Matrix<Tw>(data2, (Ti)2, (Ti)3);

    auto matcopy = Matrix<Tw>(W, (Ti)2, (Ti)3);
    matcopy.SetValue((Tw)0);
    mat.CopyTo(matcopy);
    ASSERT_EQ(matcopy.Get((Ti)1, (Ti)2), (Tw)6);
    matcopy.Set((Ti)1, (Ti)1, (Tw)100);
    ASSERT_EQ(mat.Get((Ti)1, (Ti)1), (Tw)4);
    ASSERT_EQ(matcopy.Get((Ti)1, (Ti)1), (Tw)100);
  }

  template <class Tw> static void sub0() {
    auto data = std::vector<Tw>(25);
    std::iota(data.begin(), data.end(), (Tw)0);
    auto mat = Matrix<Tw>(&data[0], (Ti)5, (Ti)5);

    Tw data2[] = {-1, -2, -3, -4, -5, -6};
    auto mat1 = Matrix<Tw>(data2, (Ti)2, (Ti)3);
    mat.SetSub(1, 1, mat1, 0, 0, 2, 3);
    Tw W[40];
    auto getmat = Matrix<Tw>(W, (Ti)3, (Ti)3);
    mat.GetSub(1, 1, 3, 3, getmat);
    ASSERT_EQ(getmat.ToString0(',', ';', 0),
              std::string("-1,-3,-5;-2,-4,-6;8,13,18"));

    // source transposed
    data = std::vector<Tw>(25);
    std::iota(data.begin(), data.end(), (Tw)0);
    auto mat0 = Matrix<Tw>(&data[0], (Ti)5, (Ti)5);
    mat0.SetSub_t(1, 1, mat1, 0, 0, 3, 2);
    ASSERT_EQ(mat0.ToString0(',', ';', 0),
              std::string("0,5,10,15,20;1,-1,-2,16,21;2,-3,-4,17,22;3,-5,-6,18,"
                          "23;4,9,14,19,24"));

    // start from (1,1)
    getmat = Matrix<Tw>(W, (Ti)4, (Ti)4);
    getmat.SetValue(0);
    mat0.GetSub(1, 1, 3, 3, getmat, 1, 1);
    ASSERT_EQ(getmat.ToString0(',', ';', 0),
              std::string("0,0,0,0;0,-1,-2,16;0,-3,-4,17;0,-5,-6,18"));

    // vector
    std::iota(data.begin(), data.end(), (Tw)0);
    auto mat4 = Matrix<Tw>(&data[0], (Ti)16);
    Tw D[] = {-1, -2, -3, -4};
    mat1 = Matrix<Tw>(D, (Ti)4, (Ti)1);
    mat4.SetSubVector(1, mat1, 1, 3);
    Tw K[6];
    getmat = Matrix<Tw>(K, (Ti)6);
    mat4.GetSubVector(1, 6, getmat, 0);
    ASSERT_EQ(getmat.ToString0(',', ';', 0), std::string("-2;-3;-4;4;5;6"));

    // row
    Tw D0[] = {1, 2, 3, 4, 5, 6};
    auto mat5 = Matrix<Tw>(D0, (Ti)2, (Ti)3);
    Tw H[] = {10, 11, 12};
    mat5.SetRow(1, Matrix<Tw>(H, (Ti)3));
    Tw P[3] = {0, 0, 0};
    auto row = Matrix<Tw>(P, (Ti)3, (Ti)1);
    mat5.GetRow(1, row);
    Tw J[] = {10, 11, 12};
    ASSERT_EQ(row.Equals(Matrix<Tw>(J, (Ti)3, (Ti)1)), true);

    //  row   from row
    Tw HH[] = {10, 11, 12, 13, 14, 15};
    auto rowmat = Matrix<Tw>(HH, (Ti)2, (Ti)3);
    mat5.SetRowFromRow(1, rowmat, 1);
    mat5.GetRow(1, row);
    ASSERT_EQ(row.Equals(Matrix<Tw>(new Tw[]{11, 13, 15}, (Ti)3, (Ti)1)), true);

    // row    from diag
    Tw HHH[] = {10, 11, 12, 13, 14, 15, 16, 17, 18};
    auto diagmat = Matrix<Tw>(HHH, (Ti)3, (Ti)3);
    mat5.SetRowFromDiag(1, diagmat);
    mat5.GetRow(1, row);
    ASSERT_EQ(row.Equals(Matrix<Tw>(new Tw[]{10, 14, 18}, (Ti)3, (Ti)1)), true);

    // row    sub
    diagmat.SetSubRow(1, 1, rowmat, 2);

    // column
    Tw D1[] = {1, 2, 3, 4, 5, 6};
    auto mat6 = Matrix<Tw>(D1, (Ti)2, (Ti)3);
    Tw F[] = {10, 11};
    mat6.SetColumn(1, Matrix<Tw>(F, (Ti)2));
    Tw Df[] = {0, 0};
    auto col = Matrix<Tw>(Df, (Ti)2, (Ti)1);
    mat6.GetColumn(1, col);
    Tw FF[] = {10, 11};
    ASSERT_EQ(col.Equals(Matrix<Tw>(FF, (Ti)2, (Ti)1)), true);

    // column     from column
    mat6.SetColumnFromColumn(1, rowmat, 1);
    mat6.GetColumn(1, col);
    ASSERT_EQ(col.Equals(Matrix<Tw>(new Tw[]{12, 13}, (Ti)2, (Ti)1)), true);

    // column     from diag
    Tw GGG[] = {10, 11, 12, 13};
    auto diagcmat = Matrix<Tw>(GGG, (Ti)2, (Ti)2);
    mat6.SetColumnFromDiag(1, diagcmat);
    mat6.GetColumn(1, col);
    ASSERT_EQ(col.Equals(Matrix<Tw>(new Tw[]{10, 13}, (Ti)2, (Ti)1)), true);

    // diag
    diagmat.GetDiag(row);
    ASSERT_EQ(row.Equals(Matrix<Tw>(new Tw[]{10, 10, 18}, (Ti)3, (Ti)1)), true);
  }

  template <class Tw> static void sub_inds0() {

    Tw D1[20] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
                 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    auto mat = Matrix<Tw>(D1, (Ti)4, (Ti)5);

    Tw S1[6];
    Matrix<Tw> storage(S1, (Ti)2, (Ti)3);
    mat.GetSub(1, 2, std::vector<Ti>({0, 4, 3}), true, storage, (Ti)0, 0);
    ASSERT_EQ(true, storage.Equals(Matrix<Tw>(new Tw[6]{2, 3, 18, 19, 14, 15},
                                              (Ti)2, (Ti)3)));

    mat.GetSub(1, 2, std::vector<Ti>({2, 1}), true, storage, (Ti)0, 0,
               true); // exclude indexes
    ASSERT_EQ(true, storage.Equals(
                        Matrix<Tw>(new Tw[6]{2, 3, 14, 15, 18, 19}, (Ti)2,
                                   (Ti)3))); // the order might be different

    mat.GetSub(1, 3, std::vector<Ti>({2, 1}), false, storage, 0, 0);
    ASSERT_EQ(true, storage.Equals(Matrix<Tw>(new Tw[6]{7, 6, 11, 10, 15, 14},
                                              (Ti)2, (Ti)3)));

    mat.GetSub(1, 3, std::vector<Ti>({3, 0}), false, storage, 0, 0,
               true); // exclude
    ASSERT_EQ(true, storage.Equals(Matrix<Tw>(new Tw[6]{6, 7, 10, 11, 14, 15},
                                              (Ti)2, (Ti)3)));

    auto v = std::vector<Ti>({1, 3, 2});
    auto w = std::vector<Ti>({2, 1});
    mat.GetSub(w, v, storage, 0, 0);
    ASSERT_EQ(true, storage.Equals(Matrix<Tw>(new Tw[6]{7, 6, 15, 14, 11, 10},
                                              (Ti)2, (Ti)3)));

    // source start
    Tw S2[12];
    Matrix<Tw> storage1(0, S2, (Ti)3, (Ti)4);
    mat.GetSub(1, 2, std::vector<Ti>({0, 3, 4}), true, storage1, 1, 1);
    ASSERT_EQ(true, storage1.Equals(Matrix<Tw>(
                        new Tw[12]{0, 0, 0, 0, 2, 3, 0, 14, 15, 0, 18, 19},
                        (Ti)3, (Ti)4)));

    mat.GetSub(1, 3, std::vector<Ti>({1, 2}), false, storage1, 1, 1);
    ASSERT_EQ(true, storage1.Equals(Matrix<Tw>(
                        new Tw[12]{0, 0, 0, 0, 6, 7, 0, 10, 11, 0, 14, 15},
                        (Ti)3, (Ti)4)));

    v = std::vector<Ti>({1, 2, 3});
    w = std::vector<Ti>({1, 2});
    mat.GetSub(w, v, storage1, 1, 1);
    ASSERT_EQ(true, storage1.Equals(Matrix<Tw>(
                        new Tw[12]{0, 0, 0, 0, 6, 7, 0, 10, 11, 0, 14, 15},
                        (Ti)3, (Ti)4)));
  }

  template <class Tw> static void translateindexes0() {
    Tw D1[20] = {2, 5, 6, 4, 1, 1, 8, 1, 7, 1, 8, 2, 6, 2, 5, 8, 5, 2, 5, 8};
    auto mat = Matrix<Tw>(D1, (Ti)4, (Ti)5);
    Ti c, r;

    // I test translate here
    mat.TranslateIndex(0, r, c);
    ASSERT_EQ(r, (Ti)0);
    ASSERT_EQ(c, (Ti)0);
    mat.TranslateIndex(1, r, c);
    ASSERT_EQ(r, (Ti)1);
    ASSERT_EQ(c, (Ti)0);
    mat.TranslateIndex(2, r, c);
    ASSERT_EQ(r, (Ti)2);
    ASSERT_EQ(c, (Ti)0);
    mat.TranslateIndex(3, r, c);
    ASSERT_EQ(r, (Ti)3);
    ASSERT_EQ(c, (Ti)0);
    mat.TranslateIndex(4, r, c);
    ASSERT_EQ(r, (Ti)0);
    ASSERT_EQ(c, (Ti)1);
    mat.TranslateIndex(5, r, c);
    ASSERT_EQ(r, (Ti)1);
    ASSERT_EQ(c, (Ti)1);
  }

  template <class Tw> static void maxmin_ind0() {

    Tw D1[20] = {2, 5, 6, 4, 1, 1, 8, 1, 7, 1, 8, 2, 6, 2, 5, 8, 5, 2, 5, 8};
    auto mat = Matrix<Tw>(D1, (Ti)4, (Ti)5);
    Ti c, r;
    ASSERT_EQ(mat.MaximumInRow(2, c), (Tw)8);
    ASSERT_EQ(c, (Ti)1);
    ASSERT_EQ(mat.MinimumInRow(3, c), (Tw)1);
    ASSERT_EQ(c, (Ti)1);
    ASSERT_EQ(mat.MaximumInColumn(2, c), (Tw)8);
    ASSERT_EQ(c, (Ti)2);
    ASSERT_EQ(mat.MinimumInColumn(3, c), (Tw)2);
    ASSERT_EQ(c, (Ti)1);

    // global max/min
    ASSERT_EQ(mat.max(r, c), (Tw)8);
    ASSERT_EQ(r, (Ti)2);
    ASSERT_EQ(c, (Ti)1);
    ASSERT_EQ(mat.min(r, c), (Tw)1);
    ASSERT_EQ(r, (Ti)0);
    ASSERT_EQ(c, (Ti)1);
  }

  template <class Tw> static void helpers0() {
    Tw W[10];
    // diagonal
    auto diagm = Matrix<Tw>(W, (Ti)3, (Ti)3);
    Matrix<Tw>::Diagonal(diagm, 2, 3);
    ASSERT_EQ(diagm.ToString0(',', ';', 0), std::string("2,3,3;3,2,3;3,3,2"));

    // triangular

    // with diagonal
    Ti n = 3;
    Tw h[]{1, 2, 3, 4, 5, 6};
    auto elem = Matrix<Tw>(h, (Ti)6);
    auto storage = Matrix<Tw>(W, (Ti)3, (Ti)3);

    // up -- by row
    storage.SetValue(0);
    Matrix<Tw>::MakeTriangular(storage, elem, 1, true, true);
    ASSERT_EQ(storage.ToString0(',', ';', 0), std::string("1,2,3;0,4,5;0,0,6"));
    Matrix<Tw>::MakeTriangular(storage, elem, 0, true, true);
    ASSERT_EQ(storage.ToString0(',', ';', 0), std::string("1,2,3;2,4,5;3,5,6"));

    // up -- by column
    // res = Matrix.Triangular(n, elem, true, true, false);
    // Assert.IsTrue(Matrix.Parse("1,2,4;0,3,5;0,0,6").IsEqualTo(res));
    // res = Matrix.Triangular(n, elem, null, true, false);
    // Assert.IsTrue(Matrix.Parse("1,2,4;2,3,5;4,5,6").IsEqualTo(res));

    // low -- by column
    storage.SetValue(0);
    Matrix<Tw>::MakeTriangular(storage, elem, -1, true, false);
    ASSERT_EQ(storage.ToString0(',', ';', 0), std::string("1,0,0;2,4,0;3,5,6"));

    // no diagonal
    n = 3;
    Tw m[] = {2, 3, 5};
    auto elem1 = Matrix<Tw>(m, (Ti)3);
    storage.SetValue(0);
    Matrix<Tw>::MakeTriangular(storage, elem1, 1, false, true);
    ASSERT_EQ(storage.ToString0(',', ';', 0), std::string("0,2,3;0,0,5;0,0,0"));
    storage.SetValue(0);
    Matrix<Tw>::MakeTriangular(storage, elem1, 0, false, true);
    ASSERT_EQ(storage.ToString0(',', ';', 0), std::string("0,2,3;2,0,5;3,5,0"));
    storage.SetValue(0);
    Matrix<Tw>::MakeTriangular(storage, elem1, -1, false, false);
    ASSERT_EQ(storage.ToString0(',', ';', 0), std::string("0,0,0;2,0,0;3,5,0"));
  }

  template <class Tw> static void remove_nan0() {

    Tw data1[] = {NAN, 1,   2, 3, 4,   5, 6,   NAN, NAN,
                  NAN, NAN, 2, 3, NAN, 5, NAN, NAN, NAN};
    auto mat = Matrix<Tw>(data1, (Ti)9, (Ti)2);
    mat.RemoveNanVector_in();
    ASSERT_EQ((Ti)9, mat.length());

    // remove nan any column

    Tw data20[] = {NAN, 2,   3,   4,  NAN, 6,  7,  8,  9,
                   NAN, NAN, NAN, 13, 14,  15, 16, 17, NAN};
    auto mat3 = Matrix<Tw>(data20, (Ti)3, (Ti)6);
    // Tw data[18];
    // Matrix<Tw> mat = Matrix<Tw>(data, 6, 3);
    // mat1.Transpose(mat);
    mat3.RemoveColumnsAnyNan_in(false);
    ASSERT_EQ((Ti)2, mat3.ColsCount);
    ASSERT_EQ((Tw)7, mat3.Get(0, 0));
    ASSERT_EQ((Tw)8, mat3.Get(1, 0));
    ASSERT_EQ((Tw)9, mat3.Get(2, 0));
    ASSERT_EQ((Tw)13, mat3.Get(0, 1));

    // locate NAN
    Tw data5[] = {NAN, 2,   3,   4,  NAN, 6,  7,  8,  9,
                  NAN, NAN, NAN, 13, 14,  15, 16, 17, NAN};
    auto mat8 = Matrix<Tw>(data5, (Ti)3, (Ti)6);
    // TODO: column

    // row
    Tw data0[18];
    Matrix<Tw> mat0 = Matrix<Tw>(data0, (Ti)6, (Ti)3);
    mat8.Transpose(mat0);

    auto v = std::vector<Ti>();
    mat0.GetAnyNanRow(v);
    ASSERT_EQ(v.at(0), (Ti)2);
    ASSERT_EQ(v.at(1), (Ti)4);
    auto w = std::vector<Ti>({1, 2});
    v.clear();
    mat0.GetAnyNanRow(v, false, &w);
    ASSERT_EQ(v.at(0), (Ti)0);
    ASSERT_EQ(v.at(1), (Ti)2);
    ASSERT_EQ(v.at(2), (Ti)4);
  }

  template <class Tw> static void range_nan0() {
    // column
    Tw data1[] = {NAN, 1,   2, 3, 4,   5, 6,   NAN, NAN,
                  NAN, NAN, 2, 3, NAN, 5, NAN, NAN, NAN};
    bool hasmissing;
    auto mat = Matrix<Tw>(data1, (Ti)9, (Ti)2);
    auto range = mat.GetRangeColumn(hasmissing, 1);
    ASSERT_EQ(range.StartIndex, (Ti)2);
    ASSERT_EQ(range.EndIndex, (Ti)5);
    ASSERT_EQ(hasmissing, true);

    range = mat.GetRangeColumn(hasmissing, 0);
    ASSERT_EQ(range.StartIndex, (Ti)1);
    ASSERT_EQ(range.EndIndex, (Ti)6);
    ASSERT_EQ(hasmissing, false);

    // whole Matrix
    Tw data2[] = {NAN, NAN, NAN, 2,   NAN, 0,  (Tw)-0.1, 0, NAN,
                  NAN, NAN, NAN, 1,   0,   -1, 1,        1, NAN,
                  NAN, NAN, NAN, 1,   0,   -1, 0,        0, NAN,
                  NAN, 9,   8,   NAN, 7,   6,  7,        3, 0};
    mat = Matrix<Tw>(data2, (Ti)9, (Ti)4);
    range = mat.GetRange(hasmissing);
    ASSERT_EQ(range.StartIndex, (Ti)3);
    ASSERT_EQ(range.EndIndex, (Ti)7);

    Tw data3[] = {NAN, NAN, NAN, 2, 1, -2, NAN, NAN, NAN, 1,   1, -1,
                  1,   NAN, NAN, 0, 1, 1,  NAN, 1,   0,   NAN, 1, 0};
    mat = Matrix<Tw>(data3, (Ti)6, (Ti)4);
    range = mat.GetRange(hasmissing);
    ASSERT_EQ(range.StartIndex, (Ti)3);
    ASSERT_EQ(range.EndIndex, (Ti)5);
  }

  template <class Tw> static void la_elementwise0() {

    Tw A0[] = {1, 2, 3, 4, 5, 6};
    auto mat1 = Matrix<Tw>(A0, (Ti)3, (Ti)2);
    Tw B0[] = {2, 3, 4, 5, 6, 7};
    auto mat2 = Matrix<Tw>(B0, (Ti)3, (Ti)2);
    Tw Z[] = {0, 0, 0, 0, 0, 0};
    auto storage = Matrix<Tw>(Z, (Ti)3, (Ti)2);

    // add
    mat1.Add(2, storage);
    Tw A[] = {3, 4, 5, 6, 7, 8};
    ASSERT_EQ(storage.Equals(Matrix<Tw>(A, (Ti)3, (Ti)2)), true);
    storage.Add_in((Tw)1);
    Tw B[] = {4, 5, 6, 7, 8, 9};
    ASSERT_EQ(storage.Equals(Matrix<Tw>(B, (Ti)3, (Ti)2)), true);
    mat1.Add(mat2, storage);
    Tw C[] = {3, 5, 7, 9, 11, 13};
    ASSERT_EQ(storage.Equals(Matrix<Tw>(C, (Ti)3, (Ti)2)), true);
    storage.Add_in(mat2);
    Tw D[]{5, 8, 11, 14, 17, 20};
    ASSERT_EQ(storage.Equals(Matrix<Tw>(D, (Ti)3, (Ti)2)), true);

    // multiply
    mat1.Multiply(2, storage);
    Tw E[] = {2, 4, 6, 8, 10, 12};
    ASSERT_EQ(storage.Equals(Matrix<Tw>(E, (Ti)3, (Ti)2)), true);
    storage.Multiply_in(2);
    Tw F[] = {4, 8, 12, 16, 20, 24};
    ASSERT_EQ(storage.Equals(Matrix<Tw>(F, (Ti)3, (Ti)2)), true);
    mat1.Multiply(mat2, storage);
    Tw H[] = {2, 6, 12, 20, 30, 42};
    ASSERT_EQ(storage.Equals(Matrix<Tw>(H, (Ti)3, (Ti)2)), true);
    storage.Multiply_in(mat1);
    Tw I[] = {2, 12, 36, 80, 150, 6 * 42};
    ASSERT_EQ(storage.Equals(Matrix<Tw>(I, (Ti)3, (Ti)2)), true);

    // ...
  }

  template <class Tw> static void Dot0() {
    // vector-vector
    Tw s1[] = {1, 2, 3, 4, 5, 6};
    auto mat1 = Matrix<Tw>(s1, (Ti)6, (Ti)1);
    Tw s2[] = {2, 3, 4, 5, 6, 7};
    auto mat2 = Matrix<Tw>(s2, (Ti)6, (Ti)1);
    ASSERT_NEAR(mat1.VectorDotVector(mat2), (Tw)(2 + 6 + 12 + 20 + 30 + 42),
                (Tw)1e-10);

    // Matrix-vector
    Tw s3[] = {1, 2, 3};
    auto mat12 = Matrix<Tw>(s3, (Ti)3);
    Tw s4[] = {4, 5, 6, 7, 8, 9};
    auto mat22 = Matrix<Tw>(s4, (Ti)2, (Ti)3);
    Tw ss2[] = {1, 1};
    auto storage = Matrix<Tw>(ss2, (Ti)2, (Ti)1);
    mat22.DotVector(mat12, storage);
    Tw s5[] = {40, 19 + 27};
    auto m1 = Matrix<Tw>(s5, (Ti)2, (Ti)1);
    ASSERT_EQ(storage.Equals(m1), true);

    Tw s6[] = {1, 2, 3};
    auto mat13 = Matrix<Tw>(s6, (Ti)3);
    Tw s7[] = {4, 6, 8, 5, 7, 9};
    auto mat23 = Matrix<Tw>(s7, (Ti)3, (Ti)2);
    Tw s8[] = {1, 1};
    storage = Matrix<Tw>(s8, (Ti)2, (Ti)1);
    mat23.tDotVector(mat13, storage);
    Tw s9[] = {40, 19 + 27};
    auto m2 = Matrix<Tw>(s9, (Ti)2, (Ti)1);
    ASSERT_EQ(storage.Equals(m2), true);

    // Matrix-Matrix
    Tw s10[] = {1, 2, 3, 4, 5, 6};
    auto mat14 = Matrix<Tw>(s10, (Ti)2, (Ti)3);
    Tw s11[] = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    auto mat24 = Matrix<Tw>(s11, (Ti)3, (Ti)4);
    Tw s12[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    storage = Matrix<Tw>(s12, (Ti)2, (Ti)4);
    mat14.Dot(mat24, storage);
    Tw s13[] = {49, 64, 76, 100, 103, 136, 130, 172};
    auto m3 = Matrix<Tw>(s13, (Ti)2, (Ti)4);
    ASSERT_EQ(storage.Equals(m3), true);

    Tw s14[] = {1, 2, 3, 4, 5, 6};
    auto mat15 = Matrix<Tw>(s14, (Ti)2, (Ti)3);
    Tw s15[] = {4, 7, 10, 13, 5, 8, 11, 14, 6, 9, 12, 15};
    auto mat25 = Matrix<Tw>(s15, (Ti)4, (Ti)3);
    Tw s16[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    storage = Matrix<Tw>(s16, (Ti)2, (Ti)4);
    mat15.DotTr(mat25, storage);
    Tw s17[] = {49, 64, 76, 100, 103, 136, 130, 172};
    auto m4 = Matrix<Tw>(s17, (Ti)2, (Ti)4);
    ASSERT_EQ(storage.Equals(m4), true);

    Tw s18[] = {1, 3, 5, 2, 4, 6};
    auto mat16 = Matrix<Tw>(s18, (Ti)3, (Ti)2);
    Tw s19[] = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    auto mat26 = Matrix<Tw>(s19, (Ti)3, (Ti)4);
    Tw s20[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    storage = Matrix<Tw>(s20, (Ti)2, (Ti)4);
    mat16.TrDot(mat26, storage);
    Tw s21[] = {49, 64, 76, 100, 103, 136, 130, 172};
    auto m5 = Matrix<Tw>(s21, (Ti)2, (Ti)4);
    ASSERT_EQ(storage.Equals(m5), true);

    Tw s22[] = {1, 3, 5, 2, 4, 6};
    auto mat17 = Matrix<Tw>(s22, (Ti)3, (Ti)2);
    Tw s23[] = {4, 7, 10, 13, 5, 8, 11, 14, 6, 9, 12, 15};
    auto mat27 = Matrix<Tw>(s23, (Ti)4, (Ti)3);
    Tw s24[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    storage = Matrix<Tw>(s24, (Ti)2, (Ti)4);
    mat17.TrDotTr(mat27, storage);
    Tw s27[] = {49, 64, 76, 100, 103, 136, 130, 172};
    auto m6 = Matrix<Tw>(s27, (Ti)2, (Ti)4);
    ASSERT_EQ(storage.Equals(m6), true);

    // Matrix-Matrix<Tw> (symmetric)
    Tw s37[] = {1, NAN, 2, 1};
    auto mat18 = Matrix<Tw>(s37, (Ti)2, (Ti)2);
    Tw s38[] = {4, 5, 6, 7, 8, 9};
    auto mat28 = Matrix<Tw>(s38, (Ti)2, (Ti)3);
    Tw s39[] = {1, 1, 1, 1, 1, 1};
    storage = Matrix<Tw>(s39, (Ti)2, (Ti)3);
    mat18.SymDot(mat28, storage, true, true);
    Tw s40[] = {14, 13, 20, 19, 26, 25};
    auto m7 = Matrix<Tw>(s40, (Ti)2, (Ti)3);
    ASSERT_EQ(storage.Equals(m7), true);

    Tw s41[] = {4, 5, 6, 7, 8, 9};
    auto mat29 = Matrix<Tw>(s41, (Ti)3, (Ti)2);
    Tw s42[] = {1, 1, 1, 1, 1, 1};
    auto storage1 = Matrix<Tw>(s42, (Ti)3, (Ti)2);
    mat29.DotSym(mat18, storage1, false, true);
    Tw s43[] = {18, 21, 24, 15, 18, 21};
    auto m8 = Matrix<Tw>(s43, (Ti)3, (Ti)2);
    ASSERT_EQ(storage1.Equals(m8), true);
  }

  template <class Tw> static void dot_diag0() {

    Tw M1[]{1, 3, 5, 2, 4, 6};
    auto mat1 = Matrix<Tw>(M1, (Ti)3, (Ti)2);
    Tw M2[]{10, 12};
    auto dmat = Matrix<Tw>(M2, (Ti)2);
    Tw M3[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    auto S = Matrix<Tw>(M3, (Ti)3, (Ti)2);
    mat1.DotDiag(dmat, S);

    Tw M4[]{10, 0, 0, 12};
    auto dmat0 = Matrix<Tw>(M4, (Ti)2, (Ti)2);
    Tw M5[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    auto S0 = Matrix<Tw>(M5, (Ti)3, (Ti)2);
    mat1.Dot(dmat0, S0);
    ASSERT_EQ(S.Equals(S0), true);

    // Left
    Tw M6[]{10, 12, 87};
    auto dmat2 = Matrix<Tw>(M6, (Ti)3);
    dmat2.DiagDot(mat1, S);
    Tw M7[]{10, 0, 0, 0, 12, 0, 0, 0, 87};
    auto dmat3 = Matrix<Tw>(M7, (Ti)3, (Ti)3);
    dmat3.Dot(mat1, S0);
    ASSERT_EQ(S0.Equals(S), true);
  }

  template <class Tw> static void dot_kron0() {
    auto mat1 = Matrix<Tw>(new std::vector<Tw>{1, 3, 5, 2, 4, 6}, (Ti)2, (Ti)3);
    auto mat2 = Matrix<Tw>(new std::vector<Tw>{10, 100}, (Ti)1, (Ti)2);
    auto kmat = Matrix<Tw>(new std::vector<Tw>(12), (Ti)2, (Ti)6);
    mat1.Kron(mat2, kmat);
    auto S = Matrix<Tw>(new std::vector<Tw>{10, 30, 100, 300, 50, 20, 500, 200,
                                            40, 60, 400, 600},
                        (Ti)2, (Ti)6);
    ASSERT_EQ(S.Equals(S), true);

    //  I
    auto kmatI = Matrix<Tw>(new std::vector<Tw>(8), (Ti)2, (Ti)4);
    mat2.KronIden(2, kmatI);
    auto S1 = Matrix<Tw>(new std::vector<Tw>{10, 0, 0, 10, 100, 0, 0, 100},
                         (Ti)2, (Ti)4);
    ASSERT_EQ(kmatI.Equals(S1), true);

    // I first
    mat2.IdenKron(2, kmatI);
    auto S11 = Matrix<Tw>(new std::vector<Tw>{10, 0, 100, 0, 0, 10, 0, 100},
                          (Ti)2, (Ti)4);
    ASSERT_EQ(kmatI.Equals(S11), true);

    // t I
    auto tkmatI = Matrix<Tw>(new std::vector<Tw>(8), (Ti)4, (Ti)2);
    mat2.TrKronIden(2, tkmatI);
    auto S2 = Matrix<Tw>(new std::vector<Tw>{10, 0, 100, 0, 0, 10, 0, 100},
                         (Ti)4, (Ti)2);
    ASSERT_EQ(tkmatI.Equals(S2), true);
  }

  template <class Tw> static void dot_AA0() {
    Tw s[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    auto mat = Matrix<Tw>(s, (Ti)2, (Ti)6);

    // AA'

    auto storage1 = Matrix<Tw>(new Tw[4], (Ti)2, (Ti)2);
    mat.Dot_AAt(storage1);
    Tw r1[] = {286, 322, 322, 364};
    ASSERT_EQ(storage1.Equals(Matrix<Tw>(r1, (Ti)2, (Ti)2)), true);

    // A'A

    auto storage2 = Matrix<Tw>(new Tw[36], (Ti)6, (Ti)6);
    mat.Dot_AtA(storage2);
    Tw r2[] = {5,  11, 17,  23,  29,  35,  11, 25, 39,  53,  67,  81,
               17, 39, 61,  83,  105, 127, 23, 53, 83,  113, 143, 173,
               29, 67, 105, 143, 181, 219, 35, 81, 127, 173, 219, 265};
    ASSERT_EQ(storage2.Equals(Matrix<Tw>(r2, (Ti)6, (Ti)6)), true);
  }

  template <class Tw> static void Inv0() {
    auto mat = Matrix<Tw>(new std::vector<Tw>{1, 2, 3, 4}, (Ti)2, (Ti)2);
    auto storage = Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0}, (Ti)2, (Ti)2);
    mat.Inv(storage);
    // multiply
    auto mul = Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0}, (Ti)2, (Ti)2);
    mat.Dot(storage, mul);
    ASSERT_EQ(
        mul.Equals(Matrix<Tw>(new std::vector<Tw>{1, 0, 0, 1}, (Ti)2, (Ti)2)),
        true);
  }

  template <class Tw> static void Inv2x20() {

    auto mat = Matrix<Tw>(new std::vector<Tw>{1, 2, 3, 4}, (Ti)2, (Ti)2);
    auto storage = Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0}, (Ti)2, (Ti)2);
    mat.CopyTo00(storage);
    storage.Inv2x2();
    // multiply
    auto mul = Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0}, (Ti)2, (Ti)2);
    mat.Dot(storage, mul);
    ASSERT_EQ(
        mul.Equals(Matrix<Tw>(new std::vector<Tw>{1, 0, 0, 1}, (Ti)2, (Ti)2)),
        true);
  }

  template <class Tw> static void Chol0() {
    auto mat = Matrix<Tw>(new std::vector<Tw>{5, 1, 1, 3}, (Ti)2, (Ti)2);
    auto storage = Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0}, (Ti)2, (Ti)2);
    auto info = mat.Chol(storage, true);
    auto mulstorage = Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0}, (Ti)2, (Ti)2);
    storage.TrDot(storage, mulstorage);
    ASSERT_EQ(mat.Equals(mulstorage, (Tw)1e-12), true);
    // lower
    info = mat.Chol(storage, false);
    storage.DotTr(storage, mulstorage);
    ASSERT_EQ(mat.Equals(mulstorage, (Tw)1e-12), true);
  }

  template <class Tw> static void SVD0() {
    Tw eps = (Tw)1e-10;
    if constexpr (std::is_same<float, Tw>()) {
      eps = (Tw)1e-5;
    }

    auto mat1 = Matrix<Tw>(new std::vector<Tw>{1, 3, 5, 2, 4, 6}, (Ti)2, (Ti)3);
    auto svd1 = MatrixSvd<Tv>(2, 3);
    auto W1 = new Tw[svd1.WorkSize];
    auto S1 = new Tw[svd1.StorageSize];
    svd1.Calculate(mat1, S1, W1);

    auto S0 = Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0, 0, 0}, (Ti)2, (Ti)3);
    for (Ti i = 0; i < svd1.S.length(); i++)
      S0.Set(i, i, svd1.S.Data[i]);
    auto US = Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0, 0, 0}, (Ti)2, (Ti)3);
    svd1.U.Dot(S0, US);
    auto USVT = Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0, 0, 0}, (Ti)2, (Ti)3);
    US.Dot(svd1.VT, USVT);
    ASSERT_EQ(mat1.Equals(USVT, (Tw)eps), true);

    //      check jobs VT:'A'
    auto svd2 = MatrixSvd<Tv>(2, 3, 'N', 'A');
    auto W2 = new Tw[svd2.WorkSize];
    auto S2 = new Tw[svd2.StorageSize];
    svd2.Calculate(mat1, S2, W2);
    ASSERT_EQ(svd2.VT.Equals(svd1.VT, (Tw)eps), true);

    //      check job VT:'S'
    auto svd3 = MatrixSvd<Tv>(2, 3, 'N', 'S');
    auto W3 = new Tw[svd3.WorkSize];
    auto S3 = new Tw[svd3.StorageSize];
    svd3.Calculate(mat1, S3, W3);
    ASSERT_NEAR(svd3.VT.Get(0, 0), (Tw)-0.3180747629905496, (Tw)eps);
    ASSERT_NEAR(svd3.VT.Get(1, 1), (Tw)-0.8147843370349828, (Tw)eps);

    //      check jobs U:'A'
    auto svd4 = MatrixSvd<Tv>(2, 3, 'A', 'N');
    auto W4 = new Tw[svd4.WorkSize];
    auto S4 = new Tw[svd4.StorageSize];
    svd4.Calculate(mat1, S4, W4);
    ASSERT_EQ(svd4.U.Equals(svd1.U, (Tw)eps), true);

    //      check job U:'S'
    auto svd5 = MatrixSvd<Tv>(2, 3, 'S', 'S');
    auto W5 = new Tw[svd5.WorkSize];
    auto S5 = new Tw[svd5.StorageSize];
    svd5.Calculate(mat1, S5, W5);
    ASSERT_NEAR(svd5.VT.Get(0, 0), (Tw)-0.3180747629905496, (Tw)eps);
    ASSERT_NEAR(svd5.VT.Get(1, 1), (Tw)-0.8147843370349828, (Tw)eps);
    ASSERT_EQ(svd5.U.Equals(svd1.U, (Tw)eps), true); // it is minimum

    delete[] W1;
    delete[] S1;
    delete[] W2;
    delete[] S2;
    delete[] W3;
    delete[] S3;
    delete[] W4;
    delete[] S4;
    delete[] W5;
    delete[] S5;
  }

  template <class Tw> static void svd00() { // same but with larger rows
    Tw eps = (Tw)1e-10;
    if constexpr (std::is_same<float, Tw>()) {
      eps = (Tw)1e-5;
    }

    auto mat1 = Matrix<Tw>(new std::vector<Tw>{1, 3, 5, 2, 4, 6}, (Ti)3, (Ti)2);
    auto svd1 = MatrixSvd<Tv>(3, 2);
    auto W1 = new Tw[svd1.WorkSize];
    auto S1 = new Tw[svd1.StorageSize];
    svd1.Calculate(mat1, S1, W1);

    auto S0 = Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0, 0, 0}, (Ti)3, (Ti)2);
    for (Ti i = 0; i < svd1.S.length(); i++)
      S0.Set(i, i, svd1.S.Data[i]);
    auto US = Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0, 0, 0}, (Ti)3, (Ti)2);
    svd1.U.Dot(S0, US);
    auto USVT = Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0, 0, 0}, (Ti)3, (Ti)2);
    US.Dot(svd1.VT, USVT);
    ASSERT_EQ(mat1.Equals(USVT, (Tw)eps), true);

    //      check jobs VT:'A'
    auto svd2 = MatrixSvd<Tv>(3, 2, 'N', 'A');
    auto W2 = new Tw[svd2.WorkSize];
    auto S2 = new Tw[svd2.StorageSize];
    svd2.Calculate(mat1, S2, W2);
    ASSERT_EQ(svd2.VT.Equals(svd1.VT, (Tw)eps), true);

    //      check job VT:'S'
    auto svd3 = MatrixSvd<Tv>(3, 2, 'N', 'S');
    auto W3 = new Tw[svd3.WorkSize];
    auto S3 = new Tw[svd3.StorageSize];
    svd3.Calculate(mat1, S3, W3);
    ASSERT_NEAR(std::abs(svd3.VT.Get(0, 0)), (Tw)0.6196294838293402, (Tw)eps);
    ASSERT_NEAR(std::abs(svd3.VT.Get(1, 1)), (Tw)0.6196294838293402, (Tw)eps);

    //      check jobs U:'A'
    auto svd4 = MatrixSvd<Tv>(3, 2, 'A', 'N');
    auto W4 = new Tw[svd4.WorkSize];
    auto S4 = new Tw[svd4.StorageSize];
    svd4.Calculate(mat1, S4, W4);
    ASSERT_EQ(svd4.U.Equals(svd1.U, (Tw)eps), true);

    //      check job U:'S'
    auto svd5 = MatrixSvd<Tv>(3, 2, 'S', 'S');
    auto W5 = new Tw[svd5.WorkSize];
    auto S5 = new Tw[svd5.StorageSize];
    svd5.Calculate(mat1, S5, W5);
    ASSERT_NEAR(std::abs(svd5.U.Get(0, 0)), (Tw)0.2298476964000714, (Tw)eps);
    ASSERT_NEAR(std::abs(svd5.U.Get(2, 1)), (Tw)0.4018960334334317, (Tw)eps);
    ASSERT_EQ(svd5.VT.Equals(svd1.VT, (Tw)eps), true); // it is minimum

    delete[] W1;
    delete[] S1;
    delete[] W2;
    delete[] S2;
    delete[] W3;
    delete[] S3;
    delete[] W4;
    delete[] S4;
    delete[] W5;
    delete[] S5;
  }

  template <class Tw> static void norm_T0() {
    auto mat = Matrix<Tw>(new std::vector<Tw>{1, 2, 3, 4, 5, 6}, (Ti)2, (Ti)3);
    auto n = mat.Norm('I');
    ASSERT_EQ((Tw)12, n);
    n = mat.Norm('1');
    ASSERT_EQ((Tw)11, n);
    n = mat.Norm('F');
    ASSERT_NEAR((Tw)9.539392, n, (Tw)1e-4);
    n = mat.Norm('M');
    ASSERT_EQ((Tw)6, n);
  }

  template <class Tw> static void determinants0() {

    auto mat = Matrix<Tw>(new std::vector<Tw>{2, -1, 0, -1, 2, -1, 0, -1, 2},
                          (Ti)3, (Ti)3);
    auto d = mat.Det_pd0();
    ASSERT_NEAR(d, (Tw)4, (Tw)1e-4);

    auto mat1 = Matrix<Tw>(new std::vector<Tw>{2, -1, 0, -1, 2, -1, 0, -1, 2},
                           (Ti)3, (Ti)3);
    d = mat1.Det();
    ASSERT_NEAR(d, (Tw)4, (Tw)1e-4);

    auto mat2 = Matrix<Tw>(new std::vector<Tw>{1, 4, 2, 6, 3, 6, 7, 4, 3, 6},
                           (Ti)3, (Ti)3);
    d = mat2.Det();
    ASSERT_NEAR(d, (Tw)87, (Tw)1e-4);
  }

  /*
  TEST(Matrix_T,QR_dec) {


          auto mat = Matrix<Tw>(new std::vector<Tw>{ 1, 2, 3, 4, 5, 6, 7, 8 },
  4, 2); auto Q = Matrix<Tw>(new std::vector<Tw>{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 4, 4); auto R = Matrix<Tw>(new
  std::vector<Tw>{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 2, 2); auto d =
  mat.QR(Q, R); auto res = Matrix<Tw>(new std::vector<Tw>{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },4, 2); Q.Dot(R, res);
          ASSERT_EQ(res.Equals(mat, 1e-15), true);

          mat = Matrix<Tw>(new std::vector<Tw> {1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
  11, 12}, 3, 4 ); Q = Matrix<Tw>(new std::vector<Tw>{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },3, 3); R = Matrix<Tw>(new
  std::vector<Tw>{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 4, 4); d =
  mat.QR(Q, R); res = Matrix<Tw>(new std::vector<Tw>{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 3, 4); Q.Dot(R, res);
          ASSERT_EQ(res.Equals(mat, 1e-15), true);

          mat = Matrix<Tw>(new std::vector<Tw> {1, 2, 3, 4, 5, 6, 7, 8, 9}, 3,
  3); Q = Matrix<Tw>(new std::vector<Tw>{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 3, 3); R = Matrix<Tw>(new
  std::vector<Tw>{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 3, 3); d =
  mat.QR(Q, R); res = Matrix<Tw>(new std::vector<Tw>{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 3, 3); Q.Dot(R, res);
          ASSERT_EQ(res.Equals(mat, 1e-15), true);
  }
  */

  template <class Tw> static void SolveTrian0() {
    auto A = Matrix<Tw>(new std::vector<Tw>{1, 2, 3, 0, 2, 1, 0, 0, 3}, (Ti)3,
                        (Ti)3);
    auto b = Matrix<Tw>(new std::vector<Tw>{1, 2, 3, 3, 2, 1}, (Ti)3, (Ti)2);
    auto d = A.SolveTrian(b, false, false, false);
    ASSERT_EQ(b.Equals(Matrix<Tw>(new std::vector<Tw>{1, 0, 0, 3, -2, -2},
                                  (Ti)3, (Ti)2),
                       (Tw)1e-2),
              true);

    auto b1 = Matrix<Tw>(new std::vector<Tw>{1, 2, 3, 3, 2, 1}, (Ti)3, (Ti)2);
    d = A.SolveTrian(b1, false, true, false);
    ASSERT_EQ(b1.Equals(Matrix<Tw>(new std::vector<Tw>{-3, (Tw)0.5, 1, (Tw)0.33,
                                                       (Tw)0.83, (Tw)0.33},
                                   (Ti)3, (Ti)2),
                        (Tw)1e-2),
              true);
  }

  template <class Tw> static void solve_pos0() {
    auto A = Matrix<Tw>(new std::vector<Tw>{2, -1, 0, -1, 2, -1, 0, -1, 2},
                        (Ti)3, (Ti)3);
    auto b = Matrix<Tw>(new std::vector<Tw>{3, 2, 1}, (Ti)3, (Ti)1);
    auto info = A.SolvePos(b, false);

    auto A0 = Matrix<Tw>(new std::vector<Tw>{2, -1, 0, -1, 2, -1, 0, -1, 2},
                         (Ti)3, (Ti)3);
    auto b0 = Matrix<Tw>(new std::vector<Tw>{3, 2, 1}, (Ti)3, (Ti)1);
    auto b1 = Matrix<Tw>(new Tw[3], (Ti)3, (Ti)1);
    auto Ainv = Matrix<Tw>(new Tw[9], (Ti)3, (Ti)3);
    A0.Inv(Ainv);
    Ainv.Dot(b0, b1);

    ASSERT_EQ(b.Equals(b1, (Tw)1e-5), true);
  }

  template <class Tw> static void transpose0() {
    auto A = Matrix<Tw>(new std::vector<Tw>{1, 2, 3, 4, 5, 6}, (Ti)3, (Ti)2);
    auto storage =
        Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0, 0, 0, 0}, (Ti)2, (Ti)3);
    A.Transpose(storage);
    ASSERT_EQ(storage.Equals(Matrix<Tw>(new std::vector<Tw>{1, 4, 2, 5, 3, 6},
                                        (Ti)2, (Ti)3)),
              true);

    // in-place
    // square
    auto A1 = Matrix<Tw>(new std::vector<Tw>{1, 2, 3, 4, 5, 6, 7, 8, 9}, (Ti)3,
                         (Ti)3);
    A1.Transpose();
    ASSERT_EQ(
        A1.Equals(Matrix<Tw>(new std::vector<Tw>{1, 4, 7, 2, 5, 8, 3, 6, 9},
                             (Ti)3, (Ti)3)),
        true);

    // not square
    auto v0 = std::vector<Tw>{1, 2, 3, 4, 5, 6, 7, 8};
    A = Matrix<Tw>(&v0, 2, 4);
    A.Transpose();
    ASSERT_EQ(
        A.Equals(Matrix<Tw>(new std::vector<Tw>{1, 3, 5, 7, 2, 4, 6, 8}, 4, 2)),
        true);

    auto v1 = std::vector<Tw>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    A = Matrix<Tw>(&v1, 3, 4);
    A.Transpose();
    ASSERT_EQ(
        A.Equals(Matrix<Tw>(
            new std::vector<Tw>{1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12}, 4, 3)),
        true);
  }

  template <class Tw> static void statistics_cols0() {
    Tw eps = (Tw)1e-8;
    if constexpr (std::is_same<float, Tw>()) {
      eps = (Tw)1e-4;
    }

    Tw data[] = {1, 3, 5, 2, 4, 6, 7, 18};
    auto mat = Matrix<Tw>(data, (Ti)4, (Ti)2);

    // means
    auto me = Matrix<Tw>(new std::vector<Tw>{0, 0}, (Ti)2);
    auto v = std::vector<Ti>();
    mat.ColumnsMean(me, v);
    ASSERT_NEAR(me.Get(0, 0), (Tw)2.75, eps);

    // var Matrix
    auto me1 = Matrix<Tw>(new std::vector<Tw>{0, 0, 0, 0}, (Ti)2, (Ti)2);
    v = std::vector<Ti>();
    mat.ColumnsVariance(me1, v, true);

    ASSERT_EQ(
        me1.Equals(
            Matrix<Tw>(new std::vector<Tw>{(Tw)2.91666667, (Tw)-1.0833333333,
                                           (Tw)-1.0833333333, (Tw)39.583333333},
                       (Ti)2, (Ti)2),
            eps),
        true);
  }

  template <class Tw> static void sort_v_T0() {

    Tw data[] = {3, 2, 4, 5, 1, 6, 7, 9};
    auto mat = Matrix<Tw>(data, (Ti)8, (Ti)1);

    auto storage = Matrix<Tw>(new Tw[8], (Ti)8, (Ti)1);
    std::vector<Ti> indexes = std::vector<Ti>();
    mat.SortIndicesVector(indexes, true);

    mat.SortByVector(storage, indexes);
    ASSERT_EQ((Tw)1, storage.Data[0]);
  }
};

TEST(Matrix_T, initialize) {
  matrix_t::initialize0<Tv>();
  matrix_t::initialize0<Ti>();
}

TEST(Matrix_T, iterators) {
  matrix_t::iterators0<Tv>();
  matrix_t::iterators0<Ti>();
}

TEST(Matrix_T, getdata) {
  matrix_t::getdata0<Tv>();
  matrix_t::getdata0<Ti>();
}

TEST(Matrix_T, issomething) {
  matrix_t::issomething0<Tv>();
  matrix_t::issomething0<Ti>();
}

TEST(Matrix_T, equality) {
  matrix_t::equality0<Tv>();
  matrix_t::equality0<Ti>();
}

TEST(Matrix_T, applycopy) {
  matrix_t::applycopy0<Tv>();
  matrix_t::applycopy0<Ti>();
}

TEST(Matrix_T, sub) {
  matrix_t::sub0<Tv>();
  matrix_t::sub0<Ti>();
}

TEST(Matrix_T, sub_inds) {

  matrix_t::sub_inds0<Tv>();
  matrix_t::sub_inds0<Ti>();
}

TEST(Matrix_T, translateindexes) {
  matrix_t::translateindexes0<Tv>();
  matrix_t::translateindexes0<Ti>();
}

TEST(Matrix_T, maxmin_ind) {
  matrix_t::maxmin_ind0<Tv>();
  matrix_t::maxmin_ind0<Ti>();
}

TEST(Matrix_T, helpers) {
  matrix_t::helpers0<Tv>();
  matrix_t::helpers0<Ti>();
}

TEST(Matrix_T, remove_nan) { matrix_t::remove_nan0<Tv>(); };

TEST(Matrix_T, range_nan) { matrix_t::range_nan0<Tv>(); }

TEST(Matrix_T, la_elementwise) { matrix_t::la_elementwise0<Tv>(); }

TEST(Matrix_T, dot) { matrix_t::Dot0<Tv>(); }

TEST(Matrix_T, dot_diag) { matrix_t::dot_diag0<Tv>(); }

TEST(Matrix_T, dot_kron) {
  matrix_t::dot_kron0<Tv>();
  matrix_t::dot_kron0<Ti>();
}

TEST(Matrix_T, inv) { matrix_t::Inv0<Tv>(); }

TEST(Matrix_T, dot_AA) { matrix_t::dot_AA0<Tv>(); }

TEST(Matrix_T, Inv2x2) { matrix_t::Inv2x20<Tv>(); }

TEST(Matrix_T, Chol) { matrix_t::Chol0<Tv>(); }

TEST(Matrix_T, svd) {

  matrix_t::SVD0<Tv>();

  matrix_t::svd00<Tv>();
}

TEST(Matrix_T, norm_T) { matrix_t::norm_T0<Tv>(); }

TEST(Matrix_T, determinants) { matrix_t::determinants0<Tv>(); }

TEST(Matrix_T, SolveTrian) { matrix_t::SolveTrian0<Tv>(); }

TEST(Matrix_T, solve_pos) { matrix_t::solve_pos0<Tv>(); }

TEST(Matrix_T, transpose) {
  matrix_t::transpose0<Tv>();
  matrix_t::transpose0<Ti>();
}

TEST(Matrix_T, variance_T) {
  Tv data[] = {1, 3, 5, 2, 4, 6, 7, 18};
  auto mat = Matrix<Tv>(data, (Ti)8);
  Tv mean = 0;
  Tv var = mat.Variance(mean, true, false);
  ASSERT_NEAR(28.5, var, 1e-15);

  Tv data2[] = {1, 3, NAN, 5, 2, 4, 6, 7, NAN, 18};
  auto mat2 = Matrix<Tv>(data2, (Ti)10);
  var = mat2.Variance(mean, true, true);
  ASSERT_NEAR(28.5, var, 1e-15);

  // column
  Tv data3[] = {2, 2, 2,   2, 2, 2, 2, 2, 2,   11,
                1, 3, NAN, 5, 2, 4, 6, 7, NAN, 18};
  auto mat3 = Matrix<Tv>(data3, (Ti)10, (Ti)2);
  Ti count = 0;
  var = mat3.VarianceColumn(1, mean, count, true, true);
  ASSERT_NEAR(28.5, var, 1e-15);
}

TEST(Matrix_T, covariance_T) {
  Tv data[] = {1, 3, 5, 2, 4, 6, 7, 18, 3, 4, 5, 7};
  auto mat = Matrix<Tv>(data, (Ti)4, (Ti)3);
  Tv mean1, mean2;
  Ti count = 0;
  auto covar = mat.CovarianceColumn(0, 1, mean1, mean2, count, true, false);
  ASSERT_NEAR(-1.083333333333333, covar, 1e-15);

  Tv data2[] = {1, NAN, 3, 5, 1, 2, 4, 4, 6, 7, NAN, 18, 3, 3, 4, 5, NAN, 7};
  auto mat2 = Matrix<Tv>(data2, (Ti)6, (Ti)3);
  covar = mat2.CovarianceColumn(0, 1, mean1, mean2, count, true, true);
  ASSERT_NEAR(-1.083333333333333, covar, 1e-15);
}

TEST(Matrix_T, correlation_T) {
  Tv data[] = {1, 3, 5, 2, 4, 6, 7, 18, 3, 4, 5, 7};
  auto mat = Matrix<Tv>(data, (Ti)4, (Ti)3);
  Tv mean1, mean2, var1, var2;
  Ti count = 0;
  auto cor =
      mat.CorrelationColumn(0, 1, mean1, mean2, var1, var2, count, true, false);
  ASSERT_NEAR(-0.1008236754628326, cor, 1e-15);

  Tv data2[] = {1, NAN, 3, 5, 1, 2, 4, 4, 6, 7, NAN, 18, 3, 3, 4, 5, NAN, 7};
  auto mat2 = Matrix<Tv>(data2, (Ti)6, (Ti)3);
  cor =
      mat2.CorrelationColumn(0, 1, mean1, mean2, var1, var2, count, true, true);
  ASSERT_NEAR(-0.1008236754628326, cor, 1e-15);

  // one constant,  both constant
  Tv data3[] = {1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 3, 4, 5, 6};
  auto mat3 = Matrix<Tv>(data3, (Ti)6, (Ti)3);
  cor =
      mat3.CorrelationColumn(0, 1, mean1, mean2, var1, var2, count, true, true);
  ASSERT_EQ(1.0, cor);

  cor =
      mat3.CorrelationColumn(1, 2, mean1, mean2, var1, var2, count, true, true);
  ASSERT_EQ(true, std::isnan(cor));
}

TEST(Matrix_T, statistics_cols) { matrix_t::statistics_cols0<Tv>(); }

TEST(Matrix_T, sort_v_T) { matrix_t::sort_v_T0<Tv>(); }

TEST(Matrix_T, sym_set) {
  auto mat_diag = MatrixSym<true>(new Tv[15], (Ti)5);

  mat_diag.Set(1, 1, 100.0);
  ASSERT_EQ(100.0, mat_diag.Get(1, 1));
  mat_diag.Set(1, 2, 200.0);
  ASSERT_EQ(200.0, mat_diag.Get(1, 2));
  ASSERT_EQ(200.0, mat_diag.Get(2, 1));
  mat_diag.Set(2, 1, 300.0);
  ASSERT_EQ(300.0, mat_diag.Get(2, 1));
  ASSERT_EQ(300.0, mat_diag.Get(1, 2));
  mat_diag.Set(3, 3, 400.0);
  ASSERT_EQ(400.0, mat_diag.Get(3, 3));
  mat_diag.Set(3, 4, 400.0);
  ASSERT_EQ(400.0, mat_diag.Get(3, 4));
  ASSERT_EQ(400.0, mat_diag.Get(4, 3));
  mat_diag.Set(4, 3, 500.0);
  ASSERT_EQ(500.0, mat_diag.Get(3, 4));
  ASSERT_EQ(500.0, mat_diag.Get(4, 3));

  auto mat_nodiag = MatrixSym<false>(new Tv[10], (Ti)5);

  mat_diag.Set(1, 2, 200.0);
  ASSERT_EQ(200.0, mat_diag.Get(1, 2));
  ASSERT_EQ(200.0, mat_diag.Get(2, 1));
  mat_diag.Set(2, 1, 300.0);
  ASSERT_EQ(300.0, mat_diag.Get(2, 1));
  ASSERT_EQ(300.0, mat_diag.Get(1, 2));
  mat_diag.Set(3, 4, 400.0);
  ASSERT_EQ(400.0, mat_diag.Get(3, 4));
  ASSERT_EQ(400.0, mat_diag.Get(4, 3));
  mat_diag.Set(4, 3, 500.0);
  ASSERT_EQ(500.0, mat_diag.Get(3, 4));
  ASSERT_EQ(500.0, mat_diag.Get(4, 3));
}

TEST(Matrix_T, sym_get) {
  auto mat_diag = MatrixSym<true>(
      new Tv[15]{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, (Ti)5);

  ASSERT_EQ(1.0, mat_diag.Get(0, 0));
  ASSERT_EQ(2.0, mat_diag.Get(1, 0));
  ASSERT_EQ(2.0, mat_diag.Get(0, 1));
  ASSERT_EQ(7.0, mat_diag.Get(2, 1));
  ASSERT_EQ(7.0, mat_diag.Get(1, 2));
  ASSERT_EQ(14.0, mat_diag.Get(3, 4));
  ASSERT_EQ(14.0, mat_diag.Get(4, 3));

  auto mat_nodiag =
      MatrixSym<false>(new Tv[10]{2, 3, 4, 5, 7, 8, 9, 11, 12, 14}, (Ti)5);

  // same must be true. but we cannot have diagonal

  ASSERT_EQ(2.0, mat_nodiag.Get(1, 0));
  ASSERT_EQ(2.0, mat_nodiag.Get(0, 1));
  ASSERT_EQ(7.0, mat_nodiag.Get(2, 1));
  ASSERT_EQ(7.0, mat_nodiag.Get(1, 2));
  ASSERT_EQ(14.0, mat_diag.Get(3, 4));
  ASSERT_EQ(14.0, mat_diag.Get(4, 3));
}

TEST(Matrix_T, dataset) {
  Tv data1[] = {NAN, 1,   2, 3, 4,   5, 6,   NAN, NAN,
                NAN, NAN, 2, 3, NAN, 5, NAN, NAN, NAN};
  auto mat = Matrix<Tv>(data1, 9, 2);

  auto dataset = Dataset<Tv>(9, 2, true);
  auto S = new Tv[dataset.StorageSize];
  dataset.Calculate(mat, new std::vector<Ti>({1}), S);
  ASSERT_EQ(dataset.Result.RowsCount, 3);

  dataset.Calculate(mat, new std::vector<Ti>({0, 1}), S);
  ASSERT_EQ(dataset.Result.RowsCount, 3);
}

TEST(Matrix_T, datasetTs) {

  Tv data1[] = {1,  4,  7,  11, 14, 17, 20, 2,  5,  8, 12,
                15, 18, 21, 3,  6,  9,  13, 16, 19, 22};
  auto mat1 = Matrix<Tv>(data1, 7, 3);
  mat1.Transpose();
  auto ds1 = ldt::DatasetTs<true>(3, 7, false, true);
  auto S = new Tv[ds1.StorageSize];
  ds1.Data(mat1);
  ASSERT_EQ(false, ds1.HasMissingData);
  auto v = std::vector<Ti>({1, 2});
  ds1.Update(&v, S);
  ASSERT_EQ(ds1.Result.Get(0, 1), 5);
  delete[] S;

  auto ds2 = ldt::DatasetTs<false>(3, 7, false, true);
  S = new Tv[ds2.StorageSize];
  ds2.Data(mat1);
  ASSERT_EQ(false, ds2.HasMissingData);
  v = std::vector<Ti>({1, 2});
  ds2.Update(&v, S);
  ASSERT_EQ(ds2.Result.Get(0, 1), 7);
  delete[] S;

  Tv data2[] = {NAN, 4,  7,  NAN, 14, 17, NAN, 2,  5,   8,  12,
                15,  18, 21, NAN, 6,  9,  13,  16, NAN, NAN};

  auto mat2 = Matrix<Tv>(data2, 7, 3);
  mat2.Transpose();
  auto ds11 = ldt::DatasetTs<true>(3, 7, true, true);
  S = new Tv[ds11.StorageSize];
  ds11.Data(mat2);
  ASSERT_EQ(true, ds11.HasMissingData);
  v = std::vector<Ti>({1, 2});
  ds11.Update(&v, S);
  ASSERT_EQ(ds11.Result.Get(0, 1), 8);
  delete[] S;

  auto ds22 = ldt::DatasetTs<false>(3, 7, true, true);
  S = new Tv[ds22.StorageSize];
  ds22.Data(mat2);
  ASSERT_EQ(false, ds22.HasMissingData);
  v = std::vector<Ti>({1, 2});
  ds22.Update(&v, S);
  ASSERT_EQ(ds22.Result.Get(0, 1), 7);
  delete[] S;
}

TEST(Matrix_T, Interpolate) {

  // by column

  Tv data[] = {NAN, NAN, 2,   3,   NAN, 4, 5, NAN, NAN, NAN, 9, 7,  NAN,
               NAN, NAN, 3,   NAN, 2,   1, 2, 3,   2,   4,   5, 1,  7,
               5,   6,   7,   4,   3,   3, 3, 3,   NAN, NAN, 2, 3,  NAN,
               4,   5,   NAN, NAN, NAN, 9, 7, NAN, NAN, NAN, 3, NAN};
  auto mat = Matrix<Tv>(data, 17, 3);
  Ti count = 0;
  auto range = mat.InterpolateColumn(count, 0);
  ASSERT_NEAR(mat.Get(7, 0), 6.0, 1e-16);
  ASSERT_NEAR(mat.Get(8, 0), 7.0, 1e-16);
  ASSERT_NEAR(mat.Get(9, 0), 8.0, 1e-16);

  range = mat.InterpolateColumn(count, 1);
  ASSERT_EQ(count, 0);

  // by row
  mat.Transpose();
  range = mat.InterpolateRow(count, 2);
  ASSERT_NEAR(mat.Get(2, 7), 6.0, 1e-16);
  ASSERT_NEAR(mat.Get(2, 8), 7.0, 1e-16);
  ASSERT_NEAR(mat.Get(2, 9), 8.0, 1e-16);
}

TEST(Matrix_T, standardized) {
  Tv data1[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1};
  auto mat = Matrix<Tv>(data1, 9, 2);

  auto mats1 = MatrixStandardized<Tv>(9, 2, false, true, false);
  auto S1 = new Tv[mats1.StorageSize];
  mats1.Calculate(mat, S1);
  ASSERT_NEAR(mats1.Result.Get(0, 0), -4, 1e-16);
  ASSERT_NEAR(mats1.Result.Get(0, 1), 4, 1e-16);

  // scale
  auto mats2 = MatrixStandardized<Tv>(9, 2, false, true, true);
  auto S2 = new Tv[mats2.StorageSize];
  mats2.Calculate(mat, S2);
  ASSERT_NEAR(mats2.Result.Get(0, 0), -1.4605934866804429, 1e-16);
  ASSERT_NEAR(mats2.Result.Get(0, 1), 1.4605934866804429, 1e-16);

  // scale no demean
  auto mats21 = MatrixStandardized<Tv>(9, 2, false, false, true);
  auto S21 = new Tv[mats2.StorageSize];
  mats21.Calculate(mat, S21);
  ASSERT_NEAR(mats21.Result.Get(0, 0), 0.3651483716701107, 1e-16);
  ASSERT_NEAR(mats21.Result.Get(0, 1), 3.2863353450309964, 1e-16);

  // zero variance
  Tv data2[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6,
                7, 8, 9, 1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 8, 7,
                6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  auto mat1 = Matrix<Tv>(data2, 9, 5);

  auto mats3 = MatrixStandardized<Tv>(9, 5, true, true, true);
  auto S3 = new Tv[mats3.StorageSize];
  mats3.Calculate(mat1, S3);
  ASSERT_EQ(0, mats3.RemovedZeroVar.at(0));
  ASSERT_EQ(2, mats3.RemovedZeroVar.at(1));
  ASSERT_EQ(4, mats3.RemovedZeroVar.at(2));
  ASSERT_NEAR(mats3.Result.Get(0, 0), -1.4605934866804429, 1e-16);
  ASSERT_NEAR(mats3.Result.Get(0, 1), 1.4605934866804429, 1e-16);

  delete[] S1;
  delete[] S2;
  delete[] S21;
  delete[] S3;
}