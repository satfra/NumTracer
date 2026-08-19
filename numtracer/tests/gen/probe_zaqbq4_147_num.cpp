#define KOKKOS_INLINE_FUNCTION inline
#define KOKKOS_FORCEINLINE_FUNCTION inline
#define KOKKOS_FUNCTION
#define __host__
#define __device__
#include <complex>
#include <cmath>
#include <random>
#include <cstdio>
#include <cstring>
#include "ZAqbq4_147_num_kernels.hh"
template<int N, class T> static inline T powr(T x){ T r=T(1); int n=N<0?-N:N; for(int i=0;i<n;++i) r*=x; return N<0?T(1)/r:r; }
using std::pow; using std::sqrt; using std::sin; using std::cos; using std::tan; using std::exp; using std::log; using std::fma; using std::fabs;
static inline std::complex<double> fma(const std::complex<double>&a,const std::complex<double>&b,const std::complex<double>&c){return a*b+c;}
template<class T> using complex = std::complex<T>;
static inline double ntStub(double seed, double x){ double h = std::sin(seed*0.1031 + x*0.3127 + 1.7)*43758.5453; return 0.4 + 0.5*(h - std::floor(h)); }
static inline double ntRe(double x) { return x; }
template <class T> static inline double ntRe(const T &z) { return z.real(); }
static inline double ntIm(double) { return 0.0; }
template <class T> static inline double ntIm(const T &z) { return z.imag(); }
static inline auto probe_full(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k)
{
  const double cosl1p1 = cos1;
  const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  double fenv[(DiFfRG::zaqbq4_147_num::nenv) > 0 ? (DiFfRG::zaqbq4_147_num::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * l1);
  const double dr_1 = powr<-1>(l1);
  const double dr_2 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_3 = ntStub(85453., 1. * k);
  const double dr_4 = ntStub(85453., 1. * l1);
  const double dr_5 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_7 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_8 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_9 = ntStub(30415., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_10 = sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_11 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_12 = ntStub(85453., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_13 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_14 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_15 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_16 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_17 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_18 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_19 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_20 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  DiFfRG::zaqbq4_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20);
  const auto _interp1 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp2 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp3 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp4 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp5 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _interp6 = ntStub(85453., 1. * k);
  const auto _interp7 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp8 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp9 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
  const auto _interp10 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp11 = ntStub(29441., 1. * l1);
  const auto _interp12 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp13 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp14 = ntStub(13492., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp15 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp16 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp17 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp18 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp20 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp21 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp22 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp23 = ntStub(12206., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp24 = ntStub(62323., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp25 = ntStub(29441., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp26 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _interp27 = ntStub(14883., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp28 = ntStub(38926., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp29 = ntStub(12206., 0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp30 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp31 = ntStub(30415., 1. * sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp32 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _interp33 = ntStub(85453., 1. * sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp34 = ntStub(14883., 0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp35 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp36 = ntStub(38926., 0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp37 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp38 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp39 = ntStub(27191., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp40 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp41 = ntStub(67246.00000000001, 1. * k);
  const auto _interp42 = ntStub(85453., 1.02 * k);
  const auto _interp43 = ntStub(30415., 1. * l1);
  const auto _interp44 = ntStub(85453., 1. * l1);
  const auto _interp45 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp46 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _interp47 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp48 = ntStub(29441., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp49 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
  const auto _interp50 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp51 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp55 = ntStub(29441., 1. * sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp56 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(fma(_interp12, _interp13, fma(_interp11, powr<2>(l1), 0.)));
  const auto _den3 = powr<-1>(fma(_interp12, _interp26, fma(_interp25, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _den5 = powr<-1>(fma(_interp12, _interp56, fma(_interp55, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _den6 = powr<-2>(fma(powr<-2>(l1) * powr<2>(_interp40 * _interp6 + _interp44 * l1), powr<2>(l1), powr<2>(_interp43)));
  const auto _den7 = powr<-1>(fma(_interp12, _interp49, fma(_interp48, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
  const auto _den8 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), powr<2>(_interp47 + _interp46 * _interp6 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))), powr<2>(_interp45)));
  const auto _den9 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), powr<2>(_interp7 + _interp5 * _interp6 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)))), powr<2>(_interp4)));
  const auto _den10 = powr<-1>(fma(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), powr<2>(_interp33 + _interp32 * _interp6 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)))), powr<2>(_interp31)));
  const auto _den11 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), powr<2>(_interp10 + _interp6 * _interp9 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))), powr<2>(_interp8)));
  const auto _interp19 = ntStub(14883., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp52 = ntStub(62323., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp53 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp54 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  // clang-format off
  using _T = decltype(complex<double>(0.,-0.02083333333333333) + _den1 + _den10 + _den11 + _den2 + _den3 + _den5 + _den6 + _den7 + _den8 + _den9 + DiFfRG::zaqbq4_147_num::tr0(fenv) + DiFfRG::zaqbq4_147_num::tr100(fenv) + DiFfRG::zaqbq4_147_num::tr101(fenv) + DiFfRG::zaqbq4_147_num::tr102(fenv) + DiFfRG::zaqbq4_147_num::tr103(fenv) + DiFfRG::zaqbq4_147_num::tr104(fenv) + DiFfRG::zaqbq4_147_num::tr105(fenv) + DiFfRG::zaqbq4_147_num::tr106(fenv) + DiFfRG::zaqbq4_147_num::tr107(fenv) + DiFfRG::zaqbq4_147_num::tr10(fenv) + DiFfRG::zaqbq4_147_num::tr11(fenv) + DiFfRG::zaqbq4_147_num::tr12(fenv) + DiFfRG::zaqbq4_147_num::tr13(fenv) + DiFfRG::zaqbq4_147_num::tr14(fenv) + DiFfRG::zaqbq4_147_num::tr15(fenv) + DiFfRG::zaqbq4_147_num::tr16(fenv) + DiFfRG::zaqbq4_147_num::tr17(fenv) + DiFfRG::zaqbq4_147_num::tr18(fenv) + DiFfRG::zaqbq4_147_num::tr19(fenv) + DiFfRG::zaqbq4_147_num::tr1(fenv) + DiFfRG::zaqbq4_147_num::tr20(fenv) + DiFfRG::zaqbq4_147_num::tr21(fenv) + DiFfRG::zaqbq4_147_num::tr22(fenv) + DiFfRG::zaqbq4_147_num::tr23(fenv) + DiFfRG::zaqbq4_147_num::tr24(fenv) + DiFfRG::zaqbq4_147_num::tr25(fenv) + DiFfRG::zaqbq4_147_num::tr26(fenv) + DiFfRG::zaqbq4_147_num::tr27(fenv) + DiFfRG::zaqbq4_147_num::tr28(fenv) + DiFfRG::zaqbq4_147_num::tr29(fenv) + DiFfRG::zaqbq4_147_num::tr2(fenv) + DiFfRG::zaqbq4_147_num::tr30(fenv) + DiFfRG::zaqbq4_147_num::tr31(fenv) + DiFfRG::zaqbq4_147_num::tr32(fenv) + DiFfRG::zaqbq4_147_num::tr33(fenv) + DiFfRG::zaqbq4_147_num::tr34(fenv) + DiFfRG::zaqbq4_147_num::tr35(fenv) + DiFfRG::zaqbq4_147_num::tr36(fenv) + DiFfRG::zaqbq4_147_num::tr37(fenv) + DiFfRG::zaqbq4_147_num::tr38(fenv) + DiFfRG::zaqbq4_147_num::tr39(fenv) + DiFfRG::zaqbq4_147_num::tr3(fenv) + DiFfRG::zaqbq4_147_num::tr40(fenv) + DiFfRG::zaqbq4_147_num::tr41(fenv) + DiFfRG::zaqbq4_147_num::tr42(fenv) + DiFfRG::zaqbq4_147_num::tr43(fenv) + DiFfRG::zaqbq4_147_num::tr44(fenv) + DiFfRG::zaqbq4_147_num::tr45(fenv) + DiFfRG::zaqbq4_147_num::tr46(fenv) + DiFfRG::zaqbq4_147_num::tr47(fenv) + DiFfRG::zaqbq4_147_num::tr48(fenv) + DiFfRG::zaqbq4_147_num::tr49(fenv) + DiFfRG::zaqbq4_147_num::tr4(fenv) + DiFfRG::zaqbq4_147_num::tr50(fenv) + DiFfRG::zaqbq4_147_num::tr51(fenv) + DiFfRG::zaqbq4_147_num::tr52(fenv) + DiFfRG::zaqbq4_147_num::tr53(fenv) + DiFfRG::zaqbq4_147_num::tr54(fenv) + DiFfRG::zaqbq4_147_num::tr55(fenv) + DiFfRG::zaqbq4_147_num::tr56(fenv) + DiFfRG::zaqbq4_147_num::tr57(fenv) + DiFfRG::zaqbq4_147_num::tr58(fenv) + DiFfRG::zaqbq4_147_num::tr59(fenv) + DiFfRG::zaqbq4_147_num::tr5(fenv) + DiFfRG::zaqbq4_147_num::tr60(fenv) + DiFfRG::zaqbq4_147_num::tr61(fenv) + DiFfRG::zaqbq4_147_num::tr62(fenv) + DiFfRG::zaqbq4_147_num::tr63(fenv) + DiFfRG::zaqbq4_147_num::tr64(fenv) + DiFfRG::zaqbq4_147_num::tr65(fenv) + DiFfRG::zaqbq4_147_num::tr66(fenv) + DiFfRG::zaqbq4_147_num::tr67(fenv) + DiFfRG::zaqbq4_147_num::tr68(fenv) + DiFfRG::zaqbq4_147_num::tr69(fenv) + DiFfRG::zaqbq4_147_num::tr6(fenv) + DiFfRG::zaqbq4_147_num::tr70(fenv) + DiFfRG::zaqbq4_147_num::tr71(fenv) + DiFfRG::zaqbq4_147_num::tr72(fenv) + DiFfRG::zaqbq4_147_num::tr73(fenv) + DiFfRG::zaqbq4_147_num::tr74(fenv) + DiFfRG::zaqbq4_147_num::tr75(fenv) + DiFfRG::zaqbq4_147_num::tr76(fenv) + DiFfRG::zaqbq4_147_num::tr77(fenv) + DiFfRG::zaqbq4_147_num::tr78(fenv) + DiFfRG::zaqbq4_147_num::tr79(fenv) + DiFfRG::zaqbq4_147_num::tr7(fenv) + DiFfRG::zaqbq4_147_num::tr80(fenv) + DiFfRG::zaqbq4_147_num::tr81(fenv) + DiFfRG::zaqbq4_147_num::tr82(fenv) + DiFfRG::zaqbq4_147_num::tr83(fenv) + DiFfRG::zaqbq4_147_num::tr84(fenv) + DiFfRG::zaqbq4_147_num::tr85(fenv) + DiFfRG::zaqbq4_147_num::tr86(fenv) + DiFfRG::zaqbq4_147_num::tr87(fenv) + DiFfRG::zaqbq4_147_num::tr88(fenv) + DiFfRG::zaqbq4_147_num::tr89(fenv) + DiFfRG::zaqbq4_147_num::tr8(fenv) + DiFfRG::zaqbq4_147_num::tr90(fenv) + DiFfRG::zaqbq4_147_num::tr91(fenv) + DiFfRG::zaqbq4_147_num::tr92(fenv) + DiFfRG::zaqbq4_147_num::tr93(fenv) + DiFfRG::zaqbq4_147_num::tr94(fenv) + DiFfRG::zaqbq4_147_num::tr95(fenv) + DiFfRG::zaqbq4_147_num::tr96(fenv) + DiFfRG::zaqbq4_147_num::tr97(fenv) + DiFfRG::zaqbq4_147_num::tr98(fenv) + DiFfRG::zaqbq4_147_num::tr99(fenv) + DiFfRG::zaqbq4_147_num::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp6 + _interp7 + _interp8 + _interp9 + cosl1p2 + k + l1 + p);
  // clang-format on
  _T _acc{};
  { // subkernel 1
    const auto _cse1_k1 = powr<-2>(p);
    // clang-format off
    _acc += fma(complex<double>(0.,-0.005208333333333332), _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr105(fenv) * _interp1 * _interp18 * _interp19 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(complex<double>(0.,-0.002604166666666666), _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr107(fenv) * _interp17 * _interp18 * _interp19 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(complex<double>(0.,-0.01041666666666667), _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr93(fenv) * _interp1 * _interp19 * _interp2 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(complex<double>(0.,-0.005208333333333332), _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr97(fenv) * _interp17 * _interp19 * _interp2 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(0.005208333333333332, _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr106(fenv) * _interp18 * _interp19 * _interp20 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(0.01041666666666667, _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr94(fenv) * _interp19 * _interp2 * _interp20 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(-0.01041666666666667, _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr95(fenv) * _interp1 * _interp19 * _interp21 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(-0.005208333333333332, _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr98(fenv) * _interp17 * _interp19 * _interp21 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(complex<double>(0.,-0.01041666666666667), _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr96(fenv) * _interp19 * _interp20 * _interp21 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(-0.01041666666666667, _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr101(fenv) * _interp1 * _interp18 * _interp22 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(-0.005208333333333332, _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr104(fenv) * _interp17 * _interp18 * _interp22 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(-0.02083333333333333, _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr83(fenv) * _interp1 * _interp2 * _interp22 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(-0.01041666666666667, _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr90(fenv) * _interp17 * _interp2 * _interp22 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(complex<double>(0.,-0.01041666666666667), _cse1_k1 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr102(fenv) * _interp18 * _interp20 * _interp22 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), 0.))))))))))))));
    // clang-format on
  }
  { // subkernel 2
    const auto _cse1_k2 = powr<-2>(p);
    const auto _cse2_k2 = -_interp12;
    // clang-format off
    _acc += fma(complex<double>(0.,-0.02083333333333333), _cse1_k2 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr84(fenv) * _interp2 * _interp20 * _interp22 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), fma(complex<double>(0.,0.02083333333333333), _cse1_k2 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr87(fenv) * _interp1 * _interp21 * _interp22 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), fma(complex<double>(0.,0.01041666666666667), _cse1_k2 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr92(fenv) * _interp17 * _interp21 * _interp22 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), fma(-0.02083333333333333, _cse1_k2 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr88(fenv) * _interp20 * _interp21 * _interp22 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), fma(complex<double>(0.,0.01041666666666667), _cse1_k2 * _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr78(fenv) * _interp18 * _interp23 * _interp24 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), fma(complex<double>(0.,0.02083333333333333), _cse1_k2 * _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr72(fenv) * _interp2 * _interp23 * _interp24 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), fma(0.02083333333333333, _cse1_k2 * _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr74(fenv) * _interp21 * _interp23 * _interp24 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), fma(complex<double>(0.,0.005208333333333332), _cse1_k2 * _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr80(fenv) * _interp18 * _interp24 * _interp27 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), fma(complex<double>(0.,0.01041666666666667), _cse1_k2 * _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr76(fenv) * _interp2 * _interp24 * _interp27 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), fma(0.01041666666666667, _cse1_k2 * _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr77(fenv) * _interp21 * _interp24 * _interp27 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), fma(-0.01041666666666667, _cse1_k2 * _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr79(fenv) * _interp18 * _interp24 * _interp28 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), fma(-0.02083333333333333, _cse1_k2 * _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr73(fenv) * _interp2 * _interp24 * _interp28 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), fma(complex<double>(0.,0.02083333333333333), _cse1_k2 * _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr75(fenv) * _interp21 * _interp24 * _interp28 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), fma(complex<double>(0.,-0.01041666666666667), _cse1_k2 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr99(fenv) * _interp1 * _interp18 * _interp3 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k2 + _interp16) * powr<6>(k))), 0.))))))))))))));
    // clang-format on
  }
  { // subkernel 3
    const auto _cse1_k3 = powr<-2>(p);
    const auto _cse2_k3 = -_interp12;
    // clang-format off
    _acc += fma(complex<double>(0.,-0.005208333333333332), _cse1_k3 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr103(fenv) * _interp17 * _interp18 * _interp3 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), fma(complex<double>(0.,-0.02083333333333333), _cse1_k3 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr81(fenv) * _interp1 * _interp2 * _interp3 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), fma(complex<double>(0.,-0.01041666666666667), _cse1_k3 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr89(fenv) * _interp17 * _interp2 * _interp3 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), fma(0.01041666666666667, _cse1_k3 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr100(fenv) * _interp18 * _interp20 * _interp3 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), fma(0.02083333333333333, _cse1_k3 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr82(fenv) * _interp2 * _interp20 * _interp3 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), fma(-0.02083333333333333, _cse1_k3 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr85(fenv) * _interp1 * _interp21 * _interp3 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), fma(-0.01041666666666667, _cse1_k3 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr91(fenv) * _interp17 * _interp21 * _interp3 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), fma(complex<double>(0.,-0.02083333333333333), _cse1_k3 * _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr86(fenv) * _interp20 * _interp21 * _interp3 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), fma(complex<double>(0.,0.02083333333333333), _cse1_k3 * _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr63(fenv) * _interp24 * _interp29 * _interp30 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), fma(complex<double>(0.,0.01041666666666667), _cse1_k3 * _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr67(fenv) * _interp24 * _interp30 * _interp34 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), fma(complex<double>(0.,0.01041666666666667), _cse1_k3 * _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr69(fenv) * _interp24 * _interp29 * _interp35 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), fma(complex<double>(0.,0.005208333333333332), _cse1_k3 * _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr71(fenv) * _interp24 * _interp34 * _interp35 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), fma(0.02083333333333333, _cse1_k3 * _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr64(fenv) * _interp24 * _interp30 * _interp36 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), fma(0.01041666666666667, _cse1_k3 * _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr70(fenv) * _interp24 * _interp35 * _interp36 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (_cse2_k3 + _interp16) * powr<6>(k))), 0.))))))))))))));
    // clang-format on
  }
  { // subkernel 4
    const auto _cse1_k4 = powr<-1>(l1);
    const auto _cse2_k4 = powr<-2>(p);
    // clang-format off
    _acc += fma(0.01041666666666667, _cse2_k4 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr27(fenv) * _interp18 * _interp23 * _interp38 * (-_cse1_k4 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k4 * _interp39 * _interp6), fma(0.02083333333333333, _cse2_k4 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr9(fenv) * _interp2 * _interp23 * _interp38 * (-_cse1_k4 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k4 * _interp39 * _interp6), fma(complex<double>(0.,0.02083333333333333), _cse2_k4 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr13(fenv) * _interp21 * _interp23 * _interp38 * (-_cse1_k4 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k4 * _interp39 * _interp6), fma(0.005208333333333332, _cse2_k4 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr33(fenv) * _interp18 * _interp27 * _interp38 * (-_cse1_k4 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k4 * _interp39 * _interp6), fma(0.01041666666666667, _cse2_k4 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr21(fenv) * _interp2 * _interp27 * _interp38 * (-_cse1_k4 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k4 * _interp39 * _interp6), fma(complex<double>(0.,0.01041666666666667), _cse2_k4 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr23(fenv) * _interp21 * _interp27 * _interp38 * (-_cse1_k4 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k4 * _interp39 * _interp6), fma(complex<double>(0.,0.01041666666666667), _cse2_k4 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr29(fenv) * _interp18 * _interp28 * _interp38 * (-_cse1_k4 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k4 * _interp39 * _interp6), fma(complex<double>(0.,0.02083333333333333), _cse2_k4 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr11(fenv) * _interp2 * _interp28 * _interp38 * (-_cse1_k4 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k4 * _interp39 * _interp6), fma(-0.02083333333333333, _cse2_k4 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr15(fenv) * _interp21 * _interp28 * _interp38 * (-_cse1_k4 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k4 * _interp39 * _interp6), fma(-0.02083333333333333, _cse2_k4 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr36(fenv) * _interp29 * _interp30 * _interp38 * (-_cse1_k4 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k4 * _interp39 * _interp6), fma(-0.01041666666666667, _cse2_k4 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr48(fenv) * _interp30 * _interp34 * _interp38 * (-_cse1_k4 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k4 * _interp39 * _interp6), fma(-0.02083333333333333, _cse2_k4 * _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr65(fenv) * _interp24 * _interp29 * _interp37 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(-0.01041666666666667, _cse2_k4 * _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr68(fenv) * _interp24 * _interp34 * _interp37 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), fma(complex<double>(0.,0.02083333333333333), _cse2_k4 * _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr66(fenv) * _interp24 * _interp36 * _interp37 * (_interp12 * _interp14 + _interp13 * (_interp15 + 50. * _den1 * (-_interp12 + _interp16) * powr<6>(k))), 0.))))))))))))));
    // clang-format on
  }
  { // subkernel 5
    const auto _cse1_k5 = powr<-1>(l1);
    const auto _cse2_k5 = powr<-2>(p);
    // clang-format off
    _acc += fma(-0.01041666666666667, _cse2_k5 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr54(fenv) * _interp29 * _interp35 * _interp38 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), fma(-0.005208333333333332, _cse2_k5 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr60(fenv) * _interp34 * _interp35 * _interp38 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), fma(complex<double>(0.,-0.02083333333333333), _cse2_k5 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr38(fenv) * _interp30 * _interp36 * _interp38 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), fma(complex<double>(0.,-0.01041666666666667), _cse2_k5 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr56(fenv) * _interp35 * _interp36 * _interp38 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), fma(complex<double>(0.,-0.02083333333333333), _cse2_k5 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr40(fenv) * _interp29 * _interp37 * _interp38 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), fma(complex<double>(0.,-0.01041666666666667), _cse2_k5 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr50(fenv) * _interp34 * _interp37 * _interp38 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), fma(0.02083333333333333, _cse2_k5 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr42(fenv) * _interp36 * _interp37 * _interp38 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), fma(0.005208333333333332, _cse2_k5 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr31(fenv) * _interp18 * _interp23 * _interp50 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), fma(0.01041666666666667, _cse2_k5 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr17(fenv) * _interp2 * _interp23 * _interp50 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), fma(complex<double>(0.,0.01041666666666667), _cse2_k5 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr19(fenv) * _interp21 * _interp23 * _interp50 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), fma(0.002604166666666666, _cse2_k5 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr35(fenv) * _interp18 * _interp27 * _interp50 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), fma(0.005208333333333332, _cse2_k5 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr25(fenv) * _interp2 * _interp27 * _interp50 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), fma(complex<double>(0.,0.005208333333333332), _cse2_k5 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr26(fenv) * _interp21 * _interp27 * _interp50 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), fma(complex<double>(0.,0.005208333333333332), _cse2_k5 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr32(fenv) * _interp18 * _interp28 * _interp50 * (-_cse1_k5 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k5 * _interp39 * _interp6), 0.))))))))))))));
    // clang-format on
  }
  { // subkernel 6
    const auto _cse1_k6 = powr<-1>(l1);
    const auto _cse2_k6 = powr<-2>(p);
    // clang-format off
    _acc += fma(complex<double>(0.,0.01041666666666667), _cse2_k6 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr18(fenv) * _interp2 * _interp28 * _interp50 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), fma(-0.01041666666666667, _cse2_k6 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr20(fenv) * _interp21 * _interp28 * _interp50 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), fma(-0.01041666666666667, _cse2_k6 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr44(fenv) * _interp29 * _interp30 * _interp50 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), fma(-0.005208333333333332, _cse2_k6 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr52(fenv) * _interp30 * _interp34 * _interp50 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), fma(-0.005208333333333332, _cse2_k6 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr58(fenv) * _interp29 * _interp35 * _interp50 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), fma(-0.002604166666666666, _cse2_k6 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr62(fenv) * _interp34 * _interp35 * _interp50 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), fma(complex<double>(0.,-0.01041666666666667), _cse2_k6 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr45(fenv) * _interp30 * _interp36 * _interp50 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), fma(complex<double>(0.,-0.005208333333333332), _cse2_k6 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr59(fenv) * _interp35 * _interp36 * _interp50 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), fma(complex<double>(0.,-0.01041666666666667), _cse2_k6 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr46(fenv) * _interp29 * _interp37 * _interp50 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), fma(complex<double>(0.,-0.005208333333333332), _cse2_k6 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr53(fenv) * _interp34 * _interp37 * _interp50 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), fma(0.01041666666666667, _cse2_k6 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr47(fenv) * _interp36 * _interp37 * _interp50 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), fma(complex<double>(0.,-0.01041666666666667), _cse2_k6 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr28(fenv) * _interp18 * _interp23 * _interp51 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), fma(complex<double>(0.,-0.02083333333333333), _cse2_k6 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr10(fenv) * _interp2 * _interp23 * _interp51 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), fma(0.02083333333333333, _cse2_k6 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr14(fenv) * _interp21 * _interp23 * _interp51 * (-_cse1_k6 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k6 * _interp39 * _interp6), 0.))))))))))))));
    // clang-format on
  }
  { // subkernel 7
    const auto _cse1_k7 = powr<-1>(l1);
    const auto _cse2_k7 = powr<-2>(p);
    // clang-format off
    _acc += fma(complex<double>(0.,-0.005208333333333332), _cse2_k7 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr34(fenv) * _interp18 * _interp27 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), fma(complex<double>(0.,-0.01041666666666667), _cse2_k7 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr22(fenv) * _interp2 * _interp27 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), fma(0.01041666666666667, _cse2_k7 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr24(fenv) * _interp21 * _interp27 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), fma(0.01041666666666667, _cse2_k7 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr30(fenv) * _interp18 * _interp28 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), fma(0.02083333333333333, _cse2_k7 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr12(fenv) * _interp2 * _interp28 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), fma(complex<double>(0.,0.02083333333333333), _cse2_k7 * _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr16(fenv) * _interp21 * _interp28 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), fma(complex<double>(0.,0.02083333333333333), _cse2_k7 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr37(fenv) * _interp29 * _interp30 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), fma(complex<double>(0.,0.01041666666666667), _cse2_k7 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr49(fenv) * _interp30 * _interp34 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), fma(complex<double>(0.,0.01041666666666667), _cse2_k7 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr55(fenv) * _interp29 * _interp35 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), fma(complex<double>(0.,0.005208333333333332), _cse2_k7 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr61(fenv) * _interp34 * _interp35 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), fma(-0.02083333333333333, _cse2_k7 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr39(fenv) * _interp30 * _interp36 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), fma(-0.01041666666666667, _cse2_k7 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr57(fenv) * _interp35 * _interp36 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), fma(-0.02083333333333333, _cse2_k7 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr41(fenv) * _interp29 * _interp37 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), fma(-0.01041666666666667, _cse2_k7 * _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr51(fenv) * _interp34 * _interp37 * _interp51 * (-_cse1_k7 * _interp40 * (_interp41 + 50. * (_interp42 - _interp6)) - _cse1_k7 * _interp39 * _interp6), 0.))))))))))))));
    // clang-format on
  }
  { // subkernel 8
    const auto _den4 = powr<-1>(_interp12 * _interp54 + _interp53 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
    // clang-format off
    _acc += fma(complex<double>(0.,-0.02083333333333333), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr43(fenv) * _interp36 * _interp37 * _interp51 * (-_interp40 * (_interp41 + 50. * (_interp42 - _interp6)) * powr<-1>(l1) - _interp39 * _interp6 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr6(fenv) * _interp1 * _interp18 * _interp52 * (-_interp40 * (_interp41 + 50. * (_interp42 - _interp6)) * powr<-1>(l1) - _interp39 * _interp6 * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr8(fenv) * _interp17 * _interp18 * _interp52 * (-_interp40 * (_interp41 + 50. * (_interp42 - _interp6)) * powr<-1>(l1) - _interp39 * _interp6 * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr0(fenv) * _interp1 * _interp2 * _interp52 * (-_interp40 * (_interp41 + 50. * (_interp42 - _interp6)) * powr<-1>(l1) - _interp39 * _interp6 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr4(fenv) * _interp17 * _interp2 * _interp52 * (-_interp40 * (_interp41 + 50. * (_interp42 - _interp6)) * powr<-1>(l1) - _interp39 * _interp6 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr7(fenv) * _interp18 * _interp20 * _interp52 * (-_interp40 * (_interp41 + 50. * (_interp42 - _interp6)) * powr<-1>(l1) - _interp39 * _interp6 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.02083333333333333), _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr1(fenv) * _interp2 * _interp20 * _interp52 * (-_interp40 * (_interp41 + 50. * (_interp42 - _interp6)) * powr<-1>(l1) - _interp39 * _interp6 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.02083333333333333), _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr2(fenv) * _interp1 * _interp21 * _interp52 * (-_interp40 * (_interp41 + 50. * (_interp42 - _interp6)) * powr<-1>(l1) - _interp39 * _interp6 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr5(fenv) * _interp17 * _interp21 * _interp52 * (-_interp40 * (_interp41 + 50. * (_interp42 - _interp6)) * powr<-1>(l1) - _interp39 * _interp6 * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr3(fenv) * _interp20 * _interp21 * _interp52 * (-_interp40 * (_interp41 + 50. * (_interp42 - _interp6)) * powr<-1>(l1) - _interp39 * _interp6 * powr<-1>(l1)) * powr<-2>(p), 0.))))))))));
    // clang-format on
  }
  return _acc;
}
static inline auto probe_proj(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k)
{
  const double cosl1p1 = cos1;
  const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  double fenv[(DiFfRG::zaqbq4_147_num::nenv) > 0 ? (DiFfRG::zaqbq4_147_num::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * l1);
  const double dr_1 = powr<-1>(l1);
  const double dr_2 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_3 = ntStub(85453., 1. * k);
  const double dr_4 = ntStub(85453., 1. * l1);
  const double dr_5 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_7 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_8 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_9 = ntStub(30415., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_10 = sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_11 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_12 = ntStub(85453., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_13 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_14 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_15 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_16 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_17 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_18 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_19 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_20 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  DiFfRG::zaqbq4_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20);
  const auto _interp2 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp3 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp4 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp5 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp6 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _interp7 = ntStub(85453., 1. * k);
  const auto _interp8 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp9 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp10 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
  const auto _interp11 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp12 = ntStub(29441., 1. * l1);
  const auto _interp13 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp14 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp15 = ntStub(13492., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp16 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp17 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp19 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp24 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp25 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp27 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp31 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp37 = ntStub(12206., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp38 = ntStub(62323., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp39 = ntStub(29441., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp40 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _interp42 = ntStub(14883., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp44 = ntStub(38926., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp47 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp48 = ntStub(38926., 0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp53 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp55 = ntStub(12206., 0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp56 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp58 = ntStub(14883., 0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp60 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp61 = ntStub(27191., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp62 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp63 = ntStub(67246.00000000001, 1. * k);
  const auto _interp64 = ntStub(85453., 1.02 * k);
  const auto _interp65 = ntStub(30415., 1. * l1);
  const auto _interp66 = ntStub(85453., 1. * l1);
  const auto _interp67 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp68 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _interp69 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp70 = ntStub(29441., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp71 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
  const auto _interp73 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp96 = ntStub(29441., 1. * sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp97 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(fma(_interp13, _interp14, fma(_interp12, powr<2>(l1), 0.)));
  const auto _den3 = powr<-1>(fma(_interp13, _interp40, fma(_interp39, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _den5 = powr<-1>(fma(_interp13, _interp97, fma(_interp96, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _den6 = powr<-2>(fma(powr<-2>(l1) * powr<2>(_interp62 * _interp7 + _interp66 * l1), powr<2>(l1), powr<2>(_interp65)));
  const auto _den7 = powr<-1>(fma(_interp13, _interp71, fma(_interp70, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
  const auto _den8 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), powr<2>(_interp69 + _interp68 * _interp7 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))), powr<2>(_interp67)));
  const auto _den9 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), powr<2>(_interp8 + _interp6 * _interp7 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)))), powr<2>(_interp5)));
  const auto _den11 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), powr<2>(_interp11 + _interp10 * _interp7 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))), powr<2>(_interp9)));
  const auto _interp20 = ntRe(DiFfRG::zaqbq4_147_num::tr94(fenv));
  const auto _interp21 = ntStub(14883., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp22 = ntRe(DiFfRG::zaqbq4_147_num::tr106(fenv));
  const auto _interp23 = ntRe(DiFfRG::zaqbq4_147_num::tr85(fenv));
  const auto _interp26 = ntRe(DiFfRG::zaqbq4_147_num::tr91(fenv));
  const auto _interp28 = ntRe(DiFfRG::zaqbq4_147_num::tr95(fenv));
  const auto _interp29 = ntRe(DiFfRG::zaqbq4_147_num::tr98(fenv));
  const auto _interp30 = ntRe(DiFfRG::zaqbq4_147_num::tr83(fenv));
  const auto _interp32 = ntRe(DiFfRG::zaqbq4_147_num::tr90(fenv));
  const auto _interp33 = ntRe(DiFfRG::zaqbq4_147_num::tr101(fenv));
  const auto _interp34 = ntRe(DiFfRG::zaqbq4_147_num::tr104(fenv));
  const auto _interp36 = ntRe(DiFfRG::zaqbq4_147_num::tr74(fenv));
  const auto _interp1 = ntRe(DiFfRG::zaqbq4_147_num::tr82(fenv));
  const auto _interp18 = ntRe(DiFfRG::zaqbq4_147_num::tr100(fenv));
  const auto _interp35 = ntRe(DiFfRG::zaqbq4_147_num::tr88(fenv));
  const auto _interp41 = ntRe(DiFfRG::zaqbq4_147_num::tr77(fenv));
  const auto _interp43 = ntRe(DiFfRG::zaqbq4_147_num::tr73(fenv));
  const auto _interp45 = ntRe(DiFfRG::zaqbq4_147_num::tr79(fenv));
  const auto _interp46 = ntRe(DiFfRG::zaqbq4_147_num::tr64(fenv));
  const auto _interp49 = ntStub(30415., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp50 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp51 = ntStub(85453., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp52 = ntRe(DiFfRG::zaqbq4_147_num::tr70(fenv));
  const auto _interp54 = ntRe(DiFfRG::zaqbq4_147_num::tr65(fenv));
  const auto _interp57 = ntRe(DiFfRG::zaqbq4_147_num::tr68(fenv));
  const auto _interp101 = ntRe(DiFfRG::zaqbq4_147_num::tr54(fenv));
  const auto _den10 = powr<-1>(powr<2>(_interp49) + powr<2>(_interp50 * _interp7 + _interp51 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
  const auto _interp59 = ntRe(DiFfRG::zaqbq4_147_num::tr9(fenv));
  const auto _interp72 = ntRe(DiFfRG::zaqbq4_147_num::tr17(fenv));
  const auto _interp74 = ntRe(DiFfRG::zaqbq4_147_num::tr27(fenv));
  const auto _interp75 = ntRe(DiFfRG::zaqbq4_147_num::tr31(fenv));
  const auto _interp76 = ntRe(DiFfRG::zaqbq4_147_num::tr21(fenv));
  const auto _interp100 = ntRe(DiFfRG::zaqbq4_147_num::tr52(fenv));
  const auto _interp102 = ntRe(DiFfRG::zaqbq4_147_num::tr58(fenv));
  const auto _interp103 = ntRe(DiFfRG::zaqbq4_147_num::tr60(fenv));
  const auto _interp104 = ntRe(DiFfRG::zaqbq4_147_num::tr62(fenv));
  const auto _interp109 = ntRe(DiFfRG::zaqbq4_147_num::tr42(fenv));
  const auto _interp110 = ntRe(DiFfRG::zaqbq4_147_num::tr47(fenv));
  const auto _interp77 = ntRe(DiFfRG::zaqbq4_147_num::tr25(fenv));
  const auto _interp78 = ntRe(DiFfRG::zaqbq4_147_num::tr33(fenv));
  const auto _interp79 = ntRe(DiFfRG::zaqbq4_147_num::tr35(fenv));
  const auto _interp80 = ntRe(DiFfRG::zaqbq4_147_num::tr14(fenv));
  const auto _interp81 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp82 = ntRe(DiFfRG::zaqbq4_147_num::tr24(fenv));
  const auto _interp83 = ntRe(DiFfRG::zaqbq4_147_num::tr12(fenv));
  const auto _interp84 = ntRe(DiFfRG::zaqbq4_147_num::tr30(fenv));
  const auto _interp105 = ntRe(DiFfRG::zaqbq4_147_num::tr39(fenv));
  const auto _interp106 = ntRe(DiFfRG::zaqbq4_147_num::tr57(fenv));
  const auto _interp107 = ntRe(DiFfRG::zaqbq4_147_num::tr41(fenv));
  const auto _interp108 = ntRe(DiFfRG::zaqbq4_147_num::tr51(fenv));
  const auto _interp85 = ntRe(DiFfRG::zaqbq4_147_num::tr15(fenv));
  const auto _interp86 = ntRe(DiFfRG::zaqbq4_147_num::tr20(fenv));
  const auto _interp87 = ntRe(DiFfRG::zaqbq4_147_num::tr0(fenv));
  const auto _interp88 = ntStub(62323., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp89 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp90 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp91 = ntRe(DiFfRG::zaqbq4_147_num::tr4(fenv));
  const auto _interp92 = ntRe(DiFfRG::zaqbq4_147_num::tr6(fenv));
  const auto _interp93 = ntRe(DiFfRG::zaqbq4_147_num::tr8(fenv));
  const auto _interp94 = ntRe(DiFfRG::zaqbq4_147_num::tr3(fenv));
  const auto _interp95 = ntRe(DiFfRG::zaqbq4_147_num::tr36(fenv));
  const auto _interp98 = ntRe(DiFfRG::zaqbq4_147_num::tr44(fenv));
  const auto _interp99 = ntRe(DiFfRG::zaqbq4_147_num::tr48(fenv));
  // clang-format off
  using _T = decltype(_den1 + _den10 + _den11 + _den2 + _den3 + _den5 + _den6 + _den7 + _den8 + _den9 + _interp1 + _interp10 + _interp100 + _interp101 + _interp102 + _interp103 + _interp104 + _interp105 + _interp106 + _interp107 + _interp108 + _interp109 + _interp11 + _interp110 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + _interp98 + _interp99 + cosl1p2 + k + l1 + p);
  // clang-format on
  _T _acc{};
  { // subkernel 1
  // clang-format off
  _acc += fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp21 * _interp24 * _interp25 * _interp28 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den11 * _den2 * _den9 * _interp21 * _interp25 * _interp27 * _interp29 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp23 * _interp24 * _interp25 * _interp3 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp25 * _interp26 * _interp27 * _interp3 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp2 * _interp24 * _interp30 * _interp31 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp2 * _interp27 * _interp31 * _interp32 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp19 * _interp24 * _interp31 * _interp33 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den11 * _den2 * _den9 * _interp19 * _interp27 * _interp31 * _interp34 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den11 * _den2 * _den3 * _interp25 * _interp36 * _interp37 * _interp38 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp2 * _interp20 * _interp21 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.005208333333333332, _den11 * _den2 * _den9 * _interp19 * _interp21 * _interp22 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), 0.)))))))))));
  // clang-format on
  }
  { // subkernel 2
  // clang-format off
  _acc += fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp18 * _interp19 * _interp3 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den11 * _den2 * _den9 * _interp1 * _interp2 * _interp3 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp25 * _interp31 * _interp35 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den3 * _interp25 * _interp38 * _interp41 * _interp42 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den3 * _interp2 * _interp38 * _interp43 * _interp44 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den3 * _interp19 * _interp38 * _interp44 * _interp45 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den10 * _den2 * _den3 * _interp38 * _interp46 * _interp47 * _interp48 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den10 * _den2 * _den3 * _interp38 * _interp48 * _interp52 * _interp53 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den10 * _den2 * _den3 * _interp38 * _interp54 * _interp55 * _interp56 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den10 * _den2 * _den3 * _interp38 * _interp56 * _interp57 * _interp58 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp101 * _interp53 * _interp55 * _interp60 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), 0.)))))))))));
  // clang-format on
  }
  { // subkernel 3
  // clang-format off
  _acc += fma(0.02083333333333333, _den5 * _den6 * _den8 * _interp109 * _interp48 * _interp56 * _interp60 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp103 * _interp53 * _interp58 * _interp60 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp2 * _interp37 * _interp59 * _interp60 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp102 * _interp53 * _interp55 * _interp73 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den5 * _den6 * _den8 * _interp110 * _interp48 * _interp56 * _interp73 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp100 * _interp47 * _interp58 * _interp73 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.002604166666666666, _den5 * _den6 * _den8 * _interp104 * _interp53 * _interp58 * _interp73 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp2 * _interp37 * _interp72 * _interp73 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp19 * _interp37 * _interp60 * _interp74 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp19 * _interp37 * _interp73 * _interp75 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp2 * _interp42 * _interp60 * _interp76 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), 0.)))))))))));
  // clang-format on
  }
  { // subkernel 4
  // clang-format off
  _acc += fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp2 * _interp42 * _interp73 * _interp77 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp19 * _interp42 * _interp60 * _interp78 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(0.002604166666666666, _den6 * _den7 * _den8 * _interp19 * _interp42 * _interp73 * _interp79 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp105 * _interp47 * _interp48 * _interp81 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp106 * _interp48 * _interp53 * _interp81 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp107 * _interp55 * _interp56 * _interp81 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp108 * _interp56 * _interp58 * _interp81 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp25 * _interp37 * _interp80 * _interp81 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp25 * _interp42 * _interp81 * _interp82 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp2 * _interp44 * _interp81 * _interp83 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp19 * _interp44 * _interp81 * _interp84 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), 0.)))))))))));
  // clang-format on
  }
  { // subkernel 5
    const auto _den4 = powr<-1>(_interp13 * _interp90 + _interp89 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
    // clang-format off
    _acc += fma(-0.02083333333333333, _den6 * _den7 * _den8 * _interp25 * _interp44 * _interp60 * _interp85 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den6 * _den7 * _den8 * _interp25 * _interp44 * _interp73 * _interp86 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den4 * _den6 * _den7 * _interp2 * _interp24 * _interp87 * _interp88 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den4 * _den6 * _den7 * _interp2 * _interp27 * _interp88 * _interp91 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den4 * _den6 * _den7 * _interp19 * _interp24 * _interp88 * _interp92 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den4 * _den6 * _den7 * _interp19 * _interp27 * _interp88 * _interp93 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den4 * _den6 * _den7 * _interp25 * _interp4 * _interp88 * _interp94 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp47 * _interp55 * _interp60 * _interp95 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp47 * _interp55 * _interp73 * _interp98 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp47 * _interp58 * _interp60 * _interp99 * (-_interp62 * (_interp63 + 50. * (_interp64 - _interp7)) * powr<-1>(l1) - _interp61 * _interp7 * powr<-1>(l1)) * powr<-2>(p), 0.))))))))));
    // clang-format on
  }
  return _acc;
}
static inline auto probe_repart(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k)
{
  const double cosl1p1 = cos1;
  const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  double fenv[(DiFfRG::zaqbq4_147_num::nenv) > 0 ? (DiFfRG::zaqbq4_147_num::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * l1);
  const double dr_1 = powr<-1>(l1);
  const double dr_2 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_3 = ntStub(85453., 1. * k);
  const double dr_4 = ntStub(85453., 1. * l1);
  const double dr_5 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_7 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_8 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_9 = ntStub(30415., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_10 = sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_11 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_12 = ntStub(85453., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_13 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_14 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_15 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_16 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_17 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_18 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_19 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_20 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  DiFfRG::zaqbq4_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20);
  const auto _interp2 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp3 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp4 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp5 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp6 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _interp7 = ntStub(85453., 1. * k);
  const auto _interp8 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp9 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp10 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
  const auto _interp11 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp12 = ntStub(29441., 1. * l1);
  const auto _interp13 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp14 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp15 = ntStub(13492., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp16 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp17 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp19 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp21 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp29 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp34 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp51 = ntStub(12206., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp52 = ntStub(62323., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp53 = ntStub(29441., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp54 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _interp57 = ntStub(14883., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp62 = ntStub(38926., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp66 = ntStub(12206., 0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp67 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp68 = ntStub(30415., 1. * sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp69 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _interp70 = ntStub(85453., 1. * sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp72 = ntStub(14883., 0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp74 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp77 = ntStub(38926., 0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp80 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
  const auto _interp84 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp85 = ntStub(27191., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp86 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp87 = ntStub(67246.00000000001, 1. * k);
  const auto _interp88 = ntStub(85453., 1.02 * k);
  const auto _interp89 = ntStub(30415., 1. * l1);
  const auto _interp90 = ntStub(85453., 1. * l1);
  const auto _interp91 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp92 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _interp93 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp94 = ntStub(29441., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp95 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
  const auto _interp97 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp105 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp137 = ntStub(29441., 1. * sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp138 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(fma(_interp13, _interp14, fma(_interp12, powr<2>(l1), 0.)));
  const auto _den3 = powr<-1>(fma(_interp13, _interp54, fma(_interp53, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _den5 = powr<-1>(fma(_interp13, _interp138, fma(_interp137, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _den6 = powr<-2>(fma(powr<-2>(l1) * powr<2>(_interp7 * _interp86 + _interp90 * l1), powr<2>(l1), powr<2>(_interp89)));
  const auto _den7 = powr<-1>(fma(_interp13, _interp95, fma(_interp94, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
  const auto _den8 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), powr<2>(_interp93 + _interp7 * _interp92 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))), powr<2>(_interp91)));
  const auto _den9 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), powr<2>(_interp8 + _interp6 * _interp7 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)))), powr<2>(_interp5)));
  const auto _den10 = powr<-1>(fma(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), powr<2>(_interp70 + _interp69 * _interp7 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)))), powr<2>(_interp68)));
  const auto _den11 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), powr<2>(_interp11 + _interp10 * _interp7 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))), powr<2>(_interp9)));
  const auto _interp1 = ntIm(DiFfRG::zaqbq4_147_num::tr81(fenv));
  const auto _interp18 = ntIm(DiFfRG::zaqbq4_147_num::tr89(fenv));
  const auto _interp20 = ntIm(DiFfRG::zaqbq4_147_num::tr99(fenv));
  const auto _interp22 = ntIm(DiFfRG::zaqbq4_147_num::tr103(fenv));
  const auto _interp23 = ntIm(DiFfRG::zaqbq4_147_num::tr93(fenv));
  const auto _interp24 = ntStub(14883., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp25 = ntIm(DiFfRG::zaqbq4_147_num::tr97(fenv));
  const auto _interp26 = ntIm(DiFfRG::zaqbq4_147_num::tr105(fenv));
  const auto _interp27 = ntIm(DiFfRG::zaqbq4_147_num::tr107(fenv));
  const auto _interp28 = ntRe(DiFfRG::zaqbq4_147_num::tr82(fenv));
  const auto _interp31 = ntRe(DiFfRG::zaqbq4_147_num::tr94(fenv));
  const auto _interp32 = ntRe(DiFfRG::zaqbq4_147_num::tr106(fenv));
  const auto _interp36 = ntRe(DiFfRG::zaqbq4_147_num::tr95(fenv));
  const auto _interp37 = ntRe(DiFfRG::zaqbq4_147_num::tr98(fenv));
  const auto _interp39 = ntIm(DiFfRG::zaqbq4_147_num::tr96(fenv));
  const auto _interp30 = ntRe(DiFfRG::zaqbq4_147_num::tr100(fenv));
  const auto _interp33 = ntRe(DiFfRG::zaqbq4_147_num::tr85(fenv));
  const auto _interp35 = ntRe(DiFfRG::zaqbq4_147_num::tr91(fenv));
  const auto _interp38 = ntIm(DiFfRG::zaqbq4_147_num::tr86(fenv));
  const auto _interp40 = ntRe(DiFfRG::zaqbq4_147_num::tr83(fenv));
  const auto _interp41 = ntStub(38926., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp42 = ntRe(DiFfRG::zaqbq4_147_num::tr90(fenv));
  const auto _interp43 = ntRe(DiFfRG::zaqbq4_147_num::tr101(fenv));
  const auto _interp44 = ntRe(DiFfRG::zaqbq4_147_num::tr104(fenv));
  const auto _interp45 = ntIm(DiFfRG::zaqbq4_147_num::tr84(fenv));
  const auto _interp46 = ntIm(DiFfRG::zaqbq4_147_num::tr102(fenv));
  const auto _interp47 = ntIm(DiFfRG::zaqbq4_147_num::tr87(fenv));
  const auto _interp48 = ntIm(DiFfRG::zaqbq4_147_num::tr92(fenv));
  const auto _interp49 = ntRe(DiFfRG::zaqbq4_147_num::tr88(fenv));
  const auto _interp50 = ntIm(DiFfRG::zaqbq4_147_num::tr72(fenv));
  const auto _interp55 = ntIm(DiFfRG::zaqbq4_147_num::tr78(fenv));
  const auto _interp56 = ntIm(DiFfRG::zaqbq4_147_num::tr76(fenv));
  const auto _interp58 = ntIm(DiFfRG::zaqbq4_147_num::tr80(fenv));
  const auto _interp59 = ntRe(DiFfRG::zaqbq4_147_num::tr74(fenv));
  const auto _interp60 = ntRe(DiFfRG::zaqbq4_147_num::tr77(fenv));
  const auto _interp61 = ntRe(DiFfRG::zaqbq4_147_num::tr73(fenv));
  const auto _interp63 = ntRe(DiFfRG::zaqbq4_147_num::tr79(fenv));
  const auto _interp64 = ntIm(DiFfRG::zaqbq4_147_num::tr75(fenv));
  const auto _interp65 = ntIm(DiFfRG::zaqbq4_147_num::tr63(fenv));
  const auto _interp71 = ntIm(DiFfRG::zaqbq4_147_num::tr67(fenv));
  const auto _interp73 = ntIm(DiFfRG::zaqbq4_147_num::tr69(fenv));
  const auto _interp75 = ntIm(DiFfRG::zaqbq4_147_num::tr71(fenv));
  const auto _interp76 = ntRe(DiFfRG::zaqbq4_147_num::tr64(fenv));
  const auto _interp78 = ntRe(DiFfRG::zaqbq4_147_num::tr70(fenv));
  const auto _interp79 = ntRe(DiFfRG::zaqbq4_147_num::tr65(fenv));
  const auto _interp81 = ntRe(DiFfRG::zaqbq4_147_num::tr68(fenv));
  const auto _interp82 = ntIm(DiFfRG::zaqbq4_147_num::tr66(fenv));
  const auto _interp104 = ntIm(DiFfRG::zaqbq4_147_num::tr10(fenv));
  const auto _interp106 = ntIm(DiFfRG::zaqbq4_147_num::tr28(fenv));
  const auto _interp124 = ntRe(DiFfRG::zaqbq4_147_num::tr0(fenv));
  const auto _interp125 = ntStub(62323., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp126 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp127 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp128 = ntRe(DiFfRG::zaqbq4_147_num::tr4(fenv));
  const auto _interp129 = ntRe(DiFfRG::zaqbq4_147_num::tr6(fenv));
  const auto _interp130 = ntRe(DiFfRG::zaqbq4_147_num::tr8(fenv));
  const auto _interp131 = ntIm(DiFfRG::zaqbq4_147_num::tr1(fenv));
  const auto _interp132 = ntIm(DiFfRG::zaqbq4_147_num::tr7(fenv));
  const auto _interp133 = ntIm(DiFfRG::zaqbq4_147_num::tr2(fenv));
  const auto _interp134 = ntIm(DiFfRG::zaqbq4_147_num::tr5(fenv));
  const auto _interp135 = ntRe(DiFfRG::zaqbq4_147_num::tr3(fenv));
  const auto _interp107 = ntIm(DiFfRG::zaqbq4_147_num::tr22(fenv));
  const auto _interp108 = ntIm(DiFfRG::zaqbq4_147_num::tr34(fenv));
  const auto _interp113 = ntRe(DiFfRG::zaqbq4_147_num::tr14(fenv));
  const auto _interp114 = ntRe(DiFfRG::zaqbq4_147_num::tr24(fenv));
  const auto _interp119 = ntRe(DiFfRG::zaqbq4_147_num::tr12(fenv));
  const auto _interp120 = ntRe(DiFfRG::zaqbq4_147_num::tr30(fenv));
  const auto _interp123 = ntIm(DiFfRG::zaqbq4_147_num::tr16(fenv));
  const auto _interp146 = ntIm(DiFfRG::zaqbq4_147_num::tr37(fenv));
  const auto _interp147 = ntIm(DiFfRG::zaqbq4_147_num::tr49(fenv));
  const auto _interp148 = ntIm(DiFfRG::zaqbq4_147_num::tr55(fenv));
  const auto _interp149 = ntIm(DiFfRG::zaqbq4_147_num::tr61(fenv));
  const auto _interp154 = ntRe(DiFfRG::zaqbq4_147_num::tr39(fenv));
  const auto _interp155 = ntRe(DiFfRG::zaqbq4_147_num::tr57(fenv));
  const auto _interp160 = ntRe(DiFfRG::zaqbq4_147_num::tr41(fenv));
  const auto _interp100 = ntRe(DiFfRG::zaqbq4_147_num::tr21(fenv));
  const auto _interp102 = ntRe(DiFfRG::zaqbq4_147_num::tr33(fenv));
  const auto _interp109 = ntIm(DiFfRG::zaqbq4_147_num::tr13(fenv));
  const auto _interp111 = ntIm(DiFfRG::zaqbq4_147_num::tr23(fenv));
  const auto _interp115 = ntIm(DiFfRG::zaqbq4_147_num::tr11(fenv));
  const auto _interp117 = ntIm(DiFfRG::zaqbq4_147_num::tr29(fenv));
  const auto _interp121 = ntRe(DiFfRG::zaqbq4_147_num::tr15(fenv));
  const auto _interp136 = ntRe(DiFfRG::zaqbq4_147_num::tr36(fenv));
  const auto _interp140 = ntRe(DiFfRG::zaqbq4_147_num::tr48(fenv));
  const auto _interp142 = ntRe(DiFfRG::zaqbq4_147_num::tr54(fenv));
  const auto _interp144 = ntRe(DiFfRG::zaqbq4_147_num::tr60(fenv));
  const auto _interp150 = ntIm(DiFfRG::zaqbq4_147_num::tr38(fenv));
  const auto _interp161 = ntRe(DiFfRG::zaqbq4_147_num::tr51(fenv));
  const auto _interp164 = ntIm(DiFfRG::zaqbq4_147_num::tr43(fenv));
  const auto _interp83 = ntRe(DiFfRG::zaqbq4_147_num::tr9(fenv));
  const auto _interp101 = ntRe(DiFfRG::zaqbq4_147_num::tr25(fenv));
  const auto _interp103 = ntRe(DiFfRG::zaqbq4_147_num::tr35(fenv));
  const auto _interp110 = ntIm(DiFfRG::zaqbq4_147_num::tr19(fenv));
  const auto _interp112 = ntIm(DiFfRG::zaqbq4_147_num::tr26(fenv));
  const auto _interp116 = ntIm(DiFfRG::zaqbq4_147_num::tr18(fenv));
  const auto _interp118 = ntIm(DiFfRG::zaqbq4_147_num::tr32(fenv));
  const auto _interp122 = ntRe(DiFfRG::zaqbq4_147_num::tr20(fenv));
  const auto _interp139 = ntRe(DiFfRG::zaqbq4_147_num::tr44(fenv));
  const auto _interp141 = ntRe(DiFfRG::zaqbq4_147_num::tr52(fenv));
  const auto _interp152 = ntIm(DiFfRG::zaqbq4_147_num::tr56(fenv));
  const auto _interp156 = ntIm(DiFfRG::zaqbq4_147_num::tr40(fenv));
  const auto _interp158 = ntIm(DiFfRG::zaqbq4_147_num::tr50(fenv));
  const auto _interp162 = ntRe(DiFfRG::zaqbq4_147_num::tr42(fenv));
  const auto _interp96 = ntRe(DiFfRG::zaqbq4_147_num::tr17(fenv));
  const auto _interp98 = ntRe(DiFfRG::zaqbq4_147_num::tr27(fenv));
  const auto _interp99 = ntRe(DiFfRG::zaqbq4_147_num::tr31(fenv));
  const auto _interp143 = ntRe(DiFfRG::zaqbq4_147_num::tr58(fenv));
  const auto _interp145 = ntRe(DiFfRG::zaqbq4_147_num::tr62(fenv));
  const auto _interp151 = ntIm(DiFfRG::zaqbq4_147_num::tr45(fenv));
  const auto _interp153 = ntIm(DiFfRG::zaqbq4_147_num::tr59(fenv));
  const auto _interp157 = ntIm(DiFfRG::zaqbq4_147_num::tr46(fenv));
  const auto _interp159 = ntIm(DiFfRG::zaqbq4_147_num::tr53(fenv));
  const auto _interp163 = ntRe(DiFfRG::zaqbq4_147_num::tr47(fenv));
  // clang-format off
  using _T = decltype(_den1 + _den10 + _den11 + _den2 + _den3 + _den5 + _den6 + _den7 + _den8 + _den9 + _interp1 + _interp10 + _interp100 + _interp101 + _interp102 + _interp103 + _interp104 + _interp105 + _interp106 + _interp107 + _interp108 + _interp109 + _interp11 + _interp110 + _interp111 + _interp112 + _interp113 + _interp114 + _interp115 + _interp116 + _interp117 + _interp118 + _interp119 + _interp12 + _interp120 + _interp121 + _interp122 + _interp123 + _interp124 + _interp125 + _interp126 + _interp127 + _interp128 + _interp129 + _interp13 + _interp130 + _interp131 + _interp132 + _interp133 + _interp134 + _interp135 + _interp136 + _interp137 + _interp138 + _interp139 + _interp14 + _interp140 + _interp141 + _interp142 + _interp143 + _interp144 + _interp145 + _interp146 + _interp147 + _interp148 + _interp149 + _interp15 + _interp150 + _interp151 + _interp152 + _interp153 + _interp154 + _interp155 + _interp156 + _interp157 + _interp158 + _interp159 + _interp16 + _interp160 + _interp161 + _interp162 + _interp163 + _interp164 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + _interp98 + _interp99 + cosl1p2 + k + l1 + p);
  // clang-format on
  _T _acc{};
  { // subkernel 1
  // clang-format off
  _acc += fma(0.005208333333333332, _den11 * _den2 * _den9 * _interp2 * _interp21 * _interp24 * _interp26 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.002604166666666666, _den11 * _den2 * _den9 * _interp19 * _interp21 * _interp24 * _interp27 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp2 * _interp23 * _interp24 * _interp3 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.005208333333333332, _den11 * _den2 * _den9 * _interp19 * _interp24 * _interp25 * _interp3 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp24 * _interp29 * _interp3 * _interp31 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.005208333333333332, _den11 * _den2 * _den9 * _interp21 * _interp24 * _interp29 * _interp32 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp2 * _interp24 * _interp34 * _interp36 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den11 * _den2 * _den9 * _interp19 * _interp24 * _interp34 * _interp37 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp24 * _interp29 * _interp34 * _interp39 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp2 * _interp20 * _interp21 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.005208333333333332, _den11 * _den2 * _den9 * _interp19 * _interp21 * _interp22 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp18 * _interp19 * _interp3 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den11 * _den2 * _den9 * _interp1 * _interp2 * _interp3 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den11 * _den2 * _den9 * _interp28 * _interp29 * _interp3 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), 0.))))))))))))));
  // clang-format on
  }
  { // subkernel 2
  // clang-format off
  _acc += fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp21 * _interp29 * _interp30 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp2 * _interp33 * _interp34 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp19 * _interp34 * _interp35 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den11 * _den2 * _den9 * _interp29 * _interp34 * _interp38 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp2 * _interp3 * _interp40 * _interp41 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp19 * _interp3 * _interp41 * _interp42 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp2 * _interp21 * _interp41 * _interp43 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den11 * _den2 * _den9 * _interp19 * _interp21 * _interp41 * _interp44 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den11 * _den2 * _den9 * _interp29 * _interp3 * _interp41 * _interp45 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp21 * _interp29 * _interp41 * _interp46 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp2 * _interp34 * _interp41 * _interp47 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp19 * _interp34 * _interp41 * _interp48 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp29 * _interp34 * _interp41 * _interp49 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den3 * _interp3 * _interp50 * _interp51 * _interp52 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), 0.))))))))))))));
  // clang-format on
  }
  { // subkernel 3
  // clang-format off
  _acc += fma(-0.01041666666666667, _den11 * _den2 * _den3 * _interp21 * _interp51 * _interp52 * _interp55 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den3 * _interp3 * _interp52 * _interp56 * _interp57 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den11 * _den2 * _den3 * _interp21 * _interp52 * _interp57 * _interp58 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den11 * _den2 * _den3 * _interp34 * _interp51 * _interp52 * _interp59 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den3 * _interp34 * _interp52 * _interp57 * _interp60 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den3 * _interp3 * _interp52 * _interp61 * _interp62 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den3 * _interp21 * _interp52 * _interp62 * _interp63 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den3 * _interp34 * _interp52 * _interp62 * _interp64 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den10 * _den2 * _den3 * _interp52 * _interp65 * _interp66 * _interp67 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den10 * _den2 * _den3 * _interp52 * _interp67 * _interp71 * _interp72 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den10 * _den2 * _den3 * _interp52 * _interp66 * _interp73 * _interp74 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den10 * _den2 * _den3 * _interp52 * _interp72 * _interp74 * _interp75 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den10 * _den2 * _den3 * _interp52 * _interp67 * _interp76 * _interp77 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den10 * _den2 * _den3 * _interp52 * _interp74 * _interp77 * _interp78 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), 0.))))))))))))));
  // clang-format on
  }
  { // subkernel 4
    const auto _den4 = powr<-1>(_interp127 * _interp13 + _interp126 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
    // clang-format off
    _acc += fma(-0.02083333333333333, _den10 * _den2 * _den3 * _interp52 * _interp66 * _interp79 * _interp80 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den10 * _den2 * _den3 * _interp52 * _interp72 * _interp80 * _interp81 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den10 * _den2 * _den3 * _interp52 * _interp77 * _interp80 * _interp82 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den4 * _den6 * _den7 * _interp125 * _interp130 * _interp19 * _interp21 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den4 * _den6 * _den7 * _interp125 * _interp129 * _interp2 * _interp21 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den4 * _den6 * _den7 * _interp125 * _interp132 * _interp21 * _interp29 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den4 * _den6 * _den7 * _interp125 * _interp128 * _interp19 * _interp3 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den4 * _den6 * _den7 * _interp124 * _interp125 * _interp2 * _interp3 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den4 * _den6 * _den7 * _interp125 * _interp131 * _interp29 * _interp3 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den4 * _den6 * _den7 * _interp125 * _interp134 * _interp19 * _interp34 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den4 * _den6 * _den7 * _interp125 * _interp133 * _interp2 * _interp34 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den4 * _den6 * _den7 * _interp125 * _interp135 * _interp29 * _interp34 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp105 * _interp106 * _interp21 * _interp51 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp104 * _interp105 * _interp3 * _interp51 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), 0.))))))))))))));
    // clang-format on
  }
  { // subkernel 5
  // clang-format off
  _acc += fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp105 * _interp113 * _interp34 * _interp51 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp105 * _interp108 * _interp21 * _interp57 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp105 * _interp107 * _interp3 * _interp57 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp105 * _interp114 * _interp34 * _interp57 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp105 * _interp120 * _interp21 * _interp62 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp105 * _interp119 * _interp3 * _interp62 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den6 * _den7 * _den8 * _interp105 * _interp123 * _interp34 * _interp62 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp105 * _interp146 * _interp66 * _interp67 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp105 * _interp147 * _interp67 * _interp72 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp105 * _interp148 * _interp66 * _interp74 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp105 * _interp149 * _interp72 * _interp74 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp105 * _interp154 * _interp67 * _interp77 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp105 * _interp155 * _interp74 * _interp77 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp105 * _interp160 * _interp66 * _interp80 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), 0.))))))))))))));
  // clang-format on
  }
  { // subkernel 6
  // clang-format off
  _acc += fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp105 * _interp161 * _interp72 * _interp80 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den5 * _den6 * _den8 * _interp105 * _interp164 * _interp77 * _interp80 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den6 * _den7 * _den8 * _interp109 * _interp34 * _interp51 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp102 * _interp21 * _interp57 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp100 * _interp3 * _interp57 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den6 * _den7 * _den8 * _interp111 * _interp34 * _interp57 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den6 * _den7 * _den8 * _interp117 * _interp21 * _interp62 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den6 * _den7 * _den8 * _interp115 * _interp3 * _interp62 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den6 * _den7 * _den8 * _interp121 * _interp34 * _interp62 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp136 * _interp66 * _interp67 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp140 * _interp67 * _interp72 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp142 * _interp66 * _interp74 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp144 * _interp72 * _interp74 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den5 * _den6 * _den8 * _interp150 * _interp67 * _interp77 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), 0.))))))))))))));
  // clang-format on
  }
  { // subkernel 7
  // clang-format off
  _acc += fma(0.01041666666666667, _den5 * _den6 * _den8 * _interp152 * _interp74 * _interp77 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den5 * _den6 * _den8 * _interp156 * _interp66 * _interp80 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den5 * _den6 * _den8 * _interp158 * _interp72 * _interp80 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den5 * _den6 * _den8 * _interp162 * _interp77 * _interp80 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp3 * _interp51 * _interp83 * _interp84 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den6 * _den7 * _den8 * _interp110 * _interp34 * _interp51 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.002604166666666666, _den6 * _den7 * _den8 * _interp103 * _interp21 * _interp57 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp101 * _interp3 * _interp57 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den6 * _den7 * _den8 * _interp112 * _interp34 * _interp57 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den6 * _den7 * _den8 * _interp118 * _interp21 * _interp62 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den6 * _den7 * _den8 * _interp116 * _interp3 * _interp62 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den6 * _den7 * _den8 * _interp122 * _interp34 * _interp62 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp139 * _interp66 * _interp67 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp141 * _interp67 * _interp72 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), 0.))))))))))))));
  // clang-format on
  }
  { // subkernel 8
  // clang-format off
  _acc += fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp143 * _interp66 * _interp74 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.002604166666666666, _den5 * _den6 * _den8 * _interp145 * _interp72 * _interp74 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den5 * _den6 * _den8 * _interp151 * _interp67 * _interp77 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den5 * _den6 * _den8 * _interp153 * _interp74 * _interp77 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den5 * _den6 * _den8 * _interp157 * _interp66 * _interp80 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den5 * _den6 * _den8 * _interp159 * _interp72 * _interp80 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den5 * _den6 * _den8 * _interp163 * _interp77 * _interp80 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp3 * _interp51 * _interp96 * _interp97 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp21 * _interp51 * _interp84 * _interp98 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp21 * _interp51 * _interp97 * _interp99 * (-_interp7 * _interp85 * powr<-1>(l1) - _interp86 * (_interp87 + 50. * (-_interp7 + _interp88)) * powr<-1>(l1)) * powr<-2>(p), 0.))))))))));
  // clang-format on
  }
  return _acc;
}
int main(int argc, char** argv){
  const char* outf=nullptr; const char* macro=nullptr;
  for(int i=1;i<argc;++i){ if(!std::strcmp(argv[i],"-o") && i+1<argc) outf=argv[++i];
                           else if(!std::strcmp(argv[i],"-m") && i+1<argc) macro=argv[++i]; }
  std::mt19937_64 rng(12345); std::uniform_real_distribution<double> U(0.25,3.0),Uc(-0.9,0.9),Uph(0.1,6.2);
  double mim=0,mdiff=0,mre=0,mrim=0,mrdiff=0,mrrep=0; long ok=0, nrep=0;
  for(int n=0;n<4000;++n){ double l1 = U(rng); double cos1 = Uc(rng); double cos2 = Uc(rng); double p = U(rng); double k = U(rng);
    std::complex<double> f = probe_full(l1, cos1, cos2, p, k); double pj = probe_proj(l1, cos1, cos2, p, k);
    double rp = probe_repart(l1, cos1, cos2, p, k);
    double im=std::imag(f), re=std::real(f), df=std::abs(f-pj);
    { double rr = std::fabs(rp-re)/(std::abs(f)+1.0);
      if(std::isfinite(rr)){ mrrep=std::max(mrrep, rr); if(rr > 1.e-9) ++nrep; } }
    if(std::isfinite(im)&&std::isfinite(re)&&std::isfinite(df)){ mim=std::max(mim,std::fabs(im)); mdiff=std::max(mdiff,df); mre=std::max(mre,std::fabs(re));
      mrim=std::max(mrim, std::fabs(im)/(std::fabs(re)+1.0)); mrdiff=std::max(mrdiff, df/(std::abs(f)+1.0)); ++ok; } }
  if(ok < 1){ std::fprintf(stderr, "[probe] no usable points\n"); return 2; }
  const int verdict = (mrim > 1.e-9) ? 0 : ((mrdiff <= 1.e-9) ? 2 : 1);
  std::printf("%.10e %.10e %.10e %.10e %.10e %.10e %ld %ld %d\n", mim, mdiff, mre, mrim, mrdiff, mrrep, nrep, ok, verdict);
  if(nrep * 100 > ok){ std::fprintf(stderr, "[probe] the RePart projection does not reproduce Re(integrand): %ld of %ld points disagree by more than 1.e-9 (worst rel=%.3e).\n"
      "[probe] A few isolated outliers would be catastrophic cancellation; this many is a NumTracer emitter bug (ntRePartIntegrand).\n", nrep, ok, mrrep); return 4; }
  if(outf && macro){ std::FILE* f = std::fopen(outf, "w");
    if(!f){ std::fprintf(stderr, "[probe] cannot write %s\n", outf); return 3; }
    std::fprintf(f, "// GENERATED by the numtrace step — do not edit.\n");
    std::fprintf(f, "// 2 = Pure (imaginary coefficients dropped), 1 = RePart (re/im split), 0 = complex.\n");
    std::fprintf(f, "#pragma once\n#define %s %d\n", macro, verdict);
    std::fclose(f); }
  return 0; }
