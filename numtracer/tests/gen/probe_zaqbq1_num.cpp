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
#include "ZAqbq1_num_kernels.hh"
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
  double fenv[(DiFfRG::zaqbq1_num::nenv) > 0 ? (DiFfRG::zaqbq1_num::nenv) : 1];
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
  DiFfRG::zaqbq1_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20);
  const auto _interp1 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp2 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp3 = ntStub(12206., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp4 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp5 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp6 = ntStub(85453., 1. * k);
  const auto _interp7 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp8 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp9 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp10 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp11 = ntStub(29441., 1. * l1);
  const auto _interp12 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp13 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp14 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp15 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp16 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp17 = ntStub(12206., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (2. * cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp18 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp19 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp20 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp21 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp22 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 1.5 * powr<2>(p)));
  const auto _interp23 = ntStub(30415., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp24 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp25 = ntStub(85453., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp26 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp27 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp28 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp29 = ntStub(67246.00000000001, 1. * k);
  const auto _interp30 = ntStub(85453., 1.02 * k);
  const auto _interp31 = ntStub(30415., 1. * l1);
  const auto _interp32 = ntStub(85453., 1. * l1);
  const auto _interp33 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp34 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp35 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp36 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp37 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp38 = ntStub(62323., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp39 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp40 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp41 = ntStub(29441., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp42 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(_interp12 * _interp13 + _interp11 * powr<2>(l1));
  const auto _den3 = powr<-1>(_interp12 * _interp20 + _interp19 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _den4 = powr<-1>(_interp12 * _interp40 + _interp39 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den5 = powr<-1>(_interp12 * _interp42 + _interp41 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den6 = powr<-2>(powr<2>(_interp31) + powr<2>(_interp28 * _interp6 + _interp32 * l1));
  const auto _den7 = powr<-1>(_interp12 * _interp37 + _interp36 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _den8 = powr<-1>(powr<2>(_interp33) + powr<2>(_interp34 * _interp6 + _interp35 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
  const auto _den9 = powr<-1>(powr<2>(_interp4) + powr<2>(_interp5 * _interp6 + _interp7 * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
  const auto _den10 = powr<-1>(powr<2>(_interp23) + powr<2>(_interp24 * _interp6 + _interp25 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
  const auto _den11 = powr<-1>(powr<2>(_interp8) + powr<2>(_interp10) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp6 * _interp9 * (_interp6 * _interp9 + 2. * _interp10 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
  const auto _cse1 = powr<-1>(l1);
  const auto _cse2 = -_interp6;
  const auto _cse3 = _cse2 + _interp30;
  const auto _cse4 = 50. * _cse3;
  const auto _cse5 = _cse4 + _interp29;
  const auto _cse6 = -_cse1 * _cse5 * _interp28;
  const auto _cse7 = -_cse1 * _interp27 * _interp6;
  const auto _cse8 = _cse6 + _cse7;
  const auto _cse9 = -_interp12;
  const auto _cse10 = _cse9 + _interp16;
  const auto _cse11 = 50. * _cse10 * _den1 * powr<6>(k);
  // clang-format off
  return 0.02083333333333333 * fma(-1., _cse11 * _den11 * _den2 * _den3 * DiFfRG::zaqbq1_num::tr4(fenv) * _interp13 * _interp17 * _interp18 * _interp2, fma(-1., _den11 * _den2 * _den3 * DiFfRG::zaqbq1_num::tr4(fenv) * _interp12 * _interp14 * _interp17 * _interp18 * _interp2, fma(-1., _den11 * _den2 * _den3 * DiFfRG::zaqbq1_num::tr4(fenv) * _interp13 * _interp15 * _interp17 * _interp18 * _interp2, fma(-1., _cse11 * _den10 * _den2 * _den3 * DiFfRG::zaqbq1_num::tr3(fenv) * _interp13 * _interp18 * _interp21 * _interp22, fma(-1., _den10 * _den2 * _den3 * DiFfRG::zaqbq1_num::tr3(fenv) * _interp12 * _interp14 * _interp18 * _interp21 * _interp22, fma(-1., _den10 * _den2 * _den3 * DiFfRG::zaqbq1_num::tr3(fenv) * _interp13 * _interp15 * _interp18 * _interp21 * _interp22, fma(complex<double>(0.,1.), _cse8 * _den6 * _den7 * _den8 * DiFfRG::zaqbq1_num::tr1(fenv) * _interp17 * _interp2 * _interp26, fma(complex<double>(0.,-1.), _cse8 * _den5 * _den6 * _den8 * DiFfRG::zaqbq1_num::tr2(fenv) * _interp21 * _interp22 * _interp26, fma(_cse11, _den11 * _den2 * _den9 * DiFfRG::zaqbq1_num::tr5(fenv) * _interp1 * _interp13 * _interp2 * _interp3, fma(_den11, _den2 * _den9 * DiFfRG::zaqbq1_num::tr5(fenv) * _interp1 * _interp12 * _interp14 * _interp2 * _interp3, fma(_den11, _den2 * _den9 * DiFfRG::zaqbq1_num::tr5(fenv) * _interp1 * _interp13 * _interp15 * _interp2 * _interp3, fma(complex<double>(0.,-1.), _cse8 * _den4 * _den6 * _den7 * DiFfRG::zaqbq1_num::tr0(fenv) * _interp1 * _interp2 * _interp38, 0.))))))))))));
  // clang-format on
}
static inline auto probe_proj(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k)
{
  const double cosl1p1 = cos1;
  const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  double fenv[(DiFfRG::zaqbq1_num::nenv) > 0 ? (DiFfRG::zaqbq1_num::nenv) : 1];
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
  DiFfRG::zaqbq1_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20);
  const auto _interp1 = ntRe(DiFfRG::zaqbq1_num::tr5(fenv));
  const auto _interp2 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp3 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp4 = ntStub(12206., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp5 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp6 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp7 = ntStub(85453., 1. * k);
  const auto _interp8 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp9 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp10 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp11 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp12 = ntStub(29441., 1. * l1);
  const auto _interp13 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp14 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp15 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp16 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp17 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp18 = ntRe(DiFfRG::zaqbq1_num::tr4(fenv));
  const auto _interp19 = ntStub(12206., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (2. * cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp20 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp21 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp22 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp23 = ntRe(DiFfRG::zaqbq1_num::tr3(fenv));
  const auto _interp24 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp25 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 1.5 * powr<2>(p)));
  const auto _interp26 = ntStub(30415., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp27 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp28 = ntStub(85453., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(_interp13 * _interp14 + _interp12 * powr<2>(l1));
  const auto _den3 = powr<-1>(_interp13 * _interp22 + _interp21 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _den4 = powr<-1>(powr<2>(_interp5) + powr<2>(_interp6 * _interp7 + _interp8 * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
  const auto _den5 = powr<-1>(powr<2>(_interp26) + powr<2>(_interp27 * _interp7 + _interp28 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
  const auto _den6 = powr<-1>(powr<2>(_interp9) + powr<2>(_interp11) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp10 * _interp7 * (_interp10 * _interp7 + 2. * _interp11 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
  const auto _cse1 = -_interp13;
  const auto _cse2 = _cse1 + _interp17;
  const auto _cse3 = 50. * _cse2 * _den1 * powr<6>(k);
  const auto _cse4 = _cse3 + _interp16;
  const auto _cse5 = _cse4 * _interp14;
  const auto _cse6 = _interp13 * _interp15;
  const auto _cse7 = _cse5 + _cse6;
  return 0.02083333333333333 * fma(-1., _cse7 * _den2 * _den3 * _den5 * _interp20 * _interp23 * _interp24 * _interp25, fma(-1., _cse7 * _den2 * _den3 * _den6 * _interp18 * _interp19 * _interp20 * _interp3, fma(_cse7, _den2 * _den4 * _den6 * _interp1 * _interp2 * _interp3 * _interp4, 0.)));
}
static inline auto probe_repart(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k)
{
  const double cosl1p1 = cos1;
  const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  double fenv[(DiFfRG::zaqbq1_num::nenv) > 0 ? (DiFfRG::zaqbq1_num::nenv) : 1];
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
  DiFfRG::zaqbq1_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20);
  const auto _interp1 = ntRe(DiFfRG::zaqbq1_num::tr5(fenv));
  const auto _interp2 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp3 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp4 = ntStub(12206., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp5 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp6 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp7 = ntStub(85453., 1. * k);
  const auto _interp8 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp9 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp10 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp11 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp12 = ntStub(29441., 1. * l1);
  const auto _interp13 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp14 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp15 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp16 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp17 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp18 = ntRe(DiFfRG::zaqbq1_num::tr4(fenv));
  const auto _interp19 = ntStub(12206., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (2. * cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp20 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp21 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp22 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp23 = ntRe(DiFfRG::zaqbq1_num::tr3(fenv));
  const auto _interp24 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp25 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 1.5 * powr<2>(p)));
  const auto _interp26 = ntStub(30415., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp27 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp28 = ntStub(85453., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp29 = ntIm(DiFfRG::zaqbq1_num::tr1(fenv));
  const auto _interp30 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp31 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp32 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp33 = ntStub(67246.00000000001, 1. * k);
  const auto _interp34 = ntStub(85453., 1.02 * k);
  const auto _interp35 = ntStub(30415., 1. * l1);
  const auto _interp36 = ntStub(85453., 1. * l1);
  const auto _interp37 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp38 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp39 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp40 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp41 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp42 = ntIm(DiFfRG::zaqbq1_num::tr0(fenv));
  const auto _interp43 = ntStub(62323., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp44 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp45 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp46 = ntIm(DiFfRG::zaqbq1_num::tr2(fenv));
  const auto _interp47 = ntStub(29441., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp48 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(_interp13 * _interp14 + _interp12 * powr<2>(l1));
  const auto _den3 = powr<-1>(_interp13 * _interp22 + _interp21 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _den4 = powr<-1>(_interp13 * _interp45 + _interp44 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den5 = powr<-1>(_interp13 * _interp48 + _interp47 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den6 = powr<-2>(powr<2>(_interp35) + powr<2>(_interp32 * _interp7 + _interp36 * l1));
  const auto _den7 = powr<-1>(_interp13 * _interp41 + _interp40 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _den8 = powr<-1>(powr<2>(_interp37) + powr<2>(_interp38 * _interp7 + _interp39 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
  const auto _den9 = powr<-1>(powr<2>(_interp5) + powr<2>(_interp6 * _interp7 + _interp8 * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
  const auto _den10 = powr<-1>(powr<2>(_interp26) + powr<2>(_interp27 * _interp7 + _interp28 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
  const auto _den11 = powr<-1>(powr<2>(_interp9) + powr<2>(_interp11) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp10 * _interp7 * (_interp10 * _interp7 + 2. * _interp11 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
  const auto _cse1 = powr<-1>(l1);
  const auto _cse2 = -_interp7;
  const auto _cse3 = _cse2 + _interp34;
  const auto _cse4 = 50. * _cse3;
  const auto _cse5 = _cse4 + _interp33;
  // clang-format off
  return fma(-0.02083333333333333, _den6 * _den7 * _den8 * _interp19 * _interp29 * _interp3 * _interp30 * (-_cse1 * _cse5 * _interp32 - _cse1 * _interp31 * _interp7), fma(0.02083333333333333, _den4 * _den6 * _den7 * _interp2 * _interp3 * _interp42 * _interp43 * (-_cse1 * _cse5 * _interp32 - _cse1 * _interp31 * _interp7), fma(0.02083333333333333, _den5 * _den6 * _den8 * _interp24 * _interp25 * _interp30 * _interp46 * (-_cse1 * _cse5 * _interp32 - _cse1 * _interp31 * _interp7), fma(-0.02083333333333333, _den10 * _den2 * _den3 * _interp20 * _interp23 * _interp24 * _interp25 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))), fma(-0.02083333333333333, _den11 * _den2 * _den3 * _interp18 * _interp19 * _interp20 * _interp3 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))), fma(0.02083333333333333, _den11 * _den2 * _den9 * _interp1 * _interp2 * _interp3 * _interp4 * (_interp13 * _interp15 + _interp14 * (_interp16 + 50. * _den1 * (-_interp13 + _interp17) * powr<6>(k))), 0.))))));
  // clang-format on
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
