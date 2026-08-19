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
#include "ZAAqbq1SmallRef_num_kernels.hh"
template<int N, class T> static inline T powr(T x){ T r=T(1); int n=N<0?-N:N; for(int i=0;i<n;++i) r*=x; return N<0?T(1)/r:r; }
using std::pow; using std::sqrt; using std::sin; using std::cos; using std::tan; using std::exp; using std::log; using std::fma; using std::fabs;
static inline std::complex<double> fma(const std::complex<double>&a,const std::complex<double>&b,const std::complex<double>&c){return a*b+c;}
template<class T> using complex = std::complex<T>;
static inline double ntStub(double seed, double x){ double h = std::sin(seed*0.1031 + x*0.3127 + 1.7)*43758.5453; return 0.4 + 0.5*(h - std::floor(h)); }
static inline double ntRe(double x) { return x; }
template <class T> static inline double ntRe(const T &z) { return z.real(); }
static inline double ntIm(double) { return 0.0; }
template <class T> static inline double ntIm(const T &z) { return z.imag(); }
static inline auto probe_full(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k)
{
  const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
  const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
  const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
  const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
  double fenv[(DiFfRG::zaaqbq1_small_ref::nenv) > 0 ? (DiFfRG::zaaqbq1_small_ref::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * l1);
  const double dr_1 = powr<-1>(l1);
  const double dr_2 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_3 = ntStub(85453., 1. * k);
  const double dr_4 = ntStub(85453., 1. * l1);
  DiFfRG::zaaqbq1_small_ref::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4);
  const auto _interp1 = ntStub(13076., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp2 = ntStub(29441., 1. * l1);
  const auto _interp3 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp4 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp5 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp6 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp7 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp8 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp9 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p1 - 3.23606797749979 * cosl1p2) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp10 = ntStub(98454., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp11 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp12 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp13 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp14 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp15 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp16 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + 2. * (cosl1p2 + cosl1p3)) * l1 * p + 5. * powr<2>(p)));
  const auto _interp17 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp18 = ntStub(85453., 1. * k);
  const auto _interp19 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp20 = ntStub(67246.00000000001, 1. * k);
  const auto _interp21 = ntStub(85453., 1.02 * k);
  const auto _interp22 = ntStub(30415., 1. * l1);
  const auto _interp23 = ntStub(85453., 1. * l1);
  const auto _interp24 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp25 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp26 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _interp27 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p2 - 3.23606797749979 * cosl1p3) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp28 = ntStub(98454., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp29 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp30 = ntStub(29441., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp31 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp32 = ntStub(98454., 0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp2 * powr<2>(l1));
  const auto _den3 = powr<-1>(_interp14 * _interp3 + _interp13 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _den4 = powr<-1>(_interp3 * _interp31 + _interp30 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den5 = powr<-2>(powr<2>(_interp22) + powr<2>(_interp18 * _interp19 + _interp23 * l1));
  const auto _den6 = 3. * powr<-1>(3. * _interp3 * _interp9 + _interp8 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _den7 = 3. * powr<-1>(3. * _interp27 * _interp3 + _interp26 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _den8 = powr<-1>(_interp25 * _interp3 + _interp24 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _cse1 = powr<-1>(l1);
  const auto _cse2 = powr<-2>(p);
  const auto _cse3 = _interp3 * _interp5;
  const auto _cse4 = -_interp3;
  const auto _cse5 = _cse4 + _interp7;
  const auto _cse6 = 50. * _cse5 * _den1 * powr<6>(k);
  const auto _cse7 = _cse6 + _interp6;
  const auto _cse8 = _cse7 * _interp4;
  const auto _cse9 = _cse3 + _cse8;
  return fma(-0.0003662109375, _cse2 * _cse9 * _den2 * _den6 * DiFfRG::zaaqbq1_small_ref::tr2(fenv) * _interp1 * _interp10, fma(-0.000732421875, _cse2 * _cse9 * _den2 * _den3 * _den6 * DiFfRG::zaaqbq1_small_ref::tr1(fenv) * _interp10 * _interp11 * _interp12, fma(complex<double>(0.,0.000732421875), _cse2 * _den5 * _den7 * _den8 * DiFfRG::zaaqbq1_small_ref::tr3(fenv) * _interp15 * _interp16 * (-_cse1 * _interp17 * _interp18 - _cse1 * _interp19 * (_interp20 + 50. * (-_interp18 + _interp21))) * _interp28, fma(-0.0003662109375, _cse2 * _cse9 * _den2 * _den3 * _den4 * DiFfRG::zaaqbq1_small_ref::tr0(fenv) * _interp11 * _interp29 * _interp32, 0.))));
}
static inline auto probe_proj(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k)
{
  const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
  const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
  const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
  const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
  double fenv[(DiFfRG::zaaqbq1_small_ref::nenv) > 0 ? (DiFfRG::zaaqbq1_small_ref::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * l1);
  const double dr_1 = powr<-1>(l1);
  const double dr_2 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_3 = ntStub(85453., 1. * k);
  const double dr_4 = ntStub(85453., 1. * l1);
  DiFfRG::zaaqbq1_small_ref::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4);
  const auto _interp1 = ntRe(DiFfRG::zaaqbq1_small_ref::tr2(fenv));
  const auto _interp2 = ntStub(13076., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp3 = ntStub(29441., 1. * l1);
  const auto _interp4 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp5 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp6 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp7 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp8 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp9 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp10 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p1 - 3.23606797749979 * cosl1p2) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp11 = ntStub(98454., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp12 = ntRe(DiFfRG::zaaqbq1_small_ref::tr1(fenv));
  const auto _interp13 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp14 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp15 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp16 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp17 = ntRe(DiFfRG::zaaqbq1_small_ref::tr0(fenv));
  const auto _interp18 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp19 = ntStub(29441., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp20 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp21 = ntStub(98454., 0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(_interp4 * _interp5 + _interp3 * powr<2>(l1));
  const auto _den3 = powr<-1>(_interp16 * _interp4 + _interp15 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _den4 = powr<-1>(_interp20 * _interp4 + _interp19 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den5 = 3. * powr<-1>(3. * _interp10 * _interp4 + _interp9 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _cse1 = powr<-2>(p);
  const auto _cse2 = -_interp4;
  const auto _cse3 = _cse2 + _interp8;
  const auto _cse4 = 50. * _cse3 * _den1 * powr<6>(k);
  const auto _cse5 = _cse4 + _interp7;
  const auto _cse6 = _cse5 * _interp5;
  const auto _cse7 = _interp4 * _interp6;
  const auto _cse8 = _cse6 + _cse7;
  return -0.0003662109375 * fma(2., _cse1 * _cse8 * _den2 * _den3 * _den5 * _interp11 * _interp12 * _interp13 * _interp14, fma(_cse1, _cse8 * _den2 * _den5 * _interp1 * _interp11 * _interp2, fma(_cse1, _cse8 * _den2 * _den3 * _den4 * _interp13 * _interp17 * _interp18 * _interp21, 0.)));
}
static inline auto probe_repart(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k)
{
  const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
  const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
  const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
  const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
  double fenv[(DiFfRG::zaaqbq1_small_ref::nenv) > 0 ? (DiFfRG::zaaqbq1_small_ref::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * l1);
  const double dr_1 = powr<-1>(l1);
  const double dr_2 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_3 = ntStub(85453., 1. * k);
  const double dr_4 = ntStub(85453., 1. * l1);
  DiFfRG::zaaqbq1_small_ref::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4);
  const auto _interp1 = ntRe(DiFfRG::zaaqbq1_small_ref::tr2(fenv));
  const auto _interp2 = ntStub(13076., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp3 = ntStub(29441., 1. * l1);
  const auto _interp4 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp5 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp6 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp7 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp8 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp9 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp10 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p1 - 3.23606797749979 * cosl1p2) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp11 = ntStub(98454., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp12 = ntRe(DiFfRG::zaaqbq1_small_ref::tr1(fenv));
  const auto _interp13 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp14 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp15 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp16 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp17 = ntIm(DiFfRG::zaaqbq1_small_ref::tr3(fenv));
  const auto _interp18 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp19 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + 2. * (cosl1p2 + cosl1p3)) * l1 * p + 5. * powr<2>(p)));
  const auto _interp20 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp21 = ntStub(85453., 1. * k);
  const auto _interp22 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp23 = ntStub(67246.00000000001, 1. * k);
  const auto _interp24 = ntStub(85453., 1.02 * k);
  const auto _interp25 = ntStub(30415., 1. * l1);
  const auto _interp26 = ntStub(85453., 1. * l1);
  const auto _interp27 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp28 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp29 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _interp30 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p2 - 3.23606797749979 * cosl1p3) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp31 = ntStub(98454., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp32 = ntRe(DiFfRG::zaaqbq1_small_ref::tr0(fenv));
  const auto _interp33 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp34 = ntStub(29441., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp35 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp36 = ntStub(98454., 0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(_interp4 * _interp5 + _interp3 * powr<2>(l1));
  const auto _den3 = powr<-1>(_interp16 * _interp4 + _interp15 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _den4 = powr<-1>(_interp35 * _interp4 + _interp34 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den5 = powr<-2>(powr<2>(_interp25) + powr<2>(_interp21 * _interp22 + _interp26 * l1));
  const auto _den6 = 3. * powr<-1>(3. * _interp10 * _interp4 + _interp9 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _den7 = 3. * powr<-1>(3. * _interp30 * _interp4 + _interp29 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _den8 = powr<-1>(_interp28 * _interp4 + _interp27 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _cse1 = powr<-1>(l1);
  const auto _cse2 = powr<-2>(p);
  const auto _cse3 = -_interp4;
  const auto _cse4 = _cse3 + _interp8;
  const auto _cse5 = 50. * _cse4 * _den1 * powr<6>(k);
  const auto _cse6 = _cse5 + _interp7;
  const auto _cse7 = _cse6 * _interp5;
  const auto _cse8 = _interp4 * _interp6;
  const auto _cse9 = _cse7 + _cse8;
  return fma(-0.000732421875, _cse2 * _cse9 * _den2 * _den3 * _den6 * _interp11 * _interp12 * _interp13 * _interp14, fma(-0.0003662109375, _cse2 * _cse9 * _den2 * _den6 * _interp1 * _interp11 * _interp2, fma(-0.000732421875, _cse2 * _den5 * _den7 * _den8 * _interp17 * _interp18 * _interp19 * (-_cse1 * _interp20 * _interp21 - _cse1 * _interp22 * (_interp23 + 50. * (-_interp21 + _interp24))) * _interp31, fma(-0.0003662109375, _cse2 * _cse9 * _den2 * _den3 * _den4 * _interp13 * _interp32 * _interp33 * _interp36, 0.))));
}
int main(int argc, char** argv){
  const char* outf=nullptr; const char* macro=nullptr;
  for(int i=1;i<argc;++i){ if(!std::strcmp(argv[i],"-o") && i+1<argc) outf=argv[++i];
                           else if(!std::strcmp(argv[i],"-m") && i+1<argc) macro=argv[++i]; }
  std::mt19937_64 rng(12345); std::uniform_real_distribution<double> U(0.25,3.0),Uc(-0.9,0.9),Uph(0.1,6.2);
  double mim=0,mdiff=0,mre=0,mrim=0,mrdiff=0,mrrep=0; long ok=0, nrep=0;
  for(int n=0;n<4000;++n){ double l1 = U(rng); double cos1 = Uc(rng); double cos2 = Uc(rng); double phi = Uph(rng); double p = U(rng); double k = U(rng);
    std::complex<double> f = probe_full(l1, cos1, cos2, phi, p, k); double pj = probe_proj(l1, cos1, cos2, phi, p, k);
    double rp = probe_repart(l1, cos1, cos2, phi, p, k);
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
