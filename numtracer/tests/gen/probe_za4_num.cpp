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
#include "ZA4_num_kernels.hh"
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
  double fenv[(DiFfRG::za4_num::nenv) > 0 ? (DiFfRG::za4_num::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_1 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_2 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_3 = ntStub(85453., 1. * k);
  const double dr_4 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_5 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_7 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_8 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_9 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_10 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_11 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_12 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_13 = ntStub(30415., 1. * l1);
  const double dr_14 = powr<-1>(l1);
  const double dr_15 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_16 = ntStub(85453., 1. * l1);
  DiFfRG::za4_num::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16);
  const auto _interp1 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp2 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp3 = ntStub(12206., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp4 = ntStub(12206., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp5 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp6 = ntStub(85453., 1. * k);
  const auto _interp7 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp8 = ntStub(67246.00000000001, 1. * k);
  const auto _interp9 = ntStub(85453., 1.02 * k);
  const auto _interp10 = ntStub(30415., 1. * l1);
  const auto _interp11 = ntStub(85453., 1. * l1);
  const auto _interp12 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp13 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp14 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp15 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp16 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp17 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp18 = ntStub(30415., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp19 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p1 - 3.23606797749979 * cosl1p2) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp20 = ntStub(85453., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp21 = ntStub(13076., 0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
  const auto _interp22 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp23 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp24 = ntStub(29441., 1. * l1);
  const auto _interp25 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp26 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp27 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp28 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp29 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp30 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp31 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp32 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp33 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp34 = ntStub(13076., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp35 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp36 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp37 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p1 - 3.23606797749979 * cosl1p2) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp38 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp39 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp40 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp41 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp42 = ntStub(13076., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp43 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _interp44 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p2 - 3.23606797749979 * cosl1p3) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp45 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp46 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp47 = ntStub(27736., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp48 = ntStub(27736., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp49 = ntStub(20621., 1. * l1);
  const auto _interp50 = ntStub(20621., 1. * k);
  const auto _interp51 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp52 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp53 = ntStub(20621., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp54 = ntStub(20621., 1.02 * k);
  const auto _interp55 = ntStub(95661., 1. * k);
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(_interp25 * _interp26 + _interp24 * powr<2>(l1));
  const auto _den3 = powr<-2>(_interp26 * _interp50 + _interp49 * powr<2>(l1));
  const auto _den4 = powr<-1>(_interp25 * _interp41 + _interp40 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _den5 = powr<-1>(_interp25 * _interp31 + _interp30 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den6 = -powr<-1>(_interp41 * _interp50 + _interp51 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _den7 = powr<-2>(powr<2>(_interp10) + powr<2>(_interp6 * _interp7 + _interp11 * l1));
  const auto _den8 = 3. * powr<-1>(3. * _interp25 * _interp37 + _interp36 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _den9 = 3. * powr<-1>(3. * _interp25 * _interp44 + _interp43 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _den10 = powr<-1>(_interp25 * _interp33 + _interp32 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _den11 = -3. * powr<-1>(3. * _interp37 * _interp50 + _interp53 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _den12 = -powr<-1>(_interp33 * _interp50 + _interp52 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _den13 = powr<-1>(powr<2>(_interp12) + powr<2>(_interp13 * _interp6 + _interp14 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
  const auto _den14 = powr<-1>(powr<2>(_interp18) + powr<2>(_interp19 * _interp6 + _interp20 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p))));
  const auto _den15 = powr<-1>(powr<2>(_interp15) + powr<2>(_interp17) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)) + _interp16 * _interp6 * (_interp16 * _interp6 + 2. * _interp17 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p))));
  // clang-format off
  return fma(0.01020408163265306, _den11 * _den12 * _den3 * _den6 * DiFfRG::za4_num::tr4(fenv) * _interp45 * _interp46 * _interp47 * _interp48 * (_interp27 * _interp50 + _interp26 * (50. * (-_interp50 + _interp54) + _interp55)), fma(0.002551020408163265, _den10 * _den2 * _den5 * DiFfRG::za4_num::tr2(fenv) * _interp21 * _interp22 * _interp23 * (_interp25 * _interp27 + _interp26 * (_interp28 + 50. * _den1 * (-_interp25 + _interp29) * powr<6>(k))), fma(0.00510204081632653, _den10 * _den2 * _den8 * DiFfRG::za4_num::tr1(fenv) * _interp23 * _interp34 * _interp35 * (_interp25 * _interp27 + _interp26 * (_interp28 + 50. * _den1 * (-_interp25 + _interp29) * powr<6>(k))), fma(0.00510204081632653, _den10 * _den2 * _den4 * _den8 * DiFfRG::za4_num::tr0(fenv) * _interp23 * _interp35 * _interp38 * _interp39 * (_interp25 * _interp27 + _interp26 * (_interp28 + 50. * _den1 * (-_interp25 + _interp29) * powr<6>(k))), fma(0.001275510204081632, _den2 * _den9 * DiFfRG::za4_num::tr3(fenv) * powr<2>(_interp42) * (_interp25 * _interp27 + _interp26 * (_interp28 + 50. * _den1 * (-_interp25 + _interp29) * powr<6>(k))), fma(complex<double>(0.,-0.02040816326530612), _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_num::tr5(fenv) * _interp1 * _interp2 * _interp3 * _interp4 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), 0.))))));
  // clang-format on
}
static inline auto probe_proj(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k)
{
  const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
  const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
  const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
  const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
  double fenv[(DiFfRG::za4_num::nenv) > 0 ? (DiFfRG::za4_num::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_1 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_2 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_3 = ntStub(85453., 1. * k);
  const double dr_4 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_5 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_7 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_8 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_9 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_10 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_11 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_12 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_13 = ntStub(30415., 1. * l1);
  const double dr_14 = powr<-1>(l1);
  const double dr_15 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_16 = ntStub(85453., 1. * l1);
  DiFfRG::za4_num::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16);
  const auto _interp1 = ntRe(DiFfRG::za4_num::tr2(fenv));
  const auto _interp2 = ntStub(13076., 0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
  const auto _interp3 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp4 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp5 = ntStub(29441., 1. * l1);
  const auto _interp6 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp7 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp8 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp9 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp10 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp11 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp12 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp13 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp14 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp15 = ntRe(DiFfRG::za4_num::tr1(fenv));
  const auto _interp16 = ntStub(13076., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp17 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp18 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp19 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p1 - 3.23606797749979 * cosl1p2) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp20 = ntRe(DiFfRG::za4_num::tr0(fenv));
  const auto _interp21 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp22 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp23 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp24 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp25 = ntRe(DiFfRG::za4_num::tr3(fenv));
  const auto _interp26 = ntStub(13076., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp27 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _interp28 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p2 - 3.23606797749979 * cosl1p3) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp29 = ntRe(DiFfRG::za4_num::tr4(fenv));
  const auto _interp30 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp31 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp32 = ntStub(27736., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp33 = ntStub(27736., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp34 = ntStub(20621., 1. * l1);
  const auto _interp35 = ntStub(20621., 1. * k);
  const auto _interp36 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp37 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp38 = ntStub(20621., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp39 = ntStub(20621., 1.02 * k);
  const auto _interp40 = ntStub(95661., 1. * k);
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(_interp6 * _interp7 + _interp5 * powr<2>(l1));
  const auto _den3 = powr<-2>(_interp35 * _interp7 + _interp34 * powr<2>(l1));
  const auto _den4 = powr<-1>(_interp24 * _interp6 + _interp23 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _den5 = powr<-1>(_interp12 * _interp6 + _interp11 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den6 = -powr<-1>(_interp24 * _interp35 + _interp36 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _den7 = 3. * powr<-1>(3. * _interp19 * _interp6 + _interp18 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _den8 = 3. * powr<-1>(3. * _interp28 * _interp6 + _interp27 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _den9 = powr<-1>(_interp14 * _interp6 + _interp13 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _den10 = -3. * powr<-1>(3. * _interp19 * _interp35 + _interp38 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _den11 = -powr<-1>(_interp14 * _interp35 + _interp37 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _cse1 = _interp6 * _interp8;
  const auto _cse2 = -_interp6;
  const auto _cse3 = _cse2 + _interp10;
  const auto _cse4 = 50. * _cse3 * _den1 * powr<6>(k);
  const auto _cse5 = _cse4 + _interp9;
  const auto _cse6 = _cse5 * _interp7;
  const auto _cse7 = _cse1 + _cse6;
  return fma(0.001275510204081632, _cse7 * _den2 * _den8 * _interp25 * powr<2>(_interp26), fma(0.00510204081632653, _cse7 * _den2 * _den7 * _den9 * _interp15 * _interp16 * _interp17 * _interp4, fma(0.00510204081632653, _cse7 * _den2 * _den4 * _den7 * _den9 * _interp17 * _interp20 * _interp21 * _interp22 * _interp4, fma(0.002551020408163265, _cse7 * _den2 * _den5 * _den9 * _interp1 * _interp2 * _interp3 * _interp4, fma(0.01020408163265306, _den10 * _den11 * _den3 * _den6 * _interp29 * _interp30 * _interp31 * _interp32 * _interp33 * ((50. * (-_interp35 + _interp39) + _interp40) * _interp7 + _interp35 * _interp8), 0.)))));
}
static inline auto probe_repart(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k)
{
  const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
  const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
  const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
  const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
  double fenv[(DiFfRG::za4_num::nenv) > 0 ? (DiFfRG::za4_num::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_1 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_2 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_3 = ntStub(85453., 1. * k);
  const double dr_4 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_5 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_7 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_8 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_9 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_10 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_11 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_12 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_13 = ntStub(30415., 1. * l1);
  const double dr_14 = powr<-1>(l1);
  const double dr_15 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_16 = ntStub(85453., 1. * l1);
  DiFfRG::za4_num::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16);
  const auto _interp1 = ntIm(DiFfRG::za4_num::tr5(fenv));
  const auto _interp2 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp3 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp4 = ntStub(12206., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp5 = ntStub(12206., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp6 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp7 = ntStub(85453., 1. * k);
  const auto _interp8 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp9 = ntStub(67246.00000000001, 1. * k);
  const auto _interp10 = ntStub(85453., 1.02 * k);
  const auto _interp11 = ntStub(30415., 1. * l1);
  const auto _interp12 = ntStub(85453., 1. * l1);
  const auto _interp13 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp14 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp15 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp16 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp17 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp18 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp19 = ntStub(30415., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp20 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p1 - 3.23606797749979 * cosl1p2) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp21 = ntStub(85453., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp22 = ntRe(DiFfRG::za4_num::tr2(fenv));
  const auto _interp23 = ntStub(13076., 0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
  const auto _interp24 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp25 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp26 = ntStub(29441., 1. * l1);
  const auto _interp27 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp28 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp29 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp30 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp31 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp32 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp33 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp34 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp35 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp36 = ntRe(DiFfRG::za4_num::tr1(fenv));
  const auto _interp37 = ntStub(13076., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp38 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp39 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp40 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p1 - 3.23606797749979 * cosl1p2) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp41 = ntRe(DiFfRG::za4_num::tr0(fenv));
  const auto _interp42 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp43 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp44 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp45 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp46 = ntRe(DiFfRG::za4_num::tr3(fenv));
  const auto _interp47 = ntStub(13076., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp48 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _interp49 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p2 - 3.23606797749979 * cosl1p3) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp50 = ntRe(DiFfRG::za4_num::tr4(fenv));
  const auto _interp51 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp52 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp53 = ntStub(27736., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp54 = ntStub(27736., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp55 = ntStub(20621., 1. * l1);
  const auto _interp56 = ntStub(20621., 1. * k);
  const auto _interp57 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp58 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp59 = ntStub(20621., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp60 = ntStub(20621., 1.02 * k);
  const auto _interp61 = ntStub(95661., 1. * k);
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(_interp27 * _interp28 + _interp26 * powr<2>(l1));
  const auto _den3 = powr<-2>(_interp28 * _interp56 + _interp55 * powr<2>(l1));
  const auto _den4 = powr<-1>(_interp27 * _interp45 + _interp44 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _den5 = powr<-1>(_interp27 * _interp33 + _interp32 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den6 = -powr<-1>(_interp45 * _interp56 + _interp57 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _den7 = powr<-2>(powr<2>(_interp11) + powr<2>(_interp7 * _interp8 + _interp12 * l1));
  const auto _den8 = 3. * powr<-1>(3. * _interp27 * _interp40 + _interp39 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _den9 = 3. * powr<-1>(3. * _interp27 * _interp49 + _interp48 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _den10 = powr<-1>(_interp27 * _interp35 + _interp34 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _den11 = -3. * powr<-1>(3. * _interp40 * _interp56 + _interp59 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _den12 = -powr<-1>(_interp35 * _interp56 + _interp58 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _den13 = powr<-1>(powr<2>(_interp13) + powr<2>(_interp14 * _interp7 + _interp15 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
  const auto _den14 = powr<-1>(powr<2>(_interp19) + powr<2>(_interp20 * _interp7 + _interp21 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p))));
  const auto _den15 = powr<-1>(powr<2>(_interp16) + powr<2>(_interp18) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)) + _interp17 * _interp7 * (_interp17 * _interp7 + 2. * _interp18 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p))));
  // clang-format off
  return fma(0.01020408163265306, _den11 * _den12 * _den3 * _den6 * _interp50 * _interp51 * _interp52 * _interp53 * _interp54 * (_interp29 * _interp56 + _interp28 * (50. * (-_interp56 + _interp60) + _interp61)), fma(0.002551020408163265, _den10 * _den2 * _den5 * _interp22 * _interp23 * _interp24 * _interp25 * (_interp27 * _interp29 + _interp28 * (_interp30 + 50. * _den1 * (-_interp27 + _interp31) * powr<6>(k))), fma(0.00510204081632653, _den10 * _den2 * _den8 * _interp25 * _interp36 * _interp37 * _interp38 * (_interp27 * _interp29 + _interp28 * (_interp30 + 50. * _den1 * (-_interp27 + _interp31) * powr<6>(k))), fma(0.00510204081632653, _den10 * _den2 * _den4 * _den8 * _interp25 * _interp38 * _interp41 * _interp42 * _interp43 * (_interp27 * _interp29 + _interp28 * (_interp30 + 50. * _den1 * (-_interp27 + _interp31) * powr<6>(k))), fma(0.001275510204081632, _den2 * _den9 * _interp46 * powr<2>(_interp47) * (_interp27 * _interp29 + _interp28 * (_interp30 + 50. * _den1 * (-_interp27 + _interp31) * powr<6>(k))), fma(0.02040816326530612, _den13 * _den14 * _den15 * _den7 * _interp1 * _interp2 * _interp3 * _interp4 * _interp5 * (-_interp6 * _interp7 * powr<-1>(l1) - _interp8 * (50. * (_interp10 - _interp7) + _interp9) * powr<-1>(l1)), 0.))))));
  // clang-format on
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
