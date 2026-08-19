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
#include "Zq_collect_kernels.hh"
template<int N, class T> static inline T powr(T x){ T r=T(1); int n=N<0?-N:N; for(int i=0;i<n;++i) r*=x; return N<0?T(1)/r:r; }
using std::pow; using std::sqrt; using std::sin; using std::cos; using std::tan; using std::exp; using std::log; using std::fma; using std::fabs;
static inline std::complex<double> fma(const std::complex<double>&a,const std::complex<double>&b,const std::complex<double>&c){return a*b+c;}
template<class T> using complex = std::complex<T>;
static inline double ntStub(double seed, double x){ double h = std::sin(seed*0.1031 + x*0.3127 + 1.7)*43758.5453; return 0.4 + 0.5*(h - std::floor(h)); }
static inline double ntRe(double x) { return x; }
template <class T> static inline double ntRe(const T &z) { return z.real(); }
static inline double ntIm(double) { return 0.0; }
template <class T> static inline double ntIm(const T &z) { return z.imag(); }
static inline auto probe_full(const double& l1, const double& cos1, const double& p, const double& k)
{
  double fenv[(DiFfRG::zq_collect::nenv) > 0 ? (DiFfRG::zq_collect::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * l1);
  const double dr_1 = sqrt(powr<-1>(powr<2>(l1)));
  const double dr_2 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_3 = ntStub(85453., 1. * k);
  const double dr_4 = ntStub(85453., 1. * l1);
  const double dr_5 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const double dr_7 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const double dr_8 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  DiFfRG::zq_collect::fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8);
  const auto _interp1 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
  const auto _interp2 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp3 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp4 = ntStub(85453., 1. * k);
  const auto _interp5 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp6 = ntStub(29441., 1. * l1);
  const auto _interp7 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp8 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp9 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp10 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp11 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp12 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp13 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp14 = ntStub(67246.00000000001, 1. * k);
  const auto _interp15 = ntStub(85453., 1.02 * k);
  const auto _interp16 = ntStub(30415., 1. * l1);
  const auto _interp17 = ntStub(85453., 1. * l1);
  const auto _interp18 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp19 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(_interp7 * _interp8 + _interp6 * powr<2>(l1));
  const auto _den3 = powr<-1>(_interp19 * _interp7 + _interp18 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _den4 = powr<-2>(powr<2>(_interp16) + powr<2>(_interp13 * _interp4 + _interp17 * l1));
  const auto _den5 = powr<-1>(powr<2>(_interp2) + powr<2>(_interp3 * _interp4 + _interp5 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
  const auto _cse1 = powr<-1>(l1);
  const auto _cse2 = powr<-2>(p);
  return fma(-0.0833333333333333, _cse2 * _den3 * _den4 * DiFfRG::zq_collect::tr0(fenv) * powr<2>(_interp1) * (-_cse1 * _interp13 * (_interp14 + 50. * (_interp15 - _interp4)) - _cse1 * _interp12 * _interp4), fma(complex<double>(0.,0.0833333333333333), _cse2 * _den2 * _den5 * DiFfRG::zq_collect::tr1(fenv) * powr<2>(_interp1) * (_interp7 * _interp9 + _interp8 * (_interp10 + 50. * _den1 * (_interp11 - _interp7) * powr<6>(k))), 0.));
}
static inline auto probe_proj(const double& l1, const double& cos1, const double& p, const double& k)
{
  double fenv[(DiFfRG::zq_collect::nenv) > 0 ? (DiFfRG::zq_collect::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * l1);
  const double dr_1 = sqrt(powr<-1>(powr<2>(l1)));
  const double dr_2 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_3 = ntStub(85453., 1. * k);
  const double dr_4 = ntStub(85453., 1. * l1);
  const double dr_5 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const double dr_7 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const double dr_8 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  DiFfRG::zq_collect::fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8);
  const auto _interp1 = ntRe(DiFfRG::zq_collect::tr0(fenv));
  const auto _interp2 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
  const auto _interp3 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp4 = ntStub(85453., 1. * k);
  const auto _interp5 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp6 = ntStub(67246.00000000001, 1. * k);
  const auto _interp7 = ntStub(85453., 1.02 * k);
  const auto _interp8 = ntStub(30415., 1. * l1);
  const auto _interp9 = ntStub(85453., 1. * l1);
  const auto _interp10 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp11 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp12 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _den1 = powr<-1>(_interp11 * _interp12 + _interp10 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _den2 = powr<-2>(powr<2>(_interp8) + powr<2>(_interp4 * _interp5 + _interp9 * l1));
  const auto _cse1 = powr<-1>(l1);
  return -0.0833333333333333 * _den1 * _den2 * _interp1 * powr<2>(_interp2) * powr<-2>(p) * fma(-1., _cse1 * _interp3 * _interp4, fma(-1., _cse1 * _interp5 * (_interp6 + 50. * (-_interp4 + _interp7)), 0.));
}
static inline auto probe_repart(const double& l1, const double& cos1, const double& p, const double& k)
{
  double fenv[(DiFfRG::zq_collect::nenv) > 0 ? (DiFfRG::zq_collect::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * l1);
  const double dr_1 = sqrt(powr<-1>(powr<2>(l1)));
  const double dr_2 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_3 = ntStub(85453., 1. * k);
  const double dr_4 = ntStub(85453., 1. * l1);
  const double dr_5 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const double dr_7 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const double dr_8 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  DiFfRG::zq_collect::fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8);
  const auto _interp1 = ntIm(DiFfRG::zq_collect::tr1(fenv));
  const auto _interp2 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
  const auto _interp3 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp4 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp5 = ntStub(85453., 1. * k);
  const auto _interp6 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp7 = ntStub(29441., 1. * l1);
  const auto _interp8 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp9 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp10 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp11 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp12 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp13 = ntRe(DiFfRG::zq_collect::tr0(fenv));
  const auto _interp14 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp15 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp16 = ntStub(67246.00000000001, 1. * k);
  const auto _interp17 = ntStub(85453., 1.02 * k);
  const auto _interp18 = ntStub(30415., 1. * l1);
  const auto _interp19 = ntStub(85453., 1. * l1);
  const auto _interp20 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp21 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(_interp8 * _interp9 + _interp7 * powr<2>(l1));
  const auto _den3 = powr<-1>(_interp21 * _interp8 + _interp20 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _den4 = powr<-2>(powr<2>(_interp18) + powr<2>(_interp15 * _interp5 + _interp19 * l1));
  const auto _den5 = powr<-1>(powr<2>(_interp3) + powr<2>(_interp4 * _interp5 + _interp6 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
  const auto _cse1 = powr<-1>(l1);
  const auto _cse2 = powr<-2>(p);
  return fma(-0.0833333333333333, _cse2 * _den3 * _den4 * _interp13 * powr<2>(_interp2) * (-_cse1 * _interp15 * (_interp16 + 50. * (_interp17 - _interp5)) - _cse1 * _interp14 * _interp5), fma(-0.0833333333333333, _cse2 * _den2 * _den5 * _interp1 * powr<2>(_interp2) * (_interp10 * _interp8 + _interp9 * (_interp11 + 50. * _den1 * (_interp12 - _interp8) * powr<6>(k))), 0.));
}
int main(int argc, char** argv){
  const char* outf=nullptr; const char* macro=nullptr;
  for(int i=1;i<argc;++i){ if(!std::strcmp(argv[i],"-o") && i+1<argc) outf=argv[++i];
                           else if(!std::strcmp(argv[i],"-m") && i+1<argc) macro=argv[++i]; }
  std::mt19937_64 rng(12345); std::uniform_real_distribution<double> U(0.25,3.0),Uc(-0.9,0.9),Uph(0.1,6.2);
  double mim=0,mdiff=0,mre=0,mrim=0,mrdiff=0,mrrep=0; long ok=0, nrep=0;
  for(int n=0;n<4000;++n){ double l1 = U(rng); double cos1 = Uc(rng); double p = U(rng); double k = U(rng);
    std::complex<double> f = probe_full(l1, cos1, p, k); double pj = probe_proj(l1, cos1, p, k);
    double rp = probe_repart(l1, cos1, p, k);
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
