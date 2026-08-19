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
#include "EtaPiL_num_kernels.hh"
template<int N, class T> static inline T powr(T x){ T r=T(1); int n=N<0?-N:N; for(int i=0;i<n;++i) r*=x; return N<0?T(1)/r:r; }
using std::pow; using std::sqrt; using std::sin; using std::cos; using std::tan; using std::exp; using std::log; using std::fma; using std::fabs;
static inline std::complex<double> fma(const std::complex<double>&a,const std::complex<double>&b,const std::complex<double>&c){return a*b+c;}
template<class T> using complex = std::complex<T>;
static inline double ntStub(double seed, double x){ double h = std::sin(seed*0.1031 + x*0.3127 + 1.7)*43758.5453; return 0.4 + 0.5*(h - std::floor(h)); }
static inline double ntRe(double x) { return x; }
template <class T> static inline double ntRe(const T &z) { return z.real(); }
static inline double ntIm(double) { return 0.0; }
template <class T> static inline double ntIm(const T &z) { return z.imag(); }
static inline auto probe_full(const double& l1, const double& cos1, const double& k)
{
  double fenv[(DiFfRG::etapil_num::nenv) > 0 ? (DiFfRG::etapil_num::nenv) : 1];
  const double dr_0 = 0.7071067811865475;
  const double dr_1 = sqrt(ntStub(19694., 0.));
  const double dr_2 = ntStub(69647.99999999999, 1. * l1);
  const double dr_3 = powr<-1>(l1);
  const double dr_4 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_5 = ntStub(85453., 1. * k);
  const double dr_6 = ntStub(85453., 1. * l1);
  const double dr_7 = sqrt(powr<3>(ntStub(19694., 0.)));
  const double dr_8 = ntStub(75561.00000000001, 0.816496580927726 * sqrt(1.e-6 - 0.001 * cos1 * l1 + powr<2>(l1)));
  const double dr_9 = ntStub(69647.99999999999, 1. * sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
  const double dr_10 = ntStub(75561.00000000001, 0.816496580927726 * l1);
  const double dr_11 = ntStub(19694., 0.);
  const double dr_12 = sqrt(powr<-1>(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
  const double dr_13 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
  const double dr_14 = powr<-1>(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
  const double dr_15 = cos1;
  const double dr_16 = l1;
  const double dr_17 = 1.414213562373095;
  const double dr_18 = ntStub(85453., 1. * sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
  DiFfRG::etapil_num::fill(fenv, l1, cos1, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18);
  const auto _interp1 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp2 = ntStub(85453., 1. * k);
  const auto _interp3 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp4 = ntStub(67246.00000000001, 1. * k);
  const auto _interp5 = ntStub(85453., 1.02 * k);
  const auto _interp6 = ntStub(19694., 0.);
  const auto _interp7 = ntStub(69647.99999999999, 1. * l1);
  const auto _interp8 = ntStub(85453., 1. * l1);
  const auto _interp9 = ntStub(69647.99999999999, 0.001 * sqrt(1. - 2000. * (cos1 - 500. * l1) * l1));
  const auto _interp10 = ntStub(80376., 1.618033988749895e-6 + 1. * powr<2>(k) - 0.00323606797749979 * cos1 * l1 + 1.618033988749895 * powr<2>(l1));
  const auto _interp11 = ntStub(85453., 0.001 * sqrt(1. - 2000. * (cos1 - 500. * l1) * l1));
  const auto _den1 = powr<-3>(0.5 * _interp6 * powr<2>(_interp7) + powr<2>(_interp2 * _interp3 + _interp8 * l1));
  const auto _den2 = powr<-1>(0.5 * _interp6 * powr<2>(_interp9) + powr<2>(_interp10 * _interp2 + _interp11 * sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1))));
  const auto _cse1 = powr<-1>(l1);
  return -1.e6 * _den1 * _den2 * DiFfRG::etapil_num::tr0(fenv) * fma(-1., _cse1 * _interp1 * _interp2, fma(-1., _cse1 * _interp3 * (_interp4 + 50. * (-_interp2 + _interp5)), 0.));
}
static inline auto probe_proj(const double& l1, const double& cos1, const double& k)
{
  double fenv[(DiFfRG::etapil_num::nenv) > 0 ? (DiFfRG::etapil_num::nenv) : 1];
  const double dr_0 = 0.7071067811865475;
  const double dr_1 = sqrt(ntStub(19694., 0.));
  const double dr_2 = ntStub(69647.99999999999, 1. * l1);
  const double dr_3 = powr<-1>(l1);
  const double dr_4 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_5 = ntStub(85453., 1. * k);
  const double dr_6 = ntStub(85453., 1. * l1);
  const double dr_7 = sqrt(powr<3>(ntStub(19694., 0.)));
  const double dr_8 = ntStub(75561.00000000001, 0.816496580927726 * sqrt(1.e-6 - 0.001 * cos1 * l1 + powr<2>(l1)));
  const double dr_9 = ntStub(69647.99999999999, 1. * sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
  const double dr_10 = ntStub(75561.00000000001, 0.816496580927726 * l1);
  const double dr_11 = ntStub(19694., 0.);
  const double dr_12 = sqrt(powr<-1>(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
  const double dr_13 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
  const double dr_14 = powr<-1>(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
  const double dr_15 = cos1;
  const double dr_16 = l1;
  const double dr_17 = 1.414213562373095;
  const double dr_18 = ntStub(85453., 1. * sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
  DiFfRG::etapil_num::fill(fenv, l1, cos1, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18);
  const auto _interp1 = ntRe(DiFfRG::etapil_num::tr0(fenv));
  const auto _interp2 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp3 = ntStub(85453., 1. * k);
  const auto _interp4 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp5 = ntStub(67246.00000000001, 1. * k);
  const auto _interp6 = ntStub(85453., 1.02 * k);
  const auto _interp7 = ntStub(19694., 0.);
  const auto _interp8 = ntStub(69647.99999999999, 1. * l1);
  const auto _interp9 = ntStub(85453., 1. * l1);
  const auto _interp10 = ntStub(69647.99999999999, 0.001 * sqrt(1. - 2000. * (cos1 - 500. * l1) * l1));
  const auto _interp11 = ntStub(80376., 1.618033988749895e-6 + 1. * powr<2>(k) - 0.00323606797749979 * cos1 * l1 + 1.618033988749895 * powr<2>(l1));
  const auto _interp12 = ntStub(85453., 0.001 * sqrt(1. - 2000. * (cos1 - 500. * l1) * l1));
  const auto _den1 = powr<-3>(0.5 * _interp7 * powr<2>(_interp8) + powr<2>(_interp3 * _interp4 + _interp9 * l1));
  const auto _den2 = powr<-1>(0.5 * powr<2>(_interp10) * _interp7 + powr<2>(_interp11 * _interp3 + _interp12 * sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1))));
  const auto _cse1 = powr<-1>(l1);
  return -1.e6 * _den1 * _den2 * _interp1 * fma(-1., _cse1 * _interp2 * _interp3, fma(-1., _cse1 * _interp4 * (_interp5 + 50. * (-_interp3 + _interp6)), 0.));
}
static inline auto probe_repart(const double& l1, const double& cos1, const double& k)
{
  double fenv[(DiFfRG::etapil_num::nenv) > 0 ? (DiFfRG::etapil_num::nenv) : 1];
  const double dr_0 = 0.7071067811865475;
  const double dr_1 = sqrt(ntStub(19694., 0.));
  const double dr_2 = ntStub(69647.99999999999, 1. * l1);
  const double dr_3 = powr<-1>(l1);
  const double dr_4 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const double dr_5 = ntStub(85453., 1. * k);
  const double dr_6 = ntStub(85453., 1. * l1);
  const double dr_7 = sqrt(powr<3>(ntStub(19694., 0.)));
  const double dr_8 = ntStub(75561.00000000001, 0.816496580927726 * sqrt(1.e-6 - 0.001 * cos1 * l1 + powr<2>(l1)));
  const double dr_9 = ntStub(69647.99999999999, 1. * sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
  const double dr_10 = ntStub(75561.00000000001, 0.816496580927726 * l1);
  const double dr_11 = ntStub(19694., 0.);
  const double dr_12 = sqrt(powr<-1>(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
  const double dr_13 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
  const double dr_14 = powr<-1>(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
  const double dr_15 = cos1;
  const double dr_16 = l1;
  const double dr_17 = 1.414213562373095;
  const double dr_18 = ntStub(85453., 1. * sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
  DiFfRG::etapil_num::fill(fenv, l1, cos1, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18);
  const auto _interp1 = ntRe(DiFfRG::etapil_num::tr0(fenv));
  const auto _interp2 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp3 = ntStub(85453., 1. * k);
  const auto _interp4 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp5 = ntStub(67246.00000000001, 1. * k);
  const auto _interp6 = ntStub(85453., 1.02 * k);
  const auto _interp7 = ntStub(19694., 0.);
  const auto _interp8 = ntStub(69647.99999999999, 1. * l1);
  const auto _interp9 = ntStub(85453., 1. * l1);
  const auto _interp10 = ntStub(69647.99999999999, 0.001 * sqrt(1. - 2000. * (cos1 - 500. * l1) * l1));
  const auto _interp11 = ntStub(80376., 1.618033988749895e-6 + 1. * powr<2>(k) - 0.00323606797749979 * cos1 * l1 + 1.618033988749895 * powr<2>(l1));
  const auto _interp12 = ntStub(85453., 0.001 * sqrt(1. - 2000. * (cos1 - 500. * l1) * l1));
  const auto _den1 = powr<-3>(0.5 * _interp7 * powr<2>(_interp8) + powr<2>(_interp3 * _interp4 + _interp9 * l1));
  const auto _den2 = powr<-1>(0.5 * powr<2>(_interp10) * _interp7 + powr<2>(_interp11 * _interp3 + _interp12 * sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1))));
  const auto _cse1 = powr<-1>(l1);
  return -1.e6 * _den1 * _den2 * _interp1 * fma(-1., _cse1 * _interp2 * _interp3, fma(-1., _cse1 * _interp4 * (_interp5 + 50. * (-_interp3 + _interp6)), 0.));
}
int main(int argc, char** argv){
  const char* outf=nullptr; const char* macro=nullptr;
  for(int i=1;i<argc;++i){ if(!std::strcmp(argv[i],"-o") && i+1<argc) outf=argv[++i];
                           else if(!std::strcmp(argv[i],"-m") && i+1<argc) macro=argv[++i]; }
  std::mt19937_64 rng(12345); std::uniform_real_distribution<double> U(0.25,3.0),Uc(-0.9,0.9),Uph(0.1,6.2);
  double mim=0,mdiff=0,mre=0,mrim=0,mrdiff=0,mrrep=0; long ok=0, nrep=0;
  for(int n=0;n<4000;++n){ double l1 = U(rng); double cos1 = Uc(rng); double k = U(rng);
    std::complex<double> f = probe_full(l1, cos1, k); double pj = probe_proj(l1, cos1, k);
    double rp = probe_repart(l1, cos1, k);
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
