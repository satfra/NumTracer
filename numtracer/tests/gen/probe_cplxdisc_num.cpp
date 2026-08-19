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
#include "Cplxdisc_num_kernels.hh"
template<int N, class T> static inline T powr(T x){ T r=T(1); int n=N<0?-N:N; for(int i=0;i<n;++i) r*=x; return N<0?T(1)/r:r; }
using std::pow; using std::sqrt; using std::sin; using std::cos; using std::tan; using std::exp; using std::log; using std::fma; using std::fabs;
static inline std::complex<double> fma(const std::complex<double>&a,const std::complex<double>&b,const std::complex<double>&c){return a*b+c;}
template<class T> using complex = std::complex<T>;
static inline double ntStub(double seed, double x){ double h = std::sin(seed*0.1031 + x*0.3127 + 1.7)*43758.5453; return 0.4 + 0.5*(h - std::floor(h)); }
static inline double ntRe(double x) { return x; }
template <class T> static inline double ntRe(const T &z) { return z.real(); }
static inline double ntIm(double) { return 0.0; }
template <class T> static inline double ntIm(const T &z) { return z.imag(); }
static inline auto probe_full(const double& l1, const double& cos1, const double& cos2, const double& p)
{
  double fenv[(numtracer_kernels::cplxdisc_num::nenv) > 0 ? (numtracer_kernels::cplxdisc_num::nenv) : 1];
  numtracer_kernels::cplxdisc_num::fill(fenv, l1, cos1, cos2, p);
  return complex<double>(0.,1.) * fma(numtracer_kernels::cplxdisc_num::tr1(fenv), numtracer_kernels::cplxdisc_num::tr3(fenv) * numtracer_kernels::cplxdisc_num::tr4(fenv), fma(numtracer_kernels::cplxdisc_num::tr2(fenv), numtracer_kernels::cplxdisc_num::tr5(fenv) * numtracer_kernels::cplxdisc_num::tr6(fenv) * numtracer_kernels::cplxdisc_num::tr7(fenv), numtracer_kernels::cplxdisc_num::tr0(fenv)));
}
static inline auto probe_proj(const double& l1, const double& cos1, const double& cos2, const double& p)
{
  double fenv[(numtracer_kernels::cplxdisc_num::nenv) > 0 ? (numtracer_kernels::cplxdisc_num::nenv) : 1];
  numtracer_kernels::cplxdisc_num::fill(fenv, l1, cos1, cos2, p);
  return 0.;
}
static inline auto probe_repart(const double& l1, const double& cos1, const double& cos2, const double& p)
{
  double fenv[(numtracer_kernels::cplxdisc_num::nenv) > 0 ? (numtracer_kernels::cplxdisc_num::nenv) : 1];
  numtracer_kernels::cplxdisc_num::fill(fenv, l1, cos1, cos2, p);
  const auto _interp1 = ntIm(numtracer_kernels::cplxdisc_num::tr0(fenv));
  const auto _interp2 = ntIm(numtracer_kernels::cplxdisc_num::tr1(fenv));
  const auto _interp3 = ntIm(numtracer_kernels::cplxdisc_num::tr3(fenv));
  const auto _interp4 = ntIm(numtracer_kernels::cplxdisc_num::tr4(fenv));
  const auto _interp5 = ntIm(numtracer_kernels::cplxdisc_num::tr5(fenv));
  const auto _interp6 = ntIm(numtracer_kernels::cplxdisc_num::tr6(fenv));
  const auto _interp7 = ntIm(numtracer_kernels::cplxdisc_num::tr7(fenv));
  const auto _interp8 = ntRe(numtracer_kernels::cplxdisc_num::tr2(fenv));
  const auto _interp9 = ntRe(numtracer_kernels::cplxdisc_num::tr1(fenv));
  const auto _interp10 = ntRe(numtracer_kernels::cplxdisc_num::tr3(fenv));
  const auto _interp11 = ntRe(numtracer_kernels::cplxdisc_num::tr4(fenv));
  const auto _interp12 = ntIm(numtracer_kernels::cplxdisc_num::tr2(fenv));
  const auto _interp13 = ntRe(numtracer_kernels::cplxdisc_num::tr5(fenv));
  const auto _interp14 = ntRe(numtracer_kernels::cplxdisc_num::tr6(fenv));
  const auto _interp15 = ntRe(numtracer_kernels::cplxdisc_num::tr7(fenv));
  return fma(-1., _interp1, fma(-1., _interp12 * _interp13 * _interp14 * _interp15, fma(-1., _interp10 * _interp11 * _interp2, fma(_interp2, _interp3 * _interp4, fma(_interp12, _interp15 * _interp5 * _interp6, fma(_interp12, _interp14 * _interp5 * _interp7, fma(_interp12, _interp13 * _interp6 * _interp7, fma(-1., _interp14 * _interp15 * _interp5 * _interp8, fma(-1., _interp13 * _interp15 * _interp6 * _interp8, fma(-1., _interp13 * _interp14 * _interp7 * _interp8, fma(_interp5, _interp6 * _interp7 * _interp8, fma(-1., _interp11 * _interp3 * _interp9, fma(-1., _interp10 * _interp4 * _interp9, 0.)))))))))))));
}
int main(int argc, char** argv){
  const char* outf=nullptr; const char* macro=nullptr;
  for(int i=1;i<argc;++i){ if(!std::strcmp(argv[i],"-o") && i+1<argc) outf=argv[++i];
                           else if(!std::strcmp(argv[i],"-m") && i+1<argc) macro=argv[++i]; }
  std::mt19937_64 rng(12345); std::uniform_real_distribution<double> U(0.25,3.0),Uc(-0.9,0.9),Uph(0.1,6.2);
  double mim=0,mdiff=0,mre=0,mrim=0,mrdiff=0,mrrep=0; long ok=0, nrep=0;
  for(int n=0;n<4000;++n){ double l1 = U(rng); double cos1 = Uc(rng); double cos2 = Uc(rng); double p = U(rng);
    std::complex<double> f = probe_full(l1, cos1, cos2, p); double pj = probe_proj(l1, cos1, cos2, p);
    double rp = probe_repart(l1, cos1, cos2, p);
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
