#define KOKKOS_INLINE_FUNCTION inline
#define KOKKOS_FORCEINLINE_FUNCTION inline
#define KOKKOS_FUNCTION
#define __host__
#define __device__
#include <complex>
#include <cmath>
#include <algorithm>
#include <random>
#include <cstdio>
#include <cstring>
#include "Cplxrt_num_kernels.hh"
template<int N, class T> static inline T powr(T x){ T r=T(1); int n=N<0?-N:N; for(int i=0;i<n;++i) r*=x; return N<0?T(1)/r:r; }
using std::pow; using std::sqrt; using std::sin; using std::cos; using std::tan; using std::exp; using std::log; using std::fma; using std::fabs; using std::min; using std::max;
static inline std::complex<double> fma(const std::complex<double>&a,const std::complex<double>&b,const std::complex<double>&c){return a*b+c;}
template<class T> using complex = std::complex<T>;
static inline double ntStub(double seed, double x){ double h = std::sin(seed*0.1031 + x*0.3127 + 1.7)*43758.5453; return 0.4 + 0.5*(h - std::floor(h)); }
static inline double ntRe(double x) { return x; }
static inline double ntIm(double) { return 0.0; }
template <class T> static inline auto ntRe(const T &z) -> decltype(z.real()) { return z.real(); }
template <class T> static inline auto ntRe(const T &z) -> decltype(real(z)) requires (!requires { z.real(); }) { return real(z); }
template <class T> static inline auto ntIm(const T &z) -> decltype(z.imag()) { return z.imag(); }
template <class T> static inline auto ntIm(const T &z) -> decltype(imag(z)) requires (!requires { z.imag(); }) { return imag(z); }
static inline auto probe_full(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const double& muq, const double& Ep)
{
  double fenv[(numtracer_kernels::cplxrt_num::nenv) > 0 ? (numtracer_kernels::cplxrt_num::nenv) : 1];
  numtracer_kernels::cplxrt_num::fill(fenv, l0, l1, cos1, p0, p, muq, Ep);
  const auto _den1 = powr<-2>(Ep + l0 + complex<double>(0.,1.) * muq);
  return fma(complex<double>(0.,1.), numtracer_kernels::cplxrt_num::tr0(fenv), fma(_den1, numtracer_kernels::cplxrt_num::tr1(fenv), 0.));
}
static inline auto probe_proj(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const double& muq, const double& Ep)
{
  double fenv[(numtracer_kernels::cplxrt_num::nenv) > 0 ? (numtracer_kernels::cplxrt_num::nenv) : 1];
  numtracer_kernels::cplxrt_num::fill(fenv, l0, l1, cos1, p0, p, muq, Ep);
  const auto _interp1 = ntRe(numtracer_kernels::cplxrt_num::tr1(fenv));
  const auto _interp2 = ntRe(powr<-2>(Ep + l0 + complex<double>(0.,1.) * muq));
  return _interp1 * _interp2;
}
static inline auto probe_repart(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const double& muq, const double& Ep)
{
  double fenv[(numtracer_kernels::cplxrt_num::nenv) > 0 ? (numtracer_kernels::cplxrt_num::nenv) : 1];
  numtracer_kernels::cplxrt_num::fill(fenv, l0, l1, cos1, p0, p, muq, Ep);
  const auto _interp1 = ntIm(numtracer_kernels::cplxrt_num::tr0(fenv));
  const auto _interp2 = ntIm(numtracer_kernels::cplxrt_num::tr1(fenv));
  const auto _interp3 = ntIm(powr<-2>(Ep + l0 + complex<double>(0.,1.) * muq));
  const auto _interp4 = ntRe(numtracer_kernels::cplxrt_num::tr1(fenv));
  const auto _interp5 = ntRe(powr<-2>(Ep + l0 + complex<double>(0.,1.) * muq));
  return fma(-1., _interp1, fma(-1., _interp2 * _interp3, fma(_interp4, _interp5, 0.)));
}
int main(int argc, char** argv){
  const char* outf=nullptr; const char* macro=nullptr;
  for(int i=1;i<argc;++i){ if(!std::strcmp(argv[i],"-o") && i+1<argc) outf=argv[++i];
                           else if(!std::strcmp(argv[i],"-m") && i+1<argc) macro=argv[++i]; }
  std::mt19937_64 rng(12345); std::uniform_real_distribution<double> U(0.25,3.0),Uc(-0.9,0.9),Uph(0.1,6.2);
  double mim=0,mdiff=0,mre=0,mrim=0,mrdiff=0,mrrep=0; long ok=0, nrep=0;
  for(int n=0;n<4000;++n){ double l0 = U(rng); double l1 = U(rng); double cos1 = Uc(rng); double p0 = U(rng); double p = U(rng); double muq = U(rng); double Ep = U(rng);
    std::complex<double> f = probe_full(l0, l1, cos1, p0, p, muq, Ep); double pj = probe_proj(l0, l1, cos1, p0, p, muq, Ep);
    double rp = probe_repart(l0, l1, cos1, p0, p, muq, Ep);
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
