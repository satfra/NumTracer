#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Cplxrt_num_kernels.hh"
#include "numtrace_verdict.hh"

namespace numtracer_kernels
{
  class Cplxrt_num_kernel
  {
    public:
    #if NT_CPLXRT_NUM_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const double& muq, const double& Ep)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::cplxrt_num::nenv) > 0 ? (numtracer_kernels::cplxrt_num::nenv) : 1];
      numtracer_kernels::cplxrt_num::fill(fenv, l0, l1, cos1, p0, p, muq, Ep);
      const auto _interp1 = ntRe(numtracer_kernels::cplxrt_num::tr1(fenv));
      const auto _interp2 = ntRe(powr<-2>(Ep + l0 + complex<double>(0.,1.) * muq));
      return _interp1 * _interp2;
    }
    #elif NT_CPLXRT_NUM_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const double& muq, const double& Ep)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::cplxrt_num::nenv) > 0 ? (numtracer_kernels::cplxrt_num::nenv) : 1];
      numtracer_kernels::cplxrt_num::fill(fenv, l0, l1, cos1, p0, p, muq, Ep);
      const auto _interp1 = ntIm(numtracer_kernels::cplxrt_num::tr0(fenv));
      const auto _interp2 = ntIm(numtracer_kernels::cplxrt_num::tr1(fenv));
      const auto _interp3 = ntIm(powr<-2>(Ep + l0 + complex<double>(0.,1.) * muq));
      const auto _interp4 = ntRe(numtracer_kernels::cplxrt_num::tr1(fenv));
      const auto _interp5 = ntRe(powr<-2>(Ep + l0 + complex<double>(0.,1.) * muq));
      return fma(-1., _interp1, fma(-1., _interp2 * _interp3, fma(_interp4, _interp5, 0.)));
    }
    #else                              // the imaginary part survives: genuinely complex
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const double& muq, const double& Ep)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::cplxrt_num::nenv) > 0 ? (numtracer_kernels::cplxrt_num::nenv) : 1];
      numtracer_kernels::cplxrt_num::fill(fenv, l0, l1, cos1, p0, p, muq, Ep);
      const auto _den1 = powr<-2>(Ep + l0 + complex<double>(0.,1.) * muq);
      return fma(complex<double>(0.,1.), numtracer_kernels::cplxrt_num::tr0(fenv), fma(_den1, numtracer_kernels::cplxrt_num::tr1(fenv), 0.));
    }
    #endif

    static inline auto constant(const double& p)
    {
      return 0.;
    }
    private:
    static inline double ntRe(double x) { return x; }
    static inline double ntIm(double) { return 0.0; }
    template <class T> static inline auto ntRe(const T &z) -> decltype(z.real()) { return z.real(); }
    template <class T> static inline auto ntRe(const T &z) -> decltype(real(z)) requires (!requires { z.real(); }) { return real(z); }
    template <class T> static inline auto ntIm(const T &z) -> decltype(z.imag()) { return z.imag(); }
    template <class T> static inline auto ntIm(const T &z) -> decltype(imag(z)) requires (!requires { z.imag(); }) { return imag(z); }
  };
}
using numtracer_kernels::Cplxrt_num_kernel;