#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Cplxrtro_num_kernels.hh"
#include "numtrace_verdict.hh"

namespace numtracer_kernels
{
  class Cplxrtro_num_kernel
  {
    public:
    #if NT_CPLXRTRO_NUM_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const double& muq, const double& Ep)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::cplxrtro_num::nenv) > 0 ? (numtracer_kernels::cplxrtro_num::nenv) : 1];
      numtracer_kernels::cplxrtro_num::fill(fenv, l0, l1, cos1, p0, p, muq, Ep);
      const auto _interp1 = ntRe(numtracer_kernels::cplxrtro_num::tr1(fenv));
      const auto _interp2 = ntRe(powr<-2>(Ep + l0 + complex<double>(0.,1.) * muq));
      return _interp1 * _interp2;
    }
    #else                              // 1 = RePart; 0 = complex, truncated by RealOutput
    #  if NT_CPLXRTRO_NUM_VERDICT == 0
    #    warning "NumTracer: flow 'cplxrtro_num' probed GENUINELY COMPLEX (verdict 0) but was generated with RealOutput -> True. The kernel returns only the real part; the imaginary part of the integrand is discarded. That is a truncation of the flow equation, not an identity. If it is not what you intended, regenerate without RealOutput and give the consumer a complex integrator."
    #  endif
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const double& muq, const double& Ep)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::cplxrtro_num::nenv) > 0 ? (numtracer_kernels::cplxrtro_num::nenv) : 1];
      numtracer_kernels::cplxrtro_num::fill(fenv, l0, l1, cos1, p0, p, muq, Ep);
      const auto _interp1 = ntIm(numtracer_kernels::cplxrtro_num::tr0(fenv));
      const auto _interp2 = ntIm(numtracer_kernels::cplxrtro_num::tr1(fenv));
      const auto _interp3 = ntIm(powr<-2>(Ep + l0 + complex<double>(0.,1.) * muq));
      const auto _interp4 = ntRe(numtracer_kernels::cplxrtro_num::tr1(fenv));
      const auto _interp5 = ntRe(powr<-2>(Ep + l0 + complex<double>(0.,1.) * muq));
      return fma(-1., _interp1, fma(-1., _interp2 * _interp3, fma(_interp4, _interp5, 0.)));
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
using numtracer_kernels::Cplxrtro_num_kernel;