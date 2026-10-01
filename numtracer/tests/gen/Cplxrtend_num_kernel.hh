#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Cplxrtend_num_kernels.hh"

namespace numtracer_kernels
{
  class Cplxrtend_num_kernel
  {
    public:
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const double& muq, const double& Ep)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::cplxrtend_num::nenv) > 0 ? (numtracer_kernels::cplxrtend_num::nenv) : 1];
      numtracer_kernels::cplxrtend_num::fill(fenv, l0, l1, cos1, p0, p, muq, Ep);
      const auto _den1 = powr<-2>(Ep + l0 + complex<double>(0.,1.) * muq);
      return ntRe(fma(complex<double>(0.,1.), numtracer_kernels::cplxrtend_num::tr0(fenv), fma(_den1, numtracer_kernels::cplxrtend_num::tr1(fenv), 0.)));
    }

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
using numtracer_kernels::Cplxrtend_num_kernel;