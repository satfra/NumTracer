#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Cplxcoef_num_kernels.hh"

namespace numtracer_kernels
{
  class Cplxcoef_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::cplxcoef_num::nenv) > 0 ? (numtracer_kernels::cplxcoef_num::nenv) : 1];
      numtracer_kernels::cplxcoef_num::fill(fenv, l1, cos1, cos2, p);
      return ntRe(complex<double>(0.,1.) * numtracer_kernels::cplxcoef_num::tr0(fenv) * numtracer_kernels::cplxcoef_num::tr1(fenv) * numtracer_kernels::cplxcoef_num::tr2(fenv));
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
using numtracer_kernels::Cplxcoef_num_kernel;