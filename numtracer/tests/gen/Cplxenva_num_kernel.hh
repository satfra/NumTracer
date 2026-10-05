#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Cplxenva_num_kernels.hh"

namespace numtracer_kernels
{
  class Cplxenva_num_kernel
  {
    public:
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const double& muq)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::cplxenva_num::nenv) > 0 ? (numtracer_kernels::cplxenva_num::nenv) : 1];
      numtracer_kernels::cplxenva_num::fill(fenv, l0, l1, cos1, p0, p, muq);
      return ntRe(numtracer_kernels::cplxenva_num::tr0(fenv));
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
using numtracer_kernels::Cplxenva_num_kernel;