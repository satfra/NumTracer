#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Cplxenvb_num_kernels.hh"

namespace numtracer_kernels
{
  class Cplxenvb_num_kernel
  {
    public:
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const double& muq, const auto& Zq)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::cplxenvb_num::nenv) > 0 ? (numtracer_kernels::cplxenvb_num::nenv) : 1];
      const double dr_0 = l0;
      const double dr_0_im = muq;
      numtracer_kernels::cplxenvb_num::fill(fenv, l0, l1, cos1, p0, p, muq, dr_0, dr_0_im);
      const auto _interp1 = Zq(l1);
      return ntRe(complex<double>(0.,1.) * _interp1 * numtracer_kernels::cplxenvb_num::tr0(fenv));
    }

    static inline auto constant(const double& p, const auto& Zq)
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
using numtracer_kernels::Cplxenvb_num_kernel;