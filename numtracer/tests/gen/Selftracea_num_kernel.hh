#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Selftracea_num_kernels.hh"

namespace numtracer_kernels
{
  class Selftracea_num_kernel
  {
    public:
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const auto& Mq, const auto& Zq)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::selftracea_num::nenv) > 0 ? (numtracer_kernels::selftracea_num::nenv) : 1];
      numtracer_kernels::selftracea_num::fill(fenv, l0, l1, cos1, p0, p);
      const auto _interp1 = Mq(l1);
      return _interp1 * numtracer_kernels::selftracea_num::tr0(fenv);
    }

    static inline auto constant(const double& p, const auto& Mq, const auto& Zq)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::Selftracea_num_kernel;