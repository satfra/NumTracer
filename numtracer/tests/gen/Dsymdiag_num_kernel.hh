#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Dsymdiag_num_kernels.hh"

namespace numtracer_kernels
{
  class Dsymdiag_num_kernel
  {
    public:
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const auto& ZA8, const auto& ZA4)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::dsymdiag_num::nenv) > 0 ? (numtracer_kernels::dsymdiag_num::nenv) : 1];
      numtracer_kernels::dsymdiag_num::fill(fenv, l0, l1, cos1, p0, p);
      const auto _interp1 = ZA4(l1);
      const auto _interp2 = ZA8(l1);
      return 0.0833333333333333 * _interp1 * _interp2 * numtracer_kernels::dsymdiag_num::tr0(fenv);
    }

    static inline auto constant(const double& p, const auto& ZA8, const auto& ZA4)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::Dsymdiag_num_kernel;