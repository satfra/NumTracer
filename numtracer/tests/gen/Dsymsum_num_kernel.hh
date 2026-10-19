#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Dsymsum_num_kernels.hh"

namespace numtracer_kernels
{
  class Dsymsum_num_kernel
  {
    public:
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::dsymsum_num::nenv) > 0 ? (numtracer_kernels::dsymsum_num::nenv) : 1];
      numtracer_kernels::dsymsum_num::fill(fenv, l0, l1, cos1, p0, p);
      return numtracer_kernels::dsymsum_num::tr0(fenv);
    }

    static inline auto constant(const double& p)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::Dsymsum_num_kernel;