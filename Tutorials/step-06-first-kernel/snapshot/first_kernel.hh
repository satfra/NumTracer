#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "first_kernels.hh"

namespace numtracer_kernels
{
  class first_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::step06::nenv) > 0 ? (numtracer_kernels::step06::nenv) : 1];
      numtracer_kernels::step06::fill(fenv, l1, cos1, p);
      return numtracer_kernels::step06::tr0(fenv) * powr<-2>(p);
    }

    static inline auto constant(const double& p)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::first_kernel;