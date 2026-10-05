#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Diagfusepin_num_kernels.hh"

namespace numtracer_kernels
{
  class Diagfusepin_num_kernel
  {
    public:
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const auto& Z)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::diagfusepin_num::nenv) > 0 ? (numtracer_kernels::diagfusepin_num::nenv) : 1];
      numtracer_kernels::diagfusepin_num::fill(fenv, l0, l1, cos1, p0, p);
      const auto _interp1 = Z(l1);
      return 0.5 * _interp1 * numtracer_kernels::diagfusepin_num::tr0(fenv);
    }

    static inline auto constant(const double& p, const auto& Z)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::Diagfusepin_num_kernel;