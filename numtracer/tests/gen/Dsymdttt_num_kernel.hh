#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Dsymdttt_num_kernels.hh"

namespace numtracer_kernels
{
  class Dsymdttt_num_kernel
  {
    public:
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const auto& Mq)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::dsymdttt_num::nenv) > 0 ? (numtracer_kernels::dsymdttt_num::nenv) : 1];
      numtracer_kernels::dsymdttt_num::fill(fenv, l0, l1, cos1, p0, p);
      const auto _interp1 = Mq(l1);
      return _interp1 * numtracer_kernels::dsymdttt_num::tr0(fenv);
    }

    static inline auto constant(const double& p, const auto& Mq)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::Dsymdttt_num_kernel;