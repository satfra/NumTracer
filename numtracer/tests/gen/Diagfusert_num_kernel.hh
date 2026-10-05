#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Diagfusert_num_kernels.hh"

namespace numtracer_kernels
{
  class Diagfusert_num_kernel
  {
    public:
    static inline auto kernel(const double& l0, const double& l1, const double& cos1, const double& p0, const double& p, const auto& Z, const auto& Gu, const auto& Gd)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::diagfusert_num::nenv) > 0 ? (numtracer_kernels::diagfusert_num::nenv) : 1];
      numtracer_kernels::diagfusert_num::fill(fenv, l0, l1, cos1, p0, p);
      const auto _interp1 = Gd(l1);
      const auto _interp2 = Gu(l1);
      const auto _interp3 = Z(l1);
      return fma(powr<2>(_interp1), _interp3 * numtracer_kernels::diagfusert_num::tr0(fenv), fma(powr<2>(_interp2), _interp3 * numtracer_kernels::diagfusert_num::tr0(fenv), 0.));
    }

    static inline auto constant(const double& p, const auto& Z, const auto& Gu, const auto& Gd)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::Diagfusert_num_kernel;