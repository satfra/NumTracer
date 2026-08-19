#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Kaon_tadpole_fund_kernels.hh"

namespace DiFfRG
{
  class Kaon_tadpole_fund_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GK, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GKdot)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::kaon_tadpole_fund::nenv) > 0 ? (DiFfRG::kaon_tadpole_fund::nenv) : 1];
      DiFfRG::kaon_tadpole_fund::fill(fenv, l1, cos1, p);
      const auto _interp1 = GK(sqrt(powr<2>(l1)));
      const auto _interp2 = GKdot(sqrt(powr<2>(l1)));
      return -0.5 * DiFfRG::kaon_tadpole_fund::tr0(fenv) * powr<2>(_interp1) * _interp2;
    }

    static inline auto constant(const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GK, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GKdot)
    {
      return 0.;
    }
  };
}
using DiFfRG::Kaon_tadpole_fund_kernel;