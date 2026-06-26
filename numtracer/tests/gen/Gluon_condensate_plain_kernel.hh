#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "numtracer/expr/real_cse.hpp"
#include "Gluon_condensate_plain_kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class Gluon_condensate_plain_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAflat, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAflatDot)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[DiFfRG::gluon_condensate_plain::nenv];
      DiFfRG::gluon_condensate_plain::fill(fenv, l1, cos1, p);
      const auto _interp1 = ZAflat(sqrt(powr<2>(l1)));
      const auto _interp2 = ZAflat(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp3 = ZAflatDot(sqrt(powr<2>(l1)));
      return 2. * DiFfRG::gluon_condensate_plain::tr0(fenv) * powr<2>(_interp1) * _interp2 * _interp3;
    }

    static inline auto constant(const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAflat, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAflatDot)
    {
      return 0.;
    }private: static inline auto RB(const auto &k2, const auto &p2) { return REG::RB(k2, p2); }
    static inline auto RF(const auto &k2, const auto &p2) { return REG::RF(k2, p2); }
    static inline auto RBdot(const auto &k2, const auto &p2) { return REG::RBdot(k2, p2); }
    static inline auto RFdot(const auto &k2, const auto &p2) { return REG::RFdot(k2, p2); }
    static inline auto dq2RB(const auto &k2, const auto &p2) { return REG::dq2RB(k2, p2); }
    static inline auto dq2RF(const auto &k2, const auto &p2) { return REG::dq2RF(k2, p2); }
  };
}
using DiFfRG::Gluon_condensate_plain_kernel;