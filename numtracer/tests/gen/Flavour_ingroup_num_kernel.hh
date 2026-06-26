#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "numtracer/expr/real_cse.hpp"
#include "Flavour_ingroup_num_kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class Flavour_ingroup_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const std::array<SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>, 2>& G, const std::array<SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>, 2>& Gdot)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[DiFfRG::flavour_ingroup_num::nenv];
      DiFfRG::flavour_ingroup_num::fill(fenv, l1, cos1, p);
      const auto _interp1 = ntDiagDress(G, 0., sqrt(powr<2>(l1)));
      const auto _interp2 = ntDiagDress(G, 0., sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp3 = ntDiagDress(Gdot, 0., sqrt(powr<2>(l1)));
      const auto _interp4 = ntDiagDress(G, 1., sqrt(powr<2>(l1)));
      const auto _interp5 = ntDiagDress(G, 1., sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp6 = ntDiagDress(Gdot, 1., sqrt(powr<2>(l1)));
      return 6. * fma(DiFfRG::flavour_ingroup_num::tr0(fenv), powr<2>(_interp1) * _interp2 * _interp3, fma(DiFfRG::flavour_ingroup_num::tr0(fenv), powr<2>(_interp4) * _interp5 * _interp6, 0.));
    }

    static inline auto constant(const double& p, const std::array<SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>, 2>& G, const std::array<SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>, 2>& Gdot)
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
using DiFfRG::Flavour_ingroup_num_kernel;