#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "numtracer/expr/real_cse.hpp"
#include "Gluon_condensate_diag_kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class Gluon_condensate_diag_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const std::array<SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>, 8>& ZA, const std::array<SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>, 8>& ZAdot)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[DiFfRG::gluon_condensate_diag::nenv];
      DiFfRG::gluon_condensate_diag::fill(fenv, l1, cos1, p);
      const auto _interp1 = ntDiagDress(ZA, 0., sqrt(powr<2>(l1)));
      const auto _interp2 = ntDiagDress(ZA, 0., sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp3 = ntDiagDress(ZAdot, 0., sqrt(powr<2>(l1)));
      const auto _interp4 = ntDiagDress(ZA, 1., sqrt(powr<2>(l1)));
      const auto _interp5 = ntDiagDress(ZA, 1., sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp6 = ntDiagDress(ZAdot, 1., sqrt(powr<2>(l1)));
      const auto _interp7 = ntDiagDress(ZA, 2., sqrt(powr<2>(l1)));
      const auto _interp8 = ntDiagDress(ZA, 2., sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp9 = ntDiagDress(ZAdot, 2., sqrt(powr<2>(l1)));
      const auto _interp10 = ntDiagDress(ZA, 3., sqrt(powr<2>(l1)));
      const auto _interp11 = ntDiagDress(ZA, 3., sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp12 = ntDiagDress(ZAdot, 3., sqrt(powr<2>(l1)));
      const auto _interp13 = ntDiagDress(ZA, 4., sqrt(powr<2>(l1)));
      const auto _interp14 = ntDiagDress(ZA, 4., sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp15 = ntDiagDress(ZAdot, 4., sqrt(powr<2>(l1)));
      const auto _interp16 = ntDiagDress(ZA, 5., sqrt(powr<2>(l1)));
      const auto _interp17 = ntDiagDress(ZA, 5., sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp18 = ntDiagDress(ZAdot, 5., sqrt(powr<2>(l1)));
      const auto _interp19 = ntDiagDress(ZA, 6., sqrt(powr<2>(l1)));
      const auto _interp20 = ntDiagDress(ZA, 6., sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp21 = ntDiagDress(ZAdot, 6., sqrt(powr<2>(l1)));
      const auto _interp22 = ntDiagDress(ZA, 7., sqrt(powr<2>(l1)));
      const auto _interp23 = ntDiagDress(ZA, 7., sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp24 = ntDiagDress(ZAdot, 7., sqrt(powr<2>(l1)));
      return 6. * fma(DiFfRG::gluon_condensate_diag::tr0(fenv), powr<2>(_interp10) * _interp11 * _interp12, fma(DiFfRG::gluon_condensate_diag::tr0(fenv), powr<2>(_interp13) * _interp14 * _interp15, fma(DiFfRG::gluon_condensate_diag::tr0(fenv), powr<2>(_interp16) * _interp17 * _interp18, fma(DiFfRG::gluon_condensate_diag::tr0(fenv), powr<2>(_interp19) * _interp20 * _interp21, fma(DiFfRG::gluon_condensate_diag::tr0(fenv), powr<2>(_interp22) * _interp23 * _interp24, fma(DiFfRG::gluon_condensate_diag::tr0(fenv), powr<2>(_interp1) * _interp2 * _interp3, fma(DiFfRG::gluon_condensate_diag::tr0(fenv), powr<2>(_interp4) * _interp5 * _interp6, fma(DiFfRG::gluon_condensate_diag::tr0(fenv), powr<2>(_interp7) * _interp8 * _interp9, 0.))))))));
    }

    static inline auto constant(const double& p, const std::array<SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>, 8>& ZA, const std::array<SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>, 8>& ZAdot)
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
using DiFfRG::Gluon_condensate_diag_kernel;