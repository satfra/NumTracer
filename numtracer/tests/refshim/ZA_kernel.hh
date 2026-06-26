#pragma once

#include "shim.hpp"

namespace DiFfRG
{
  template <typename _Regulator> class ZA_kernel
  {
  public:
    using Regulator = _Regulator;

    // clang-format off
static KOKKOS_FORCEINLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const double& Nf, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
{
using namespace DiFfRG;
      // clang-format on
      using namespace DiFfRG::compute;
      const auto _interp1 = dtZA(pow(1. + powr<6>(k), 0.16666666666666666667));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp3 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp4 = ZA(pow(1. + powr<6>(k), 0.16666666666666666667));
      const auto _interp5 = ZA(1.02 * pow(1. + powr<6>(k), 0.16666666666666666667));
      const auto _interp6 = ZA(l1);
      const auto _interp7 = RB(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp8 = ZA(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp9 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp10 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + powr<2>(p)));
      const auto _interp11 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp12 = dtZc(k);
      const auto _interp13 = Zc(k);
      const auto _interp14 = Zc(1.02 * k);
      const auto _interp15 = Zc(l1);
      const auto _interp16 = Zc(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp17 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp18 = dtZq(k);
      const auto _interp19 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp20 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp21 = Zq(k);
      const auto _interp22 = Zq(1.02 * k);
      const auto _interp23 = Mq(l1);
      const auto _interp24 = Zq(l1);
      const auto _interp25 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp26 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp27 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _cse1 = powr<2>(l1);
      const auto _cse2 = -2. * cos1 * l1 * p;
      const auto _cse3 = powr<2>(p);
      const auto _cse4 = _cse1 + _cse2 + _cse3;
      const auto _cse5 = _interp19 * _interp21;
      const auto _cse6 = _interp24 * l1;
      const auto _cse7 = _cse5 + _cse6;
      const auto _cse8 = powr<2>(_interp23);
      const auto _cse9 = powr<2>(_cse7);
      const auto _cse10 = sqrt(_cse4);
      const auto _cse11 = powr<-1>(_cse3);
      const auto _cse12 = powr<2>(cos1);
      const auto _cse13 = -6. * _cse12;
      const auto _cse14 = 6. + _cse13;
      const auto _cse15 = 3. * _cse12;
      const auto _cse16 = -3. + _cse15;
      const auto _cse17 = -2. * cos1 * l1 * p;
      const auto _cse18 = _cse1 + _cse17 + _cse3;
      const auto _cse19 = powr<6>(k);
      const auto _cse20 = 1. + _cse19;
      const auto _cse21 = _interp2 * _interp4;
      const auto _cse22 = _cse1 * _interp6;
      const auto _cse23 = _cse21 + _cse22;
      const auto _cse24 = powr<-2>(_cse23); // clang-format off

return fma(-1.,
_cse11 * (7. - _cse12) * _cse24 * _interp10 *
(_interp3 * _interp4 +
_interp2 * (_interp1 + 50. * _cse19 * powr<-1>(_cse20) * (-_interp4 + _interp5))),
fma(2.,
_cse1 * _cse11 * (-1. + _cse12) * powr<2>(_interp11) *
powr<-2>(_cse1 * _interp15 + _interp13 * _interp2) *
(_interp12 * _interp2 + (-50. * _interp13 + 50. * _interp14) * _interp2 +
_interp13 * _interp3) *
powr<-1>(_cse18 * _interp16 + _interp13 * _interp7),
fma(-4.,
_cse11 *powr<-1>(_cse18) * powr<-1>(1. + _cse19) * _cse24 *
((1. + 1. * _cse19) * _interp1 * _interp2 + _cse20 * _interp3 * _interp4 +
_cse19 * _interp2 * (-50. * _interp4 + 50. * _interp5)) *
powr<-1>(_interp4 * _interp7 + _cse18 * _interp8) * powr<2>(_interp9) *
(powr<2>(_cse1) * _cse16 + _cse1 * (-8. + 7. * _cse12 + powr<2>(_cse12)) * _cse3 +
_cse16 * powr<2>(_cse3) + _cse14 * cos1 * powr<3>(l1) * p +
_cse14 * cos1 * l1 * powr<3>(p)),
fma(-1.333333333333333,
_cse11 *powr<-2>(_cse8 + _cse9) * powr<2>(_interp17) *
(-_interp18 * _interp19 - _interp20 * _interp21 +
_interp19 * (50. * _interp21 - 50. * _interp22)) *
powr<-1>(powr<2>(_interp25) + powr<2>(_interp21) * powr<2>(_interp26) +
2. * _cse10 * _interp21 * _interp26 * _interp27 +
_cse4 * powr<2>(_interp27)) *
powr<-1>(l1) * Nf *
(6. * _cse7 * _interp23 * _interp25 * l1 +
sqrt(powr<-1>(_cse4)) * (-_cse8 + _cse9) *
(_interp21 * _interp26 + _cse10 * _interp27) * l1 *
(l1 + 2. * _cse12 * l1 - 3. * cos1 * p)),
0.))));
      // clang-format on
    }

    // clang-format off
static KOKKOS_FORCEINLINE_FUNCTION auto constant(const double& p, const double& k, const double& Nf, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
{
using namespace DiFfRG;
      // clang-format on
      using namespace DiFfRG::compute;
      return 0.;
    }

  private:
    static KOKKOS_FORCEINLINE_FUNCTION auto RB(const auto &k2, const auto &p2) { return Regulator::RB(k2, p2); }

    static KOKKOS_FORCEINLINE_FUNCTION auto RF(const auto &k2, const auto &p2) { return Regulator::RF(k2, p2); }

    static KOKKOS_FORCEINLINE_FUNCTION auto RBdot(const auto &k2, const auto &p2) { return Regulator::RBdot(k2, p2); }

    static KOKKOS_FORCEINLINE_FUNCTION auto RFdot(const auto &k2, const auto &p2) { return Regulator::RFdot(k2, p2); }

    static KOKKOS_FORCEINLINE_FUNCTION auto dq2RB(const auto &k2, const auto &p2) { return Regulator::dq2RB(k2, p2); }

    static KOKKOS_FORCEINLINE_FUNCTION auto dq2RF(const auto &k2, const auto &p2) { return Regulator::dq2RF(k2, p2); }
  };
} // namespace DiFfRG
using DiFfRG::ZA_kernel;