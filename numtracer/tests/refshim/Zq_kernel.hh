#pragma once

// Copied FormTracer-generated quark wave-function (Zq) kernel, adapted to the
// standalone NumTracer shim (the only change vs the DiFfRG source is the include
// line: the two DiFfRG physics headers are replaced by shim.hpp, which provides
// powr/pow/sqrt/fma, the interpolator/coordinate/memory type stubs, and the
// Regulator wrappers). This is the numeric ORACLE for compare_zq.
#include "shim.hpp"

namespace DiFfRG
{
  template <typename _Regulator> class Zq_kernel
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
      const auto _interp7 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp8 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp9 = Zq(k);
      const auto _interp10 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp11 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp12 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cos1 * l1 * p + powr<2>(p)));
      const auto _interp13 = dtZq(k);
      const auto _interp14 = RF(powr<2>(k), powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp15 = RFdot(powr<2>(k), powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp16 = Zq(1.02 * k);
      const auto _interp17 = Mq(sqrt(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp18 = Zq(sqrt(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _cse1 = powr<2>(l1);
      const auto _cse2 = powr<2>(p);
      const auto _cse3 = 2. * cos1 * l1 * p;
      const auto _cse4 = _cse1 + _cse2 + _cse3;
      const auto _cse5 = sqrt(_cse4);
      const auto _cse6 = powr<2>(_interp17);
      const auto _cse7 = 0.5 * _cse1;
      const auto _cse8 = 1. * cos1 * l1 * p;
      const auto _cse9 = 0.5 * _cse2;
      const auto _cse10 = _cse7 + _cse8 + _cse9;
      const auto _cse11 = powr<2>(_interp18);
      const auto _cse12 = _cse10 * _cse11;
      const auto _cse13 = 1. * _cse5 * _interp14 * _interp18 * _interp9;
      const auto _cse14 = powr<2>(_interp14);
      const auto _cse15 = powr<2>(_interp9);
      const auto _cse16 = 0.5 * _cse14 * _cse15;
      const auto _cse17 = powr<-1>(p);
      const auto _cse18 = powr<2>(cos1);
      const auto _cse19 = powr<6>(k);
      const auto _cse20 = _interp2 * _interp4;
      const auto _cse21 = _cse1 * _interp6;
      const auto _cse22 = _cse20 + _cse21;
      const auto _cse23 = -2. * cos1 * l1 * p;
      const auto _cse24 = _cse1 + _cse2 + _cse23;
      const auto _cse25 = sqrt(_cse24);
      const auto _cse26 = -2. * cos1 * l1 * p;
      const auto _cse27 = _cse1 + _cse2 + _cse26; // clang-format off

return fma(1.333333333333333,
_cse17 *powr<-1>(1. + _cse19) * powr<-2>(_cse22) * _cse25 * powr<-1>(_cse27) *
((1. + 1. * _cse19) * _interp1 * _interp2 + (1. + _cse19) * _interp3 * _interp4 +
_cse19 * _interp2 * (-50. * _interp4 + 50. * _interp5)) *
powr<2>(_interp7) * (_cse25 * _interp10 + _interp8 * _interp9) *
powr<-1>(_cse27 * powr<2>(_interp10) + powr<2>(_interp11) + _cse15 * powr<2>(_interp8) +
2. * _cse25 * _interp10 * _interp8 * _interp9) *
(-3. * cos1 * l1 + p + 2. * _cse18 * p),
fma(2.,
_cse17 *powr<-1>(_cse22) * _cse5 * (_cse12 + _cse13 + _cse16 - 0.5 * _cse6) *
powr<-2>(_cse12 + _cse13 + _cse16 + 0.5 * _cse6) * powr<2>(_interp12) *
(_interp13 * _interp14 + _interp14 * (50. * _interp16 - 50. * _interp9) +
_interp15 * _interp9) *
(1. * cos1 * l1 + 0.3333333333333333 * p + 0.6666666666666666 * _cse18 * p) *
powr<-1>(_cse1 + _cse2 + 2. * cos1 * l1 * p),
0.));
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
using DiFfRG::Zq_kernel;
