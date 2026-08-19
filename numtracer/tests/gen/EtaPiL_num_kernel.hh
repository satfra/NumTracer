#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#ifndef NT_TRACE_COMPLEX
#define NT_TRACE_COMPLEX DiFfRG::complex<double>
#endif
#include "EtaPiL_num_kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  class EtaPiL_num_kernel
  {
    public:
    #if NT_ETAPIL_NUM_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& k, const double& etaPiL, const double& etaSigL, const double& d1V, const double& d2V, const double& rhoL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hSigL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zSigL)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::etapil_num::nenv) > 0 ? (DiFfRG::etapil_num::nenv) : 1];
      const double dr_0 = 0.7071067811865475;
      const double dr_1 = sqrt(rhoL);
      const double dr_2 = hSigL(l1);
      const double dr_3 = powr<-1>(l1);
      const double dr_4 = RF(powr<2>(k), powr<2>(l1));
      const double dr_5 = Zq(k);
      const double dr_6 = Zq(l1);
      const double dr_7 = sqrt(powr<3>(rhoL));
      const double dr_8 = hPiL(0.816496580927726 * sqrt(1.e-6 - 0.001 * cos1 * l1 + powr<2>(l1)));
      const double dr_9 = hSigL(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      const double dr_10 = hPiL(0.816496580927726 * l1);
      const double dr_11 = rhoL;
      const double dr_12 = sqrt(powr<-1>(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      const double dr_13 = RF(powr<2>(k), 1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
      const double dr_14 = powr<-1>(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
      const double dr_15 = cos1;
      const double dr_16 = l1;
      const double dr_17 = 1.414213562373095;
      const double dr_18 = Zq(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      DiFfRG::etapil_num::fill(fenv, l1, cos1, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18);
      const auto _interp1 = ntRe(DiFfRG::etapil_num::tr0(fenv));
      const auto _interp2 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = Zq(k);
      const auto _interp4 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp5 = dtZq(k);
      const auto _interp6 = Zq(1.02 * k);
      const auto _interp7 = hSigL(l1);
      const auto _interp8 = Zq(l1);
      const auto _interp9 = hSigL(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      const auto _interp10 = RF(powr<2>(k), 1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
      const auto _interp11 = Zq(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      const auto _den1 = powr<-3>(powr<2>(_interp3 * _interp4 + _interp8 * l1) + 0.5 * powr<2>(_interp7) * rhoL);
      const auto _den2 = powr<-1>(powr<2>(_interp10 * _interp3 + _interp11 * sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1))) + 0.5 * powr<2>(_interp9) * rhoL);
      const auto _cse1 = powr<-1>(l1);
      return -1.e6 * _den1 * _den2 * _interp1 * fma(-1., _cse1 * _interp2 * _interp3, fma(-1., _cse1 * _interp4 * (_interp5 + 50. * (-_interp3 + _interp6)), 0.));
    }
    #elif NT_ETAPIL_NUM_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& k, const double& etaPiL, const double& etaSigL, const double& d1V, const double& d2V, const double& rhoL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hSigL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zSigL)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::etapil_num::nenv) > 0 ? (DiFfRG::etapil_num::nenv) : 1];
      const double dr_0 = 0.7071067811865475;
      const double dr_1 = sqrt(rhoL);
      const double dr_2 = hSigL(l1);
      const double dr_3 = powr<-1>(l1);
      const double dr_4 = RF(powr<2>(k), powr<2>(l1));
      const double dr_5 = Zq(k);
      const double dr_6 = Zq(l1);
      const double dr_7 = sqrt(powr<3>(rhoL));
      const double dr_8 = hPiL(0.816496580927726 * sqrt(1.e-6 - 0.001 * cos1 * l1 + powr<2>(l1)));
      const double dr_9 = hSigL(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      const double dr_10 = hPiL(0.816496580927726 * l1);
      const double dr_11 = rhoL;
      const double dr_12 = sqrt(powr<-1>(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      const double dr_13 = RF(powr<2>(k), 1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
      const double dr_14 = powr<-1>(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
      const double dr_15 = cos1;
      const double dr_16 = l1;
      const double dr_17 = 1.414213562373095;
      const double dr_18 = Zq(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      DiFfRG::etapil_num::fill(fenv, l1, cos1, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18);
      const auto _interp1 = ntRe(DiFfRG::etapil_num::tr0(fenv));
      const auto _interp2 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = Zq(k);
      const auto _interp4 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp5 = dtZq(k);
      const auto _interp6 = Zq(1.02 * k);
      const auto _interp7 = hSigL(l1);
      const auto _interp8 = Zq(l1);
      const auto _interp9 = hSigL(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      const auto _interp10 = RF(powr<2>(k), 1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
      const auto _interp11 = Zq(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      const auto _den1 = powr<-3>(powr<2>(_interp3 * _interp4 + _interp8 * l1) + 0.5 * powr<2>(_interp7) * rhoL);
      const auto _den2 = powr<-1>(powr<2>(_interp10 * _interp3 + _interp11 * sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1))) + 0.5 * powr<2>(_interp9) * rhoL);
      const auto _cse1 = powr<-1>(l1);
      return -1.e6 * _den1 * _den2 * _interp1 * fma(-1., _cse1 * _interp2 * _interp3, fma(-1., _cse1 * _interp4 * (_interp5 + 50. * (-_interp3 + _interp6)), 0.));
    }
    #else                              // the imaginary part survives: genuinely complex
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& k, const double& etaPiL, const double& etaSigL, const double& d1V, const double& d2V, const double& rhoL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hSigL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zSigL)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::etapil_num::nenv) > 0 ? (DiFfRG::etapil_num::nenv) : 1];
      const double dr_0 = 0.7071067811865475;
      const double dr_1 = sqrt(rhoL);
      const double dr_2 = hSigL(l1);
      const double dr_3 = powr<-1>(l1);
      const double dr_4 = RF(powr<2>(k), powr<2>(l1));
      const double dr_5 = Zq(k);
      const double dr_6 = Zq(l1);
      const double dr_7 = sqrt(powr<3>(rhoL));
      const double dr_8 = hPiL(0.816496580927726 * sqrt(1.e-6 - 0.001 * cos1 * l1 + powr<2>(l1)));
      const double dr_9 = hSigL(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      const double dr_10 = hPiL(0.816496580927726 * l1);
      const double dr_11 = rhoL;
      const double dr_12 = sqrt(powr<-1>(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      const double dr_13 = RF(powr<2>(k), 1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
      const double dr_14 = powr<-1>(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
      const double dr_15 = cos1;
      const double dr_16 = l1;
      const double dr_17 = 1.414213562373095;
      const double dr_18 = Zq(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      DiFfRG::etapil_num::fill(fenv, l1, cos1, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18);
      const auto _interp1 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp2 = Zq(k);
      const auto _interp3 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp4 = dtZq(k);
      const auto _interp5 = Zq(1.02 * k);
      const auto _interp6 = hSigL(l1);
      const auto _interp7 = Zq(l1);
      const auto _interp8 = hSigL(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      const auto _interp9 = RF(powr<2>(k), 1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
      const auto _interp10 = Zq(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
      const auto _den1 = powr<-3>(powr<2>(_interp2 * _interp3 + _interp7 * l1) + 0.5 * powr<2>(_interp6) * rhoL);
      const auto _den2 = powr<-1>(powr<2>(_interp2 * _interp9 + _interp10 * sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1))) + 0.5 * powr<2>(_interp8) * rhoL);
      const auto _cse1 = powr<-1>(l1);
      return -1.e6 * _den1 * _den2 * DiFfRG::etapil_num::tr0(fenv) * fma(-1., _cse1 * _interp1 * _interp2, fma(-1., _cse1 * _interp3 * (_interp4 + 50. * (-_interp2 + _interp5)), 0.));
    }
    #endif

    // clang-format off
    static inline auto constant(const double& k, const double& etaPiL, const double& etaSigL, const double& d1V, const double& d2V, const double& rhoL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hSigL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zSigL)
    // clang-format on
    {
      return 0.;
    }
    private:
    static inline double ntRe(double x) { return x; }
    template <class T> static inline double ntRe(const T &z) { return z.real(); }
    static inline double ntIm(double) { return 0.0; }
    template <class T> static inline double ntIm(const T &z) { return z.imag(); }
  };
}
using DiFfRG::EtaPiL_num_kernel;