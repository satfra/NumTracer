#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#ifndef NT_TRACE_COMPLEX
#define NT_TRACE_COMPLEX DiFfRG::complex<double>
#endif
#include "Zq_num_kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  class Zq_num_kernel
  {
    public:
    #if NT_ZQ_NUM_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zq_num::nenv) > 0 ? (DiFfRG::zq_num::nenv) : 1];
      const double dr_0 = Mq(l1);
      const double dr_1 = sqrt(powr<-1>(powr<2>(l1)));
      const double dr_2 = RF(powr<2>(k), powr<2>(l1));
      const double dr_3 = Zq(k);
      const double dr_4 = Zq(l1);
      const double dr_5 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_7 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const double dr_8 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      DiFfRG::zq_num::fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8);
      const auto _interp1 = ntRe(DiFfRG::zq_num::tr0(fenv));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = ZA(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp5 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp6 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp7 = Zq(k);
      const auto _interp8 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp9 = dtZq(k);
      const auto _interp10 = Zq(1.02 * k);
      const auto _interp11 = Mq(l1);
      const auto _interp12 = Zq(l1);
      const auto _den1 = powr<-1>(_interp2 * _interp3 + _interp4 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den2 = powr<-2>(powr<2>(_interp11) + powr<2>(_interp7 * _interp8 + _interp12 * l1));
      const auto _cse1 = powr<-1>(l1);
      return -0.0833333333333333 * _den1 * _den2 * _interp1 * powr<2>(_interp5) * powr<-2>(p) * fma(-1., _cse1 * _interp6 * _interp7, fma(-1., _cse1 * _interp8 * (50. * (_interp10 - _interp7) + _interp9), 0.));
    }
    #elif NT_ZQ_NUM_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zq_num::nenv) > 0 ? (DiFfRG::zq_num::nenv) : 1];
      const double dr_0 = Mq(l1);
      const double dr_1 = sqrt(powr<-1>(powr<2>(l1)));
      const double dr_2 = RF(powr<2>(k), powr<2>(l1));
      const double dr_3 = Zq(k);
      const double dr_4 = Zq(l1);
      const double dr_5 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_7 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const double dr_8 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      DiFfRG::zq_num::fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8);
      const auto _interp1 = ntRe(DiFfRG::zq_num::tr0(fenv));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = ZA(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp5 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp6 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp7 = Zq(k);
      const auto _interp8 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp9 = dtZq(k);
      const auto _interp10 = Zq(1.02 * k);
      const auto _interp11 = Mq(l1);
      const auto _interp12 = Zq(l1);
      const auto _interp13 = ntIm(DiFfRG::zq_num::tr1(fenv));
      const auto _interp14 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp15 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp16 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp17 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp18 = ZA(l1);
      const auto _interp19 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp20 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp21 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp15 * _interp3 + _interp18 * powr<2>(l1));
      const auto _den3 = powr<-1>(_interp2 * _interp3 + _interp4 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den4 = powr<-2>(powr<2>(_interp11) + powr<2>(_interp7 * _interp8 + _interp12 * l1));
      const auto _den5 = powr<-1>(powr<2>(_interp19) + powr<2>(_interp20 * _interp7 + _interp21 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _cse1 = powr<-1>(l1);
      const auto _cse2 = powr<-2>(p);
      return fma(-0.0833333333333333, _cse2 * _den3 * _den4 * _interp1 * powr<2>(_interp5) * (-_cse1 * _interp6 * _interp7 - _cse1 * _interp8 * (50. * (_interp10 - _interp7) + _interp9)), fma(-0.0833333333333333, _cse2 * _den2 * _den5 * _interp13 * powr<2>(_interp5) * (_interp14 * _interp3 + _interp15 * (_interp16 + 50. * _den1 * (_interp17 - _interp3) * powr<6>(k))), 0.));
    }
    #else                              // the imaginary part survives: genuinely complex
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zq_num::nenv) > 0 ? (DiFfRG::zq_num::nenv) : 1];
      const double dr_0 = Mq(l1);
      const double dr_1 = sqrt(powr<-1>(powr<2>(l1)));
      const double dr_2 = RF(powr<2>(k), powr<2>(l1));
      const double dr_3 = Zq(k);
      const double dr_4 = Zq(l1);
      const double dr_5 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_7 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const double dr_8 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      DiFfRG::zq_num::fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8);
      const auto _interp1 = RB(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = ZA(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp4 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp5 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp6 = Zq(k);
      const auto _interp7 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp8 = dtZq(k);
      const auto _interp9 = Zq(1.02 * k);
      const auto _interp10 = Mq(l1);
      const auto _interp11 = Zq(l1);
      const auto _interp12 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp13 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp14 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp15 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp16 = ZA(l1);
      const auto _interp17 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp18 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp19 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp13 * _interp2 + _interp16 * powr<2>(l1));
      const auto _den3 = powr<-1>(_interp1 * _interp2 + _interp3 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den4 = powr<-2>(powr<2>(_interp10) + powr<2>(_interp6 * _interp7 + _interp11 * l1));
      const auto _den5 = powr<-1>(powr<2>(_interp17) + powr<2>(_interp18 * _interp6 + _interp19 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _cse1 = powr<-1>(l1);
      const auto _cse2 = powr<-2>(p);
      return fma(-0.0833333333333333, _cse2 * _den3 * _den4 * DiFfRG::zq_num::tr0(fenv) * powr<2>(_interp4) * (-_cse1 * _interp5 * _interp6 - _cse1 * _interp7 * (_interp8 + 50. * (-_interp6 + _interp9))), fma(complex<double>(0.,0.0833333333333333), _cse2 * _den2 * _den5 * DiFfRG::zq_num::tr1(fenv) * powr<2>(_interp4) * (_interp12 * _interp2 + _interp13 * (_interp14 + 50. * _den1 * (_interp15 - _interp2) * powr<6>(k))), 0.));
    }
    #endif

    static inline auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq)
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
using DiFfRG::Zq_num_kernel;