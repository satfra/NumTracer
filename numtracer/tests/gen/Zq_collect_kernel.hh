#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "numtracer/expr/real_cse.hpp"
#include "Zq_collect_kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class Zq_collect_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[DiFfRG::zq_collect::nenv];
      const double dr_0 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_1 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_2 = Mq(sqrt(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p)));
      DiFfRG::zq_collect::fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3);
      const auto _interp1 = ntIm(DiFfRG::zq_collect::tr0(fenv));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp5 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp7 = ZA(l1);
      const auto _interp8 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp9 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp10 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp11 = Zq(k);
      const auto _interp12 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp13 = ntRe(DiFfRG::zq_collect::tr1(fenv));
      const auto _interp14 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cos1 * l1 * p + powr<2>(p)));
      const auto _interp15 = RFdot(powr<2>(k), powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp16 = RF(powr<2>(k), powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp17 = dtZq(k);
      const auto _interp18 = Zq(1.02 * k);
      const auto _interp19 = Mq(sqrt(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp20 = Zq(sqrt(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
      const auto _den3 = powr<-1>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
      const auto _den4 = powr<-1>(powr<2>(_interp9) + powr<2>(_interp10 * _interp11 + _interp12 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _den5 = powr<-2>(powr<2>(_interp19) + powr<2>(_interp11 * _interp16 + _interp20 * sqrt(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _cse1 = 2. * cos1 * l1 * p;
      const auto _cse4 = powr<-2>(p);
      const auto _cse2 = _cse1 + powr<-1>(_cse4) + powr<2>(l1);
      const auto _cse3 = sqrt(powr<-1>(_cse2));
      return fma(-0.0833333333333333, _cse4 * _den3 * _den5 * _interp13 * powr<2>(_interp14) * (-_cse3 * _interp11 * _interp15 - _cse3 * _interp16 * (_interp17 + 50. * (-_interp11 + _interp18))), fma(-0.0833333333333333, _cse4 * _den2 * _den4 * _interp1 * powr<2>(_interp8) * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), 0.));
    }

    static inline auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq)
    {
      return 0.;
    }private: static inline auto RB(const auto &k2, const auto &p2) { return REG::RB(k2, p2); }
    static inline auto RF(const auto &k2, const auto &p2) { return REG::RF(k2, p2); }
    static inline auto RBdot(const auto &k2, const auto &p2) { return REG::RBdot(k2, p2); }
    static inline auto RFdot(const auto &k2, const auto &p2) { return REG::RFdot(k2, p2); }
    static inline auto dq2RB(const auto &k2, const auto &p2) { return REG::dq2RB(k2, p2); }
    static inline auto dq2RF(const auto &k2, const auto &p2) { return REG::dq2RF(k2, p2); }

    static inline double ntRe(double x) { return x; }
    template <class T> static inline double ntRe(const T &z) { return z.real(); }
    static inline double ntIm(double) { return 0.0; }
    template <class T> static inline double ntIm(const T &z) { return z.imag(); }
  };
}
using DiFfRG::Zq_collect_kernel;