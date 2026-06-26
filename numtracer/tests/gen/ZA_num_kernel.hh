#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "numtracer/expr/real_cse.hpp"
#include "ZA_num_kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class ZA_num_kernel
  {
    public:
    // clang-format off
static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
{
using namespace DiFfRG;
    // clang-format on

    using namespace DiFfRG::compute;
    using namespace numtracer;
    double fenv[DiFfRG::za_num::nenv];
    const double dr_0 = Mq(l1);
    const double dr_1 = -sqrt(powr<-1>(powr<2>(l1))) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
    const double dr_2 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
    const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
    DiFfRG::za_num::fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3);
    const auto _interp1 = ntRe(DiFfRG::za_num::tr0(fenv));
    const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
    const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
    const auto _interp4 = RB(powr<2>(k), powr<2>(l1));
    const auto _interp5 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
    const auto _interp6 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
    const auto _interp7 = ZA(l1);
    const auto _interp8 = RB(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
    const auto _interp9 = ZA(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
    const auto _interp10 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
    const auto _interp11 = ntRe(DiFfRG::za_num::tr1(fenv));
    const auto _interp12 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + powr<2>(p)));
    const auto _interp13 = ntRe(DiFfRG::za_num::tr2(fenv));
    const auto _interp14 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
    const auto _interp15 = Zc(k);
    const auto _interp16 = dtZc(k);
    const auto _interp17 = Zc(1.02 * k);
    const auto _interp18 = Zc(l1);
    const auto _interp19 = Zc(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
    const auto _interp20 = ntIm(DiFfRG::za_num::tr3(fenv));
    const auto _interp21 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
    const auto _interp22 = RFdot(powr<2>(k), powr<2>(l1));
    const auto _interp23 = Zq(k);
    const auto _interp24 = RF(powr<2>(k), powr<2>(l1));
    const auto _interp25 = dtZq(k);
    const auto _interp26 = Zq(1.02 * k);
    const auto _interp27 = Mq(l1);
    const auto _interp28 = Zq(l1);
    const auto _interp29 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
    const auto _interp30 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
    const auto _interp31 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
    const auto _den1 = powr<-1>(1. + powr<6>(k));
    const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
    const auto _den3 = powr<-2>(_interp15 * _interp4 + _interp18 * powr<2>(l1));
    const auto _den4 = powr<-1>(_interp3 * _interp8 + _interp9 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
    const auto _den5 = -powr<-1>(_interp15 * _interp8 + _interp19 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
    const auto _den6 = powr<-2>(powr<2>(_interp27) + powr<2>(_interp23 * _interp24 + _interp28 * sqrt(powr<2>(l1))));
    const auto _den7 = powr<-1>(powr<2>(_interp29) + powr<2>(_interp23 * _interp30 + _interp31 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
    const auto _cse1 = sqrt(powr<-1>(powr<2>(l1)));
    const auto _cse2 = powr<-2>(p);
    const auto _cse3 = _interp2 * _interp3;
    const auto _cse4 = -_interp3;
    const auto _cse5 = _cse4 + _interp6;
    const auto _cse6 = 50. * _cse5 * _den1 * powr<6>(k);
    const auto _cse7 = _cse6 + _interp5;
    const auto _cse8 = _cse7 * _interp4;
    const auto _cse9 = _cse3 + _cse8;
    return fma(-0.04166666666666666, _cse2 * _cse9 * _den2 * _den4 * _interp1 * powr<2>(_interp10), fma(-0.02083333333333333, _cse2 * _cse9 * _den2 * _interp11 * _interp12, fma(0.1666666666666667, _cse2 * _den6 * _den7 * _interp20 * powr<2>(_interp21) * (-_cse1 * _interp22 * _interp23 - _cse1 * _interp24 * (_interp25 + 50. * (-_interp23 + _interp26))), fma(0.0833333333333333, _cse2 * _den3 * _den5 * _interp13 * powr<2>(_interp14) * (_interp15 * _interp2 + (_interp16 + 50. * (-_interp15 + _interp17)) * _interp4), 0.))));
  }

  // clang-format off
static inline auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
{
return 0.;
  // clang-format on

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
using DiFfRG::ZA_num_kernel;