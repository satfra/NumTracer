#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "numtracer/expr/real_cse.hpp"
#include "ZA4_num_kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class ZA4_num_kernel
  {
    public:
    // clang-format off
static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
{
using namespace DiFfRG;
    // clang-format on

    using namespace DiFfRG::compute;
    using namespace numtracer;
    const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
    const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
    const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
    const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
    double fenv[DiFfRG::za4_num::nenv];
    const double dr_0 = Mq(l1);
    const double dr_1 = -powr<-1>(l1) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
    const double dr_2 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const double dr_4 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
    const double dr_5 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
    const double dr_6 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
    const double dr_7 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
    DiFfRG::za4_num::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7);
    const auto _interp1 = ntRe(DiFfRG::za4_num::tr0(fenv));
    const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
    const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
    const auto _interp4 = RB(powr<2>(k), powr<2>(l1));
    const auto _interp5 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
    const auto _interp6 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
    const auto _interp7 = ZA(l1);
    const auto _interp8 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
    const auto _interp9 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp10 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
    const auto _interp11 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
    const auto _interp12 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
    const auto _interp13 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
    const auto _interp14 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp15 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
    const auto _interp16 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
    const auto _interp17 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
    const auto _interp18 = ntRe(DiFfRG::za4_num::tr1(fenv));
    const auto _interp19 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
    const auto _interp20 = ntRe(DiFfRG::za4_num::tr3(fenv));
    const auto _interp21 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
    const auto _interp22 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
    const auto _interp23 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
    const auto _interp24 = ntRe(DiFfRG::za4_num::tr2(fenv));
    const auto _interp25 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
    const auto _interp26 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
    const auto _interp27 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
    const auto _interp28 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
    const auto _interp29 = ntRe(DiFfRG::za4_num::tr4(fenv));
    const auto _interp30 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp31 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
    const auto _interp32 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
    const auto _interp33 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
    const auto _interp34 = Zc(k);
    const auto _interp35 = dtZc(k);
    const auto _interp36 = Zc(1.02 * k);
    const auto _interp37 = Zc(l1);
    const auto _interp38 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp39 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
    const auto _interp40 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
    const auto _interp41 = ntIm(DiFfRG::za4_num::tr5(fenv));
    const auto _interp42 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp43 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
    const auto _interp44 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
    const auto _interp45 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
    const auto _interp46 = RFdot(powr<2>(k), powr<2>(l1));
    const auto _interp47 = Zq(k);
    const auto _interp48 = RF(powr<2>(k), powr<2>(l1));
    const auto _interp49 = dtZq(k);
    const auto _interp50 = Zq(1.02 * k);
    const auto _interp51 = Mq(l1);
    const auto _interp52 = Zq(l1);
    const auto _interp53 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp54 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
    const auto _interp55 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp56 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
    const auto _interp57 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
    const auto _interp58 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
    const auto _interp59 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
    const auto _interp60 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
    const auto _interp61 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
    const auto _den1 = powr<-1>(1. + powr<6>(k));
    const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
    const auto _den3 = powr<-2>(_interp34 * _interp4 + _interp37 * powr<2>(l1));
    const auto _den4 = powr<-1>(_interp3 * _interp8 + _interp9 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _den5 = powr<-1>(_interp25 * _interp3 + _interp26 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
    const auto _den6 = -powr<-1>(_interp34 * _interp8 + _interp38 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _den7 = powr<-2>(powr<2>(_interp51) + powr<2>(_interp47 * _interp48 + _interp52 * l1));
    const auto _den8 = 3. * powr<-1>(3. * _interp12 * _interp3 + _interp13 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
    const auto _den9 = 3. * powr<-1>(3. * _interp21 * _interp3 + _interp22 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
    const auto _den10 = powr<-1>(_interp10 * _interp3 + _interp11 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
    const auto _den11 = -3. * powr<-1>(3. * _interp12 * _interp34 + _interp40 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
    const auto _den12 = -powr<-1>(_interp10 * _interp34 + _interp39 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
    const auto _den13 = powr<-1>(powr<2>(_interp53) + powr<2>(_interp47 * _interp54 + _interp55 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
    const auto _den14 = powr<-1>(powr<2>(_interp59) + powr<2>(_interp47 * _interp60 + _interp61 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p))));
    const auto _den15 = powr<-1>(powr<2>(_interp56) + powr<2>(_interp58) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)) + _interp47 * _interp57 * (_interp47 * _interp57 + 2. * _interp58 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p))));// clang-format off

    return fma(-0.01020408163265306, _den11 * _den12 * _den3 * _den6 * _interp29 * _interp30 * _interp31 * _interp32 * _interp33 * (_interp2 * _interp34 + (_interp35 + 50. * (-_interp34 + _interp36)) * _interp4), fma(0.00510204081632653, _den10 * _den2 * _den4 * _den8 * _interp1 * _interp14 * _interp15 * _interp16 * _interp17 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.00510204081632653, _den10 * _den2 * _den8 * _interp15 * _interp17 * _interp18 * _interp19 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.001275510204081632, _den2 * _den9 * _interp20 * powr<2>(_interp23) * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.002551020408163265, _den10 * _den2 * _den5 * _interp15 * _interp24 * _interp27 * _interp28 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(-0.02040816326530612, _den13 * _den14 * _den15 * _den7 * _interp41 * _interp42 * _interp43 * _interp44 * _interp45 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), 0.))))));
    // clang-format on

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
using DiFfRG::ZA4_num_kernel;