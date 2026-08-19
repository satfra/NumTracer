#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#ifndef NT_TRACE_COMPLEX
#define NT_TRACE_COMPLEX DiFfRG::complex<double>
#endif
#include "ZAqbq1_num_kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  class ZAqbq1_num_kernel
  {
    public:
    #if NT_ZAQBQ1_NUM_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      double fenv[(DiFfRG::zaqbq1_num::nenv) > 0 ? (DiFfRG::zaqbq1_num::nenv) : 1];
      const double dr_0 = Mq(l1);
      const double dr_1 = powr<-1>(l1);
      const double dr_2 = RF(powr<2>(k), powr<2>(l1));
      const double dr_3 = Zq(k);
      const double dr_4 = Zq(l1);
      const double dr_5 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_7 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const double dr_8 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_9 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_10 = sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_11 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const double dr_12 = Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_13 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_14 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_15 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const double dr_16 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_17 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_18 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_19 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const double dr_20 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      DiFfRG::zaqbq1_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20);
      const auto _interp1 = ntRe(DiFfRG::zaqbq1_num::tr4(fenv));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp5 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp7 = ZA(l1);
      const auto _interp8 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp9 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp10 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp11 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp12 = ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp13 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp14 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp15 = Zq(k);
      const auto _interp16 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp17 = ntRe(DiFfRG::zaqbq1_num::tr5(fenv));
      const auto _interp18 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp19 = ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp20 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp21 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp22 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp23 = ntRe(DiFfRG::zaqbq1_num::tr3(fenv));
      const auto _interp24 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp25 = ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (-cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp26 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp27 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp28 = Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
      const auto _den3 = powr<-1>(_interp3 * _interp8 + _interp9 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den4 = powr<-1>(powr<2>(_interp20) + powr<2>(_interp15 * _interp21 + _interp22 * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den5 = powr<-1>(powr<2>(_interp26) + powr<2>(_interp15 * _interp27 + _interp28 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den6 = powr<-1>(powr<2>(_interp13) + powr<2>(_interp16) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp14 * _interp15 * (_interp14 * _interp15 + 2. * _interp16 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
      const auto _cse1 = -_interp3;
      const auto _cse2 = _cse1 + _interp6;
      const auto _cse3 = 50. * _cse2 * _den1 * powr<6>(k);
      const auto _cse4 = _cse3 + _interp5;
      const auto _cse5 = _cse4 * _interp4;
      const auto _cse6 = _interp2 * _interp3;
      const auto _cse7 = _cse5 + _cse6;
      return 0.02083333333333333 * fma(-1., _cse7 * _den2 * _den3 * _den6 * _interp1 * _interp10 * _interp11 * _interp12, fma(_cse7, _den2 * _den4 * _den6 * _interp11 * _interp17 * _interp18 * _interp19, fma(-1., _cse7 * _den2 * _den3 * _den5 * _interp10 * _interp23 * _interp24 * _interp25, 0.)));
    }
    #elif NT_ZAQBQ1_NUM_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      double fenv[(DiFfRG::zaqbq1_num::nenv) > 0 ? (DiFfRG::zaqbq1_num::nenv) : 1];
      const double dr_0 = Mq(l1);
      const double dr_1 = powr<-1>(l1);
      const double dr_2 = RF(powr<2>(k), powr<2>(l1));
      const double dr_3 = Zq(k);
      const double dr_4 = Zq(l1);
      const double dr_5 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_7 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const double dr_8 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_9 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_10 = sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_11 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const double dr_12 = Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_13 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_14 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_15 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const double dr_16 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_17 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_18 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_19 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const double dr_20 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      DiFfRG::zaqbq1_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20);
      const auto _interp1 = ntIm(DiFfRG::zaqbq1_num::tr0(fenv));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp5 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp6 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp7 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp8 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp9 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp10 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp11 = Zq(k);
      const auto _interp12 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp13 = dtZq(k);
      const auto _interp14 = Zq(1.02 * k);
      const auto _interp15 = Mq(l1);
      const auto _interp16 = Zq(l1);
      const auto _interp17 = ntIm(DiFfRG::zaqbq1_num::tr1(fenv));
      const auto _interp18 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp19 = ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp20 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp21 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp22 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp23 = ntIm(DiFfRG::zaqbq1_num::tr2(fenv));
      const auto _interp24 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp25 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp26 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp27 = ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (-cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp28 = ntRe(DiFfRG::zaqbq1_num::tr4(fenv));
      const auto _interp29 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp30 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp31 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp32 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp33 = ZA(l1);
      const auto _interp34 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp35 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp36 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp37 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp38 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp39 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp40 = ntRe(DiFfRG::zaqbq1_num::tr5(fenv));
      const auto _interp41 = ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp42 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp43 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp44 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp45 = ntRe(DiFfRG::zaqbq1_num::tr3(fenv));
      const auto _interp46 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp47 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp48 = Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp3 * _interp30 + _interp33 * powr<2>(l1));
      const auto _den3 = powr<-1>(_interp3 * _interp34 + _interp35 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den4 = powr<-1>(_interp2 * _interp3 + _interp4 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den5 = powr<-1>(_interp24 * _interp3 + _interp25 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den6 = powr<-2>(powr<2>(_interp15) + powr<2>(_interp11 * _interp12 + _interp16 * l1));
      const auto _den7 = powr<-1>(_interp3 * _interp5 + _interp6 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den8 = powr<-1>(powr<2>(_interp20) + powr<2>(_interp11 * _interp21 + _interp22 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
      const auto _den9 = powr<-1>(powr<2>(_interp42) + powr<2>(_interp11 * _interp43 + _interp44 * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den10 = powr<-1>(powr<2>(_interp46) + powr<2>(_interp11 * _interp47 + _interp48 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den11 = powr<-1>(powr<2>(_interp37) + powr<2>(_interp39) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp11 * _interp38 * (_interp11 * _interp38 + 2. * _interp39 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
      const auto _cse1 = powr<-1>(l1);
      const auto _cse2 = -_cse1 * _interp10 * _interp11;
      const auto _cse3 = -_interp11;
      const auto _cse4 = _cse3 + _interp14;
      const auto _cse5 = 50. * _cse4;
      // clang-format off
      return fma(0.02083333333333333, _den5 * _den6 * _den8 * (_cse2 - _cse1 * _interp12 * (_cse5 + _interp13)) * _interp18 * _interp23 * _interp26 * _interp27, fma(-0.02083333333333333, _den6 * _den7 * _den8 * (_cse2 - _cse1 * _interp12 * (_cse5 + _interp13)) * _interp17 * _interp18 * _interp19 * _interp9, fma(0.02083333333333333, _den4 * _den6 * _den7 * _interp1 * (_cse2 - _cse1 * _interp12 * (_cse5 + _interp13)) * _interp7 * _interp8 * _interp9, fma(-0.02083333333333333, _den10 * _den2 * _den3 * _interp26 * _interp27 * _interp36 * _interp45 * (_interp29 * _interp3 + _interp30 * (_interp31 + 50. * _den1 * (-_interp3 + _interp32) * powr<6>(k))), fma(-0.02083333333333333, _den11 * _den2 * _den3 * _interp19 * _interp28 * _interp36 * _interp9 * (_interp29 * _interp3 + _interp30 * (_interp31 + 50. * _den1 * (-_interp3 + _interp32) * powr<6>(k))), fma(0.02083333333333333, _den11 * _den2 * _den9 * _interp40 * _interp41 * _interp8 * _interp9 * (_interp29 * _interp3 + _interp30 * (_interp31 + 50. * _den1 * (-_interp3 + _interp32) * powr<6>(k))), 0.))))));
      // clang-format on
    }
    #else                              // the imaginary part survives: genuinely complex
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      double fenv[(DiFfRG::zaqbq1_num::nenv) > 0 ? (DiFfRG::zaqbq1_num::nenv) : 1];
      const double dr_0 = Mq(l1);
      const double dr_1 = powr<-1>(l1);
      const double dr_2 = RF(powr<2>(k), powr<2>(l1));
      const double dr_3 = Zq(k);
      const double dr_4 = Zq(l1);
      const double dr_5 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_7 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const double dr_8 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_9 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_10 = sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_11 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const double dr_12 = Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_13 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_14 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_15 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const double dr_16 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_17 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_18 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_19 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const double dr_20 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      DiFfRG::zaqbq1_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20);
      const auto _interp1 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp5 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp6 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp7 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp8 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp9 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp10 = Zq(k);
      const auto _interp11 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp12 = dtZq(k);
      const auto _interp13 = Zq(1.02 * k);
      const auto _interp14 = Mq(l1);
      const auto _interp15 = Zq(l1);
      const auto _interp16 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp17 = ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp18 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp19 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp20 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp21 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp22 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp23 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp24 = ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (-cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp25 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp26 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp27 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp28 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp29 = ZA(l1);
      const auto _interp30 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp31 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp32 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp33 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp34 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp35 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp36 = ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp37 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp38 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp39 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp40 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp41 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp42 = Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp2 * _interp26 + _interp29 * powr<2>(l1));
      const auto _den3 = powr<-1>(_interp2 * _interp30 + _interp31 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den4 = powr<-1>(_interp1 * _interp2 + _interp3 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den5 = powr<-1>(_interp2 * _interp21 + _interp22 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den6 = powr<-2>(powr<2>(_interp14) + powr<2>(_interp10 * _interp11 + _interp15 * l1));
      const auto _den7 = powr<-1>(_interp2 * _interp4 + _interp5 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den8 = powr<-1>(powr<2>(_interp18) + powr<2>(_interp10 * _interp19 + _interp20 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
      const auto _den9 = powr<-1>(powr<2>(_interp37) + powr<2>(_interp10 * _interp38 + _interp39 * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den10 = powr<-1>(powr<2>(_interp40) + powr<2>(_interp10 * _interp41 + _interp42 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den11 = powr<-1>(powr<2>(_interp33) + powr<2>(_interp35) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp10 * _interp34 * (_interp10 * _interp34 + 2. * _interp35 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
      const auto _cse1 = powr<-1>(l1);
      const auto _cse2 = -_interp10;
      const auto _cse3 = _cse2 + _interp13;
      const auto _cse4 = 50. * _cse3;
      const auto _cse5 = _cse4 + _interp12;
      const auto _cse6 = -_cse1 * _cse5 * _interp11;
      const auto _cse7 = -_cse1 * _interp10 * _interp9;
      const auto _cse8 = _cse6 + _cse7;
      const auto _cse9 = -_interp2;
      const auto _cse10 = _cse9 + _interp28;
      const auto _cse11 = 50. * _cse10 * _den1 * powr<6>(k);
      // clang-format off
      return 0.02083333333333333 * fma(complex<double>(0.,-1.), _cse8 * _den5 * _den6 * _den8 * DiFfRG::zaqbq1_num::tr2(fenv) * _interp16 * _interp23 * _interp24, fma(-1., _den10 * _den2 * _den3 * DiFfRG::zaqbq1_num::tr3(fenv) * _interp2 * _interp23 * _interp24 * _interp25 * _interp32, fma(-1., _cse11 * _den10 * _den2 * _den3 * DiFfRG::zaqbq1_num::tr3(fenv) * _interp23 * _interp24 * _interp26 * _interp32, fma(-1., _den10 * _den2 * _den3 * DiFfRG::zaqbq1_num::tr3(fenv) * _interp23 * _interp24 * _interp26 * _interp27 * _interp32, fma(complex<double>(0.,1.), _cse8 * _den6 * _den7 * _den8 * DiFfRG::zaqbq1_num::tr1(fenv) * _interp16 * _interp17 * _interp8, fma(-1., _den11 * _den2 * _den3 * DiFfRG::zaqbq1_num::tr4(fenv) * _interp17 * _interp2 * _interp25 * _interp32 * _interp8, fma(-1., _cse11 * _den11 * _den2 * _den3 * DiFfRG::zaqbq1_num::tr4(fenv) * _interp17 * _interp26 * _interp32 * _interp8, fma(-1., _den11 * _den2 * _den3 * DiFfRG::zaqbq1_num::tr4(fenv) * _interp17 * _interp26 * _interp27 * _interp32 * _interp8, fma(_den11, _den2 * _den9 * DiFfRG::zaqbq1_num::tr5(fenv) * _interp2 * _interp25 * _interp36 * _interp7 * _interp8, fma(_cse11, _den11 * _den2 * _den9 * DiFfRG::zaqbq1_num::tr5(fenv) * _interp26 * _interp36 * _interp7 * _interp8, fma(_den11, _den2 * _den9 * DiFfRG::zaqbq1_num::tr5(fenv) * _interp26 * _interp27 * _interp36 * _interp7 * _interp8, fma(complex<double>(0.,-1.), _cse8 * _den4 * _den6 * _den7 * DiFfRG::zaqbq1_num::tr0(fenv) * _interp6 * _interp7 * _interp8, 0.))))))))))));
      // clang-format on
    }
    #endif

    // clang-format off
    static inline auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
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
using DiFfRG::ZAqbq1_num_kernel;