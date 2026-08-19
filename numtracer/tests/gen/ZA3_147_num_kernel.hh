#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#ifndef NT_TRACE_COMPLEX
#define NT_TRACE_COMPLEX DiFfRG::complex<double>
#endif
#include "ZA3_147_num_kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  class ZA3_147_num_kernel
  {
    public:
    #if NT_ZA3_147_NUM_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      double fenv[(DiFfRG::za3_147_num::nenv) > 0 ? (DiFfRG::za3_147_num::nenv) : 1];
      const double dr_0 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_1 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_2 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_3 = Mq(l1);
      const double dr_4 = powr<-1>(l1);
      const double dr_5 = RF(powr<2>(k), powr<2>(l1));
      const double dr_6 = Zq(k);
      const double dr_7 = Zq(l1);
      const double dr_8 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_9 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_10 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_11 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_12 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_13 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const double dr_14 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_15 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_16 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_17 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_18 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_19 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_20 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const double dr_21 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      DiFfRG::za3_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20, dr_21);
      const auto _interp1 = ntRe(DiFfRG::za3_147_num::tr0(fenv));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp5 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp7 = ZA(l1);
      const auto _interp8 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp9 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp10 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp11 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp12 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp13 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp14 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp15 = ntRe(DiFfRG::za3_147_num::tr1(fenv));
      const auto _interp16 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
      const auto _interp17 = ntRe(DiFfRG::za3_147_num::tr2(fenv));
      const auto _interp18 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp19 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp20 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp21 = Zc(k);
      const auto _interp22 = dtZc(k);
      const auto _interp23 = Zc(1.02 * k);
      const auto _interp24 = Zc(l1);
      const auto _interp25 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp26 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp27 = ntRe(DiFfRG::za3_147_num::tr3(fenv));
      const auto _interp28 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp29 = Zq(k);
      const auto _interp30 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp31 = dtZq(k);
      const auto _interp32 = Zq(1.02 * k);
      const auto _interp33 = Mq(l1);
      const auto _interp34 = Zq(l1);
      const auto _interp35 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp36 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp37 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp38 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp39 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp40 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
      const auto _den3 = powr<-2>(_interp21 * _interp4 + _interp24 * powr<2>(l1));
      const auto _den4 = powr<-1>(_interp3 * _interp8 + _interp9 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den5 = -powr<-1>(_interp21 * _interp8 + _interp25 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den6 = powr<-2>(powr<2>(_interp33) + powr<2>(_interp29 * _interp30 + _interp34 * l1));
      const auto _den7 = powr<-1>(_interp10 * _interp3 + _interp11 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den8 = -powr<-1>(_interp10 * _interp21 + _interp26 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den9 = powr<-1>(powr<2>(_interp35) + powr<2>(_interp29 * _interp36 + _interp37 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
      const auto _den10 = powr<-1>(powr<2>(_interp38) + powr<2>(_interp40) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp29 * _interp39 * (_interp29 * _interp39 + 2. * _interp40 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
      const auto _cse1 = powr<-1>(l1);
      const auto _cse2 = powr<-2>(p);
      const auto _cse3 = -_interp3;
      const auto _cse4 = _cse3 + _interp6;
      const auto _cse5 = 50. * _cse4 * _den1 * powr<6>(k);
      const auto _cse6 = _cse5 + _interp5;
      const auto _cse7 = _cse6 * _interp4;
      const auto _cse8 = _interp2 * _interp3;
      const auto _cse9 = _cse7 + _cse8;
      return fma(-0.005050505050505051, _cse2 * _cse9 * _den2 * _den4 * _den7 * _interp1 * _interp12 * _interp13 * _interp14, fma(-0.005050505050505051, _cse2 * _cse9 * _den2 * _den7 * _interp13 * _interp15 * _interp16, fma(0.0202020202020202, _cse2 * _den10 * _den6 * _den9 * _interp27 * (-_cse1 * _interp28 * _interp29 - _cse1 * _interp30 * (_interp31 + 50. * (-_interp29 + _interp32))), fma(-0.0101010101010101, _cse2 * _den3 * _den5 * _den8 * _interp17 * _interp18 * _interp19 * _interp20 * (_interp2 * _interp21 + (_interp22 + 50. * (-_interp21 + _interp23)) * _interp4), 0.))));
    }
    #elif NT_ZA3_147_NUM_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      double fenv[(DiFfRG::za3_147_num::nenv) > 0 ? (DiFfRG::za3_147_num::nenv) : 1];
      const double dr_0 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_1 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_2 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_3 = Mq(l1);
      const double dr_4 = powr<-1>(l1);
      const double dr_5 = RF(powr<2>(k), powr<2>(l1));
      const double dr_6 = Zq(k);
      const double dr_7 = Zq(l1);
      const double dr_8 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_9 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_10 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_11 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_12 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_13 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const double dr_14 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_15 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_16 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_17 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_18 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_19 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_20 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const double dr_21 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      DiFfRG::za3_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20, dr_21);
      const auto _interp1 = ntRe(DiFfRG::za3_147_num::tr0(fenv));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp5 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp7 = ZA(l1);
      const auto _interp8 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp9 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp10 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp11 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp12 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp13 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp14 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp15 = ntRe(DiFfRG::za3_147_num::tr1(fenv));
      const auto _interp16 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
      const auto _interp17 = ntRe(DiFfRG::za3_147_num::tr2(fenv));
      const auto _interp18 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp19 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp20 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp21 = Zc(k);
      const auto _interp22 = dtZc(k);
      const auto _interp23 = Zc(1.02 * k);
      const auto _interp24 = Zc(l1);
      const auto _interp25 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp26 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp27 = ntRe(DiFfRG::za3_147_num::tr3(fenv));
      const auto _interp28 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp29 = Zq(k);
      const auto _interp30 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp31 = dtZq(k);
      const auto _interp32 = Zq(1.02 * k);
      const auto _interp33 = Mq(l1);
      const auto _interp34 = Zq(l1);
      const auto _interp35 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp36 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp37 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp38 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp39 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp40 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
      const auto _den3 = powr<-2>(_interp21 * _interp4 + _interp24 * powr<2>(l1));
      const auto _den4 = powr<-1>(_interp3 * _interp8 + _interp9 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den5 = -powr<-1>(_interp21 * _interp8 + _interp25 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den6 = powr<-2>(powr<2>(_interp33) + powr<2>(_interp29 * _interp30 + _interp34 * l1));
      const auto _den7 = powr<-1>(_interp10 * _interp3 + _interp11 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den8 = -powr<-1>(_interp10 * _interp21 + _interp26 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den9 = powr<-1>(powr<2>(_interp35) + powr<2>(_interp29 * _interp36 + _interp37 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
      const auto _den10 = powr<-1>(powr<2>(_interp38) + powr<2>(_interp40) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp29 * _interp39 * (_interp29 * _interp39 + 2. * _interp40 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
      const auto _cse1 = powr<-1>(l1);
      const auto _cse2 = powr<-2>(p);
      const auto _cse3 = -_interp3;
      const auto _cse4 = _cse3 + _interp6;
      const auto _cse5 = 50. * _cse4 * _den1 * powr<6>(k);
      const auto _cse6 = _cse5 + _interp5;
      const auto _cse7 = _cse6 * _interp4;
      const auto _cse8 = _interp2 * _interp3;
      const auto _cse9 = _cse7 + _cse8;
      return fma(-0.005050505050505051, _cse2 * _cse9 * _den2 * _den4 * _den7 * _interp1 * _interp12 * _interp13 * _interp14, fma(-0.005050505050505051, _cse2 * _cse9 * _den2 * _den7 * _interp13 * _interp15 * _interp16, fma(0.0202020202020202, _cse2 * _den10 * _den6 * _den9 * _interp27 * (-_cse1 * _interp28 * _interp29 - _cse1 * _interp30 * (_interp31 + 50. * (-_interp29 + _interp32))), fma(-0.0101010101010101, _cse2 * _den3 * _den5 * _den8 * _interp17 * _interp18 * _interp19 * _interp20 * (_interp2 * _interp21 + (_interp22 + 50. * (-_interp21 + _interp23)) * _interp4), 0.))));
    }
    #else                              // the imaginary part survives: genuinely complex
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      double fenv[(DiFfRG::za3_147_num::nenv) > 0 ? (DiFfRG::za3_147_num::nenv) : 1];
      const double dr_0 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_1 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_2 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_3 = Mq(l1);
      const double dr_4 = powr<-1>(l1);
      const double dr_5 = RF(powr<2>(k), powr<2>(l1));
      const double dr_6 = Zq(k);
      const double dr_7 = Zq(l1);
      const double dr_8 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_9 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_10 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_11 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_12 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_13 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const double dr_14 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_15 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_16 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_17 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_18 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_19 = sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_20 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const double dr_21 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      DiFfRG::za3_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20, dr_21);
      const auto _interp1 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp4 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp5 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(l1);
      const auto _interp7 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp8 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp9 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp10 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp11 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp12 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp13 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp14 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
      const auto _interp15 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp16 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp17 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp18 = Zc(k);
      const auto _interp19 = dtZc(k);
      const auto _interp20 = Zc(1.02 * k);
      const auto _interp21 = Zc(l1);
      const auto _interp22 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp23 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp24 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp25 = Zq(k);
      const auto _interp26 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp27 = dtZq(k);
      const auto _interp28 = Zq(1.02 * k);
      const auto _interp29 = Mq(l1);
      const auto _interp30 = Zq(l1);
      const auto _interp31 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp32 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp33 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp34 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp35 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp36 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp2 * _interp3 + _interp6 * powr<2>(l1));
      const auto _den3 = powr<-2>(_interp18 * _interp3 + _interp21 * powr<2>(l1));
      const auto _den4 = powr<-1>(_interp2 * _interp7 + _interp8 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den5 = -powr<-1>(_interp18 * _interp7 + _interp22 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den6 = powr<-2>(powr<2>(_interp29) + powr<2>(_interp25 * _interp26 + _interp30 * l1));
      const auto _den7 = powr<-1>(_interp2 * _interp9 + _interp10 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den8 = -powr<-1>(_interp18 * _interp9 + _interp23 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den9 = powr<-1>(powr<2>(_interp31) + powr<2>(_interp25 * _interp32 + _interp33 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
      const auto _den10 = powr<-1>(powr<2>(_interp34) + powr<2>(_interp36) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp25 * _interp35 * (_interp25 * _interp35 + 2. * _interp36 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
      const auto _cse1 = powr<-1>(l1);
      const auto _cse2 = powr<-2>(p);
      const auto _cse3 = -_interp2;
      const auto _cse4 = _cse3 + _interp5;
      const auto _cse5 = 50. * _cse4 * _den1 * powr<6>(k);
      const auto _cse6 = _cse5 + _interp4;
      const auto _cse7 = _cse6 * _interp3;
      const auto _cse8 = _interp1 * _interp2;
      const auto _cse9 = _cse7 + _cse8;
      return fma(-0.005050505050505051, _cse2 * _cse9 * _den2 * _den4 * _den7 * DiFfRG::za3_147_num::tr0(fenv) * _interp11 * _interp12 * _interp13, fma(-0.005050505050505051, _cse2 * _cse9 * _den2 * _den7 * DiFfRG::za3_147_num::tr1(fenv) * _interp12 * _interp14, fma(0.0202020202020202, _cse2 * _den10 * _den6 * _den9 * DiFfRG::za3_147_num::tr3(fenv) * (-_cse1 * _interp24 * _interp25 - _cse1 * _interp26 * (_interp27 + 50. * (-_interp25 + _interp28))), fma(-0.0101010101010101, _cse2 * _den3 * _den5 * _den8 * DiFfRG::za3_147_num::tr2(fenv) * _interp15 * _interp16 * _interp17 * (_interp1 * _interp18 + (_interp19 + 50. * (-_interp18 + _interp20)) * _interp3), 0.))));
    }
    #endif

    // clang-format off
    static inline auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
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
using DiFfRG::ZA3_147_num_kernel;