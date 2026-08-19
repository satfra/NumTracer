#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#ifndef NT_TRACE_COMPLEX
#define NT_TRACE_COMPLEX DiFfRG::complex<double>
#endif
#include "ZAAqbq1SmallRef_num_kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  class ZAAqbq1SmallRef_num_kernel
  {
    public:
    #if NT_ZAAQBQ1_SMALL_REF_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      double fenv[(DiFfRG::zaaqbq1_small_ref::nenv) > 0 ? (DiFfRG::zaaqbq1_small_ref::nenv) : 1];
      const double dr_0 = Mq(l1);
      const double dr_1 = powr<-1>(l1);
      const double dr_2 = RF(powr<2>(k), powr<2>(l1));
      const double dr_3 = Zq(k);
      const double dr_4 = Zq(l1);
      DiFfRG::zaaqbq1_small_ref::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4);
      const auto _interp1 = ntRe(DiFfRG::zaaqbq1_small_ref::tr1(fenv));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp5 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp7 = ZA(l1);
      const auto _interp8 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp9 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp10 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp11 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp12 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp13 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp14 = ZAAqbq1(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp15 = ntRe(DiFfRG::zaaqbq1_small_ref::tr2(fenv));
      const auto _interp16 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp17 = ntRe(DiFfRG::zaaqbq1_small_ref::tr0(fenv));
      const auto _interp18 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp19 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp20 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp21 = ZAAqbq1(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
      const auto _den3 = powr<-1>(_interp3 * _interp8 + _interp9 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den4 = powr<-1>(_interp18 * _interp3 + _interp19 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den5 = 3. * powr<-1>(3. * _interp10 * _interp3 + _interp11 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _cse1 = powr<-2>(p);
      const auto _cse2 = -_interp3;
      const auto _cse3 = _cse2 + _interp6;
      const auto _cse4 = 50. * _cse3 * _den1 * powr<6>(k);
      const auto _cse5 = _cse4 + _interp5;
      const auto _cse6 = _cse5 * _interp4;
      const auto _cse7 = _interp2 * _interp3;
      const auto _cse8 = _cse6 + _cse7;
      return -0.0003662109375 * fma(2., _cse1 * _cse8 * _den2 * _den3 * _den5 * _interp1 * _interp12 * _interp13 * _interp14, fma(_cse1, _cse8 * _den2 * _den5 * _interp14 * _interp15 * _interp16, fma(_cse1, _cse8 * _den2 * _den3 * _den4 * _interp12 * _interp17 * _interp20 * _interp21, 0.)));
    }
    #elif NT_ZAAQBQ1_SMALL_REF_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      double fenv[(DiFfRG::zaaqbq1_small_ref::nenv) > 0 ? (DiFfRG::zaaqbq1_small_ref::nenv) : 1];
      const double dr_0 = Mq(l1);
      const double dr_1 = powr<-1>(l1);
      const double dr_2 = RF(powr<2>(k), powr<2>(l1));
      const double dr_3 = Zq(k);
      const double dr_4 = Zq(l1);
      DiFfRG::zaaqbq1_small_ref::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4);
      const auto _interp1 = ntRe(DiFfRG::zaaqbq1_small_ref::tr1(fenv));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp5 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp7 = ZA(l1);
      const auto _interp8 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp9 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp10 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp11 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp12 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp13 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp14 = ZAAqbq1(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp15 = ntRe(DiFfRG::zaaqbq1_small_ref::tr2(fenv));
      const auto _interp16 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp17 = ntRe(DiFfRG::zaaqbq1_small_ref::tr0(fenv));
      const auto _interp18 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp19 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp20 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp21 = ZAAqbq1(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
      const auto _interp22 = ntIm(DiFfRG::zaaqbq1_small_ref::tr3(fenv));
      const auto _interp23 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _interp24 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp25 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp26 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp27 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * (cosl1p2 + cosl1p3)) * l1 * p + 10. * powr<2>(p)));
      const auto _interp28 = ZAAqbq1(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp29 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp30 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp31 = Zq(k);
      const auto _interp32 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp33 = dtZq(k);
      const auto _interp34 = Zq(1.02 * k);
      const auto _interp35 = Mq(l1);
      const auto _interp36 = Zq(l1);
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
      const auto _den3 = powr<-1>(_interp3 * _interp8 + _interp9 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den4 = powr<-1>(_interp18 * _interp3 + _interp19 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den5 = powr<-2>(powr<2>(_interp35) + powr<2>(_interp31 * _interp32 + _interp36 * l1));
      const auto _den6 = 3. * powr<-1>(3. * _interp10 * _interp3 + _interp11 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den7 = 3. * powr<-1>(3. * _interp25 * _interp3 + _interp26 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den8 = powr<-1>(_interp23 * _interp3 + _interp24 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _cse1 = powr<-1>(l1);
      const auto _cse2 = powr<-2>(p);
      const auto _cse3 = -_interp3;
      const auto _cse4 = _cse3 + _interp6;
      const auto _cse5 = 50. * _cse4 * _den1 * powr<6>(k);
      const auto _cse6 = _cse5 + _interp5;
      const auto _cse7 = _cse6 * _interp4;
      const auto _cse8 = _interp2 * _interp3;
      const auto _cse9 = _cse7 + _cse8;
      return fma(-0.000732421875, _cse2 * _cse9 * _den2 * _den3 * _den6 * _interp1 * _interp12 * _interp13 * _interp14, fma(-0.0003662109375, _cse2 * _cse9 * _den2 * _den6 * _interp14 * _interp15 * _interp16, fma(-0.0003662109375, _cse2 * _cse9 * _den2 * _den3 * _den4 * _interp12 * _interp17 * _interp20 * _interp21, fma(-0.000732421875, _cse2 * _den5 * _den7 * _den8 * _interp22 * _interp27 * _interp28 * _interp29 * (-_cse1 * _interp30 * _interp31 - _cse1 * _interp32 * (_interp33 + 50. * (-_interp31 + _interp34))), 0.))));
    }
    #else                              // the imaginary part survives: genuinely complex
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      double fenv[(DiFfRG::zaaqbq1_small_ref::nenv) > 0 ? (DiFfRG::zaaqbq1_small_ref::nenv) : 1];
      const double dr_0 = Mq(l1);
      const double dr_1 = powr<-1>(l1);
      const double dr_2 = RF(powr<2>(k), powr<2>(l1));
      const double dr_3 = Zq(k);
      const double dr_4 = Zq(l1);
      DiFfRG::zaaqbq1_small_ref::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4);
      const auto _interp1 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp4 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp5 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(l1);
      const auto _interp7 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp8 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp9 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp10 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp11 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp12 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp13 = ZAAqbq1(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp14 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp15 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp16 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp17 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp18 = ZAAqbq1(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
      const auto _interp19 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _interp20 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp21 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp22 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp23 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * (cosl1p2 + cosl1p3)) * l1 * p + 10. * powr<2>(p)));
      const auto _interp24 = ZAAqbq1(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp25 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp26 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp27 = Zq(k);
      const auto _interp28 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp29 = dtZq(k);
      const auto _interp30 = Zq(1.02 * k);
      const auto _interp31 = Mq(l1);
      const auto _interp32 = Zq(l1);
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp2 * _interp3 + _interp6 * powr<2>(l1));
      const auto _den3 = powr<-1>(_interp2 * _interp7 + _interp8 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den4 = powr<-1>(_interp15 * _interp2 + _interp16 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den5 = powr<-2>(powr<2>(_interp31) + powr<2>(_interp27 * _interp28 + _interp32 * l1));
      const auto _den6 = 3. * powr<-1>(3. * _interp2 * _interp9 + _interp10 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den7 = 3. * powr<-1>(3. * _interp2 * _interp21 + _interp22 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den8 = powr<-1>(_interp19 * _interp2 + _interp20 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _cse1 = powr<-1>(l1);
      const auto _cse2 = powr<-2>(p);
      const auto _cse3 = -_interp2;
      const auto _cse4 = _cse3 + _interp5;
      const auto _cse5 = 50. * _cse4 * _den1 * powr<6>(k);
      const auto _cse6 = _cse5 + _interp4;
      const auto _cse7 = _cse6 * _interp3;
      const auto _cse8 = _interp1 * _interp2;
      const auto _cse9 = _cse7 + _cse8;
      return fma(-0.000732421875, _cse2 * _cse9 * _den2 * _den3 * _den6 * DiFfRG::zaaqbq1_small_ref::tr1(fenv) * _interp11 * _interp12 * _interp13, fma(-0.0003662109375, _cse2 * _cse9 * _den2 * _den6 * DiFfRG::zaaqbq1_small_ref::tr2(fenv) * _interp13 * _interp14, fma(-0.0003662109375, _cse2 * _cse9 * _den2 * _den3 * _den4 * DiFfRG::zaaqbq1_small_ref::tr0(fenv) * _interp11 * _interp17 * _interp18, fma(complex<double>(0.,0.000732421875), _cse2 * _den5 * _den7 * _den8 * DiFfRG::zaaqbq1_small_ref::tr3(fenv) * _interp23 * _interp24 * _interp25 * (-_cse1 * _interp26 * _interp27 - _cse1 * _interp28 * (_interp29 + 50. * (-_interp27 + _interp30))), 0.))));
    }
    #endif

    // clang-format off
    static inline auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
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
using DiFfRG::ZAAqbq1SmallRef_num_kernel;