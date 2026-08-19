#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#ifndef NT_TRACE_COMPLEX
#define NT_TRACE_COMPLEX DiFfRG::complex<double>
#endif
#include "ZAqbq4_147_num_kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  class ZAqbq4_147_num_kernel
  {
    public:
    #if NT_ZAQBQ4_147_NUM_VERDICT == 2   // Pure: the Complex -> Re projection is exact
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
      double fenv[(DiFfRG::zaqbq4_147_num::nenv) > 0 ? (DiFfRG::zaqbq4_147_num::nenv) : 1];
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
      DiFfRG::zaqbq4_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20);
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp5 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp6 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp8 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp9 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp10 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp11 = Zq(k);
      const auto _interp12 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp13 = dtZq(k);
      const auto _interp14 = Zq(1.02 * k);
      const auto _interp15 = Mq(l1);
      const auto _interp16 = Zq(l1);
      const auto _interp18 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp19 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp21 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp23 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp26 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp27 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp28 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp29 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp30 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp32 = RB(powr<2>(k), fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp33 = ZA(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp34 = ZAqbq1(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp35 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp37 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp39 = ZAqbq4(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp41 = ZAqbq4(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp44 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp47 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp55 = ZAqbq7(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp59 = ZAqbq7(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp65 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp71 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp72 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp73 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp74 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp75 = ZA(l1);
      const auto _interp79 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp80 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp81 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp87 = Mq(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp88 = RF(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp89 = Zq(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp100 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp3, _interp72, fma(_interp75, powr<2>(l1), 0.)));
      const auto _den5 = powr<-1>(fma(_interp3, _interp32, fma(_interp33, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den6 = powr<-2>(fma(powr<-2>(l1) * powr<2>(_interp11 * _interp12 + _interp16 * l1), powr<2>(l1), powr<2>(_interp15)));
      const auto _den7 = powr<-1>(fma(_interp3, _interp5, fma(_interp6, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
      const auto _den8 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), powr<2>(_interp30 + _interp11 * _interp29 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))), powr<2>(_interp28)));
      const auto _den9 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), powr<2>(_interp89 + _interp11 * _interp88 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)))), powr<2>(_interp87)));
      const auto _den11 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), powr<2>(_interp81 + _interp11 * _interp80 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))), powr<2>(_interp79)));
      const auto _interp70 = ntRe(DiFfRG::zaqbq4_147_num::tr74(fenv));
      const auto _interp76 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp77 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp78 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp82 = ntRe(DiFfRG::zaqbq4_147_num::tr73(fenv));
      const auto _interp83 = ntRe(DiFfRG::zaqbq4_147_num::tr79(fenv));
      const auto _interp84 = ntRe(DiFfRG::zaqbq4_147_num::tr77(fenv));
      const auto _interp101 = ntRe(DiFfRG::zaqbq4_147_num::tr95(fenv));
      const auto _interp102 = ntRe(DiFfRG::zaqbq4_147_num::tr98(fenv));
      const auto _interp103 = ntRe(DiFfRG::zaqbq4_147_num::tr106(fenv));
      const auto _interp104 = ntRe(DiFfRG::zaqbq4_147_num::tr64(fenv));
      const auto _interp105 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp106 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp107 = Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp108 = ntRe(DiFfRG::zaqbq4_147_num::tr65(fenv));
      const auto _interp109 = ntRe(DiFfRG::zaqbq4_147_num::tr68(fenv));
      const auto _interp110 = ntRe(DiFfRG::zaqbq4_147_num::tr70(fenv));
      const auto _den10 = powr<-1>(powr<2>(_interp105) + powr<2>(_interp106 * _interp11 + _interp107 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _interp31 = ntRe(DiFfRG::zaqbq4_147_num::tr36(fenv));
      const auto _interp85 = ntRe(DiFfRG::zaqbq4_147_num::tr82(fenv));
      const auto _interp86 = ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp90 = ntRe(DiFfRG::zaqbq4_147_num::tr85(fenv));
      const auto _interp91 = ntRe(DiFfRG::zaqbq4_147_num::tr83(fenv));
      const auto _interp92 = ZAqbq4(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp93 = ntRe(DiFfRG::zaqbq4_147_num::tr88(fenv));
      const auto _interp94 = ntRe(DiFfRG::zaqbq4_147_num::tr91(fenv));
      const auto _interp95 = ntRe(DiFfRG::zaqbq4_147_num::tr90(fenv));
      const auto _interp96 = ntRe(DiFfRG::zaqbq4_147_num::tr100(fenv));
      const auto _interp97 = ntRe(DiFfRG::zaqbq4_147_num::tr101(fenv));
      const auto _interp98 = ntRe(DiFfRG::zaqbq4_147_num::tr104(fenv));
      const auto _interp99 = ntRe(DiFfRG::zaqbq4_147_num::tr94(fenv));
      const auto _interp36 = ntRe(DiFfRG::zaqbq4_147_num::tr14(fenv));
      const auto _interp38 = ntRe(DiFfRG::zaqbq4_147_num::tr39(fenv));
      const auto _interp42 = ntRe(DiFfRG::zaqbq4_147_num::tr15(fenv));
      const auto _interp43 = ntRe(DiFfRG::zaqbq4_147_num::tr41(fenv));
      const auto _interp45 = ntRe(DiFfRG::zaqbq4_147_num::tr42(fenv));
      const auto _interp48 = ntRe(DiFfRG::zaqbq4_147_num::tr44(fenv));
      const auto _interp49 = ntRe(DiFfRG::zaqbq4_147_num::tr20(fenv));
      const auto _interp50 = ntRe(DiFfRG::zaqbq4_147_num::tr47(fenv));
      const auto _interp51 = ntRe(DiFfRG::zaqbq4_147_num::tr27(fenv));
      const auto _interp52 = ntRe(DiFfRG::zaqbq4_147_num::tr30(fenv));
      const auto _interp53 = ntRe(DiFfRG::zaqbq4_147_num::tr31(fenv));
      const auto _interp54 = ntRe(DiFfRG::zaqbq4_147_num::tr48(fenv));
      const auto _interp56 = ntRe(DiFfRG::zaqbq4_147_num::tr51(fenv));
      const auto _interp57 = ntRe(DiFfRG::zaqbq4_147_num::tr52(fenv));
      const auto _interp60 = ntRe(DiFfRG::zaqbq4_147_num::tr24(fenv));
      const auto _interp62 = ntRe(DiFfRG::zaqbq4_147_num::tr33(fenv));
      const auto _interp63 = ntRe(DiFfRG::zaqbq4_147_num::tr35(fenv));
      const auto _interp64 = ntRe(DiFfRG::zaqbq4_147_num::tr54(fenv));
      const auto _interp66 = ntRe(DiFfRG::zaqbq4_147_num::tr57(fenv));
      const auto _interp67 = ntRe(DiFfRG::zaqbq4_147_num::tr58(fenv));
      const auto _interp68 = ntRe(DiFfRG::zaqbq4_147_num::tr60(fenv));
      const auto _interp69 = ntRe(DiFfRG::zaqbq4_147_num::tr62(fenv));
      const auto _interp1 = ntRe(DiFfRG::zaqbq4_147_num::tr0(fenv));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp4 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp7 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp17 = ntRe(DiFfRG::zaqbq4_147_num::tr3(fenv));
      const auto _interp20 = ntRe(DiFfRG::zaqbq4_147_num::tr4(fenv));
      const auto _interp22 = ntRe(DiFfRG::zaqbq4_147_num::tr6(fenv));
      const auto _interp24 = ntRe(DiFfRG::zaqbq4_147_num::tr8(fenv));
      const auto _interp25 = ntRe(DiFfRG::zaqbq4_147_num::tr9(fenv));
      const auto _interp40 = ntRe(DiFfRG::zaqbq4_147_num::tr12(fenv));
      const auto _interp46 = ntRe(DiFfRG::zaqbq4_147_num::tr17(fenv));
      const auto _interp58 = ntRe(DiFfRG::zaqbq4_147_num::tr21(fenv));
      const auto _interp61 = ntRe(DiFfRG::zaqbq4_147_num::tr25(fenv));
      // clang-format off
      using _T = decltype(_den1 + _den10 + _den11 + _den2 + _den5 + _den6 + _den7 + _den8 + _den9 + _interp1 + _interp10 + _interp100 + _interp101 + _interp102 + _interp103 + _interp104 + _interp105 + _interp106 + _interp107 + _interp108 + _interp109 + _interp11 + _interp110 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + _interp98 + _interp99 + cosl1p1 + cosl1p2 + k + l1 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den3 = powr<-1>(_interp3 * _interp76 + _interp77 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += fma(-0.005208333333333332, _den11 * _den2 * _den9 * _interp100 * _interp102 * _interp19 * _interp21 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(0.005208333333333332, _den11 * _den2 * _den9 * _interp100 * _interp103 * _interp18 * _interp23 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den10 * _den2 * _den3 * _interp104 * _interp35 * _interp39 * _interp78 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den10 * _den2 * _den3 * _interp108 * _interp34 * _interp44 * _interp78 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den10 * _den2 * _den3 * _interp109 * _interp44 * _interp55 * _interp78 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den10 * _den2 * _den3 * _interp110 * _interp39 * _interp65 * _interp78 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den11 * _den2 * _den3 * _interp19 * _interp27 * _interp70 * _interp78 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp100 * _interp101 * _interp19 * _interp8 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den3 * _interp23 * _interp41 * _interp78 * _interp83 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den3 * _interp19 * _interp59 * _interp78 * _interp84 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den3 * _interp41 * _interp78 * _interp82 * _interp9 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), 0.)))))))))));
        // clang-format on
      }
      { // subkernel 2
      // clang-format off
      _acc += fma(0.02083333333333333, _den11 * _den2 * _den9 * _interp18 * _interp85 * _interp86 * _interp9 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp19 * _interp8 * _interp86 * _interp90 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp8 * _interp9 * _interp91 * _interp92 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp18 * _interp19 * _interp92 * _interp93 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp19 * _interp21 * _interp86 * _interp94 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp21 * _interp9 * _interp92 * _interp95 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp18 * _interp23 * _interp86 * _interp96 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp23 * _interp8 * _interp92 * _interp97 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den11 * _den2 * _den9 * _interp21 * _interp23 * _interp92 * _interp98 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp100 * _interp18 * _interp9 * _interp99 * (_interp3 * _interp71 + _interp72 * (_interp73 + 50. * _den1 * (-_interp3 + _interp74) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp26 * _interp31 * _interp34 * _interp35 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), 0.)))))))))));
      // clang-format on
      }
      { // subkernel 3
      // clang-format off
      _acc += fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp19 * _interp27 * _interp36 * _interp37 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp35 * _interp37 * _interp38 * _interp39 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den6 * _den7 * _den8 * _interp19 * _interp26 * _interp41 * _interp42 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp34 * _interp37 * _interp43 * _interp44 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den5 * _den6 * _den8 * _interp26 * _interp39 * _interp44 * _interp45 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp34 * _interp35 * _interp47 * _interp48 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den6 * _den7 * _den8 * _interp19 * _interp41 * _interp47 * _interp49 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den5 * _den6 * _den8 * _interp39 * _interp44 * _interp47 * _interp50 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp23 * _interp26 * _interp27 * _interp51 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp23 * _interp37 * _interp41 * _interp52 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp23 * _interp27 * _interp47 * _interp53 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), 0.)))))))))));
      // clang-format on
      }
      { // subkernel 4
      // clang-format off
      _acc += fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp26 * _interp35 * _interp54 * _interp55 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp37 * _interp44 * _interp55 * _interp56 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp35 * _interp47 * _interp55 * _interp57 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp19 * _interp37 * _interp59 * _interp60 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp23 * _interp26 * _interp59 * _interp62 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.002604166666666666, _den6 * _den7 * _den8 * _interp23 * _interp47 * _interp59 * _interp63 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp26 * _interp34 * _interp64 * _interp65 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp37 * _interp39 * _interp65 * _interp66 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp34 * _interp47 * _interp65 * _interp67 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp26 * _interp55 * _interp65 * _interp68 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.002604166666666666, _den5 * _den6 * _den8 * _interp47 * _interp55 * _interp65 * _interp69 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), 0.)))))))))));
      // clang-format on
      }
      { // subkernel 5
        const auto _den4 = powr<-1>(_interp2 * _interp3 + _interp4 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += fma(0.02083333333333333, _den4 * _den6 * _den7 * _interp17 * _interp18 * _interp19 * _interp7 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den4 * _den6 * _den7 * _interp21 * _interp23 * _interp24 * _interp7 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den4 * _den6 * _den7 * _interp22 * _interp23 * _interp7 * _interp8 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp25 * _interp26 * _interp27 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp37 * _interp40 * _interp41 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp27 * _interp46 * _interp47 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp26 * _interp58 * _interp59 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp47 * _interp59 * _interp61 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den4 * _den6 * _den7 * _interp20 * _interp21 * _interp7 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den4 * _den6 * _den7 * _interp1 * _interp7 * _interp8 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), 0.))))))))));
        // clang-format on
      }
      return _acc;
    }
    #elif NT_ZAQBQ4_147_NUM_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
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
      double fenv[(DiFfRG::zaqbq4_147_num::nenv) > 0 ? (DiFfRG::zaqbq4_147_num::nenv) : 1];
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
      DiFfRG::zaqbq4_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20);
      const auto _interp2 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = ZA(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp5 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp6 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp7 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp8 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp9 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp10 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp11 = Zq(k);
      const auto _interp12 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp13 = dtZq(k);
      const auto _interp14 = Zq(1.02 * k);
      const auto _interp15 = Mq(l1);
      const auto _interp16 = Zq(l1);
      const auto _interp18 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp20 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp23 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp26 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp30 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp31 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp32 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp33 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp34 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp36 = RB(powr<2>(k), fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp37 = ZA(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp38 = ZAqbq1(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp39 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp41 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp46 = ZAqbq4(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp49 = ZAqbq4(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp54 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp59 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp74 = ZAqbq7(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp81 = ZAqbq7(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp91 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp101 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp102 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp103 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp104 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp105 = ZA(l1);
      const auto _interp106 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp107 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp108 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp109 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp110 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp111 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp121 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp122 = Mq(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp123 = RF(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp124 = Zq(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp129 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp144 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp154 = Mq(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp155 = RF(powr<2>(k), fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp156 = Zq(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp102, _interp3, fma(_interp105, powr<2>(l1), 0.)));
      const auto _den3 = powr<-1>(fma(_interp106, _interp3, fma(_interp107, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
      const auto _den4 = powr<-1>(fma(_interp2, _interp3, fma(_interp4, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den5 = powr<-1>(fma(_interp3, _interp36, fma(_interp37, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den6 = powr<-2>(fma(powr<-2>(l1) * powr<2>(_interp11 * _interp12 + _interp16 * l1), powr<2>(l1), powr<2>(_interp15)));
      const auto _den7 = powr<-1>(fma(_interp3, _interp5, fma(_interp6, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
      const auto _den8 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), powr<2>(_interp34 + _interp11 * _interp33 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))), powr<2>(_interp32)));
      const auto _den9 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), powr<2>(_interp124 + _interp11 * _interp123 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)))), powr<2>(_interp122)));
      const auto _den10 = powr<-1>(fma(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), powr<2>(_interp156 + _interp11 * _interp155 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)))), powr<2>(_interp154)));
      const auto _den11 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), powr<2>(_interp111 + _interp11 * _interp110 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))), powr<2>(_interp109)));
      const auto _interp112 = ntRe(DiFfRG::zaqbq4_147_num::tr74(fenv));
      const auto _interp115 = ntIm(DiFfRG::zaqbq4_147_num::tr78(fenv));
      const auto _interp127 = ntIm(DiFfRG::zaqbq4_147_num::tr86(fenv));
      const auto _interp132 = ntRe(DiFfRG::zaqbq4_147_num::tr88(fenv));
      const auto _interp134 = ntRe(DiFfRG::zaqbq4_147_num::tr91(fenv));
      const auto _interp136 = ntIm(DiFfRG::zaqbq4_147_num::tr92(fenv));
      const auto _interp138 = ntRe(DiFfRG::zaqbq4_147_num::tr100(fenv));
      const auto _interp140 = ntIm(DiFfRG::zaqbq4_147_num::tr102(fenv));
      const auto _interp141 = ntIm(DiFfRG::zaqbq4_147_num::tr103(fenv));
      const auto _interp142 = ntRe(DiFfRG::zaqbq4_147_num::tr104(fenv));
      const auto _interp147 = ntIm(DiFfRG::zaqbq4_147_num::tr96(fenv));
      const auto _interp149 = ntRe(DiFfRG::zaqbq4_147_num::tr98(fenv));
      const auto _interp151 = ntRe(DiFfRG::zaqbq4_147_num::tr106(fenv));
      const auto _interp152 = ntIm(DiFfRG::zaqbq4_147_num::tr107(fenv));
      const auto _interp114 = ntIm(DiFfRG::zaqbq4_147_num::tr75(fenv));
      const auto _interp116 = ntRe(DiFfRG::zaqbq4_147_num::tr79(fenv));
      const auto _interp126 = ntRe(DiFfRG::zaqbq4_147_num::tr85(fenv));
      const auto _interp131 = ntIm(DiFfRG::zaqbq4_147_num::tr87(fenv));
      const auto _interp137 = ntIm(DiFfRG::zaqbq4_147_num::tr99(fenv));
      const auto _interp139 = ntRe(DiFfRG::zaqbq4_147_num::tr101(fenv));
      const auto _interp146 = ntRe(DiFfRG::zaqbq4_147_num::tr95(fenv));
      const auto _interp150 = ntIm(DiFfRG::zaqbq4_147_num::tr105(fenv));
      const auto _interp153 = ntIm(DiFfRG::zaqbq4_147_num::tr63(fenv));
      const auto _interp157 = ntRe(DiFfRG::zaqbq4_147_num::tr64(fenv));
      const auto _interp158 = ntRe(DiFfRG::zaqbq4_147_num::tr65(fenv));
      const auto _interp159 = ntIm(DiFfRG::zaqbq4_147_num::tr66(fenv));
      const auto _interp160 = ntIm(DiFfRG::zaqbq4_147_num::tr67(fenv));
      const auto _interp161 = ntRe(DiFfRG::zaqbq4_147_num::tr68(fenv));
      const auto _interp100 = ntIm(DiFfRG::zaqbq4_147_num::tr72(fenv));
      const auto _interp113 = ntRe(DiFfRG::zaqbq4_147_num::tr73(fenv));
      const auto _interp117 = ntIm(DiFfRG::zaqbq4_147_num::tr76(fenv));
      const auto _interp118 = ntRe(DiFfRG::zaqbq4_147_num::tr77(fenv));
      const auto _interp119 = ntIm(DiFfRG::zaqbq4_147_num::tr80(fenv));
      const auto _interp120 = ntIm(DiFfRG::zaqbq4_147_num::tr81(fenv));
      const auto _interp125 = ntRe(DiFfRG::zaqbq4_147_num::tr82(fenv));
      const auto _interp128 = ntRe(DiFfRG::zaqbq4_147_num::tr83(fenv));
      const auto _interp130 = ntIm(DiFfRG::zaqbq4_147_num::tr84(fenv));
      const auto _interp133 = ntIm(DiFfRG::zaqbq4_147_num::tr89(fenv));
      const auto _interp135 = ntRe(DiFfRG::zaqbq4_147_num::tr90(fenv));
      const auto _interp143 = ntIm(DiFfRG::zaqbq4_147_num::tr93(fenv));
      const auto _interp145 = ntRe(DiFfRG::zaqbq4_147_num::tr94(fenv));
      const auto _interp148 = ntIm(DiFfRG::zaqbq4_147_num::tr97(fenv));
      const auto _interp35 = ntRe(DiFfRG::zaqbq4_147_num::tr36(fenv));
      const auto _interp42 = ntIm(DiFfRG::zaqbq4_147_num::tr37(fenv));
      const auto _interp43 = ntIm(DiFfRG::zaqbq4_147_num::tr13(fenv));
      const auto _interp44 = ntRe(DiFfRG::zaqbq4_147_num::tr14(fenv));
      const auto _interp45 = ntIm(DiFfRG::zaqbq4_147_num::tr38(fenv));
      const auto _interp47 = ntRe(DiFfRG::zaqbq4_147_num::tr39(fenv));
      const auto _interp51 = ntRe(DiFfRG::zaqbq4_147_num::tr15(fenv));
      const auto _interp52 = ntIm(DiFfRG::zaqbq4_147_num::tr16(fenv));
      const auto _interp53 = ntIm(DiFfRG::zaqbq4_147_num::tr40(fenv));
      const auto _interp55 = ntRe(DiFfRG::zaqbq4_147_num::tr41(fenv));
      const auto _interp56 = ntRe(DiFfRG::zaqbq4_147_num::tr42(fenv));
      const auto _interp162 = ntIm(DiFfRG::zaqbq4_147_num::tr69(fenv));
      const auto _interp163 = ntRe(DiFfRG::zaqbq4_147_num::tr70(fenv));
      const auto _interp164 = ntIm(DiFfRG::zaqbq4_147_num::tr71(fenv));
      const auto _interp21 = ntRe(DiFfRG::zaqbq4_147_num::tr3(fenv));
      const auto _interp24 = ntIm(DiFfRG::zaqbq4_147_num::tr5(fenv));
      const auto _interp27 = ntIm(DiFfRG::zaqbq4_147_num::tr7(fenv));
      const auto _interp28 = ntRe(DiFfRG::zaqbq4_147_num::tr8(fenv));
      const auto _interp57 = ntIm(DiFfRG::zaqbq4_147_num::tr43(fenv));
      const auto _interp60 = ntRe(DiFfRG::zaqbq4_147_num::tr44(fenv));
      const auto _interp61 = ntIm(DiFfRG::zaqbq4_147_num::tr19(fenv));
      const auto _interp62 = ntIm(DiFfRG::zaqbq4_147_num::tr45(fenv));
      const auto _interp64 = ntRe(DiFfRG::zaqbq4_147_num::tr20(fenv));
      const auto _interp65 = ntIm(DiFfRG::zaqbq4_147_num::tr46(fenv));
      const auto _interp66 = ntRe(DiFfRG::zaqbq4_147_num::tr47(fenv));
      const auto _interp67 = ntRe(DiFfRG::zaqbq4_147_num::tr27(fenv));
      const auto _interp68 = ntIm(DiFfRG::zaqbq4_147_num::tr28(fenv));
      const auto _interp69 = ntIm(DiFfRG::zaqbq4_147_num::tr29(fenv));
      const auto _interp19 = ntIm(DiFfRG::zaqbq4_147_num::tr2(fenv));
      const auto _interp25 = ntRe(DiFfRG::zaqbq4_147_num::tr6(fenv));
      const auto _interp70 = ntRe(DiFfRG::zaqbq4_147_num::tr30(fenv));
      const auto _interp71 = ntRe(DiFfRG::zaqbq4_147_num::tr31(fenv));
      const auto _interp72 = ntIm(DiFfRG::zaqbq4_147_num::tr32(fenv));
      const auto _interp73 = ntRe(DiFfRG::zaqbq4_147_num::tr48(fenv));
      const auto _interp75 = ntIm(DiFfRG::zaqbq4_147_num::tr49(fenv));
      const auto _interp76 = ntIm(DiFfRG::zaqbq4_147_num::tr50(fenv));
      const auto _interp77 = ntRe(DiFfRG::zaqbq4_147_num::tr51(fenv));
      const auto _interp78 = ntRe(DiFfRG::zaqbq4_147_num::tr52(fenv));
      const auto _interp79 = ntIm(DiFfRG::zaqbq4_147_num::tr53(fenv));
      const auto _interp83 = ntIm(DiFfRG::zaqbq4_147_num::tr23(fenv));
      const auto _interp84 = ntRe(DiFfRG::zaqbq4_147_num::tr24(fenv));
      const auto _interp86 = ntIm(DiFfRG::zaqbq4_147_num::tr26(fenv));
      const auto _interp1 = ntRe(DiFfRG::zaqbq4_147_num::tr0(fenv));
      const auto _interp17 = ntIm(DiFfRG::zaqbq4_147_num::tr1(fenv));
      const auto _interp22 = ntRe(DiFfRG::zaqbq4_147_num::tr4(fenv));
      const auto _interp29 = ntRe(DiFfRG::zaqbq4_147_num::tr9(fenv));
      const auto _interp40 = ntIm(DiFfRG::zaqbq4_147_num::tr10(fenv));
      const auto _interp48 = ntIm(DiFfRG::zaqbq4_147_num::tr11(fenv));
      const auto _interp50 = ntRe(DiFfRG::zaqbq4_147_num::tr12(fenv));
      const auto _interp58 = ntRe(DiFfRG::zaqbq4_147_num::tr17(fenv));
      const auto _interp63 = ntIm(DiFfRG::zaqbq4_147_num::tr18(fenv));
      const auto _interp80 = ntRe(DiFfRG::zaqbq4_147_num::tr21(fenv));
      const auto _interp82 = ntIm(DiFfRG::zaqbq4_147_num::tr22(fenv));
      const auto _interp87 = ntRe(DiFfRG::zaqbq4_147_num::tr33(fenv));
      const auto _interp88 = ntIm(DiFfRG::zaqbq4_147_num::tr34(fenv));
      const auto _interp89 = ntRe(DiFfRG::zaqbq4_147_num::tr35(fenv));
      const auto _interp85 = ntRe(DiFfRG::zaqbq4_147_num::tr25(fenv));
      const auto _interp90 = ntRe(DiFfRG::zaqbq4_147_num::tr54(fenv));
      const auto _interp92 = ntIm(DiFfRG::zaqbq4_147_num::tr55(fenv));
      const auto _interp93 = ntIm(DiFfRG::zaqbq4_147_num::tr56(fenv));
      const auto _interp94 = ntRe(DiFfRG::zaqbq4_147_num::tr57(fenv));
      const auto _interp95 = ntRe(DiFfRG::zaqbq4_147_num::tr58(fenv));
      const auto _interp96 = ntIm(DiFfRG::zaqbq4_147_num::tr59(fenv));
      const auto _interp97 = ntRe(DiFfRG::zaqbq4_147_num::tr60(fenv));
      const auto _interp98 = ntIm(DiFfRG::zaqbq4_147_num::tr61(fenv));
      const auto _interp99 = ntRe(DiFfRG::zaqbq4_147_num::tr62(fenv));
      // clang-format off
      using _T = decltype(_den1 + _den10 + _den11 + _den2 + _den3 + _den4 + _den5 + _den6 + _den7 + _den8 + _den9 + _interp1 + _interp10 + _interp100 + _interp101 + _interp102 + _interp103 + _interp104 + _interp105 + _interp106 + _interp107 + _interp108 + _interp109 + _interp11 + _interp110 + _interp111 + _interp112 + _interp113 + _interp114 + _interp115 + _interp116 + _interp117 + _interp118 + _interp119 + _interp12 + _interp120 + _interp121 + _interp122 + _interp123 + _interp124 + _interp125 + _interp126 + _interp127 + _interp128 + _interp129 + _interp13 + _interp130 + _interp131 + _interp132 + _interp133 + _interp134 + _interp135 + _interp136 + _interp137 + _interp138 + _interp139 + _interp14 + _interp140 + _interp141 + _interp142 + _interp143 + _interp144 + _interp145 + _interp146 + _interp147 + _interp148 + _interp149 + _interp15 + _interp150 + _interp151 + _interp152 + _interp153 + _interp154 + _interp155 + _interp156 + _interp157 + _interp158 + _interp159 + _interp16 + _interp160 + _interp161 + _interp162 + _interp163 + _interp164 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + _interp98 + _interp99 + k + l1 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
      // clang-format off
      _acc += fma(0.02083333333333333, _den11 * _den2 * _den9 * _interp121 * _interp127 * _interp18 * _interp20 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp129 * _interp132 * _interp18 * _interp20 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp144 * _interp147 * _interp18 * _interp20 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp121 * _interp134 * _interp20 * _interp23 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp129 * _interp136 * _interp20 * _interp23 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den11 * _den2 * _den9 * _interp144 * _interp149 * _interp20 * _interp23 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp121 * _interp138 * _interp18 * _interp26 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp129 * _interp140 * _interp18 * _interp26 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.005208333333333332, _den11 * _den2 * _den9 * _interp144 * _interp151 * _interp18 * _interp26 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.005208333333333332, _den11 * _den2 * _den9 * _interp121 * _interp141 * _interp23 * _interp26 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den11 * _den2 * _den9 * _interp129 * _interp142 * _interp23 * _interp26 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.002604166666666666, _den11 * _den2 * _den9 * _interp144 * _interp152 * _interp23 * _interp26 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den11 * _den2 * _den3 * _interp108 * _interp112 * _interp20 * _interp31 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den3 * _interp108 * _interp115 * _interp26 * _interp31 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 2
      // clang-format off
      _acc += fma(-0.02083333333333333, _den10 * _den2 * _den3 * _interp108 * _interp153 * _interp38 * _interp39 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den10 * _den2 * _den3 * _interp108 * _interp157 * _interp39 * _interp46 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den3 * _interp108 * _interp114 * _interp20 * _interp49 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den3 * _interp108 * _interp116 * _interp26 * _interp49 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den10 * _den2 * _den3 * _interp108 * _interp158 * _interp38 * _interp54 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den10 * _den2 * _den3 * _interp108 * _interp159 * _interp46 * _interp54 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den10 * _den2 * _den3 * _interp108 * _interp160 * _interp39 * _interp74 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den10 * _den2 * _den3 * _interp108 * _interp161 * _interp54 * _interp74 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp121 * _interp126 * _interp20 * _interp8 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp129 * _interp131 * _interp20 * _interp8 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp144 * _interp146 * _interp20 * _interp8 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp121 * _interp137 * _interp26 * _interp8 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp129 * _interp139 * _interp26 * _interp8 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.005208333333333332, _den11 * _den2 * _den9 * _interp144 * _interp150 * _interp26 * _interp8 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 3
      // clang-format off
      _acc += fma(0.01041666666666667, _den11 * _den2 * _den3 * _interp108 * _interp118 * _interp20 * _interp81 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den11 * _den2 * _den3 * _interp108 * _interp119 * _interp26 * _interp81 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den11 * _den2 * _den9 * _interp121 * _interp125 * _interp18 * _interp9 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den11 * _den2 * _den9 * _interp129 * _interp130 * _interp18 * _interp9 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp144 * _interp145 * _interp18 * _interp9 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp121 * _interp133 * _interp23 * _interp9 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * _interp129 * _interp135 * _interp23 * _interp9 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.005208333333333332, _den11 * _den2 * _den9 * _interp144 * _interp148 * _interp23 * _interp9 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den3 * _interp100 * _interp108 * _interp31 * _interp9 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den3 * _interp108 * _interp113 * _interp49 * _interp9 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den11 * _den2 * _den9 * _interp120 * _interp121 * _interp8 * _interp9 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * _interp128 * _interp129 * _interp8 * _interp9 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * _interp143 * _interp144 * _interp8 * _interp9 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den3 * _interp108 * _interp117 * _interp81 * _interp9 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 4
      // clang-format off
      _acc += fma(-0.01041666666666667, _den10 * _den2 * _den3 * _interp108 * _interp162 * _interp38 * _interp91 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den10 * _den2 * _den3 * _interp108 * _interp163 * _interp46 * _interp91 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den10 * _den2 * _den3 * _interp108 * _interp164 * _interp74 * _interp91 * (_interp101 * _interp3 + _interp102 * (_interp103 + 50. * _den1 * (_interp104 - _interp3) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp30 * _interp35 * _interp38 * _interp39 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp38 * _interp39 * _interp41 * _interp42 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den6 * _den7 * _den8 * _interp20 * _interp30 * _interp31 * _interp43 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp20 * _interp31 * _interp41 * _interp44 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den5 * _den6 * _den8 * _interp30 * _interp39 * _interp45 * _interp46 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp39 * _interp41 * _interp46 * _interp47 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den6 * _den7 * _den8 * _interp20 * _interp30 * _interp49 * _interp51 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den6 * _den7 * _den8 * _interp20 * _interp41 * _interp49 * _interp52 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den5 * _den6 * _den8 * _interp30 * _interp38 * _interp53 * _interp54 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * _interp38 * _interp41 * _interp54 * _interp55 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den5 * _den6 * _den8 * _interp30 * _interp46 * _interp54 * _interp56 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 5
      // clang-format off
      _acc += fma(0.02083333333333333, _den5 * _den6 * _den8 * _interp41 * _interp46 * _interp54 * _interp57 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp38 * _interp39 * _interp59 * _interp60 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den6 * _den7 * _den8 * _interp20 * _interp31 * _interp59 * _interp61 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den5 * _den6 * _den8 * _interp39 * _interp46 * _interp59 * _interp62 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den6 * _den7 * _den8 * _interp20 * _interp49 * _interp59 * _interp64 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den5 * _den6 * _den8 * _interp38 * _interp54 * _interp59 * _interp65 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den5 * _den6 * _den8 * _interp46 * _interp54 * _interp59 * _interp66 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp26 * _interp30 * _interp31 * _interp67 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp26 * _interp31 * _interp41 * _interp68 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den6 * _den7 * _den8 * _interp26 * _interp30 * _interp49 * _interp69 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den4 * _den6 * _den7 * _interp18 * _interp20 * _interp21 * _interp7 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den4 * _den6 * _den7 * _interp20 * _interp23 * _interp24 * _interp7 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den4 * _den6 * _den7 * _interp18 * _interp26 * _interp27 * _interp7 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den4 * _den6 * _den7 * _interp23 * _interp26 * _interp28 * _interp7 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 6
      // clang-format off
      _acc += fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp26 * _interp41 * _interp49 * _interp70 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp26 * _interp31 * _interp59 * _interp71 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den6 * _den7 * _den8 * _interp26 * _interp49 * _interp59 * _interp72 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp30 * _interp39 * _interp73 * _interp74 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp39 * _interp41 * _interp74 * _interp75 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den5 * _den6 * _den8 * _interp30 * _interp54 * _interp74 * _interp76 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp41 * _interp54 * _interp74 * _interp77 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp39 * _interp59 * _interp74 * _interp78 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den5 * _den6 * _den8 * _interp54 * _interp59 * _interp74 * _interp79 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den4 * _den6 * _den7 * _interp19 * _interp20 * _interp7 * _interp8 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den4 * _den6 * _den7 * _interp25 * _interp26 * _interp7 * _interp8 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den6 * _den7 * _den8 * _interp20 * _interp30 * _interp81 * _interp83 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp20 * _interp41 * _interp81 * _interp84 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den6 * _den7 * _den8 * _interp20 * _interp59 * _interp81 * _interp86 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 7
      // clang-format off
      _acc += fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp26 * _interp30 * _interp81 * _interp87 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp26 * _interp41 * _interp81 * _interp88 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.002604166666666666, _den6 * _den7 * _den8 * _interp26 * _interp59 * _interp81 * _interp89 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp29 * _interp30 * _interp31 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp31 * _interp40 * _interp41 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den6 * _den7 * _den8 * _interp30 * _interp48 * _interp49 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * _interp41 * _interp49 * _interp50 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp31 * _interp58 * _interp59 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den6 * _den7 * _den8 * _interp49 * _interp59 * _interp63 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den4 * _den6 * _den7 * _interp17 * _interp18 * _interp7 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den4 * _den6 * _den7 * _interp22 * _interp23 * _interp7 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den4 * _den6 * _den7 * _interp1 * _interp7 * _interp8 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp30 * _interp80 * _interp81 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * _interp41 * _interp81 * _interp82 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 8
      // clang-format off
      _acc += fma(0.005208333333333332, _den6 * _den7 * _den8 * _interp59 * _interp81 * _interp85 * _interp9 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp30 * _interp38 * _interp90 * _interp91 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp38 * _interp41 * _interp91 * _interp92 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den5 * _den6 * _den8 * _interp30 * _interp46 * _interp91 * _interp93 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * _interp41 * _interp46 * _interp91 * _interp94 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp38 * _interp59 * _interp91 * _interp95 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den5 * _den6 * _den8 * _interp46 * _interp59 * _interp91 * _interp96 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp30 * _interp74 * _interp91 * _interp97 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * _interp41 * _interp74 * _interp91 * _interp98 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), fma(-0.002604166666666666, _den5 * _den6 * _den8 * _interp59 * _interp74 * _interp91 * _interp99 * (-_interp10 * _interp11 * powr<-1>(l1) - _interp12 * (_interp13 + 50. * (-_interp11 + _interp14)) * powr<-1>(l1)) * powr<-2>(p), 0.))))))))));
      // clang-format on
      }
      return _acc;
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
      double fenv[(DiFfRG::zaqbq4_147_num::nenv) > 0 ? (DiFfRG::zaqbq4_147_num::nenv) : 1];
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
      DiFfRG::zaqbq4_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20);
      const auto _interp1 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = ZA(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp4 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp5 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp6 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp7 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp8 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp9 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp10 = Zq(k);
      const auto _interp11 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp12 = dtZq(k);
      const auto _interp13 = Zq(1.02 * k);
      const auto _interp14 = Mq(l1);
      const auto _interp15 = Zq(l1);
      const auto _interp16 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp17 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp18 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp19 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp20 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp21 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp22 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp23 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp24 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp25 = RB(powr<2>(k), fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp26 = ZA(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp27 = ZAqbq1(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp28 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp29 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp30 = ZAqbq4(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp31 = ZAqbq4(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp32 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp33 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp34 = ZAqbq7(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp35 = ZAqbq7(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp36 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp37 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp38 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp39 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp40 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp41 = ZA(l1);
      const auto _interp42 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp43 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp44 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp45 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp46 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp47 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp48 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp49 = Mq(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp50 = RF(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp51 = Zq(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp52 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp53 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp54 = Mq(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp55 = RF(powr<2>(k), fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp56 = Zq(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp2, _interp38, fma(_interp41, powr<2>(l1), 0.)));
      const auto _den3 = powr<-1>(fma(_interp2, _interp42, fma(_interp43, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
      const auto _den4 = powr<-1>(fma(_interp1, _interp2, fma(_interp3, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den5 = powr<-1>(fma(_interp2, _interp25, fma(_interp26, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den6 = powr<-2>(fma(powr<-2>(l1) * powr<2>(_interp10 * _interp11 + _interp15 * l1), powr<2>(l1), powr<2>(_interp14)));
      const auto _den7 = powr<-1>(fma(_interp2, _interp4, fma(_interp5, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
      const auto _den8 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), powr<2>(_interp24 + _interp10 * _interp23 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))), powr<2>(_interp22)));
      const auto _den9 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), powr<2>(_interp51 + _interp10 * _interp50 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)))), powr<2>(_interp49)));
      const auto _den10 = powr<-1>(fma(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), powr<2>(_interp56 + _interp10 * _interp55 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)))), powr<2>(_interp54)));
      const auto _den11 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), powr<2>(_interp47 + _interp10 * _interp46 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))), powr<2>(_interp45)));
      // clang-format off
      using _T = decltype(complex<double>(0.,-0.005208333333333332) + _den1 + _den10 + _den11 + _den2 + _den3 + _den4 + _den5 + _den6 + _den7 + _den8 + _den9 + DiFfRG::zaqbq4_147_num::tr0(fenv) + DiFfRG::zaqbq4_147_num::tr100(fenv) + DiFfRG::zaqbq4_147_num::tr101(fenv) + DiFfRG::zaqbq4_147_num::tr102(fenv) + DiFfRG::zaqbq4_147_num::tr103(fenv) + DiFfRG::zaqbq4_147_num::tr104(fenv) + DiFfRG::zaqbq4_147_num::tr105(fenv) + DiFfRG::zaqbq4_147_num::tr106(fenv) + DiFfRG::zaqbq4_147_num::tr107(fenv) + DiFfRG::zaqbq4_147_num::tr10(fenv) + DiFfRG::zaqbq4_147_num::tr11(fenv) + DiFfRG::zaqbq4_147_num::tr12(fenv) + DiFfRG::zaqbq4_147_num::tr13(fenv) + DiFfRG::zaqbq4_147_num::tr14(fenv) + DiFfRG::zaqbq4_147_num::tr15(fenv) + DiFfRG::zaqbq4_147_num::tr16(fenv) + DiFfRG::zaqbq4_147_num::tr17(fenv) + DiFfRG::zaqbq4_147_num::tr18(fenv) + DiFfRG::zaqbq4_147_num::tr19(fenv) + DiFfRG::zaqbq4_147_num::tr1(fenv) + DiFfRG::zaqbq4_147_num::tr20(fenv) + DiFfRG::zaqbq4_147_num::tr21(fenv) + DiFfRG::zaqbq4_147_num::tr22(fenv) + DiFfRG::zaqbq4_147_num::tr23(fenv) + DiFfRG::zaqbq4_147_num::tr24(fenv) + DiFfRG::zaqbq4_147_num::tr25(fenv) + DiFfRG::zaqbq4_147_num::tr26(fenv) + DiFfRG::zaqbq4_147_num::tr27(fenv) + DiFfRG::zaqbq4_147_num::tr28(fenv) + DiFfRG::zaqbq4_147_num::tr29(fenv) + DiFfRG::zaqbq4_147_num::tr2(fenv) + DiFfRG::zaqbq4_147_num::tr30(fenv) + DiFfRG::zaqbq4_147_num::tr31(fenv) + DiFfRG::zaqbq4_147_num::tr32(fenv) + DiFfRG::zaqbq4_147_num::tr33(fenv) + DiFfRG::zaqbq4_147_num::tr34(fenv) + DiFfRG::zaqbq4_147_num::tr35(fenv) + DiFfRG::zaqbq4_147_num::tr36(fenv) + DiFfRG::zaqbq4_147_num::tr37(fenv) + DiFfRG::zaqbq4_147_num::tr38(fenv) + DiFfRG::zaqbq4_147_num::tr39(fenv) + DiFfRG::zaqbq4_147_num::tr3(fenv) + DiFfRG::zaqbq4_147_num::tr40(fenv) + DiFfRG::zaqbq4_147_num::tr41(fenv) + DiFfRG::zaqbq4_147_num::tr42(fenv) + DiFfRG::zaqbq4_147_num::tr43(fenv) + DiFfRG::zaqbq4_147_num::tr44(fenv) + DiFfRG::zaqbq4_147_num::tr45(fenv) + DiFfRG::zaqbq4_147_num::tr46(fenv) + DiFfRG::zaqbq4_147_num::tr47(fenv) + DiFfRG::zaqbq4_147_num::tr48(fenv) + DiFfRG::zaqbq4_147_num::tr49(fenv) + DiFfRG::zaqbq4_147_num::tr4(fenv) + DiFfRG::zaqbq4_147_num::tr50(fenv) + DiFfRG::zaqbq4_147_num::tr51(fenv) + DiFfRG::zaqbq4_147_num::tr52(fenv) + DiFfRG::zaqbq4_147_num::tr53(fenv) + DiFfRG::zaqbq4_147_num::tr54(fenv) + DiFfRG::zaqbq4_147_num::tr55(fenv) + DiFfRG::zaqbq4_147_num::tr56(fenv) + DiFfRG::zaqbq4_147_num::tr57(fenv) + DiFfRG::zaqbq4_147_num::tr58(fenv) + DiFfRG::zaqbq4_147_num::tr59(fenv) + DiFfRG::zaqbq4_147_num::tr5(fenv) + DiFfRG::zaqbq4_147_num::tr60(fenv) + DiFfRG::zaqbq4_147_num::tr61(fenv) + DiFfRG::zaqbq4_147_num::tr62(fenv) + DiFfRG::zaqbq4_147_num::tr63(fenv) + DiFfRG::zaqbq4_147_num::tr64(fenv) + DiFfRG::zaqbq4_147_num::tr65(fenv) + DiFfRG::zaqbq4_147_num::tr66(fenv) + DiFfRG::zaqbq4_147_num::tr67(fenv) + DiFfRG::zaqbq4_147_num::tr68(fenv) + DiFfRG::zaqbq4_147_num::tr69(fenv) + DiFfRG::zaqbq4_147_num::tr6(fenv) + DiFfRG::zaqbq4_147_num::tr70(fenv) + DiFfRG::zaqbq4_147_num::tr71(fenv) + DiFfRG::zaqbq4_147_num::tr72(fenv) + DiFfRG::zaqbq4_147_num::tr73(fenv) + DiFfRG::zaqbq4_147_num::tr74(fenv) + DiFfRG::zaqbq4_147_num::tr75(fenv) + DiFfRG::zaqbq4_147_num::tr76(fenv) + DiFfRG::zaqbq4_147_num::tr77(fenv) + DiFfRG::zaqbq4_147_num::tr78(fenv) + DiFfRG::zaqbq4_147_num::tr79(fenv) + DiFfRG::zaqbq4_147_num::tr7(fenv) + DiFfRG::zaqbq4_147_num::tr80(fenv) + DiFfRG::zaqbq4_147_num::tr81(fenv) + DiFfRG::zaqbq4_147_num::tr82(fenv) + DiFfRG::zaqbq4_147_num::tr83(fenv) + DiFfRG::zaqbq4_147_num::tr84(fenv) + DiFfRG::zaqbq4_147_num::tr85(fenv) + DiFfRG::zaqbq4_147_num::tr86(fenv) + DiFfRG::zaqbq4_147_num::tr87(fenv) + DiFfRG::zaqbq4_147_num::tr88(fenv) + DiFfRG::zaqbq4_147_num::tr89(fenv) + DiFfRG::zaqbq4_147_num::tr8(fenv) + DiFfRG::zaqbq4_147_num::tr90(fenv) + DiFfRG::zaqbq4_147_num::tr91(fenv) + DiFfRG::zaqbq4_147_num::tr92(fenv) + DiFfRG::zaqbq4_147_num::tr93(fenv) + DiFfRG::zaqbq4_147_num::tr94(fenv) + DiFfRG::zaqbq4_147_num::tr95(fenv) + DiFfRG::zaqbq4_147_num::tr96(fenv) + DiFfRG::zaqbq4_147_num::tr97(fenv) + DiFfRG::zaqbq4_147_num::tr98(fenv) + DiFfRG::zaqbq4_147_num::tr99(fenv) + DiFfRG::zaqbq4_147_num::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp6 + _interp7 + _interp8 + _interp9 + k + l1 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
      // clang-format off
      _acc += fma(0.02083333333333333, _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr74(fenv) * _interp17 * _interp21 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,0.01041666666666667), _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr78(fenv) * _interp19 * _interp21 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,0.02083333333333333), _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr63(fenv) * _interp27 * _interp28 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr64(fenv) * _interp28 * _interp30 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,0.02083333333333333), _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr75(fenv) * _interp17 * _interp31 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr79(fenv) * _interp19 * _interp31 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr65(fenv) * _interp27 * _interp32 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,0.02083333333333333), _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr66(fenv) * _interp30 * _interp32 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,0.01041666666666667), _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr67(fenv) * _interp28 * _interp34 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr68(fenv) * _interp32 * _interp34 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr77(fenv) * _interp17 * _interp35 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,0.005208333333333332), _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr80(fenv) * _interp19 * _interp35 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,0.01041666666666667), _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr69(fenv) * _interp27 * _interp36 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr70(fenv) * _interp30 * _interp36 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 2
      // clang-format off
      _acc += fma(complex<double>(0.,0.005208333333333332), _den10 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr71(fenv) * _interp34 * _interp36 * _interp44 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,-0.02083333333333333), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr86(fenv) * _interp16 * _interp17 * _interp48 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr91(fenv) * _interp17 * _interp18 * _interp48 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr100(fenv) * _interp16 * _interp19 * _interp48 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,-0.005208333333333332), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr103(fenv) * _interp18 * _interp19 * _interp48 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr88(fenv) * _interp16 * _interp17 * _interp52 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,0.01041666666666667), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr92(fenv) * _interp17 * _interp18 * _interp52 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr102(fenv) * _interp16 * _interp19 * _interp52 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr104(fenv) * _interp18 * _interp19 * _interp52 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr96(fenv) * _interp16 * _interp17 * _interp53 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(-0.005208333333333332, _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr98(fenv) * _interp17 * _interp18 * _interp53 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(0.005208333333333332, _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr106(fenv) * _interp16 * _interp19 * _interp53 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,-0.002604166666666666), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr107(fenv) * _interp18 * _interp19 * _interp53 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr85(fenv) * _interp17 * _interp48 * _interp7 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 3
      // clang-format off
      _acc += fma(complex<double>(0.,-0.01041666666666667), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr99(fenv) * _interp19 * _interp48 * _interp7 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,0.02083333333333333), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr87(fenv) * _interp17 * _interp52 * _interp7 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr101(fenv) * _interp19 * _interp52 * _interp7 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr95(fenv) * _interp17 * _interp53 * _interp7 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,-0.005208333333333332), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr105(fenv) * _interp19 * _interp53 * _interp7 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,0.02083333333333333), _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr72(fenv) * _interp21 * _interp44 * _interp8 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr73(fenv) * _interp31 * _interp44 * _interp8 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,0.01041666666666667), _den11 * _den2 * _den3 * DiFfRG::zaqbq4_147_num::tr76(fenv) * _interp35 * _interp44 * _interp8 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(0.02083333333333333, _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr82(fenv) * _interp16 * _interp48 * _interp8 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr89(fenv) * _interp18 * _interp48 * _interp8 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,-0.02083333333333333), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr84(fenv) * _interp16 * _interp52 * _interp8 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(-0.01041666666666667, _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr90(fenv) * _interp18 * _interp52 * _interp8 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(0.01041666666666667, _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr94(fenv) * _interp16 * _interp53 * _interp8 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,-0.005208333333333332), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr97(fenv) * _interp18 * _interp53 * _interp8 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 4
      // clang-format off
      _acc += fma(complex<double>(0.,-0.02083333333333333), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr81(fenv) * _interp48 * _interp7 * _interp8 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(-0.02083333333333333, _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr83(fenv) * _interp52 * _interp7 * _interp8 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den11 * _den2 * _den9 * DiFfRG::zaqbq4_147_num::tr93(fenv) * _interp53 * _interp7 * _interp8 * (_interp2 * _interp37 + _interp38 * (_interp39 + 50. * _den1 * (-_interp2 + _interp40) * powr<6>(k))) * powr<-2>(p), fma(complex<double>(0.,0.02083333333333333), _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr13(fenv) * _interp17 * _interp20 * _interp21 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr27(fenv) * _interp19 * _interp20 * _interp21 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr36(fenv) * _interp20 * _interp27 * _interp28 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr14(fenv) * _interp17 * _interp21 * _interp29 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr28(fenv) * _interp19 * _interp21 * _interp29 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,0.02083333333333333), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr37(fenv) * _interp27 * _interp28 * _interp29 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.02083333333333333), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr38(fenv) * _interp20 * _interp28 * _interp30 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr39(fenv) * _interp28 * _interp29 * _interp30 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr15(fenv) * _interp17 * _interp20 * _interp31 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,0.01041666666666667), _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr29(fenv) * _interp19 * _interp20 * _interp31 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,0.02083333333333333), _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr16(fenv) * _interp17 * _interp29 * _interp31 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 5
      // clang-format off
      _acc += fma(0.01041666666666667, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr30(fenv) * _interp19 * _interp29 * _interp31 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.02083333333333333), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr40(fenv) * _interp20 * _interp27 * _interp32 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr41(fenv) * _interp27 * _interp29 * _interp32 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr42(fenv) * _interp20 * _interp30 * _interp32 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.02083333333333333), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr43(fenv) * _interp29 * _interp30 * _interp32 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,0.01041666666666667), _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr19(fenv) * _interp17 * _interp21 * _interp33 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr31(fenv) * _interp19 * _interp21 * _interp33 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr44(fenv) * _interp27 * _interp28 * _interp33 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr45(fenv) * _interp28 * _interp30 * _interp33 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr20(fenv) * _interp17 * _interp31 * _interp33 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,0.005208333333333332), _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr32(fenv) * _interp19 * _interp31 * _interp33 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr46(fenv) * _interp27 * _interp32 * _interp33 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr47(fenv) * _interp30 * _interp32 * _interp33 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr48(fenv) * _interp20 * _interp28 * _interp34 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 6
      // clang-format off
      _acc += fma(complex<double>(0.,0.01041666666666667), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr49(fenv) * _interp28 * _interp29 * _interp34 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr50(fenv) * _interp20 * _interp32 * _interp34 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr51(fenv) * _interp29 * _interp32 * _interp34 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr52(fenv) * _interp28 * _interp33 * _interp34 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.005208333333333332), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr53(fenv) * _interp32 * _interp33 * _interp34 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,0.01041666666666667), _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr23(fenv) * _interp17 * _interp20 * _interp35 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr33(fenv) * _interp19 * _interp20 * _interp35 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr24(fenv) * _interp17 * _interp29 * _interp35 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.005208333333333332), _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr34(fenv) * _interp19 * _interp29 * _interp35 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,0.005208333333333332), _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr26(fenv) * _interp17 * _interp33 * _interp35 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.002604166666666666, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr35(fenv) * _interp19 * _interp33 * _interp35 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr54(fenv) * _interp20 * _interp27 * _interp36 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,0.01041666666666667), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr55(fenv) * _interp27 * _interp29 * _interp36 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr56(fenv) * _interp20 * _interp30 * _interp36 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 7
      // clang-format off
      _acc += fma(-0.01041666666666667, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr57(fenv) * _interp29 * _interp30 * _interp36 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr58(fenv) * _interp27 * _interp33 * _interp36 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.005208333333333332), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr59(fenv) * _interp30 * _interp33 * _interp36 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr60(fenv) * _interp20 * _interp34 * _interp36 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,0.005208333333333332), _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr61(fenv) * _interp29 * _interp34 * _interp36 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.002604166666666666, _den5 * _den6 * _den8 * DiFfRG::zaqbq4_147_num::tr62(fenv) * _interp33 * _interp34 * _interp36 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr3(fenv) * _interp16 * _interp17 * _interp6 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr5(fenv) * _interp17 * _interp18 * _interp6 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr7(fenv) * _interp16 * _interp19 * _interp6 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.005208333333333332, _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr8(fenv) * _interp18 * _interp19 * _interp6 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.02083333333333333), _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr2(fenv) * _interp17 * _interp6 * _interp7 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr6(fenv) * _interp19 * _interp6 * _interp7 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr9(fenv) * _interp20 * _interp21 * _interp8 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.02083333333333333), _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr10(fenv) * _interp21 * _interp29 * _interp8 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 8
      // clang-format off
      _acc += fma(complex<double>(0.,0.02083333333333333), _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr11(fenv) * _interp20 * _interp31 * _interp8 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.02083333333333333, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr12(fenv) * _interp29 * _interp31 * _interp8 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr17(fenv) * _interp21 * _interp33 * _interp8 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,0.01041666666666667), _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr18(fenv) * _interp31 * _interp33 * _interp8 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.01041666666666667, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr21(fenv) * _interp20 * _interp35 * _interp8 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.01041666666666667), _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr22(fenv) * _interp29 * _interp35 * _interp8 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(0.005208333333333332, _den6 * _den7 * _den8 * DiFfRG::zaqbq4_147_num::tr25(fenv) * _interp33 * _interp35 * _interp8 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(complex<double>(0.,-0.02083333333333333), _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr1(fenv) * _interp16 * _interp6 * _interp8 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.01041666666666667, _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr4(fenv) * _interp18 * _interp6 * _interp8 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), fma(-0.02083333333333333, _den4 * _den6 * _den7 * DiFfRG::zaqbq4_147_num::tr0(fenv) * _interp6 * _interp7 * _interp8 * (-_interp11 * (_interp12 + 50. * (-_interp10 + _interp13)) * powr<-1>(l1) - _interp10 * _interp9 * powr<-1>(l1)) * powr<-2>(p), 0.))))))))));
      // clang-format on
      }
      return _acc;
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
using DiFfRG::ZAqbq4_147_num_kernel;