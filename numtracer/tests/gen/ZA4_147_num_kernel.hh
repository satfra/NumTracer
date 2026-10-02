#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#ifndef NT_TRACE_COMPLEX
#define NT_TRACE_COMPLEX DiFfRG::complex<double>
#endif
#include "ZA4_147_num_kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  class ZA4_147_num_kernel
  {
    public:
    #if NT_ZA4_147_NUM_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      double fenv[(DiFfRG::za4_147_num::nenv) > 0 ? (DiFfRG::za4_147_num::nenv) : 1];
      const double dr_0 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_1 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_2 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_4 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_5 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_6 = Mq(l1);
      const double dr_7 = -powr<-1>(l1) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
      DiFfRG::za4_147_num::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7);
      const auto _interp42 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp43 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp44 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp45 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp46 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp47 = Zq(k);
      const auto _interp48 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp49 = dtZq(k);
      const auto _interp50 = Zq(1.02 * k);
      const auto _interp51 = Mq(l1);
      const auto _interp52 = Zq(l1);
      const auto _interp53 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp54 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp55 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp56 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp57 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp58 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp59 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp60 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp61 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp63 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp64 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp66 = ZAqbq4(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp69 = ZAqbq4(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp74 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp79 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp86 = ZAqbq7(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp96 = ZAqbq7(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _den7 = powr<-2>(fma(powr<-2>(l1) * powr<2>(_interp47 * _interp48 + _interp52 * l1), powr<2>(l1), powr<2>(_interp51)));
      const auto _den13 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), powr<2>(_interp55 + _interp47 * _interp54 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))), powr<2>(_interp53)));
      const auto _den14 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), powr<2>(_interp61 + _interp47 * _interp60 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)))), powr<2>(_interp59)));
      const auto _den15 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), powr<2>(_interp58 + _interp47 * _interp57 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)))), powr<2>(_interp56)));
      const auto _interp1 = ntRe(DiFfRG::za4_147_num::tr0(fenv));
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
      const auto _interp18 = ntRe(DiFfRG::za4_147_num::tr1(fenv));
      const auto _interp19 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp20 = ntRe(DiFfRG::za4_147_num::tr3(fenv));
      const auto _interp21 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp22 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp23 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp24 = ntRe(DiFfRG::za4_147_num::tr2(fenv));
      const auto _interp25 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp26 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp27 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp28 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp29 = ntRe(DiFfRG::za4_147_num::tr4(fenv));
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
      const auto _interp41 = ntRe(DiFfRG::za4_147_num::tr6(fenv));
      const auto _interp62 = ntRe(DiFfRG::za4_147_num::tr13(fenv));
      const auto _interp65 = ntRe(DiFfRG::za4_147_num::tr7(fenv));
      const auto _interp67 = ntRe(DiFfRG::za4_147_num::tr16(fenv));
      const auto _interp68 = ntRe(DiFfRG::za4_147_num::tr9(fenv));
      const auto _interp70 = ntRe(DiFfRG::za4_147_num::tr18(fenv));
      const auto _interp71 = ntRe(DiFfRG::za4_147_num::tr12(fenv));
      const auto _den8 = 3. * powr<-1>(3. * _interp12 * _interp3 + _interp13 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den9 = 3. * powr<-1>(3. * _interp21 * _interp3 + _interp22 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den11 = -3. * powr<-1>(3. * _interp12 * _interp34 + _interp40 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _interp72 = ntRe(DiFfRG::za4_147_num::tr19(fenv));
      const auto _interp73 = ntRe(DiFfRG::za4_147_num::tr25(fenv));
      const auto _interp75 = ntRe(DiFfRG::za4_147_num::tr22(fenv));
      const auto _interp76 = ntRe(DiFfRG::za4_147_num::tr23(fenv));
      const auto _interp77 = ntRe(DiFfRG::za4_147_num::tr28(fenv));
      const auto _interp78 = ntRe(DiFfRG::za4_147_num::tr60(fenv));
      const auto _interp80 = ntRe(DiFfRG::za4_147_num::tr61(fenv));
      const auto _interp81 = ntRe(DiFfRG::za4_147_num::tr63(fenv));
      const auto _interp82 = ntRe(DiFfRG::za4_147_num::tr66(fenv));
      const auto _interp83 = ntRe(DiFfRG::za4_147_num::tr68(fenv));
      const auto _interp84 = ntRe(DiFfRG::za4_147_num::tr69(fenv));
      const auto _interp85 = ntRe(DiFfRG::za4_147_num::tr30(fenv));
      const auto _interp87 = ntRe(DiFfRG::za4_147_num::tr33(fenv));
      const auto _interp88 = ntRe(DiFfRG::za4_147_num::tr31(fenv));
      const auto _interp89 = ntRe(DiFfRG::za4_147_num::tr36(fenv));
      const auto _interp90 = ntRe(DiFfRG::za4_147_num::tr39(fenv));
      const auto _interp91 = ntRe(DiFfRG::za4_147_num::tr38(fenv));
      const auto _interp92 = ntRe(DiFfRG::za4_147_num::tr72(fenv));
      const auto _interp93 = ntRe(DiFfRG::za4_147_num::tr73(fenv));
      const auto _interp94 = ntRe(DiFfRG::za4_147_num::tr76(fenv));
      const auto _interp100 = ntRe(DiFfRG::za4_147_num::tr51(fenv));
      const auto _interp101 = ntRe(DiFfRG::za4_147_num::tr50(fenv));
      const auto _interp102 = ntRe(DiFfRG::za4_147_num::tr78(fenv));
      const auto _interp103 = ntRe(DiFfRG::za4_147_num::tr79(fenv));
      const auto _interp95 = ntRe(DiFfRG::za4_147_num::tr42(fenv));
      const auto _interp97 = ntRe(DiFfRG::za4_147_num::tr45(fenv));
      const auto _interp98 = ntRe(DiFfRG::za4_147_num::tr43(fenv));
      const auto _interp99 = ntRe(DiFfRG::za4_147_num::tr48(fenv));
      const auto _interp104 = ntRe(DiFfRG::za4_147_num::tr82(fenv));
      const auto _interp105 = ntRe(DiFfRG::za4_147_num::tr54(fenv));
      const auto _interp106 = ntRe(DiFfRG::za4_147_num::tr55(fenv));
      const auto _interp107 = ntRe(DiFfRG::za4_147_num::tr58(fenv));
      const auto _interp108 = ntRe(DiFfRG::za4_147_num::tr84(fenv));
      // clang-format off
      using _T = decltype(_den11 + _den13 + _den14 + _den15 + _den7 + _den8 + _den9 + _interp1 + _interp10 + _interp100 + _interp101 + _interp102 + _interp103 + _interp104 + _interp105 + _interp106 + _interp107 + _interp108 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + _interp98 + _interp99 + cosl1p1 + cosl1p2 + cosl1p3 + k + l1 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den1 = powr<-1>(1. + powr<6>(k));
        const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
        const auto _den3 = powr<-2>(_interp34 * _interp4 + _interp37 * powr<2>(l1));
        const auto _den4 = powr<-1>(_interp3 * _interp8 + _interp9 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den5 = powr<-1>(_interp25 * _interp3 + _interp26 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den6 = -powr<-1>(_interp34 * _interp8 + _interp38 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den10 = powr<-1>(_interp10 * _interp3 + _interp11 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
        const auto _den12 = -powr<-1>(_interp10 * _interp34 + _interp39 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += fma(0.01020408163265306, _den11 * _den12 * _den3 * _den6 * _interp29 * _interp30 * _interp31 * _interp32 * _interp33 * (_interp2 * _interp34 + (_interp35 + 50. * (-_interp34 + _interp36)) * _interp4), fma(0.00510204081632653, _den10 * _den2 * _den4 * _den8 * _interp1 * _interp14 * _interp15 * _interp16 * _interp17 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.00510204081632653, _den10 * _den2 * _den8 * _interp15 * _interp17 * _interp18 * _interp19 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.001275510204081632, _den2 * _den9 * _interp20 * powr<2>(_interp23) * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.002551020408163265, _den10 * _den2 * _den5 * _interp15 * _interp24 * _interp27 * _interp28 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(-0.02040816326530612, _den13 * _den14 * _den15 * _den7 * _interp41 * _interp42 * _interp43 * _interp44 * _interp45 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), fma(0.02040816326530612, _den13 * _den14 * _den15 * _den7 * _interp43 * _interp44 * _interp62 * _interp63 * _interp64 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), fma(-0.02040816326530612, _den13 * _den14 * _den15 * _den7 * _interp42 * _interp44 * _interp63 * _interp65 * _interp66 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), fma(-0.02040816326530612, _den13 * _den14 * _den15 * _den7 * _interp44 * _interp45 * _interp64 * _interp66 * _interp67 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), fma(-0.02040816326530612, _den13 * _den14 * _den15 * _den7 * _interp42 * _interp43 * _interp63 * _interp68 * _interp69 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), fma(-0.02040816326530612, _den13 * _den14 * _den15 * _den7 * _interp43 * _interp45 * _interp64 * _interp69 * _interp70 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), fma(0.02040816326530612, _den13 * _den14 * _den15 * _den7 * _interp42 * _interp45 * _interp66 * _interp69 * _interp71 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), 0.))))))))))));
        // clang-format on
      }
      { // subkernel 2
        const auto _cse1_k2 = powr<-1>(l1);
        const auto _cse2_k2 = -_cse1_k2 * _interp46 * _interp47;
        const auto _cse3_k2 = -_interp47;
        const auto _cse4_k2 = _cse3_k2 + _interp50;
        const auto _cse5_k2 = 50. * _cse4_k2;
        const auto _cse6_k2 = _cse5_k2 + _interp49;
        const auto _cse7_k2 = -_cse1_k2 * _cse6_k2 * _interp48;
        const auto _cse8_k2 = _cse2_k2 + _cse7_k2;
        // clang-format off
        _acc += 0.00510204081632653 * fma(-4., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp63 * _interp64 * _interp66 * _interp69 * _interp72, fma(2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp43 * _interp44 * _interp64 * _interp73 * _interp74, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp44 * _interp66 * _interp74 * _interp75, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp43 * _interp69 * _interp74 * _interp76, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp64 * _interp66 * _interp69 * _interp74 * _interp77, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp43 * _interp44 * _interp45 * _interp78 * _interp79, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp44 * _interp63 * _interp66 * _interp79 * _interp80, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp43 * _interp63 * _interp69 * _interp79 * _interp81, fma(2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp45 * _interp66 * _interp69 * _interp79 * _interp82, fma(-1., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp44 * _interp66 * _interp74 * _interp79 * _interp83, fma(-1., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp43 * _interp69 * _interp74 * _interp79 * _interp84, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp44 * _interp45 * _interp85 * _interp86, 0.))))))))))));
        // clang-format on
      }
      { // subkernel 3
        const auto _cse1_k3 = powr<-1>(l1);
        const auto _cse2_k3 = -_cse1_k3 * _interp46 * _interp47;
        const auto _cse3_k3 = -_interp47;
        const auto _cse4_k3 = _cse3_k3 + _interp50;
        const auto _cse5_k3 = 50. * _cse4_k3;
        const auto _cse6_k3 = _cse5_k3 + _interp49;
        const auto _cse7_k3 = -_cse1_k3 * _cse6_k3 * _interp48;
        const auto _cse8_k3 = _cse2_k3 + _cse7_k3;
        // clang-format off
        _acc += 0.002551020408163265 * fma(4., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp44 * _interp63 * _interp64 * _interp86 * _interp87, fma(-4., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp63 * _interp69 * _interp86 * _interp88, fma(-4., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp45 * _interp64 * _interp69 * _interp86 * _interp89, fma(2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp44 * _interp64 * _interp74 * _interp86 * _interp90, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp69 * _interp74 * _interp86 * _interp91, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp44 * _interp45 * _interp79 * _interp86 * _interp92, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp63 * _interp69 * _interp79 * _interp86 * _interp93, fma(-1., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp69 * _interp74 * _interp79 * _interp86 * _interp94, fma(2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp100 * _interp43 * _interp64 * _interp74 * _interp96, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp101 * _interp42 * _interp66 * _interp74 * _interp96, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp102 * _interp43 * _interp45 * _interp79 * _interp96, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp103 * _interp63 * _interp66 * _interp79 * _interp96, 0.))))))))))));
        // clang-format on
      }
      { // subkernel 4
        const auto _cse1_k4 = powr<-1>(l1);
        const auto _cse2_k4 = -_cse1_k4 * _interp46 * _interp47;
        const auto _cse3_k4 = -_interp47;
        const auto _cse4_k4 = _cse3_k4 + _interp50;
        const auto _cse5_k4 = 50. * _cse4_k4;
        const auto _cse6_k4 = _cse5_k4 + _interp49;
        const auto _cse7_k4 = -_cse1_k4 * _cse6_k4 * _interp48;
        const auto _cse8_k4 = _cse2_k4 + _cse7_k4;
        // clang-format off
        _acc += 0.002551020408163265 * fma(-1., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp66 * _interp74 * _interp79 * _interp96, fma(-2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp105 * _interp42 * _interp45 * _interp86 * _interp96, fma(2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp106 * _interp63 * _interp64 * _interp86 * _interp96, fma(_cse8_k4, _den13 * _den14 * _den15 * _den7 * _interp107 * _interp64 * _interp74 * _interp86 * _interp96, fma(-1., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp108 * _interp45 * _interp79 * _interp86 * _interp96, fma(-4., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp43 * _interp45 * _interp95 * _interp96, fma(4., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp43 * _interp63 * _interp64 * _interp96 * _interp97, fma(-4., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp63 * _interp66 * _interp96 * _interp98, fma(-4., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp45 * _interp64 * _interp66 * _interp96 * _interp99, 0.)))))))));
        // clang-format on
      }
      return _acc;
    }
    #elif NT_ZA4_147_NUM_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      double fenv[(DiFfRG::za4_147_num::nenv) > 0 ? (DiFfRG::za4_147_num::nenv) : 1];
      const double dr_0 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_1 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_2 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_4 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_5 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_6 = Mq(l1);
      const double dr_7 = -powr<-1>(l1) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
      DiFfRG::za4_147_num::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7);
      const auto _interp42 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp43 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp44 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp45 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp46 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp47 = Zq(k);
      const auto _interp48 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp49 = dtZq(k);
      const auto _interp50 = Zq(1.02 * k);
      const auto _interp51 = Mq(l1);
      const auto _interp52 = Zq(l1);
      const auto _interp53 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp54 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp55 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp56 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp57 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp58 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp59 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp60 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp61 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp63 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp65 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp68 = ZAqbq4(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp73 = ZAqbq4(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp82 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp91 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp104 = ZAqbq7(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp123 = ZAqbq7(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _den7 = powr<-2>(fma(powr<-2>(l1) * powr<2>(_interp47 * _interp48 + _interp52 * l1), powr<2>(l1), powr<2>(_interp51)));
      const auto _den13 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), powr<2>(_interp55 + _interp47 * _interp54 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))), powr<2>(_interp53)));
      const auto _den14 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), powr<2>(_interp61 + _interp47 * _interp60 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)))), powr<2>(_interp59)));
      const auto _den15 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), powr<2>(_interp58 + _interp47 * _interp57 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)))), powr<2>(_interp56)));
      const auto _interp1 = ntRe(DiFfRG::za4_147_num::tr0(fenv));
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
      const auto _interp18 = ntRe(DiFfRG::za4_147_num::tr1(fenv));
      const auto _interp19 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp20 = ntRe(DiFfRG::za4_147_num::tr3(fenv));
      const auto _interp21 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp22 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp23 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp24 = ntRe(DiFfRG::za4_147_num::tr2(fenv));
      const auto _interp25 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp26 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp27 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp28 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp29 = ntRe(DiFfRG::za4_147_num::tr4(fenv));
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
      const auto _interp41 = ntIm(DiFfRG::za4_147_num::tr5(fenv));
      const auto _interp62 = ntRe(DiFfRG::za4_147_num::tr6(fenv));
      const auto _interp103 = ntIm(DiFfRG::za4_147_num::tr29(fenv));
      const auto _interp105 = ntRe(DiFfRG::za4_147_num::tr30(fenv));
      const auto _interp122 = ntIm(DiFfRG::za4_147_num::tr41(fenv));
      const auto _interp124 = ntRe(DiFfRG::za4_147_num::tr42(fenv));
      const auto _interp141 = ntIm(DiFfRG::za4_147_num::tr53(fenv));
      const auto _interp142 = ntRe(DiFfRG::za4_147_num::tr54(fenv));
      const auto _den8 = 3. * powr<-1>(3. * _interp12 * _interp3 + _interp13 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den9 = 3. * powr<-1>(3. * _interp21 * _interp3 + _interp22 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den11 = -3. * powr<-1>(3. * _interp12 * _interp34 + _interp40 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _interp64 = ntRe(DiFfRG::za4_147_num::tr13(fenv));
      const auto _interp66 = ntIm(DiFfRG::za4_147_num::tr14(fenv));
      const auto _interp67 = ntRe(DiFfRG::za4_147_num::tr7(fenv));
      const auto _interp106 = ntRe(DiFfRG::za4_147_num::tr33(fenv));
      const auto _interp107 = ntIm(DiFfRG::za4_147_num::tr34(fenv));
      const auto _interp125 = ntRe(DiFfRG::za4_147_num::tr45(fenv));
      const auto _interp126 = ntIm(DiFfRG::za4_147_num::tr46(fenv));
      const auto _interp127 = ntRe(DiFfRG::za4_147_num::tr43(fenv));
      const auto _interp128 = ntIm(DiFfRG::za4_147_num::tr44(fenv));
      const auto _interp129 = ntIm(DiFfRG::za4_147_num::tr47(fenv));
      const auto _interp130 = ntRe(DiFfRG::za4_147_num::tr48(fenv));
      const auto _interp143 = ntRe(DiFfRG::za4_147_num::tr55(fenv));
      const auto _interp144 = ntIm(DiFfRG::za4_147_num::tr56(fenv));
      const auto _interp69 = ntIm(DiFfRG::za4_147_num::tr8(fenv));
      const auto _interp70 = ntIm(DiFfRG::za4_147_num::tr15(fenv));
      const auto _interp71 = ntRe(DiFfRG::za4_147_num::tr16(fenv));
      const auto _interp72 = ntRe(DiFfRG::za4_147_num::tr9(fenv));
      const auto _interp74 = ntIm(DiFfRG::za4_147_num::tr10(fenv));
      const auto _interp75 = ntIm(DiFfRG::za4_147_num::tr17(fenv));
      const auto _interp76 = ntRe(DiFfRG::za4_147_num::tr18(fenv));
      const auto _interp77 = ntIm(DiFfRG::za4_147_num::tr11(fenv));
      const auto _interp78 = ntRe(DiFfRG::za4_147_num::tr12(fenv));
      const auto _interp108 = ntRe(DiFfRG::za4_147_num::tr31(fenv));
      const auto _interp109 = ntIm(DiFfRG::za4_147_num::tr32(fenv));
      const auto _interp110 = ntIm(DiFfRG::za4_147_num::tr35(fenv));
      const auto _interp111 = ntRe(DiFfRG::za4_147_num::tr36(fenv));
      const auto _interp79 = ntRe(DiFfRG::za4_147_num::tr19(fenv));
      const auto _interp80 = ntIm(DiFfRG::za4_147_num::tr20(fenv));
      const auto _interp81 = ntIm(DiFfRG::za4_147_num::tr21(fenv));
      const auto _interp112 = ntIm(DiFfRG::za4_147_num::tr37(fenv));
      const auto _interp113 = ntRe(DiFfRG::za4_147_num::tr39(fenv));
      const auto _interp114 = ntRe(DiFfRG::za4_147_num::tr38(fenv));
      const auto _interp115 = ntIm(DiFfRG::za4_147_num::tr40(fenv));
      const auto _interp131 = ntIm(DiFfRG::za4_147_num::tr49(fenv));
      const auto _interp132 = ntRe(DiFfRG::za4_147_num::tr51(fenv));
      const auto _interp133 = ntRe(DiFfRG::za4_147_num::tr50(fenv));
      const auto _interp134 = ntIm(DiFfRG::za4_147_num::tr52(fenv));
      const auto _interp145 = ntIm(DiFfRG::za4_147_num::tr57(fenv));
      const auto _interp146 = ntRe(DiFfRG::za4_147_num::tr58(fenv));
      const auto _interp83 = ntRe(DiFfRG::za4_147_num::tr25(fenv));
      const auto _interp84 = ntRe(DiFfRG::za4_147_num::tr22(fenv));
      const auto _interp85 = ntIm(DiFfRG::za4_147_num::tr26(fenv));
      const auto _interp86 = ntRe(DiFfRG::za4_147_num::tr23(fenv));
      const auto _interp87 = ntIm(DiFfRG::za4_147_num::tr27(fenv));
      const auto _interp88 = ntIm(DiFfRG::za4_147_num::tr24(fenv));
      const auto _interp89 = ntRe(DiFfRG::za4_147_num::tr28(fenv));
      const auto _interp116 = ntIm(DiFfRG::za4_147_num::tr71(fenv));
      const auto _interp117 = ntRe(DiFfRG::za4_147_num::tr72(fenv));
      const auto _interp135 = ntIm(DiFfRG::za4_147_num::tr77(fenv));
      const auto _interp136 = ntRe(DiFfRG::za4_147_num::tr78(fenv));
      const auto _interp147 = ntIm(DiFfRG::za4_147_num::tr83(fenv));
      const auto _interp148 = ntRe(DiFfRG::za4_147_num::tr84(fenv));
      const auto _interp90 = ntIm(DiFfRG::za4_147_num::tr59(fenv));
      const auto _interp100 = ntRe(DiFfRG::za4_147_num::tr68(fenv));
      const auto _interp101 = ntRe(DiFfRG::za4_147_num::tr69(fenv));
      const auto _interp102 = ntIm(DiFfRG::za4_147_num::tr70(fenv));
      const auto _interp118 = ntRe(DiFfRG::za4_147_num::tr73(fenv));
      const auto _interp119 = ntIm(DiFfRG::za4_147_num::tr74(fenv));
      const auto _interp120 = ntIm(DiFfRG::za4_147_num::tr75(fenv));
      const auto _interp121 = ntRe(DiFfRG::za4_147_num::tr76(fenv));
      const auto _interp137 = ntRe(DiFfRG::za4_147_num::tr79(fenv));
      const auto _interp138 = ntIm(DiFfRG::za4_147_num::tr80(fenv));
      const auto _interp139 = ntIm(DiFfRG::za4_147_num::tr81(fenv));
      const auto _interp140 = ntRe(DiFfRG::za4_147_num::tr82(fenv));
      const auto _interp149 = ntIm(DiFfRG::za4_147_num::tr85(fenv));
      const auto _interp92 = ntRe(DiFfRG::za4_147_num::tr60(fenv));
      const auto _interp93 = ntRe(DiFfRG::za4_147_num::tr61(fenv));
      const auto _interp94 = ntIm(DiFfRG::za4_147_num::tr62(fenv));
      const auto _interp95 = ntRe(DiFfRG::za4_147_num::tr63(fenv));
      const auto _interp96 = ntIm(DiFfRG::za4_147_num::tr64(fenv));
      const auto _interp97 = ntIm(DiFfRG::za4_147_num::tr65(fenv));
      const auto _interp98 = ntRe(DiFfRG::za4_147_num::tr66(fenv));
      const auto _interp99 = ntIm(DiFfRG::za4_147_num::tr67(fenv));
      // clang-format off
      using _T = decltype(_den11 + _den13 + _den14 + _den15 + _den7 + _den8 + _den9 + _interp1 + _interp10 + _interp100 + _interp101 + _interp102 + _interp103 + _interp104 + _interp105 + _interp106 + _interp107 + _interp108 + _interp109 + _interp11 + _interp110 + _interp111 + _interp112 + _interp113 + _interp114 + _interp115 + _interp116 + _interp117 + _interp118 + _interp119 + _interp12 + _interp120 + _interp121 + _interp122 + _interp123 + _interp124 + _interp125 + _interp126 + _interp127 + _interp128 + _interp129 + _interp13 + _interp130 + _interp131 + _interp132 + _interp133 + _interp134 + _interp135 + _interp136 + _interp137 + _interp138 + _interp139 + _interp14 + _interp140 + _interp141 + _interp142 + _interp143 + _interp144 + _interp145 + _interp146 + _interp147 + _interp148 + _interp149 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + _interp98 + _interp99 + cosl1p1 + cosl1p2 + cosl1p3 + k + l1 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den1 = powr<-1>(1. + powr<6>(k));
        const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
        const auto _den3 = powr<-2>(_interp34 * _interp4 + _interp37 * powr<2>(l1));
        const auto _den4 = powr<-1>(_interp3 * _interp8 + _interp9 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den5 = powr<-1>(_interp25 * _interp3 + _interp26 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den6 = -powr<-1>(_interp34 * _interp8 + _interp38 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den10 = powr<-1>(_interp10 * _interp3 + _interp11 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
        const auto _den12 = -powr<-1>(_interp10 * _interp34 + _interp39 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += fma(0.01020408163265306, _den11 * _den12 * _den3 * _den6 * _interp29 * _interp30 * _interp31 * _interp32 * _interp33 * (_interp2 * _interp34 + (_interp35 + 50. * (-_interp34 + _interp36)) * _interp4), fma(0.00510204081632653, _den10 * _den2 * _den4 * _den8 * _interp1 * _interp14 * _interp15 * _interp16 * _interp17 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.00510204081632653, _den10 * _den2 * _den8 * _interp15 * _interp17 * _interp18 * _interp19 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.001275510204081632, _den2 * _den9 * _interp20 * powr<2>(_interp23) * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.002551020408163265, _den10 * _den2 * _den5 * _interp15 * _interp24 * _interp27 * _interp28 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.00510204081632653, _den13 * _den14 * _den15 * _den7 * _interp104 * _interp123 * _interp141 * _interp42 * _interp43 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), fma(0.01020408163265306, _den13 * _den14 * _den15 * _den7 * _interp122 * _interp123 * _interp42 * _interp43 * _interp44 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), fma(0.01020408163265306, _den13 * _den14 * _den15 * _den7 * _interp103 * _interp104 * _interp42 * _interp43 * _interp45 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), fma(0.02040816326530612, _den13 * _den14 * _den15 * _den7 * _interp41 * _interp42 * _interp43 * _interp44 * _interp45 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), fma(-0.00510204081632653, _den13 * _den14 * _den15 * _den7 * _interp104 * _interp123 * _interp142 * _interp43 * _interp63 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), fma(-0.01020408163265306, _den13 * _den14 * _den15 * _den7 * _interp123 * _interp124 * _interp43 * _interp44 * _interp63 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), fma(-0.01020408163265306, _den13 * _den14 * _den15 * _den7 * _interp104 * _interp105 * _interp43 * _interp45 * _interp63 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), fma(-0.02040816326530612, _den13 * _den14 * _den15 * _den7 * _interp43 * _interp44 * _interp45 * _interp62 * _interp63 * (-_interp46 * _interp47 * powr<-1>(l1) - _interp48 * (_interp49 + 50. * (-_interp47 + _interp50)) * powr<-1>(l1)), 0.)))))))))))));
        // clang-format on
      }
      { // subkernel 2
        const auto _cse1_k2 = powr<-1>(l1);
        const auto _cse2_k2 = -_cse1_k2 * _interp46 * _interp47;
        const auto _cse3_k2 = -_interp47;
        const auto _cse4_k2 = _cse3_k2 + _interp50;
        const auto _cse5_k2 = 50. * _cse4_k2;
        const auto _cse6_k2 = _cse5_k2 + _interp49;
        const auto _cse7_k2 = -_cse1_k2 * _cse6_k2 * _interp48;
        const auto _cse8_k2 = _cse2_k2 + _cse7_k2;
        // clang-format off
        _acc += 0.00510204081632653 * fma(_cse8_k2, _den13 * _den14 * _den15 * _den7 * _interp104 * _interp123 * _interp143 * _interp42 * _interp65, fma(2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp125 * _interp42 * _interp44 * _interp65, fma(2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp106 * _interp42 * _interp45 * _interp65, fma(_cse8_k2, _den13 * _den14 * _den15 * _den7 * _interp104 * _interp123 * _interp144 * _interp63 * _interp65, fma(2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp126 * _interp44 * _interp63 * _interp65, fma(2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp107 * _interp45 * _interp63 * _interp65, fma(4., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp44 * _interp45 * _interp64 * _interp65, fma(4., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp44 * _interp45 * _interp63 * _interp65 * _interp66, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp127 * _interp42 * _interp43 * _interp68, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp128 * _interp43 * _interp63 * _interp68, fma(2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp129 * _interp42 * _interp65 * _interp68, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp130 * _interp63 * _interp65 * _interp68, fma(-4., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp43 * _interp45 * _interp67 * _interp68, 0.)))))))))))));
        // clang-format on
      }
      { // subkernel 3
        const auto _cse1_k3 = powr<-1>(l1);
        const auto _cse2_k3 = -_cse1_k3 * _interp46 * _interp47;
        const auto _cse3_k3 = -_interp47;
        const auto _cse4_k3 = _cse3_k3 + _interp50;
        const auto _cse5_k3 = 50. * _cse4_k3;
        const auto _cse6_k3 = _cse5_k3 + _interp49;
        const auto _cse7_k3 = -_cse1_k3 * _cse6_k3 * _interp48;
        const auto _cse8_k3 = _cse2_k3 + _cse7_k3;
        // clang-format off
        _acc += 0.01020408163265306 * fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp43 * _interp45 * _interp63 * _interp68 * _interp69, fma(2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp45 * _interp65 * _interp68 * _interp70, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp45 * _interp63 * _interp65 * _interp68 * _interp71, fma(-1., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp108 * _interp42 * _interp43 * _interp73, fma(-1., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp109 * _interp43 * _interp63 * _interp73, fma(_cse8_k3, _den13 * _den14 * _den15 * _den7 * _interp104 * _interp110 * _interp42 * _interp65 * _interp73, fma(-1., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp111 * _interp63 * _interp65 * _interp73, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp43 * _interp44 * _interp72 * _interp73, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp43 * _interp44 * _interp63 * _interp73 * _interp74, fma(2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp44 * _interp65 * _interp73 * _interp75, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp44 * _interp63 * _interp65 * _interp73 * _interp76, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp43 * _interp68 * _interp73 * _interp77, fma(2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * _interp43 * _interp63 * _interp68 * _interp73 * _interp78, 0.)))))))))))));
        // clang-format on
      }
      { // subkernel 4
        const auto _cse1_k4 = powr<-1>(l1);
        const auto _cse2_k4 = -_cse1_k4 * _interp46 * _interp47;
        const auto _cse3_k4 = -_interp47;
        const auto _cse4_k4 = _cse3_k4 + _interp50;
        const auto _cse5_k4 = 50. * _cse4_k4;
        const auto _cse6_k4 = _cse5_k4 + _interp49;
        const auto _cse7_k4 = -_cse1_k4 * _cse6_k4 * _interp48;
        const auto _cse8_k4 = _cse2_k4 + _cse7_k4;
        // clang-format off
        _acc += 0.002551020408163265 * fma(-8., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp65 * _interp68 * _interp73 * _interp79, fma(-8., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp63 * _interp65 * _interp68 * _interp73 * _interp80, fma(_cse8_k4, _den13 * _den14 * _den15 * _den7 * _interp104 * _interp123 * _interp145 * _interp43 * _interp82, fma(2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp131 * _interp43 * _interp44 * _interp82, fma(2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp112 * _interp43 * _interp45 * _interp82, fma(_cse8_k4, _den13 * _den14 * _den15 * _den7 * _interp104 * _interp123 * _interp146 * _interp65 * _interp82, fma(2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp132 * _interp44 * _interp65 * _interp82, fma(2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp113 * _interp45 * _interp65 * _interp82, fma(-2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp133 * _interp43 * _interp68 * _interp82, fma(2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp134 * _interp65 * _interp68 * _interp82, fma(-2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp114 * _interp43 * _interp73 * _interp82, fma(2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp115 * _interp65 * _interp73 * _interp82, fma(4., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * _interp43 * _interp44 * _interp45 * _interp81 * _interp82, 0.)))))))))))));
        // clang-format on
      }
      { // subkernel 5
        const auto _cse1_k5 = powr<-1>(l1);
        const auto _cse2_k5 = -_cse1_k5 * _interp46 * _interp47;
        const auto _cse3_k5 = -_interp47;
        const auto _cse4_k5 = _cse3_k5 + _interp50;
        const auto _cse5_k5 = 50. * _cse4_k5;
        const auto _cse6_k5 = _cse5_k5 + _interp49;
        const auto _cse7_k5 = -_cse1_k5 * _cse6_k5 * _interp48;
        const auto _cse8_k5 = _cse2_k5 + _cse7_k5;
        // clang-format off
        _acc += 0.002551020408163265 * fma(4., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * _interp44 * _interp45 * _interp65 * _interp82 * _interp83, fma(-4., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * _interp43 * _interp45 * _interp68 * _interp82 * _interp84, fma(4., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * _interp45 * _interp65 * _interp68 * _interp82 * _interp85, fma(-4., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * _interp43 * _interp44 * _interp73 * _interp82 * _interp86, fma(4., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * _interp44 * _interp65 * _interp73 * _interp82 * _interp87, fma(-4., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * _interp43 * _interp68 * _interp73 * _interp82 * _interp88, fma(-4., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * _interp65 * _interp68 * _interp73 * _interp82 * _interp89, fma(_cse8_k5, _den13 * _den14 * _den15 * _den7 * _interp104 * _interp123 * _interp147 * _interp42 * _interp91, fma(2., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp135 * _interp42 * _interp44 * _interp91, fma(2., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp116 * _interp42 * _interp45 * _interp91, fma(-1., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp123 * _interp148 * _interp63 * _interp91, fma(-2., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp136 * _interp44 * _interp63 * _interp91, fma(-2., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp117 * _interp45 * _interp63 * _interp91, 0.)))))))))))));
        // clang-format on
      }
      { // subkernel 6
        const auto _cse1_k6 = powr<-1>(l1);
        const auto _cse2_k6 = -_cse1_k6 * _interp46 * _interp47;
        const auto _cse3_k6 = -_interp47;
        const auto _cse4_k6 = _cse3_k6 + _interp50;
        const auto _cse5_k6 = 50. * _cse4_k6;
        const auto _cse6_k6 = _cse5_k6 + _interp49;
        const auto _cse7_k6 = -_cse1_k6 * _cse6_k6 * _interp48;
        const auto _cse8_k6 = _cse2_k6 + _cse7_k6;
        // clang-format off
        _acc += 0.001275510204081632 * fma(-4., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp137 * _interp42 * _interp68 * _interp91, fma(-4., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp138 * _interp63 * _interp68 * _interp91, fma(-4., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp118 * _interp42 * _interp73 * _interp91, fma(-4., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp119 * _interp63 * _interp73 * _interp91, fma(_cse8_k6, _den13 * _den14 * _den15 * _den7 * _interp104 * _interp123 * _interp149 * _interp82 * _interp91, fma(2., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp139 * _interp44 * _interp82 * _interp91, fma(2., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp120 * _interp45 * _interp82 * _interp91, fma(-2., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * _interp123 * _interp140 * _interp68 * _interp82 * _interp91, fma(-4., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * _interp100 * _interp45 * _interp68 * _interp82 * _interp91, fma(-2., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * _interp104 * _interp121 * _interp73 * _interp82 * _interp91, fma(-4., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * _interp101 * _interp44 * _interp73 * _interp82 * _interp91, fma(-4., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * _interp102 * _interp68 * _interp73 * _interp82 * _interp91, fma(8., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp44 * _interp45 * _interp90 * _interp91, 0.)))))))))))));
        // clang-format on
      }
      { // subkernel 7
        const auto _cse1_k7 = powr<-1>(l1);
        const auto _cse2_k7 = -_cse1_k7 * _interp46 * _interp47;
        const auto _cse3_k7 = -_interp47;
        const auto _cse4_k7 = _cse3_k7 + _interp50;
        const auto _cse5_k7 = 50. * _cse4_k7;
        const auto _cse6_k7 = _cse5_k7 + _interp49;
        const auto _cse7_k7 = -_cse1_k7 * _cse6_k7 * _interp48;
        const auto _cse8_k7 = _cse2_k7 + _cse7_k7;
        _acc += 0.00510204081632653 * fma(-2., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * _interp44 * _interp45 * _interp63 * _interp91 * _interp92, fma(-2., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp45 * _interp68 * _interp91 * _interp93, fma(-2., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * _interp45 * _interp63 * _interp68 * _interp91 * _interp94, fma(-2., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp44 * _interp73 * _interp91 * _interp95, fma(-2., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * _interp44 * _interp63 * _interp73 * _interp91 * _interp96, fma(-2., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * _interp42 * _interp68 * _interp73 * _interp91 * _interp97, fma(2., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * _interp63 * _interp68 * _interp73 * _interp91 * _interp98, fma(_cse8_k7, _den13 * _den14 * _den15 * _den7 * _interp44 * _interp45 * _interp82 * _interp91 * _interp99, 0.))))))));
      }
      return _acc;
    }
    #else                              // the imaginary part survives: genuinely complex
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      double fenv[(DiFfRG::za4_147_num::nenv) > 0 ? (DiFfRG::za4_147_num::nenv) : 1];
      const double dr_0 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_1 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_2 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_4 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_5 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_6 = Mq(l1);
      const double dr_7 = -powr<-1>(l1) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
      DiFfRG::za4_147_num::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7);
      const auto _interp36 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp37 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp38 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp39 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp40 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp41 = Zq(k);
      const auto _interp42 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp43 = dtZq(k);
      const auto _interp44 = Zq(1.02 * k);
      const auto _interp45 = Mq(l1);
      const auto _interp46 = Zq(l1);
      const auto _interp47 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp48 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp49 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp50 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp51 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp52 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp53 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp54 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp55 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp56 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp57 = ZAqbq4(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp58 = ZAqbq4(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp59 = ZAqbq4(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp60 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp61 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp62 = ZAqbq7(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp63 = ZAqbq7(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _den7 = powr<-2>(fma(powr<-2>(l1) * powr<2>(_interp41 * _interp42 + _interp46 * l1), powr<2>(l1), powr<2>(_interp45)));
      const auto _den13 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), powr<2>(_interp49 + _interp41 * _interp48 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))), powr<2>(_interp47)));
      const auto _den14 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), powr<2>(_interp55 + _interp41 * _interp54 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)))), powr<2>(_interp53)));
      const auto _den15 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), powr<2>(_interp52 + _interp41 * _interp51 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)))), powr<2>(_interp50)));
      const auto _interp1 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp4 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp5 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(l1);
      const auto _interp7 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp8 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp9 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _interp10 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp11 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp12 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp13 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp14 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp15 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp16 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp17 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp18 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp19 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp20 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp21 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp22 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp23 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp24 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp25 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp26 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp27 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp28 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp29 = Zc(k);
      const auto _interp30 = dtZc(k);
      const auto _interp31 = Zc(1.02 * k);
      const auto _interp32 = Zc(l1);
      const auto _interp33 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp34 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp35 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _den8 = 3. * powr<-1>(3. * _interp11 * _interp2 + _interp12 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den9 = 3. * powr<-1>(3. * _interp18 * _interp2 + _interp19 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den11 = -3. * powr<-1>(3. * _interp11 * _interp29 + _interp35 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      // clang-format off
      using _T = decltype(complex<double>(0.,-4.) + _den11 + _den13 + _den14 + _den15 + _den7 + _den8 + _den9 + DiFfRG::za4_147_num::tr0(fenv) + DiFfRG::za4_147_num::tr10(fenv) + DiFfRG::za4_147_num::tr11(fenv) + DiFfRG::za4_147_num::tr12(fenv) + DiFfRG::za4_147_num::tr13(fenv) + DiFfRG::za4_147_num::tr14(fenv) + DiFfRG::za4_147_num::tr15(fenv) + DiFfRG::za4_147_num::tr16(fenv) + DiFfRG::za4_147_num::tr17(fenv) + DiFfRG::za4_147_num::tr18(fenv) + DiFfRG::za4_147_num::tr19(fenv) + DiFfRG::za4_147_num::tr1(fenv) + DiFfRG::za4_147_num::tr20(fenv) + DiFfRG::za4_147_num::tr21(fenv) + DiFfRG::za4_147_num::tr22(fenv) + DiFfRG::za4_147_num::tr23(fenv) + DiFfRG::za4_147_num::tr24(fenv) + DiFfRG::za4_147_num::tr25(fenv) + DiFfRG::za4_147_num::tr26(fenv) + DiFfRG::za4_147_num::tr27(fenv) + DiFfRG::za4_147_num::tr28(fenv) + DiFfRG::za4_147_num::tr29(fenv) + DiFfRG::za4_147_num::tr2(fenv) + DiFfRG::za4_147_num::tr30(fenv) + DiFfRG::za4_147_num::tr31(fenv) + DiFfRG::za4_147_num::tr32(fenv) + DiFfRG::za4_147_num::tr33(fenv) + DiFfRG::za4_147_num::tr34(fenv) + DiFfRG::za4_147_num::tr35(fenv) + DiFfRG::za4_147_num::tr36(fenv) + DiFfRG::za4_147_num::tr37(fenv) + DiFfRG::za4_147_num::tr38(fenv) + DiFfRG::za4_147_num::tr39(fenv) + DiFfRG::za4_147_num::tr3(fenv) + DiFfRG::za4_147_num::tr40(fenv) + DiFfRG::za4_147_num::tr41(fenv) + DiFfRG::za4_147_num::tr42(fenv) + DiFfRG::za4_147_num::tr43(fenv) + DiFfRG::za4_147_num::tr44(fenv) + DiFfRG::za4_147_num::tr45(fenv) + DiFfRG::za4_147_num::tr46(fenv) + DiFfRG::za4_147_num::tr47(fenv) + DiFfRG::za4_147_num::tr48(fenv) + DiFfRG::za4_147_num::tr49(fenv) + DiFfRG::za4_147_num::tr4(fenv) + DiFfRG::za4_147_num::tr50(fenv) + DiFfRG::za4_147_num::tr51(fenv) + DiFfRG::za4_147_num::tr52(fenv) + DiFfRG::za4_147_num::tr53(fenv) + DiFfRG::za4_147_num::tr54(fenv) + DiFfRG::za4_147_num::tr55(fenv) + DiFfRG::za4_147_num::tr56(fenv) + DiFfRG::za4_147_num::tr57(fenv) + DiFfRG::za4_147_num::tr58(fenv) + DiFfRG::za4_147_num::tr59(fenv) + DiFfRG::za4_147_num::tr5(fenv) + DiFfRG::za4_147_num::tr60(fenv) + DiFfRG::za4_147_num::tr61(fenv) + DiFfRG::za4_147_num::tr62(fenv) + DiFfRG::za4_147_num::tr63(fenv) + DiFfRG::za4_147_num::tr64(fenv) + DiFfRG::za4_147_num::tr65(fenv) + DiFfRG::za4_147_num::tr66(fenv) + DiFfRG::za4_147_num::tr67(fenv) + DiFfRG::za4_147_num::tr68(fenv) + DiFfRG::za4_147_num::tr69(fenv) + DiFfRG::za4_147_num::tr6(fenv) + DiFfRG::za4_147_num::tr70(fenv) + DiFfRG::za4_147_num::tr71(fenv) + DiFfRG::za4_147_num::tr72(fenv) + DiFfRG::za4_147_num::tr73(fenv) + DiFfRG::za4_147_num::tr74(fenv) + DiFfRG::za4_147_num::tr75(fenv) + DiFfRG::za4_147_num::tr76(fenv) + DiFfRG::za4_147_num::tr77(fenv) + DiFfRG::za4_147_num::tr78(fenv) + DiFfRG::za4_147_num::tr79(fenv) + DiFfRG::za4_147_num::tr7(fenv) + DiFfRG::za4_147_num::tr80(fenv) + DiFfRG::za4_147_num::tr81(fenv) + DiFfRG::za4_147_num::tr82(fenv) + DiFfRG::za4_147_num::tr83(fenv) + DiFfRG::za4_147_num::tr84(fenv) + DiFfRG::za4_147_num::tr85(fenv) + DiFfRG::za4_147_num::tr8(fenv) + DiFfRG::za4_147_num::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp7 + _interp8 + _interp9 + cosl1p1 + cosl1p2 + cosl1p3 + k + l1 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den1 = powr<-1>(1. + powr<6>(k));
        const auto _den2 = powr<-2>(_interp2 * _interp3 + _interp6 * powr<2>(l1));
        const auto _den3 = powr<-2>(_interp29 * _interp3 + _interp32 * powr<2>(l1));
        const auto _den4 = powr<-1>(_interp2 * _interp7 + _interp8 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den5 = powr<-1>(_interp2 * _interp21 + _interp22 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den6 = -powr<-1>(_interp29 * _interp7 + _interp33 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den10 = powr<-1>(_interp2 * _interp9 + _interp10 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
        const auto _den12 = -powr<-1>(_interp29 * _interp9 + _interp34 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += fma(0.01020408163265306, _den11 * _den12 * _den3 * _den6 * DiFfRG::za4_147_num::tr4(fenv) * _interp25 * _interp26 * _interp27 * _interp28 * (_interp1 * _interp29 + _interp3 * (_interp30 + 50. * (-_interp29 + _interp31))), fma(0.00510204081632653, _den10 * _den2 * _den4 * _den8 * DiFfRG::za4_147_num::tr0(fenv) * _interp13 * _interp14 * _interp15 * _interp16 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.00510204081632653, _den10 * _den2 * _den8 * DiFfRG::za4_147_num::tr1(fenv) * _interp14 * _interp16 * _interp17 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.001275510204081632, _den2 * _den9 * DiFfRG::za4_147_num::tr3(fenv) * powr<2>(_interp20) * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.002551020408163265, _den10 * _den2 * _den5 * DiFfRG::za4_147_num::tr2(fenv) * _interp14 * _interp23 * _interp24 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(complex<double>(0.,-0.02040816326530612), _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr5(fenv) * _interp36 * _interp37 * _interp38 * _interp39 * (-_interp40 * _interp41 * powr<-1>(l1) - _interp42 * (_interp43 + 50. * (-_interp41 + _interp44)) * powr<-1>(l1)), fma(-0.02040816326530612, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr6(fenv) * _interp37 * _interp38 * _interp39 * _interp56 * (-_interp40 * _interp41 * powr<-1>(l1) - _interp42 * (_interp43 + 50. * (-_interp41 + _interp44)) * powr<-1>(l1)), fma(0.02040816326530612, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr13(fenv) * _interp36 * _interp38 * _interp39 * _interp57 * (-_interp40 * _interp41 * powr<-1>(l1) - _interp42 * (_interp43 + 50. * (-_interp41 + _interp44)) * powr<-1>(l1)), fma(complex<double>(0.,-0.02040816326530612), _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr14(fenv) * _interp38 * _interp39 * _interp56 * _interp57 * (-_interp40 * _interp41 * powr<-1>(l1) - _interp42 * (_interp43 + 50. * (-_interp41 + _interp44)) * powr<-1>(l1)), fma(-0.02040816326530612, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr7(fenv) * _interp36 * _interp37 * _interp39 * _interp58 * (-_interp40 * _interp41 * powr<-1>(l1) - _interp42 * (_interp43 + 50. * (-_interp41 + _interp44)) * powr<-1>(l1)), fma(complex<double>(0.,0.02040816326530612), _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr8(fenv) * _interp37 * _interp39 * _interp56 * _interp58 * (-_interp40 * _interp41 * powr<-1>(l1) - _interp42 * (_interp43 + 50. * (-_interp41 + _interp44)) * powr<-1>(l1)), fma(complex<double>(0.,-0.02040816326530612), _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr15(fenv) * _interp36 * _interp39 * _interp57 * _interp58 * (-_interp40 * _interp41 * powr<-1>(l1) - _interp42 * (_interp43 + 50. * (-_interp41 + _interp44)) * powr<-1>(l1)), fma(-0.02040816326530612, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr16(fenv) * _interp39 * _interp56 * _interp57 * _interp58 * (-_interp40 * _interp41 * powr<-1>(l1) - _interp42 * (_interp43 + 50. * (-_interp41 + _interp44)) * powr<-1>(l1)), 0.)))))))))))));
        // clang-format on
      }
      { // subkernel 2
        const auto _cse1_k2 = powr<-1>(l1);
        const auto _cse2_k2 = -_cse1_k2 * _interp40 * _interp41;
        const auto _cse3_k2 = -_interp41;
        const auto _cse4_k2 = _cse3_k2 + _interp44;
        const auto _cse5_k2 = 50. * _cse4_k2;
        const auto _cse6_k2 = _cse5_k2 + _interp43;
        const auto _cse7_k2 = -_cse1_k2 * _cse6_k2 * _interp42;
        const auto _cse8_k2 = _cse2_k2 + _cse7_k2;
        // clang-format off
        _acc += complex<double>(0.,0.01020408163265306) * fma(complex<double>(0.,2.), _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr9(fenv) * _interp36 * _interp37 * _interp38 * _interp59, fma(2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr10(fenv) * _interp37 * _interp38 * _interp56 * _interp59, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr17(fenv) * _interp36 * _interp38 * _interp57 * _interp59, fma(complex<double>(0.,2.), _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr18(fenv) * _interp38 * _interp56 * _interp57 * _interp59, fma(2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr11(fenv) * _interp36 * _interp37 * _interp58 * _interp59, fma(complex<double>(0.,-2.), _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr12(fenv) * _interp37 * _interp56 * _interp58 * _interp59, fma(complex<double>(0.,2.), _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr19(fenv) * _interp36 * _interp57 * _interp58 * _interp59, fma(2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr20(fenv) * _interp56 * _interp57 * _interp58 * _interp59, fma(-1., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr21(fenv) * _interp37 * _interp38 * _interp39 * _interp60, fma(complex<double>(0.,-1.), _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr25(fenv) * _interp38 * _interp39 * _interp57 * _interp60, fma(complex<double>(0.,1.), _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr22(fenv) * _interp37 * _interp39 * _interp58 * _interp60, fma(-1., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr26(fenv) * _interp39 * _interp57 * _interp58 * _interp60, fma(complex<double>(0.,1.), _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr23(fenv) * _interp37 * _interp38 * _interp59 * _interp60, 0.)))))))))))));
        // clang-format on
      }
      { // subkernel 3
        const auto _cse1_k3 = powr<-1>(l1);
        const auto _cse2_k3 = -_cse1_k3 * _interp40 * _interp41;
        const auto _cse3_k3 = -_interp41;
        const auto _cse4_k3 = _cse3_k3 + _interp44;
        const auto _cse5_k3 = 50. * _cse4_k3;
        const auto _cse6_k3 = _cse5_k3 + _interp43;
        const auto _cse7_k3 = -_cse1_k3 * _cse6_k3 * _interp42;
        const auto _cse8_k3 = _cse2_k3 + _cse7_k3;
        // clang-format off
        _acc += complex<double>(0.,0.00510204081632653) * fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr27(fenv) * _interp38 * _interp57 * _interp59 * _interp60, fma(2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr24(fenv) * _interp37 * _interp58 * _interp59 * _interp60, fma(complex<double>(0.,2.), _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr28(fenv) * _interp57 * _interp58 * _interp59 * _interp60, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr59(fenv) * _interp36 * _interp38 * _interp39 * _interp61, fma(complex<double>(0.,2.), _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr60(fenv) * _interp38 * _interp39 * _interp56 * _interp61, fma(complex<double>(0.,2.), _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr61(fenv) * _interp36 * _interp39 * _interp58 * _interp61, fma(2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr62(fenv) * _interp39 * _interp56 * _interp58 * _interp61, fma(complex<double>(0.,2.), _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr63(fenv) * _interp36 * _interp38 * _interp59 * _interp61, fma(2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr64(fenv) * _interp38 * _interp56 * _interp59 * _interp61, fma(2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr65(fenv) * _interp36 * _interp58 * _interp59 * _interp61, fma(complex<double>(0.,-2.), _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr66(fenv) * _interp56 * _interp58 * _interp59 * _interp61, fma(-1., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr67(fenv) * _interp38 * _interp39 * _interp60 * _interp61, fma(complex<double>(0.,1.), _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr68(fenv) * _interp39 * _interp58 * _interp60 * _interp61, 0.)))))))))))));
        // clang-format on
      }
      { // subkernel 4
        const auto _cse1_k4 = powr<-1>(l1);
        const auto _cse2_k4 = -_cse1_k4 * _interp40 * _interp41;
        const auto _cse3_k4 = -_interp41;
        const auto _cse4_k4 = _cse3_k4 + _interp44;
        const auto _cse5_k4 = 50. * _cse4_k4;
        const auto _cse6_k4 = _cse5_k4 + _interp43;
        const auto _cse7_k4 = -_cse1_k4 * _cse6_k4 * _interp42;
        const auto _cse8_k4 = _cse2_k4 + _cse7_k4;
        // clang-format off
        _acc += complex<double>(0.,-0.00510204081632653) * fma(complex<double>(0.,-1.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr69(fenv) * _interp38 * _interp59 * _interp60 * _interp61, fma(-1., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr70(fenv) * _interp58 * _interp59 * _interp60 * _interp61, fma(2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr29(fenv) * _interp36 * _interp37 * _interp39 * _interp62, fma(complex<double>(0.,-2.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr30(fenv) * _interp37 * _interp39 * _interp56 * _interp62, fma(complex<double>(0.,2.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr33(fenv) * _interp36 * _interp39 * _interp57 * _interp62, fma(2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr34(fenv) * _interp39 * _interp56 * _interp57 * _interp62, fma(complex<double>(0.,-2.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr31(fenv) * _interp36 * _interp37 * _interp59 * _interp62, fma(-2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr32(fenv) * _interp37 * _interp56 * _interp59 * _interp62, fma(2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr35(fenv) * _interp36 * _interp57 * _interp59 * _interp62, fma(complex<double>(0.,-2.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr36(fenv) * _interp56 * _interp57 * _interp59 * _interp62, fma(_cse8_k4, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr37(fenv) * _interp37 * _interp39 * _interp60 * _interp62, fma(complex<double>(0.,1.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr39(fenv) * _interp39 * _interp57 * _interp60 * _interp62, fma(complex<double>(0.,-1.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr38(fenv) * _interp37 * _interp59 * _interp60 * _interp62, 0.)))))))))))));
        // clang-format on
      }
      { // subkernel 5
        const auto _cse1_k5 = powr<-1>(l1);
        const auto _cse2_k5 = -_cse1_k5 * _interp40 * _interp41;
        const auto _cse3_k5 = -_interp41;
        const auto _cse4_k5 = _cse3_k5 + _interp44;
        const auto _cse5_k5 = 50. * _cse4_k5;
        const auto _cse6_k5 = _cse5_k5 + _interp43;
        const auto _cse7_k5 = -_cse1_k5 * _cse6_k5 * _interp42;
        const auto _cse8_k5 = _cse2_k5 + _cse7_k5;
        // clang-format off
        _acc += complex<double>(0.,-0.002551020408163265) * fma(2., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr40(fenv) * _interp57 * _interp59 * _interp60 * _interp62, fma(2., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr71(fenv) * _interp36 * _interp39 * _interp61 * _interp62, fma(complex<double>(0.,-2.), _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr72(fenv) * _interp39 * _interp56 * _interp61 * _interp62, fma(complex<double>(0.,-2.), _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr73(fenv) * _interp36 * _interp59 * _interp61 * _interp62, fma(-2., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr74(fenv) * _interp56 * _interp59 * _interp61 * _interp62, fma(_cse8_k5, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr75(fenv) * _interp39 * _interp60 * _interp61 * _interp62, fma(complex<double>(0.,-1.), _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr76(fenv) * _interp59 * _interp60 * _interp61 * _interp62, fma(4., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr41(fenv) * _interp36 * _interp37 * _interp38 * _interp63, fma(complex<double>(0.,-4.), _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr42(fenv) * _interp37 * _interp38 * _interp56 * _interp63, fma(complex<double>(0.,4.), _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr45(fenv) * _interp36 * _interp38 * _interp57 * _interp63, fma(4., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr46(fenv) * _interp38 * _interp56 * _interp57 * _interp63, fma(complex<double>(0.,-4.), _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr43(fenv) * _interp36 * _interp37 * _interp58 * _interp63, fma(-4., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr44(fenv) * _interp37 * _interp56 * _interp58 * _interp63, 0.)))))))))))));
        // clang-format on
      }
      { // subkernel 6
        const auto _cse1_k6 = powr<-1>(l1);
        const auto _cse2_k6 = -_cse1_k6 * _interp40 * _interp41;
        const auto _cse3_k6 = -_interp41;
        const auto _cse4_k6 = _cse3_k6 + _interp44;
        const auto _cse5_k6 = 50. * _cse4_k6;
        const auto _cse6_k6 = _cse5_k6 + _interp43;
        const auto _cse7_k6 = -_cse1_k6 * _cse6_k6 * _interp42;
        const auto _cse8_k6 = _cse2_k6 + _cse7_k6;
        // clang-format off
        _acc += complex<double>(0.,-0.002551020408163265) * fma(4., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr47(fenv) * _interp36 * _interp57 * _interp58 * _interp63, fma(complex<double>(0.,-4.), _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr48(fenv) * _interp56 * _interp57 * _interp58 * _interp63, fma(2., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr49(fenv) * _interp37 * _interp38 * _interp60 * _interp63, fma(complex<double>(0.,2.), _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr51(fenv) * _interp38 * _interp57 * _interp60 * _interp63, fma(complex<double>(0.,-2.), _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr50(fenv) * _interp37 * _interp58 * _interp60 * _interp63, fma(2., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr52(fenv) * _interp57 * _interp58 * _interp60 * _interp63, fma(2., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr77(fenv) * _interp36 * _interp38 * _interp61 * _interp63, fma(complex<double>(0.,-2.), _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr78(fenv) * _interp38 * _interp56 * _interp61 * _interp63, fma(complex<double>(0.,-2.), _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr79(fenv) * _interp36 * _interp58 * _interp61 * _interp63, fma(-2., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr80(fenv) * _interp56 * _interp58 * _interp61 * _interp63, fma(_cse8_k6, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr81(fenv) * _interp38 * _interp60 * _interp61 * _interp63, fma(complex<double>(0.,-1.), _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr82(fenv) * _interp58 * _interp60 * _interp61 * _interp63, fma(2., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr53(fenv) * _interp36 * _interp37 * _interp62 * _interp63, 0.)))))))))))));
        // clang-format on
      }
      { // subkernel 7
        const auto _cse1_k7 = powr<-1>(l1);
        const auto _cse2_k7 = -_cse1_k7 * _interp40 * _interp41;
        const auto _cse3_k7 = -_interp41;
        const auto _cse4_k7 = _cse3_k7 + _interp44;
        const auto _cse5_k7 = 50. * _cse4_k7;
        const auto _cse6_k7 = _cse5_k7 + _interp43;
        const auto _cse7_k7 = -_cse1_k7 * _cse6_k7 * _interp42;
        const auto _cse8_k7 = _cse2_k7 + _cse7_k7;
        // clang-format off
        _acc += 0.001275510204081632 * fma(-4., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr54(fenv) * _interp37 * _interp56 * _interp62 * _interp63, fma(4., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr55(fenv) * _interp36 * _interp57 * _interp62 * _interp63, fma(complex<double>(0.,-4.), _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr56(fenv) * _interp56 * _interp57 * _interp62 * _interp63, fma(complex<double>(0.,-2.), _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr57(fenv) * _interp37 * _interp60 * _interp62 * _interp63, fma(2., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr58(fenv) * _interp57 * _interp60 * _interp62 * _interp63, fma(complex<double>(0.,-2.), _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr83(fenv) * _interp36 * _interp61 * _interp62 * _interp63, fma(-2., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr84(fenv) * _interp56 * _interp61 * _interp62 * _interp63, fma(complex<double>(0.,-1.), _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr85(fenv) * _interp60 * _interp61 * _interp62 * _interp63, 0.))))))));
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
    static inline double ntIm(double) { return 0.0; }
    template <class T> static inline auto ntRe(const T &z) -> decltype(z.real()) { return z.real(); }
    template <class T> static inline auto ntRe(const T &z) -> decltype(real(z)) requires (!requires { z.real(); }) { return real(z); }
    template <class T> static inline auto ntIm(const T &z) -> decltype(z.imag()) { return z.imag(); }
    template <class T> static inline auto ntIm(const T &z) -> decltype(imag(z)) requires (!requires { z.imag(); }) { return imag(z); }
  };
}
using DiFfRG::ZA4_147_num_kernel;