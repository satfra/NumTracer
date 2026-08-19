#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#ifndef NT_TRACE_COMPLEX
#define NT_TRACE_COMPLEX DiFfRG::complex<double>
#endif
#include "ZA_num_kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  class ZA_num_kernel
  {
    public:
    #if NT_ZA_NUM_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::za_num::nenv) > 0 ? (DiFfRG::za_num::nenv) : 1];
      const double dr_0 = Mq(l1);
      const double dr_1 = sqrt(powr<-1>(powr<2>(l1)));
      const double dr_2 = RF(powr<2>(k), powr<2>(l1));
      const double dr_3 = Zq(k);
      const double dr_4 = Zq(l1);
      const double dr_5 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_7 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const double dr_8 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      DiFfRG::za_num::fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8);
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
      const auto _interp20 = ntRe(DiFfRG::za_num::tr3(fenv));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
      const auto _den3 = powr<-2>(_interp15 * _interp4 + _interp18 * powr<2>(l1));
      const auto _den4 = powr<-1>(_interp3 * _interp8 + _interp9 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den5 = -powr<-1>(_interp15 * _interp8 + _interp19 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _cse1 = powr<-2>(p);
      const auto _cse2 = _interp2 * _interp3;
      const auto _cse3 = -_interp3;
      const auto _cse4 = _cse3 + _interp6;
      const auto _cse5 = 50. * _cse4 * _den1 * powr<6>(k);
      const auto _cse6 = _cse5 + _interp5;
      const auto _cse7 = _cse6 * _interp4;
      const auto _cse8 = _cse2 + _cse7;
      const auto _cse9 = _interp15 * _interp2;
      const auto _cse10 = -_interp15;
      const auto _cse11 = _cse10 + _interp17;
      const auto _cse12 = 50. * _cse11;
      const auto _cse13 = _cse12 + _interp16;
      const auto _cse14 = _cse13 * _interp4;
      const auto _cse15 = _cse14 + _cse9;
      return 0.02083333333333333 * fma(-2., _cse1 * _cse8 * _den2 * _den4 * _interp1 * powr<2>(_interp10), fma(-1., _cse1 * _cse8 * _den2 * _interp11 * _interp12, fma(-2., _cse1 * _cse15 * _den3 * _den5 * _interp13 * powr<2>(_interp14), fma(2., _cse1 * _cse15 * _den3 * _den5 * powr<2>(_interp14) * _interp20, 0.))));
    }
    #elif NT_ZA_NUM_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::za_num::nenv) > 0 ? (DiFfRG::za_num::nenv) : 1];
      const double dr_0 = Mq(l1);
      const double dr_1 = sqrt(powr<-1>(powr<2>(l1)));
      const double dr_2 = RF(powr<2>(k), powr<2>(l1));
      const double dr_3 = Zq(k);
      const double dr_4 = Zq(l1);
      const double dr_5 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_7 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const double dr_8 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      DiFfRG::za_num::fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8);
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
      const auto _interp20 = ntRe(DiFfRG::za_num::tr3(fenv));
      const auto _interp21 = ntIm(DiFfRG::za_num::tr4(fenv));
      const auto _interp22 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp23 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp24 = Zq(k);
      const auto _interp25 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp26 = dtZq(k);
      const auto _interp27 = Zq(1.02 * k);
      const auto _interp28 = Mq(l1);
      const auto _interp29 = Zq(l1);
      const auto _interp30 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp31 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp32 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp33 = ntIm(DiFfRG::za_num::tr5(fenv));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
      const auto _den3 = powr<-2>(_interp15 * _interp4 + _interp18 * powr<2>(l1));
      const auto _den4 = powr<-1>(_interp3 * _interp8 + _interp9 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den5 = -powr<-1>(_interp15 * _interp8 + _interp19 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den6 = powr<-2>(powr<2>(_interp28) + powr<2>(_interp24 * _interp25 + _interp29 * l1));
      const auto _den7 = powr<-1>(powr<2>(_interp30) + powr<2>(_interp24 * _interp31 + _interp32 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _cse1 = powr<-1>(l1);
      const auto _cse2 = powr<-2>(p);
      const auto _cse3 = -_cse1 * _interp23 * _interp24;
      const auto _cse4 = -_interp24;
      const auto _cse5 = _cse4 + _interp27;
      const auto _cse6 = 50. * _cse5;
      const auto _cse7 = _cse6 + _interp26;
      const auto _cse8 = -_cse1 * _cse7 * _interp25;
      const auto _cse9 = _cse3 + _cse8;
      const auto _cse10 = _interp2 * _interp3;
      const auto _cse11 = -_interp3;
      const auto _cse12 = _cse11 + _interp6;
      const auto _cse13 = 50. * _cse12 * _den1 * powr<6>(k);
      const auto _cse14 = _cse13 + _interp5;
      const auto _cse15 = _cse14 * _interp4;
      const auto _cse16 = _cse10 + _cse15;
      const auto _cse17 = _interp15 * _interp2;
      const auto _cse18 = -_interp15;
      const auto _cse19 = _cse18 + _interp17;
      const auto _cse20 = 50. * _cse19;
      const auto _cse21 = _cse20 + _interp16;
      const auto _cse22 = _cse21 * _interp4;
      const auto _cse23 = _cse17 + _cse22;
      return 0.02083333333333333 * fma(-2., _cse16 * _cse2 * _den2 * _den4 * _interp1 * powr<2>(_interp10), fma(-1., _cse16 * _cse2 * _den2 * _interp11 * _interp12, fma(-2., _cse2 * _cse23 * _den3 * _den5 * _interp13 * powr<2>(_interp14), fma(2., _cse2 * _cse23 * _den3 * _den5 * powr<2>(_interp14) * _interp20, fma(-4., _cse2 * _cse9 * _den6 * _den7 * _interp21 * powr<2>(_interp22), fma(4., _cse2 * _cse9 * _den6 * _den7 * powr<2>(_interp22) * _interp33, 0.))))));
    }
    #else                              // the imaginary part survives: genuinely complex
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::za_num::nenv) > 0 ? (DiFfRG::za_num::nenv) : 1];
      const double dr_0 = Mq(l1);
      const double dr_1 = sqrt(powr<-1>(powr<2>(l1)));
      const double dr_2 = RF(powr<2>(k), powr<2>(l1));
      const double dr_3 = Zq(k);
      const double dr_4 = Zq(l1);
      const double dr_5 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_6 = sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const double dr_7 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const double dr_8 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      DiFfRG::za_num::fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8);
      const auto _interp1 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp4 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp5 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(l1);
      const auto _interp7 = RB(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp8 = ZA(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp9 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp10 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + powr<2>(p)));
      const auto _interp11 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp12 = Zc(k);
      const auto _interp13 = dtZc(k);
      const auto _interp14 = Zc(1.02 * k);
      const auto _interp15 = Zc(l1);
      const auto _interp16 = Zc(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp17 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp18 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp19 = Zq(k);
      const auto _interp20 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp21 = dtZq(k);
      const auto _interp22 = Zq(1.02 * k);
      const auto _interp23 = Mq(l1);
      const auto _interp24 = Zq(l1);
      const auto _interp25 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp26 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp27 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp2 * _interp3 + _interp6 * powr<2>(l1));
      const auto _den3 = powr<-2>(_interp12 * _interp3 + _interp15 * powr<2>(l1));
      const auto _den4 = powr<-1>(_interp2 * _interp7 + _interp8 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den5 = -powr<-1>(_interp12 * _interp7 + _interp16 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den6 = powr<-2>(powr<2>(_interp23) + powr<2>(_interp19 * _interp20 + _interp24 * l1));
      const auto _den7 = powr<-1>(powr<2>(_interp25) + powr<2>(_interp19 * _interp26 + _interp27 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _cse1 = powr<-1>(l1);
      const auto _cse2 = powr<-2>(p);
      const auto _cse3 = -_cse1 * _interp18 * _interp19;
      const auto _cse4 = -_interp19;
      const auto _cse5 = _cse4 + _interp22;
      const auto _cse6 = 50. * _cse5;
      const auto _cse7 = _cse6 + _interp21;
      const auto _cse8 = -_cse1 * _cse7 * _interp20;
      const auto _cse9 = _cse3 + _cse8;
      const auto _cse10 = _interp1 * _interp2;
      const auto _cse11 = -_interp2;
      const auto _cse12 = _cse11 + _interp5;
      const auto _cse13 = 50. * _cse12 * _den1 * powr<6>(k);
      const auto _cse14 = _cse13 + _interp4;
      const auto _cse15 = _cse14 * _interp3;
      const auto _cse16 = _cse10 + _cse15;
      const auto _cse17 = _interp1 * _interp12;
      const auto _cse18 = -_interp12;
      const auto _cse19 = _cse18 + _interp14;
      const auto _cse20 = 50. * _cse19;
      const auto _cse21 = _cse20 + _interp13;
      const auto _cse22 = _cse21 * _interp3;
      const auto _cse23 = _cse17 + _cse22;
      return 0.02083333333333333 * fma(-1., _cse16 * _cse2 * _den2 * DiFfRG::za_num::tr1(fenv) * _interp10, fma(-2., _cse2 * _cse23 * _den3 * _den5 * DiFfRG::za_num::tr2(fenv) * powr<2>(_interp11), fma(2., _cse2 * _cse23 * _den3 * _den5 * DiFfRG::za_num::tr3(fenv) * powr<2>(_interp11), fma(complex<double>(0.,4.), _cse2 * _cse9 * _den6 * _den7 * DiFfRG::za_num::tr4(fenv) * powr<2>(_interp17), fma(complex<double>(0.,-4.), _cse2 * _cse9 * _den6 * _den7 * DiFfRG::za_num::tr5(fenv) * powr<2>(_interp17), fma(-2., _cse16 * _cse2 * _den2 * _den4 * DiFfRG::za_num::tr0(fenv) * powr<2>(_interp9), 0.))))));
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
using DiFfRG::ZA_num_kernel;