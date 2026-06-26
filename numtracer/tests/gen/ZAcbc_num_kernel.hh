#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "numtracer/expr/real_cse.hpp"
#include "ZAcbc_num_kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class ZAcbc_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[DiFfRG::zacbc_num::nenv];
      DiFfRG::zacbc_num::fill(fenv, l1, cos1, cos2, p);
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const auto _interp1 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp4 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp5 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(l1);
      const auto _interp7 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp8 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp9 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp10 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp11 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (cosl1p1 - cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp12 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp13 = Zc(k);
      const auto _interp14 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp15 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp16 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp17 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp18 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp19 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp20 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp21 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp22 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp23 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp24 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp25 = RBdot(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp26 = RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp27 = dtZc(k);
      const auto _interp28 = Zc(1.02 * k);
      const auto _interp29 = Zc(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp30 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp31 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp32 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp33 = Zc(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp34 = RBdot(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp2 * _interp3 + _interp6 * powr<2>(l1));
      const auto _den3 = powr<-1>(_interp2 * _interp3 + _interp6 * powr<2>(l1));
      const auto _den4 = powr<-1>(_interp15 * _interp2 + _interp16 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den5 = powr<-1>(_interp2 * _interp7 + _interp8 * (powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den6 = -powr<-1>(_interp12 * _interp13 + _interp14 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den7 = powr<-2>(_interp13 * _interp32 + _interp33 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den8 = -powr<-1>(_interp13 * _interp32 + _interp33 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den9 = -powr<-1>(_interp13 * _interp20 + _interp21 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den10 = powr<-2>(_interp13 * _interp26 + _interp29 * (powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den11 = -powr<-1>(_interp13 * _interp26 + _interp29 * (powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _cse1 = powr<-2>(p);
      const auto _cse2 = _interp1 * _interp2;
      const auto _cse3 = -_interp2;
      const auto _cse4 = _cse3 + _interp5;
      const auto _cse5 = 50. * _cse4 * _den1 * powr<6>(k);
      const auto _cse6 = _cse5 + _interp4;
      const auto _cse7 = _cse6 * _interp3;
      const auto _cse8 = _cse2 + _cse7;
      const auto _cse9 = -_interp13;
      const auto _cse10 = _cse9 + _interp28;
      const auto _cse11 = 50. * _cse10;
      const auto _cse12 = _cse11 + _interp27;
      const auto _cse13 = _interp13 * _interp25;
      const auto _cse14 = _cse12 * _interp26;
      const auto _cse15 = _cse13 + _cse14;
      return 0.05555555555555555 * fma(_cse1, _cse8 * _den2 * _den4 * _den9 * DiFfRG::zacbc_num::tr1(fenv) * _interp17 * _interp18 * _interp19, fma(_cse1, _cse8 * _den2 * _den6 * _den9 * DiFfRG::zacbc_num::tr2(fenv) * _interp10 * _interp18 * _interp22, fma(-1., _cse1 * _cse15 * _den10 * _den3 * _den8 * DiFfRG::zacbc_num::tr4(fenv) * _interp23 * _interp30 * _interp31, fma(-1., _cse1 * _cse12 * _den11 * _den3 * _den7 * DiFfRG::zacbc_num::tr5(fenv) * _interp23 * _interp30 * _interp31 * _interp32, fma(-1., _cse1 * _den11 * _den3 * _den7 * DiFfRG::zacbc_num::tr5(fenv) * _interp13 * _interp23 * _interp30 * _interp31 * _interp34, fma(_cse1, _cse8 * _den2 * _den5 * _den6 * DiFfRG::zacbc_num::tr0(fenv) * _interp10 * _interp11 * _interp9, fma(-1., _cse1 * _cse15 * _den10 * _den3 * _den5 * DiFfRG::zacbc_num::tr3(fenv) * _interp23 * _interp24 * _interp9, 0.)))))));
    }

    static inline auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA)
    {
      return 0.;
    }private: static inline auto RB(const auto &k2, const auto &p2) { return REG::RB(k2, p2); }
    static inline auto RF(const auto &k2, const auto &p2) { return REG::RF(k2, p2); }
    static inline auto RBdot(const auto &k2, const auto &p2) { return REG::RBdot(k2, p2); }
    static inline auto RFdot(const auto &k2, const auto &p2) { return REG::RFdot(k2, p2); }
    static inline auto dq2RB(const auto &k2, const auto &p2) { return REG::dq2RB(k2, p2); }
    static inline auto dq2RF(const auto &k2, const auto &p2) { return REG::dq2RF(k2, p2); }
  };
}
using DiFfRG::ZAcbc_num_kernel;