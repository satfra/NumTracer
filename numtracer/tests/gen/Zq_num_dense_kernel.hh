#pragma once

#include "shim.hpp"
#include "numtracer/dense/dtensor.hpp"

namespace DiFfRG
{
  template<typename REG>
  class Zq_num_dense_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double renv[7] = {};
      renv[0] = -cos1 * l1;
      renv[1] = -sqrt(1. - powr<2>(cos1)) * l1;
      renv[2] = cos1 * l1;
      renv[3] = sqrt(1. - powr<2>(cos1)) * l1;
      renv[4] = p;
      renv[5] = 1.0 / (powr<2>(l1));
      renv[6] = 1.0 / (powr<2>(l1));
      const double _tr0_0 = numtracer::dense::contract_all(numtracer::dense::gamma_axis<0, 100, 101>(), numtracer::dense::gamma_axis<2, 103, 100>(), numtracer::dense::gamma_axis<1, 102, 103>(), numtracer::dense::gamma_axis<3, 101, 102>(), numtracer::dense::transverse_projector<12, 2, 0, 3, 5>(renv), numtracer::dense::transverse_projector<12, 13, 0, 3, 5>(renv), numtracer::dense::transverse_projector<15, 3, 2, 3, 6>(renv), numtracer::dense::transverse_projector<14, 15, 0, 3, 5>(renv), numtracer::dense::transverse_projector<14, 13, 2, 3, 6>(renv), numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<0, 2, 3>(renv)), numtracer::dense::vector<0, 4, 1>(renv)), numtracer::dense::vector<1, 4, 1>(renv)).scalar_value().re;
      constexpr double _tr0_1 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 4, 5>(), numtracer::dense::delta_adj<3, 6, 5>(), numtracer::dense::delta_adj<3, 6, 7>(), numtracer::dense::T<3, 4, 11, 8>(), numtracer::dense::delta_fund<3, 8, 9>(), numtracer::dense::delta_fund<3, 10, 11>(), numtracer::dense::T<3, 7, 9, 10>()).scalar_value().re;
      const double _tr1_0 = numtracer::dense::contract_all(numtracer::dense::identity<100, 101>(), numtracer::dense::gamma_axis<1, 101, 102>(), numtracer::dense::identity<102, 103>(), numtracer::dense::gamma_axis<2, 105, 100>(), numtracer::dense::gamma_axis<0, 104, 105>(), numtracer::dense::gamma_axis<3, 103, 104>(), numtracer::dense::transverse_projector<12, 3, 0, 3, 5>(renv), numtracer::dense::transverse_projector<12, 13, 0, 3, 5>(renv), numtracer::dense::transverse_projector<13, 2, 2, 3, 6>(renv), numtracer::dense::vector<0, 4, 1>(renv), numtracer::dense::add_all(numtracer::dense::vector<1, 2, 3>(renv), numtracer::dense::vector<1, 4, 1>(renv))).scalar_value().re;
      constexpr double _tr1_1 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 4, 5>(), numtracer::dense::T<3, 5, 11, 6>(), numtracer::dense::delta_fund<3, 6, 7>(), numtracer::dense::delta_fund<3, 7, 8>(), numtracer::dense::delta_fund<3, 8, 9>(), numtracer::dense::delta_fund<3, 10, 11>(), numtracer::dense::T<3, 4, 9, 10>()).scalar_value().re;
      const double _tr2_0 = numtracer::dense::contract_all(numtracer::dense::gamma_axis<0, 100, 101>(), numtracer::dense::gamma_axis<3, 101, 102>(), numtracer::dense::gamma_axis<1, 102, 103>(), numtracer::dense::gamma_axis<4, 105, 100>(), numtracer::dense::gamma_axis<2, 104, 105>(), numtracer::dense::gamma_axis<5, 103, 104>(), numtracer::dense::transverse_projector<14, 5, 0, 3, 5>(renv), numtracer::dense::transverse_projector<14, 15, 0, 3, 5>(renv), numtracer::dense::transverse_projector<15, 4, 2, 3, 6>(renv), numtracer::dense::add_all(numtracer::dense::vector<0, 2, 3>(renv), numtracer::dense::vector<0, 4, 1>(renv)), numtracer::dense::add_all(numtracer::dense::vector<1, 2, 3>(renv), numtracer::dense::vector<1, 4, 1>(renv)), numtracer::dense::vector<2, 4, 1>(renv), numtracer::dense::add_all(numtracer::dense::vector<3, 2, 3>(renv), numtracer::dense::vector<3, 4, 1>(renv))).scalar_value().re;
      constexpr double _tr2_1 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 6, 7>(), numtracer::dense::T<3, 7, 13, 8>(), numtracer::dense::delta_fund<3, 8, 9>(), numtracer::dense::delta_fund<3, 9, 10>(), numtracer::dense::delta_fund<3, 10, 11>(), numtracer::dense::delta_fund<3, 12, 13>(), numtracer::dense::T<3, 6, 11, 12>()).scalar_value().re;
      const auto _interp1 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp4 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp5 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(l1);
      const auto _interp7 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp8 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp9 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp10 = Zq(k);
      const auto _interp11 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp12 = Mq(sqrt(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp13 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cos1 * l1 * p + powr<2>(p)));
      const auto _interp14 = RFdot(powr<2>(k), powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp15 = RF(powr<2>(k), powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp16 = dtZq(k);
      const auto _interp17 = Zq(1.02 * k);
      const auto _interp18 = Zq(sqrt(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp2 * _interp3 + _interp6 * powr<2>(l1));
      const auto _den3 = powr<-1>(_interp2 * _interp3 + _interp6 * powr<2>(l1));
      const auto _den4 = powr<-1>(powr<2>(_interp8) + powr<2>(_interp10 * _interp9 + _interp11 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _den5 = powr<-2>(powr<2>(_interp12) + powr<2>(_interp10 * _interp15 + _interp18 * sqrt(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _cse1 = 2. * cos1 * l1 * p;
      const auto _cse4 = powr<-2>(p);
      const auto _cse2 = _cse1 + powr<-1>(_cse4) + powr<2>(l1);
      const auto _cse3 = sqrt(powr<-1>(_cse2));
      const auto _cse5 = -_cse3 * _interp10 * _interp14;
      const auto _cse6 = -_interp10;
      const auto _cse7 = _cse6 + _interp17;
      const auto _cse8 = 50. * _cse7;
      const auto _cse9 = _cse8 + _interp16;
      const auto _cse10 = -_cse3 * _cse9 * _interp15;
      const auto _cse11 = _cse10 + _cse5;
      return fma(-0.0833333333333333, _cse11 * _cse4 * _den3 * _den5 * powr<2>(_interp12) * powr<2>(_interp13) * _tr1_0 * _tr1_1, fma(0.0833333333333333, _cse11 * _cse4 * _den3 * _den5 * powr<2>(_interp13) * powr<2>(-_cse3 * _interp10 * _interp15 - _interp18) * _tr2_0 * _tr2_1, fma(0.0833333333333333, _cse4 * _den2 * _den4 * powr<2>(_interp7) * _tr0_0 * _tr0_1 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))) * (-_interp11 - _interp10 * _interp9 * sqrt(powr<-1>(powr<-1>(_cse4) + powr<2>(l1) - 2. * cos1 * l1 * p))), 0.)));
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
  };
}
using DiFfRG::Zq_num_dense_kernel;