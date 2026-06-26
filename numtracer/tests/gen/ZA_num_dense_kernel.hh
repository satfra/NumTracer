#pragma once

#include "shim.hpp"
#include "numtracer/dense/dtensor.hpp"

namespace DiFfRG
{
  template<typename REG>
  class ZA_num_dense_kernel
  {
    public:
    // clang-format off
static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
{
using namespace DiFfRG;
    // clang-format on

    using namespace DiFfRG::compute;
    using namespace numtracer;
    double renv[13] = {};
    renv[0] = -cos1 * l1;
    renv[1] = -sqrt(1. - powr<2>(cos1)) * l1;
    renv[2] = cos1 * l1;
    renv[3] = sqrt(1. - powr<2>(cos1)) * l1;
    renv[4] = -p;
    renv[5] = -cos1 * l1 + p;
    renv[6] = -sqrt(1. - powr<2>(cos1)) * l1;
    renv[7] = p;
    renv[8] = 1.0 / (powr<2>(l1));
    renv[9] = 1.0 / (powr<2>(l1));
    renv[10] = 1.0 / (powr<2>(p));
    renv[11] = 1.0 / (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
    renv[12] = 1.0 / (powr<2>(p));
    constexpr double _tr0_0 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 0, 1>(), numtracer::dense::f<3, 0, 5, 2>(), numtracer::dense::delta_adj<3, 2, 3>(), numtracer::dense::delta_adj<3, 4, 5>(), numtracer::dense::delta_adj<3, 6, 3>(), numtracer::dense::delta_adj<3, 6, 7>(), numtracer::dense::f<3, 1, 7, 4>()).scalar_value().re;// clang-format off

    const double _tr0_1 = numtracer::dense::contract_all(numtracer::dense::transverse_projector<8, 9, 0, 3, 8>(renv), numtracer::dense::transverse_projector<10, 9, 2, 3, 9>(renv), numtracer::dense::transverse_projector<10, 11, 0, 3, 8>(renv), numtracer::dense::add_all(numtracer::dense::contract_all(numtracer::dense::metric<15, 8>(), numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-2., 0.0}>(numtracer::dense::add_all(numtracer::dense::vector<12, 2, 3>(renv), numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<12, 7, 1>(renv)))), numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<12, 7, 1>(renv)))), numtracer::dense::contract_all(numtracer::dense::metric<12, 15>(), numtracer::dense::add_all(numtracer::dense::vector<8, 2, 3>(renv), numtracer::dense::scale<numtracer::Cx{-2., 0.0}>(numtracer::dense::vector<8, 7, 1>(renv)))), numtracer::dense::contract_all(numtracer::dense::metric<12, 8>(), numtracer::dense::add_all(numtracer::dense::vector<15, 2, 3>(renv), numtracer::dense::vector<15, 7, 1>(renv)))), numtracer::dense::transverse_projector<12, 13, 4, 1, 10>(renv), numtracer::dense::transverse_projector<14, 15, 5, 3, 11>(renv), numtracer::dense::add_all(numtracer::dense::contract_all(numtracer::dense::metric<13, 11>(), numtracer::dense::add_all(numtracer::dense::vector<14, 2, 3>(renv), numtracer::dense::vector<14, 7, 1>(renv))), numtracer::dense::contract_all(numtracer::dense::metric<13, 14>(), numtracer::dense::add_all(numtracer::dense::vector<11, 2, 3>(renv), numtracer::dense::scale<numtracer::Cx{-2., 0.0}>(numtracer::dense::vector<11, 7, 1>(renv)))), numtracer::dense::contract_all(numtracer::dense::metric<11, 14>(), numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-2., 0.0}>(numtracer::dense::vector<13, 2, 3>(renv)), numtracer::dense::vector<13, 7, 1>(renv))))).scalar_value().re;
    // clang-format on

    constexpr double _tr1_0 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 4, 5>(), numtracer::dense::f<3, 10, 5, 4>(), numtracer::dense::f<3, 10, 9, 6>(), numtracer::dense::delta_adj<3, 6, 7>(), numtracer::dense::delta_adj<3, 8, 7>(), numtracer::dense::delta_adj<3, 8, 9>()).scalar_value().re;
    const double _tr1_1 = numtracer::dense::contract_all(numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::contract_all(numtracer::dense::metric<0, 1>(), numtracer::dense::metric<2, 3>())), numtracer::dense::contract_all(numtracer::dense::metric<0, 3>(), numtracer::dense::metric<2, 1>())), numtracer::dense::transverse_projector<3, 11, 0, 3, 8>(renv), numtracer::dense::transverse_projector<12, 1, 0, 3, 8>(renv), numtracer::dense::transverse_projector<12, 11, 2, 3, 9>(renv), numtracer::dense::transverse_projector<0, 2, 4, 1, 10>(renv)).scalar_value().re;
    constexpr double _tr2_0 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 4, 5>(), numtracer::dense::f<3, 10, 4, 6>(), numtracer::dense::delta_adj<3, 6, 7>(), numtracer::dense::delta_adj<3, 8, 7>(), numtracer::dense::delta_adj<3, 8, 9>(), numtracer::dense::f<3, 10, 5, 9>()).scalar_value().re;
    const double _tr2_1 = numtracer::dense::contract_all(numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::contract_all(numtracer::dense::metric<0, 1>(), numtracer::dense::metric<2, 3>())), numtracer::dense::contract_all(numtracer::dense::metric<2, 0>(), numtracer::dense::metric<1, 3>())), numtracer::dense::transverse_projector<3, 11, 0, 3, 8>(renv), numtracer::dense::transverse_projector<12, 1, 0, 3, 8>(renv), numtracer::dense::transverse_projector<12, 11, 2, 3, 9>(renv), numtracer::dense::transverse_projector<0, 2, 4, 1, 10>(renv)).scalar_value().re;
    constexpr double _tr3_0 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 4, 5>(), numtracer::dense::f<3, 10, 4, 9>(), numtracer::dense::delta_adj<3, 8, 9>(), numtracer::dense::delta_adj<3, 8, 7>(), numtracer::dense::delta_adj<3, 6, 7>(), numtracer::dense::f<3, 10, 5, 6>()).scalar_value().re;
    const double _tr3_1 = numtracer::dense::contract_all(numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::contract_all(numtracer::dense::metric<0, 1>(), numtracer::dense::metric<2, 3>())), numtracer::dense::contract_all(numtracer::dense::metric<2, 0>(), numtracer::dense::metric<3, 1>())), numtracer::dense::transverse_projector<1, 11, 0, 3, 8>(renv), numtracer::dense::transverse_projector<12, 3, 0, 3, 8>(renv), numtracer::dense::transverse_projector<12, 11, 2, 3, 9>(renv), numtracer::dense::transverse_projector<0, 2, 4, 1, 10>(renv)).scalar_value().re;
    constexpr double _tr4_0 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 0, 1>(), numtracer::dense::f<3, 0, 3, 4>(), numtracer::dense::delta_adj<3, 2, 3>(), numtracer::dense::delta_adj<3, 4, 5>(), numtracer::dense::delta_adj<3, 7, 2>(), numtracer::dense::delta_adj<3, 6, 7>(), numtracer::dense::f<3, 1, 5, 6>()).scalar_value().re;
    const double _tr4_1 = numtracer::dense::contract_all(numtracer::dense::transverse_projector<8, 9, 4, 1, 10>(renv), numtracer::dense::transverse_projector<9, 10, 4, 1, 10>(renv), numtracer::dense::transverse_projector<8, 11, 7, 1, 12>(renv), numtracer::dense::vector<11, 2, 3>(renv), numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<10, 2, 3>(renv)), numtracer::dense::vector<10, 7, 1>(renv))).scalar_value().re;
    const double _tr5_0 = numtracer::dense::contract_all(numtracer::dense::identity<100, 101>(), numtracer::dense::gamma_axis<2, 105, 100>(), numtracer::dense::gamma_axis<0, 104, 105>(), numtracer::dense::gamma_axis<1, 103, 104>(), numtracer::dense::identity<102, 103>(), numtracer::dense::gamma_axis<3, 101, 102>(), numtracer::dense::transverse_projector<13, 3, 4, 1, 10>(renv), numtracer::dense::transverse_projector<12, 13, 4, 1, 10>(renv), numtracer::dense::transverse_projector<12, 2, 7, 1, 12>(renv), numtracer::dense::vector<0, 2, 3>(renv), numtracer::dense::vector<1, 2, 3>(renv)).scalar_value().re;
    constexpr double _tr5_1 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 4, 5>(), numtracer::dense::T<3, 4, 7, 8>(), numtracer::dense::delta_fund<3, 6, 7>(), numtracer::dense::delta_fund<3, 8, 9>(), numtracer::dense::delta_fund<3, 11, 6>(), numtracer::dense::delta_fund<3, 10, 11>(), numtracer::dense::T<3, 5, 9, 10>()).scalar_value().re;
    const double _tr6_0 = numtracer::dense::contract_all(numtracer::dense::identity<100, 101>(), numtracer::dense::gamma_axis<1, 105, 100>(), numtracer::dense::gamma_axis<0, 104, 105>(), numtracer::dense::gamma_axis<2, 101, 102>(), numtracer::dense::identity<102, 103>(), numtracer::dense::gamma_axis<3, 103, 104>(), numtracer::dense::transverse_projector<13, 3, 4, 1, 10>(renv), numtracer::dense::transverse_projector<12, 13, 4, 1, 10>(renv), numtracer::dense::transverse_projector<12, 2, 7, 1, 12>(renv), numtracer::dense::vector<0, 2, 3>(renv), numtracer::dense::vector<1, 2, 3>(renv)).scalar_value().re;
    constexpr double _tr6_1 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 4, 5>(), numtracer::dense::T<3, 4, 7, 8>(), numtracer::dense::delta_fund<3, 6, 7>(), numtracer::dense::delta_fund<3, 8, 9>(), numtracer::dense::delta_fund<3, 11, 6>(), numtracer::dense::delta_fund<3, 10, 11>(), numtracer::dense::T<3, 5, 9, 10>()).scalar_value().re;
    const double _tr7_0 = numtracer::dense::contract_all(numtracer::dense::identity<100, 101>(), numtracer::dense::gamma_axis<1, 103, 100>(), numtracer::dense::identity<102, 103>(), numtracer::dense::gamma_axis<2, 101, 104>(), numtracer::dense::gamma_axis<0, 104, 105>(), numtracer::dense::gamma_axis<3, 105, 102>(), numtracer::dense::transverse_projector<13, 3, 4, 1, 10>(renv), numtracer::dense::transverse_projector<12, 13, 4, 1, 10>(renv), numtracer::dense::transverse_projector<12, 2, 7, 1, 12>(renv), numtracer::dense::vector<1, 2, 3>(renv), numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<0, 2, 3>(renv)), numtracer::dense::vector<0, 7, 1>(renv))).scalar_value().re;
    constexpr double _tr7_1 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 4, 5>(), numtracer::dense::T<3, 4, 7, 8>(), numtracer::dense::delta_fund<3, 6, 7>(), numtracer::dense::delta_fund<3, 8, 9>(), numtracer::dense::delta_fund<3, 11, 6>(), numtracer::dense::delta_fund<3, 10, 11>(), numtracer::dense::T<3, 5, 9, 10>()).scalar_value().re;
    const double _tr8_0 = numtracer::dense::contract_all(numtracer::dense::gamma_axis<0, 100, 101>(), numtracer::dense::gamma_axis<3, 103, 100>(), numtracer::dense::gamma_axis<1, 102, 103>(), numtracer::dense::gamma_axis<4, 101, 104>(), numtracer::dense::gamma_axis<2, 104, 105>(), numtracer::dense::gamma_axis<5, 105, 102>(), numtracer::dense::transverse_projector<15, 5, 4, 1, 10>(renv), numtracer::dense::transverse_projector<14, 15, 4, 1, 10>(renv), numtracer::dense::transverse_projector<14, 4, 7, 1, 12>(renv), numtracer::dense::vector<0, 2, 3>(renv), numtracer::dense::vector<1, 2, 3>(renv), numtracer::dense::vector<3, 2, 3>(renv), numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<2, 2, 3>(renv)), numtracer::dense::vector<2, 7, 1>(renv))).scalar_value().re;
    constexpr double _tr8_1 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 6, 7>(), numtracer::dense::T<3, 6, 9, 10>(), numtracer::dense::delta_fund<3, 8, 9>(), numtracer::dense::delta_fund<3, 10, 11>(), numtracer::dense::delta_fund<3, 13, 8>(), numtracer::dense::delta_fund<3, 12, 13>(), numtracer::dense::T<3, 7, 11, 12>()).scalar_value().re;
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
    const auto _interp17 = Mq(l1);
    const auto _interp18 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
    const auto _interp19 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
    const auto _interp20 = RFdot(powr<2>(k), powr<2>(l1));
    const auto _interp21 = Zq(k);
    const auto _interp22 = RF(powr<2>(k), powr<2>(l1));
    const auto _interp23 = dtZq(k);
    const auto _interp24 = Zq(1.02 * k);
    const auto _interp25 = Zq(l1);
    const auto _interp26 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
    const auto _interp27 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
    const auto _den1 = powr<-1>(1. + powr<6>(k));
    const auto _den2 = powr<-2>(_interp2 * _interp3 + _interp6 * powr<2>(l1));
    const auto _den3 = powr<-2>(_interp12 * _interp3 + _interp15 * powr<2>(l1));
    const auto _den4 = powr<-1>(_interp2 * _interp7 + _interp8 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
    const auto _den5 = -powr<-1>(_interp12 * _interp7 + _interp16 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
    const auto _den6 = powr<-2>(powr<2>(_interp17) + powr<2>(_interp21 * _interp22 + _interp25 * sqrt(powr<2>(l1))));
    const auto _den7 = powr<-1>(powr<2>(_interp18) + powr<2>(_interp21 * _interp26 + _interp27 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
    const auto _cse1 = sqrt(powr<-1>(powr<2>(l1)));
    const auto _cse2 = powr<-2>(p);
    const auto _cse3 = -_interp25;
    const auto _cse4 = -_cse1 * _interp21 * _interp22;
    const auto _cse5 = _cse3 + _cse4;
    const auto _cse6 = -_cse1 * _interp20 * _interp21;
    const auto _cse7 = -_interp21;
    const auto _cse8 = _cse7 + _interp24;
    const auto _cse9 = 50. * _cse8;
    const auto _cse10 = _cse9 + _interp23;
    const auto _cse11 = -_cse1 * _cse10 * _interp22;
    const auto _cse12 = _cse11 + _cse6;
    const auto _cse13 = -_interp27;
    const auto _cse14 = -2. * cos1 * l1 * p;
    const auto _cse15 = _cse14 + powr<-1>(_cse2) + powr<2>(l1);
    const auto _cse16 = sqrt(powr<-1>(_cse15));
    const auto _cse17 = -_cse16 * _interp21 * _interp26;
    const auto _cse18 = _cse13 + _cse17;
    const auto _cse19 = _interp1 * _interp2;
    const auto _cse20 = -_interp2;
    const auto _cse21 = _cse20 + _interp5;
    const auto _cse22 = 50. * _cse21 * _den1 * powr<6>(k);
    const auto _cse23 = _cse22 + _interp4;
    const auto _cse24 = _cse23 * _interp3;
    const auto _cse25 = _cse19 + _cse24;// clang-format off

    return fma(-0.04166666666666666, _cse2 * _cse25 * _den2 * _den4 * powr<2>(_interp9) * _tr0_0 * _tr0_1, fma(-0.02083333333333333, _cse2 * _cse25 * _den2 * _interp10 * _tr1_0 * _tr1_1, fma(-0.02083333333333333, _cse2 * _cse25 * _den2 * _interp10 * _tr2_0 * _tr2_1, fma(-0.02083333333333333, _cse2 * _cse25 * _den2 * _interp10 * _tr3_0 * _tr3_1, fma(0.0833333333333333, _cse2 * _den3 * _den5 * powr<2>(_interp11) * (_interp1 * _interp12 + (_interp13 + 50. * (-_interp12 + _interp14)) * _interp3) * _tr4_0 * _tr4_1, fma(0.1666666666666667, _cse12 * _cse2 * _cse5 * _den6 * _den7 * _interp17 * _interp18 * powr<2>(_interp19) * _tr5_0 * _tr5_1, fma(0.1666666666666667, _cse12 * _cse2 * _cse5 * _den6 * _den7 * _interp17 * _interp18 * powr<2>(_interp19) * _tr6_0 * _tr6_1, fma(-0.1666666666666667, _cse12 * _cse18 * _cse2 * _den6 * _den7 * powr<2>(_interp17) * powr<2>(_interp19) * _tr7_0 * _tr7_1, fma(0.1666666666666667, _cse12 * _cse18 * _cse2 * powr<2>(_cse5) * _den6 * _den7 * powr<2>(_interp19) * _tr8_0 * _tr8_1, 0.)))))))));
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
};
}
using DiFfRG::ZA_num_dense_kernel;