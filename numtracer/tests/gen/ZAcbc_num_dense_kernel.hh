#pragma once

#include "shim.hpp"
#include "numtracer/dense/dtensor.hpp"

namespace DiFfRG
{
  template<typename REG>
  class ZAcbc_num_dense_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double renv[24] = {};
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      renv[0] = -cos1 * l1;
      renv[1] = -sqrt(1. - powr<2>(cos1)) * cos2 * l1;
      renv[2] = -sqrt((-1. + powr<2>(cos1)) * (-1. + powr<2>(cos2))) * l1;
      renv[3] = cos1 * l1;
      renv[4] = sqrt(1. - powr<2>(cos1)) * cos2 * l1;
      renv[5] = sqrt((-1. + powr<2>(cos1)) * (-1. + powr<2>(cos2))) * l1;
      renv[6] = -cos1 * l1 - p;
      renv[7] = -sqrt(1. - powr<2>(cos1)) * cos2 * l1;
      renv[8] = -sqrt((-1. + powr<2>(cos1)) * (-1. + powr<2>(cos2))) * l1;
      renv[9] = p;
      renv[10] = cos1 * l1 + p;
      renv[11] = sqrt(1. - powr<2>(cos1)) * cos2 * l1;
      renv[12] = sqrt((-1. + powr<2>(cos1)) * (-1. + powr<2>(cos2))) * l1;
      renv[13] = -0.5 * p;
      renv[14] = 0.866025403784439 * p;
      renv[15] = -cos1 * l1 + p;
      renv[16] = -sqrt(1. - powr<2>(cos1)) * cos2 * l1;
      renv[17] = -sqrt((-1. + powr<2>(cos1)) * (-1. + powr<2>(cos2))) * l1;
      renv[18] = 1.0 / (powr<2>(l1));
      renv[19] = 1.0 / (powr<2>(l1));
      renv[20] = 1.0 / (powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p));
      renv[21] = 1.0 / (powr<2>(p));
      renv[22] = 1.0 / (powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p));
      renv[23] = 1.0 / (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));// clang-format off

      const double _tr0_0 = numtracer::dense::contract_all(numtracer::dense::transverse_projector<11, 12, 0, 7, 18>(renv), numtracer::dense::transverse_projector<11, 13, 0, 7, 18>(renv), numtracer::dense::transverse_projector<14, 13, 3, 7, 19>(renv), numtracer::dense::transverse_projector<14, 15, 0, 7, 18>(renv), numtracer::dense::add_all(numtracer::dense::contract_all(numtracer::dense::metric<15, 16>(), numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-2., 0.0}>(numtracer::dense::vector<18, 3, 7>(renv)), numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<18, 9, 1>(renv)))), numtracer::dense::contract_all(numtracer::dense::metric<18, 15>(), numtracer::dense::add_all(numtracer::dense::vector<16, 3, 7>(renv), numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<16, 9, 1>(renv)))), numtracer::dense::contract_all(numtracer::dense::metric<18, 16>(), numtracer::dense::add_all(numtracer::dense::vector<15, 3, 7>(renv), numtracer::dense::scale<numtracer::Cx{2., 0.0}>(numtracer::dense::vector<15, 9, 1>(renv))))), numtracer::dense::transverse_projector<16, 17, 6, 7, 20>(renv), numtracer::dense::transverse_projector<18, 19, 9, 1, 21>(renv), numtracer::dense::transverse_projector<17, 20, 10, 7, 22>(renv), numtracer::dense::vector<19, 13, 3>(renv), numtracer::dense::vector<12, 13, 3>(renv), numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<20, 3, 7>(renv)), numtracer::dense::vector<20, 13, 3>(renv))).scalar_value().re;
      // clang-format on

      constexpr double _tr0_1 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 0, 1>(), numtracer::dense::delta_adj<3, 6, 1>(), numtracer::dense::delta_adj<3, 6, 7>(), numtracer::dense::f<3, 8, 7, 4>(), numtracer::dense::delta_adj<3, 4, 5>(), numtracer::dense::f<3, 8, 9, 10>(), numtracer::dense::f<3, 0, 9, 2>(), numtracer::dense::delta_adj<3, 2, 3>(), numtracer::dense::f<3, 5, 3, 10>()).scalar_value().re;// clang-format off

      const double _tr1_0 = numtracer::dense::contract_all(numtracer::dense::transverse_projector<11, 12, 0, 7, 18>(renv), numtracer::dense::transverse_projector<13, 12, 3, 7, 19>(renv), numtracer::dense::transverse_projector<13, 14, 0, 7, 18>(renv), numtracer::dense::transverse_projector<14, 15, 3, 7, 19>(renv), numtracer::dense::add_all(numtracer::dense::contract_all(numtracer::dense::metric<20, 11>(), numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-2., 0.0}>(numtracer::dense::add_all(numtracer::dense::vector<16, 3, 7>(renv), numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<16, 9, 1>(renv)))), numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<16, 9, 1>(renv)))), numtracer::dense::contract_all(numtracer::dense::metric<16, 20>(), numtracer::dense::add_all(numtracer::dense::vector<11, 3, 7>(renv), numtracer::dense::scale<numtracer::Cx{-2., 0.0}>(numtracer::dense::vector<11, 9, 1>(renv)))), numtracer::dense::contract_all(numtracer::dense::metric<16, 11>(), numtracer::dense::add_all(numtracer::dense::vector<20, 3, 7>(renv), numtracer::dense::vector<20, 9, 1>(renv)))), numtracer::dense::transverse_projector<16, 17, 9, 1, 21>(renv), numtracer::dense::transverse_projector<18, 20, 15, 7, 23>(renv), numtracer::dense::transverse_projector<18, 19, 15, 7, 23>(renv), numtracer::dense::vector<17, 13, 3>(renv), numtracer::dense::vector<19, 13, 3>(renv), numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<15, 3, 7>(renv)), numtracer::dense::vector<15, 9, 1>(renv), numtracer::dense::vector<15, 13, 3>(renv))).scalar_value().re;
      // clang-format on

      constexpr double _tr1_1 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 0, 1>(), numtracer::dense::delta_adj<3, 6, 1>(), numtracer::dense::delta_adj<3, 6, 7>(), numtracer::dense::f<3, 8, 3, 0>(), numtracer::dense::delta_adj<3, 2, 3>(), numtracer::dense::f<3, 8, 9, 10>(), numtracer::dense::f<3, 2, 9, 4>(), numtracer::dense::delta_adj<3, 4, 5>(), numtracer::dense::f<3, 7, 5, 10>()).scalar_value().re;
      constexpr double _tr2_0 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 0, 1>(), numtracer::dense::delta_adj<3, 6, 1>(), numtracer::dense::delta_adj<3, 6, 7>(), numtracer::dense::f<3, 0, 9, 2>(), numtracer::dense::delta_adj<3, 2, 3>(), numtracer::dense::f<3, 8, 9, 10>(), numtracer::dense::f<3, 8, 3, 4>(), numtracer::dense::delta_adj<3, 4, 5>(), numtracer::dense::f<3, 7, 5, 10>()).scalar_value().re;
      const double _tr2_1 = numtracer::dense::contract_all(numtracer::dense::transverse_projector<11, 12, 0, 7, 18>(renv), numtracer::dense::transverse_projector<11, 13, 0, 7, 18>(renv), numtracer::dense::transverse_projector<14, 13, 3, 7, 19>(renv), numtracer::dense::transverse_projector<14, 15, 0, 7, 18>(renv), numtracer::dense::transverse_projector<15, 16, 3, 7, 19>(renv), numtracer::dense::vector<12, 13, 3>(renv), numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<16, 3, 7>(renv)), numtracer::dense::vector<16, 9, 1>(renv), numtracer::dense::vector<16, 13, 3>(renv))).scalar_value().re;
      const double _tr2_2 = numtracer::dense::contract_all(numtracer::dense::transverse_projector<17, 18, 9, 1, 21>(renv), numtracer::dense::transverse_projector<17, 19, 9, 1, 21>(renv), numtracer::dense::vector<18, 13, 3>(renv), numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<19, 3, 7>(renv)), numtracer::dense::vector<19, 13, 3>(renv))).scalar_value().re;// clang-format off

      const double _tr3_0 = numtracer::dense::contract_all(numtracer::dense::transverse_projector<11, 12, 0, 7, 18>(renv), numtracer::dense::transverse_projector<11, 13, 0, 7, 18>(renv), numtracer::dense::add_all(numtracer::dense::contract_all(numtracer::dense::metric<13, 14>(), numtracer::dense::add_all(numtracer::dense::scale<numtracer::Cx{-2., 0.0}>(numtracer::dense::vector<16, 3, 7>(renv)), numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<16, 9, 1>(renv)))), numtracer::dense::contract_all(numtracer::dense::metric<16, 13>(), numtracer::dense::add_all(numtracer::dense::vector<14, 3, 7>(renv), numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<14, 9, 1>(renv)))), numtracer::dense::contract_all(numtracer::dense::metric<16, 14>(), numtracer::dense::add_all(numtracer::dense::vector<13, 3, 7>(renv), numtracer::dense::scale<numtracer::Cx{2., 0.0}>(numtracer::dense::vector<13, 9, 1>(renv))))), numtracer::dense::transverse_projector<14, 15, 6, 7, 20>(renv), numtracer::dense::transverse_projector<16, 17, 9, 1, 21>(renv), numtracer::dense::transverse_projector<15, 18, 10, 7, 22>(renv), numtracer::dense::vector<17, 13, 3>(renv), numtracer::dense::vector<18, 13, 3>(renv), numtracer::dense::add_all(numtracer::dense::vector<12, 3, 7>(renv), numtracer::dense::vector<12, 9, 1>(renv), numtracer::dense::vector<12, 13, 3>(renv))).scalar_value().re;
      // clang-format on

      constexpr double _tr3_1 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 0, 1>(), numtracer::dense::delta_adj<3, 1, 2>(), numtracer::dense::delta_adj<3, 2, 7>(), numtracer::dense::f<3, 4, 9, 0>(), numtracer::dense::delta_adj<3, 3, 4>(), numtracer::dense::f<3, 8, 9, 10>(), numtracer::dense::f<3, 8, 6, 3>(), numtracer::dense::delta_adj<3, 5, 6>(), numtracer::dense::f<3, 5, 7, 10>()).scalar_value().re;
      constexpr double _tr4_0 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 0, 1>(), numtracer::dense::delta_adj<3, 1, 2>(), numtracer::dense::delta_adj<3, 2, 7>(), numtracer::dense::f<3, 8, 4, 0>(), numtracer::dense::delta_adj<3, 3, 4>(), numtracer::dense::f<3, 8, 9, 10>(), numtracer::dense::f<3, 6, 9, 3>(), numtracer::dense::delta_adj<3, 5, 6>(), numtracer::dense::f<3, 5, 7, 10>()).scalar_value().re;
      const double _tr4_1 = numtracer::dense::contract_all(numtracer::dense::transverse_projector<11, 12, 0, 7, 18>(renv), numtracer::dense::transverse_projector<11, 13, 0, 7, 18>(renv), numtracer::dense::transverse_projector<13, 14, 3, 7, 19>(renv), numtracer::dense::vector<14, 13, 3>(renv), numtracer::dense::add_all(numtracer::dense::vector<12, 3, 7>(renv), numtracer::dense::vector<12, 9, 1>(renv), numtracer::dense::vector<12, 13, 3>(renv))).scalar_value().re;
      const double _tr4_2 = numtracer::dense::contract_all(numtracer::dense::transverse_projector<15, 16, 9, 1, 21>(renv), numtracer::dense::transverse_projector<15, 17, 9, 1, 21>(renv), numtracer::dense::vector<16, 13, 3>(renv), numtracer::dense::add_all(numtracer::dense::vector<17, 3, 7>(renv), numtracer::dense::vector<17, 13, 3>(renv))).scalar_value().re;
      constexpr double _tr5_0 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 0, 1>(), numtracer::dense::delta_adj<3, 1, 2>(), numtracer::dense::delta_adj<3, 2, 7>(), numtracer::dense::f<3, 8, 7, 5>(), numtracer::dense::delta_adj<3, 5, 6>(), numtracer::dense::f<3, 8, 9, 10>(), numtracer::dense::f<3, 4, 9, 0>(), numtracer::dense::delta_adj<3, 3, 4>(), numtracer::dense::f<3, 3, 6, 10>()).scalar_value().re;
      const double _tr5_1 = numtracer::dense::contract_all(numtracer::dense::transverse_projector<11, 12, 0, 7, 18>(renv), numtracer::dense::transverse_projector<11, 13, 0, 7, 18>(renv), numtracer::dense::transverse_projector<13, 14, 3, 7, 19>(renv), numtracer::dense::vector<14, 13, 3>(renv), numtracer::dense::add_all(numtracer::dense::vector<12, 3, 7>(renv), numtracer::dense::vector<12, 9, 1>(renv), numtracer::dense::vector<12, 13, 3>(renv))).scalar_value().re;
      const double _tr5_2 = numtracer::dense::contract_all(numtracer::dense::transverse_projector<15, 16, 9, 1, 21>(renv), numtracer::dense::transverse_projector<15, 17, 9, 1, 21>(renv), numtracer::dense::vector<16, 13, 3>(renv), numtracer::dense::add_all(numtracer::dense::vector<17, 3, 7>(renv), numtracer::dense::vector<17, 13, 3>(renv))).scalar_value().re;
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
      return 0.05555555555555555 * fma(_cse1, _cse8 * _den2 * _den5 * _den6 * _interp10 * _interp11 * _interp9 * _tr0_0 * _tr0_1, fma(_cse1, _cse8 * _den2 * _den4 * _den9 * _interp17 * _interp18 * _interp19 * _tr1_0 * _tr1_1, fma(_cse1, _cse8 * _den2 * _den6 * _den9 * _interp10 * _interp18 * _interp22 * _tr2_0 * _tr2_1 * _tr2_2, fma(-1., _cse1 * _cse15 * _den10 * _den3 * _den5 * _interp23 * _interp24 * _interp9 * _tr3_0 * _tr3_1, fma(-1., _cse1 * _cse15 * _den10 * _den3 * _den8 * _interp23 * _interp30 * _interp31 * _tr4_0 * _tr4_1 * _tr4_2, fma(-1., _cse1 * _cse12 * _den11 * _den3 * _den7 * _interp23 * _interp30 * _interp31 * _interp32 * _tr5_0 * _tr5_1 * _tr5_2, fma(-1., _cse1 * _den11 * _den3 * _den7 * _interp13 * _interp23 * _interp30 * _interp31 * _interp34 * _tr5_0 * _tr5_1 * _tr5_2, 0.)))))));
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
using DiFfRG::ZAcbc_num_dense_kernel;