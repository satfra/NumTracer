#pragma once

#include "shim.hpp"
#include "numtracer/dense/dtensor.hpp"

namespace DiFfRG
{
  template<typename REG>
  class Flavour_split_num_dense_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Gu, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Gd, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GuDot, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GdDot)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double renv[3] = {};
      renv[0] = cos1 * l1;
      renv[1] = sqrt(1. - powr<2>(cos1)) * l1;
      renv[2] = p;
      const double _tr0_0 = numtracer::dense::contract_all(numtracer::dense::gamma_axis<0, 100, 101>(), numtracer::dense::gamma_axis<3, 105, 100>(), numtracer::dense::gamma_axis<2, 104, 105>(), numtracer::dense::gamma5<101, 102>(), numtracer::dense::gamma_axis<1, 102, 103>(), numtracer::dense::gamma5<103, 104>(), numtracer::dense::vector<0, 0, 3>(renv), numtracer::dense::vector<2, 0, 3>(renv), numtracer::dense::vector<3, 0, 3>(renv), numtracer::dense::add_all(numtracer::dense::vector<1, 0, 3>(renv), numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<1, 2, 1>(renv)))).scalar_value().re;
      constexpr double _tr0_1 = numtracer::dense::contract_all(numtracer::dense::delta_fund<3, 4, 5>(), numtracer::dense::delta_fund<3, 5, 6>(), numtracer::dense::delta_fund<3, 6, 7>(), numtracer::dense::delta_fund<3, 7, 8>(), numtracer::dense::delta_fund<3, 8, 9>(), numtracer::dense::delta_fund<3, 9, 4>()).scalar_value().re;
      const double _tr1_0 = numtracer::dense::contract_all(numtracer::dense::gamma_axis<0, 100, 101>(), numtracer::dense::gamma_axis<3, 105, 100>(), numtracer::dense::gamma_axis<2, 104, 105>(), numtracer::dense::gamma5<101, 102>(), numtracer::dense::gamma_axis<1, 102, 103>(), numtracer::dense::gamma5<103, 104>(), numtracer::dense::vector<0, 0, 3>(renv), numtracer::dense::vector<2, 0, 3>(renv), numtracer::dense::vector<3, 0, 3>(renv), numtracer::dense::add_all(numtracer::dense::vector<1, 0, 3>(renv), numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<1, 2, 1>(renv)))).scalar_value().re;
      constexpr double _tr1_1 = numtracer::dense::contract_all(numtracer::dense::delta_fund<3, 4, 5>(), numtracer::dense::delta_fund<3, 5, 6>(), numtracer::dense::delta_fund<3, 6, 7>(), numtracer::dense::delta_fund<3, 7, 8>(), numtracer::dense::delta_fund<3, 8, 9>(), numtracer::dense::delta_fund<3, 9, 4>()).scalar_value().re;
      const auto _interp1 = Gd(sqrt(powr<2>(l1)));
      const auto _interp2 = Gd(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp3 = GdDot(sqrt(powr<2>(l1)));
      const auto _interp4 = Gu(sqrt(powr<2>(l1)));
      const auto _interp5 = Gu(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp6 = GuDot(sqrt(powr<2>(l1)));
      return 2. * fma(powr<2>(_interp1), _interp2 * _interp3 * _tr0_0 * _tr0_1, fma(powr<2>(_interp4), _interp5 * _interp6 * _tr1_0 * _tr1_1, 0.));
    }

    static inline auto constant(const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Gu, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Gd, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GuDot, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GdDot)
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
using DiFfRG::Flavour_split_num_dense_kernel;