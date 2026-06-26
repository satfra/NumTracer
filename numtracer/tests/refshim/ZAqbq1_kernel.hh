#pragma once

#include "shim.hpp"

namespace DiFfRG
{
  template <typename _Regulator> class ZAqbq1_kernel
  {
  public:
    using Regulator = _Regulator;

    // clang-format off
static KOKKOS_FORCEINLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const double& Nf, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
{
using namespace DiFfRG;
      // clang-format on
      using namespace DiFfRG::compute;
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      // clang-format off
using _T = decltype(0.5 * powr<-1>(1. + powr<6>(k)) * powr<-1>(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)) *
sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) *
((-6. + 6. * powr<2>(cosl1p1)) * powr<4>(l1) +
(-9. * cosl1p1 + 10. * powr<3>(cosl1p1) + 6. * cosl1p2 - 4. * powr<2>(cosl1p1) * cosl1p2) *
powr<3>(l1) * p +
(-11. + 13. * powr<2>(cosl1p1) + 8. * cosl1p1 * cosl1p2 -
4. * powr<3>(cosl1p1) * cosl1p2) *
powr<2>(l1) * powr<2>(p) +
(-2. * cosl1p1 + 10. * cosl1p2 - 14. * powr<2>(cosl1p1) * cosl1p2) * l1 * powr<3>(p) +
(-3. - 6. * cosl1p1 * cosl1p2) * powr<4>(p)) *
((1. + 1. * powr<6>(k)) * dtZA(pow(1. + powr<6>(k), 0.16666666666666666667)) *
RB(powr<2>(k), powr<2>(l1)) +
(1. + powr<6>(k)) * RBdot(powr<2>(k), powr<2>(l1)) *
ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
powr<6>(k) * RB(powr<2>(k), powr<2>(l1)) *
(-50. * ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
50. * ZA(1.02 * pow(1. + powr<6>(k), 0.16666666666666666667)))) *
powr<-2>(RB(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
powr<2>(l1) * ZA(l1)) *
powr<-1>(RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)) *
ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)) *
ZA(sqrt(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)))) *
ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p1 * l1 * p + powr<2>(p))) *
ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) +
0.6666666666666666 * (cosl1p1 - cosl1p2) * l1 * p + powr<2>(p))) *
ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p))) *
(RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) +
sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) *
Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)))) *
powr<-1>(powr<2>(Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)))) +
powr<2>(RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) *
powr<2>(Zq(k)) +
2. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) *
RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) *
Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) +
(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) *
powr<2>(Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))))) +
0.5 * powr<-1>(1. + powr<6>(k)) * powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) *
powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)) *
sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
((-6. + 6. * powr<2>(cosl1p1)) * powr<4>(l1) +
(15. * cosl1p1 - 14. * powr<3>(cosl1p1) + 6. * cosl1p2 - 4. * powr<2>(cosl1p1) * cosl1p2) *
powr<3>(l1) * p +
(-11. + 5. * powr<2>(cosl1p1) + 4. * powr<4>(cosl1p1) - 8. * cosl1p1 * cosl1p2 +
4. * powr<3>(cosl1p1) * cosl1p2) *
powr<2>(l1) * powr<2>(p) +
(12. * cosl1p1 - 14. * powr<3>(cosl1p1) + 10. * cosl1p2 -
14. * powr<2>(cosl1p1) * cosl1p2) *
l1 * powr<3>(p) +
(-3. + 6. * powr<2>(cosl1p1) + 6. * cosl1p1 * cosl1p2) * powr<4>(p)) *
((1. + 1. * powr<6>(k)) * dtZA(pow(1. + powr<6>(k), 0.16666666666666666667)) *
RB(powr<2>(k), powr<2>(l1)) +
(1. + powr<6>(k)) * RBdot(powr<2>(k), powr<2>(l1)) *
ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
powr<6>(k) * RB(powr<2>(k), powr<2>(l1)) *
(-50. * ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
50. * ZA(1.02 * pow(1. + powr<6>(k), 0.16666666666666666667)))) *
powr<-2>(RB(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
powr<2>(l1) * ZA(l1)) *
powr<-1>(RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) *
ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) *
ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))) *
ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p))) *
ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) -
0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
(RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zq(k) +
sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))) *
powr<-1>(
powr<2>(Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))) +
powr<2>(RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
powr<2>(Zq(k)) +
2. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
Zq(k) * Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) +
(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)) *
powr<2>(Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))))) +
0.05555555555555555 *
(RBdot(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
RB(powr<2>(k), powr<2>(l1)) *
(dtZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
50. * powr<6>(k) * powr<-1>(1. + powr<6>(k)) *
(-ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
ZA(1.02 * pow(1. + powr<6>(k), 0.16666666666666666667))))) *
powr<-2>(RB(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
powr<2>(l1) * ZA(l1)) *
ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p))) *
ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) -
0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p))) *
powr<-1>(powr<2>(Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)))) +
powr<2>(RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) *
powr<2>(Zq(k)) +
2. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) *
RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) *
Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) +
(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) *
powr<2>(Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))))) *
((5. - 2. * powr<2>(cosl1p1)) * Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) *
Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) +
((-1. + 4. * powr<2>(cosl1p1)) * powr<2>(l1) + (cosl1p1 + 2. * cosl1p2) * l1 * p -
(-1. + 2. * powr<2>(cosl1p1) + 4. * cosl1p1 * cosl1p2 + 4. * powr<2>(cosl1p2)) *
powr<2>(p)) *
(-sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) *
RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) -
Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)))) *
(-sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
Zq(k) -
Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))))) *
powr<-1>(
powr<2>(Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))) +
powr<2>(RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
powr<2>(Zq(k)) +
2. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
Zq(k) * Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) +
(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
powr<2>(Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))))) +
1.5 * powr<-1>(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)) *
powr<-1>(powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p2 * l1 * p + powr<2>(p)) *
sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
((-1. + 1. * powr<2>(cosl1p1)) * powr<4>(l1) +
(-2.5 * cosl1p1 + 2.333333333333333 * powr<3>(cosl1p1) - cosl1p2 +
0.6666666666666666 * powr<2>(cosl1p1) * cosl1p2) *
powr<3>(l1) * p +
(-1.833333333333333 + 0.833333333333333 * powr<2>(cosl1p1) +
0.6666666666666666 * powr<4>(cosl1p1) - 1.333333333333333 * cosl1p1 * cosl1p2 +
0.6666666666666666 * powr<3>(cosl1p1) * cosl1p2) *
powr<2>(l1) * powr<2>(p) +
(-2. * cosl1p1 + 2.333333333333333 * powr<3>(cosl1p1) - 1.666666666666667 * cosl1p2 +
2.333333333333333 * powr<2>(cosl1p1) * cosl1p2) *
l1 * powr<3>(p) +
(-0.5 + 1. * powr<2>(cosl1p1) + 1. * cosl1p1 * cosl1p2) * powr<4>(p)) *
powr<-1>(RB(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
powr<2>(l1) * ZA(l1)) *
powr<-1>(RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)) *
ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)) *
ZA(sqrt(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)))) *
ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p1 * l1 * p + powr<2>(p))) *
ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) +
0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
(dtZq(k) * RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) +
RFdot(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zq(k) +
RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
(-50. * Zq(k) + 50. * Zq(1.02 * k))) *
(-0.5 * powr<2>(Mq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))) +
0.5 *
powr<2>(RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
powr<2>(Zq(k)) +
1. * sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zq(k) *
Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) +
(0.5 * powr<2>(l1) + 1. * cosl1p1 * l1 * p + 1. * cosl1p2 * l1 * p + 0.5 * powr<2>(p)) *
powr<2>(Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))))) *
powr<-2>(
0.5 * powr<2>(Mq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))) +
0.5 *
powr<2>(RF(powr<2>(k),
powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
powr<2>(Zq(k)) +
1. * sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
Zq(k) * Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) +
(0.5 * powr<2>(l1) + 1. * cosl1p1 * l1 * p + 1. * cosl1p2 * l1 * p + 0.5 * powr<2>(p)) *
powr<2>(Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))))) +
0.0833333333333333 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) *
powr<-1>(RB(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
powr<2>(l1) * ZA(l1)) *
ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p))) *
ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) +
0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p))) *
(-RFdot(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) -
RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) *
(dtZq(k) - 50. * Zq(k) + 50. * Zq(1.02 * k))) *
powr<-2>(powr<2>(Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)))) +
powr<2>(RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) *
powr<2>(Zq(k)) +
2. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) *
RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) *
Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) +
(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) *
powr<2>(Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))))) *
powr<-1>(
powr<2>(Mq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))) +
powr<2>(RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
powr<2>(Zq(k)) +
2. * sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
Zq(k) * Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) +
(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
powr<2>(Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))))) *
(1.333333333333333 * (-5. + 2. * powr<2>(cosl1p1)) *
sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) *
Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) *
Mq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
(RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) +
sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) *
Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)))) -
0.6666666666666666 *
sqrt(powr<-1>(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
((-1. + 4. * powr<2>(cosl1p1)) * powr<2>(l1) - (cosl1p1 + 2. * cosl1p2) * l1 * p -
(-1. + 2. * powr<2>(cosl1p1) + 4. * cosl1p1 * cosl1p2 + 4. * powr<2>(cosl1p2)) *
powr<2>(p)) *
(-powr<2>(Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)))) +
powr<2>(RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) *
powr<2>(Zq(k)) +
2. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) *
RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) *
Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) +
(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) *
powr<2>(Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))))) *
(RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zq(k) +
sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))))) +
0.0833333333333333 *
sqrt(powr<-1>(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
powr<-1>(RB(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k), 0.16666666666666666667)) +
powr<2>(l1) * ZA(l1)) *
ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p))) *
ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) +
0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p))) *
(-dtZq(k) * RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) -
RFdot(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zq(k) +
RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
(50. * Zq(k) - 50. * Zq(1.02 * k))) *
powr<-1>(powr<2>(Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)))) +
powr<2>(RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) *
powr<2>(Zq(k)) +
2. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) *
RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) *
Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) +
(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) *
powr<2>(Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))))) *
powr<-2>(
powr<2>(Mq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))) +
powr<2>(RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
powr<2>(Zq(k)) +
2. * sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
Zq(k) * Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) +
(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
powr<2>(Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))))) *
(1.333333333333333 * (-5. + 2. * powr<2>(cosl1p1)) *
sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) *
Mq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
(RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zq(k) +
sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))) -
0.6666666666666666 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) *
((-1. + 4. * powr<2>(cosl1p1)) * powr<2>(l1) - (cosl1p1 + 2. * cosl1p2) * l1 * p -
(-1. + 2. * powr<2>(cosl1p1) + 4. * cosl1p1 * cosl1p2 + 4. * powr<2>(cosl1p2)) *
powr<2>(p)) *
(RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) +
sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) *
Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)))) *
(-powr<2>(Mq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))) +
powr<2>(RF(powr<2>(k),
powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) *
powr<2>(Zq(k)) +
2. * sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
Zq(k) * Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) +
(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) *
powr<2>(Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))))));
      // clang-format on
      const auto _interp1 = dtZA(pow(1. + powr<6>(k), 0.16666666666666666667));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp3 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp4 = ZA(pow(1. + powr<6>(k), 0.16666666666666666667));
      const auto _interp5 = ZA(1.02 * pow(1. + powr<6>(k), 0.16666666666666666667));
      const auto _interp6 = ZA(l1);
      const auto _interp7 = RB(powr<2>(k), fma(2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp8 = ZA(sqrt(fma(2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp9 = ZA3(0.816496580927726 * sqrt(fma(cosl1p1, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp11 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp12 = RF(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp13 = Zq(k);
      const auto _interp14 = Zq(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp15 = Mq(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp19 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp21 = RF(powr<2>(k), fma(-2., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp22 = Zq(sqrt(fma(-2., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp23 = Mq(sqrt(fma(-2., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp25 = ZAqbq1(0.816496580927726 * sqrt(fma(cosl1p1 + cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp27 = dtZq(k);
      const auto _interp28 = RF(powr<2>(k), fma(2., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp29 = RFdot(powr<2>(k), fma(2., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp30 = Zq(1.02 * k);
      const auto _interp31 = Mq(sqrt(fma(2., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp32 = Zq(sqrt(fma(2., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p))));
      _T _acc{};
      { // subkernel 1
        const auto _interp10 =
          ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (cosl1p1 - cosl1p2) * l1 * p + powr<2>(p)));
        const auto _interp16 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
        const auto _interp17 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _interp18 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
        const auto _interp20 =
          ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
        const auto _cse1 = powr<2>(l1);
        const auto _cse2 = powr<2>(p);
        const auto _cse3 = powr<2>(cosl1p1);
        const auto _cse4 = powr<3>(cosl1p1);
        const auto _cse5 = -2. * cosl1p2 * l1 * p;
        const auto _cse6 = _cse1 + _cse2 + _cse5;
        const auto _cse7 = sqrt(_cse6);
        const auto _cse8 = -2. * cosl1p2 * l1 * p;
        const auto _cse9 = _cse1 + _cse2 + _cse8;
        const auto _cse10 = 2. * cosl1p1 * l1 * p;
        const auto _cse11 = _cse1 + _cse10 + _cse2;
        const auto _cse12 = powr<6>(k);
        const auto _cse13 = 1. + _cse12;
        const auto _cse14 = powr<-1>(_cse13);
        const auto _cse15 = -2. * cosl1p1 * l1 * p;
        const auto _cse16 = 6. * _cse3;
        const auto _cse17 = -6. + _cse16;
        const auto _cse18 = powr<2>(_cse1);
        const auto _cse19 = _cse17 * _cse18;
        const auto _cse20 = 6. * cosl1p2;
        const auto _cse21 = -4. * _cse3 * cosl1p2;
        const auto _cse22 = powr<3>(l1);
        const auto _cse23 = -14. * _cse4;
        const auto _cse24 = 10. * cosl1p2;
        const auto _cse25 = -14. * _cse3 * cosl1p2;
        const auto _cse26 = powr<3>(p);
        const auto _cse27 = powr<2>(_cse2);
        const auto _cse28 = cosl1p1 + cosl1p2;
        const auto _cse29 = -2. * _cse28 * l1 * p;
        const auto _cse30 = _cse1 + _cse2 + _cse29;
        const auto _cse31 = sqrt(_cse30);
        const auto _cse32 = powr<2>(_interp13);
        const auto _cse33 = _cse1 + _cse15 + _cse2 + _cse8;
        const auto _cse34 = _cse1 + _cse15 + _cse2; // clang-format off
_acc += fma(0.5,
powr<-1>(_cse11) * _cse14 * _cse7 * powr<-1>(_cse9) * _interp10 * _interp11 *
(_interp12 * _interp13 + _cse7 * _interp14) *
powr<-1>(_cse32 * powr<2>(_interp12) + 2. * _cse7 * _interp12 * _interp13 * _interp14 +
_cse9 * powr<2>(_interp14) + powr<2>(_interp15)) *
((1. + 1. * _cse12) * _interp1 * _interp2 + (1. + _cse12) * _interp3 * _interp4 +
_cse12 * _interp2 * (-50. * _interp4 + 50. * _interp5)) *
powr<-2>(_interp2 * _interp4 + _cse1 * _interp6) *
powr<-1>(_interp4 * _interp7 + _cse11 * _interp8) * _interp9 *
(_cse19 + _cse27 * (-3. - 6. * cosl1p1 * cosl1p2) +
_cse1 * _cse2 * (-11. + 13. * _cse3 - 4. * _cse4 * cosl1p2 + 8. * cosl1p1 * cosl1p2) +
_cse26 * (_cse24 + _cse25 - 2. * cosl1p1) * l1 +
_cse22 * (_cse20 + _cse21 + 10. * _cse4 - 9. * cosl1p1) * p),
fma(0.5,
_cse14 * _cse31 * powr<-1>(_cse33) * powr<-1>(_cse34) * _interp18 * _interp19 * _interp20 *
(_interp13 * _interp21 + _cse31 * _interp22) *
powr<-1>(_cse32 * powr<2>(_interp21) + 2. * _cse31 * _interp13 * _interp21 * _interp22 +
_cse33 * powr<2>(_interp22) + powr<2>(_interp23)) *
powr<-1>(_cse34 * _interp17 + _interp16 * _interp4) *
((1. + 1. * _cse12) * _interp1 * _interp2 + (1. + _cse12) * _interp3 * _interp4 +
_cse12 * _interp2 * (-50. * _interp4 + 50. * _interp5)) *
powr<-2>(_interp2 * _interp4 + _cse1 * _interp6) *
(_cse19 +
_cse1 * _cse2 *
(-11. + 5. * _cse3 + 4. * powr<2>(_cse3) + 4. * _cse4 * cosl1p2 -
8. * cosl1p1 * cosl1p2) +
_cse27 * (-3. + _cse16 + 6. * cosl1p1 * cosl1p2) +
_cse26 * (_cse23 + _cse24 + _cse25 + 12. * cosl1p1) * l1 +
_cse22 * (_cse20 + _cse21 + _cse23 + 15. * cosl1p1) * p),
0.));
        // clang-format on
      }
      { // subkernel 2
        const auto _interp24 =
          ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
        const auto _interp26 =
          ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
        const auto _cse1 = powr<2>(l1);
        const auto _cse2 = 2. * cosl1p1 * l1 * p;
        const auto _cse3 = powr<2>(p);
        const auto _cse4 = powr<2>(cosl1p1);
        const auto _cse5 = powr<3>(cosl1p1);
        const auto _cse6 = 2.333333333333333 * _cse5;
        const auto _cse7 = 1. * _cse4;
        const auto _cse8 = _cse1 + _cse2 + _cse3;
        const auto _cse9 = cosl1p1 + cosl1p2;
        const auto _cse10 = 2. * _cse9 * l1 * p;
        const auto _cse11 = _cse1 + _cse10 + _cse3;
        const auto _cse12 = sqrt(_cse11);
        const auto _cse13 = powr<2>(_interp31);
        const auto _cse14 = 0.5 * _cse1;
        const auto _cse15 = 1. * cosl1p1 * l1 * p;
        const auto _cse16 = 1. * cosl1p2 * l1 * p;
        const auto _cse17 = 0.5 * _cse3;
        const auto _cse18 = _cse14 + _cse15 + _cse16 + _cse17;
        const auto _cse19 = powr<2>(_interp32);
        const auto _cse20 = _cse18 * _cse19;
        const auto _cse21 = 1. * _cse12 * _interp13 * _interp28 * _interp32;
        const auto _cse22 = powr<2>(_interp28);
        const auto _cse23 = powr<2>(_interp13);
        const auto _cse24 = 0.5 * _cse22 * _cse23;
        const auto _cse25 = -2. * cosl1p2 * l1 * p;
        const auto _cse26 = _cse1 + _cse25 + _cse3;
        const auto _cse27 = -2. * _cse9 * l1 * p;
        const auto _cse28 = _cse1 + _cse27 + _cse3;
        const auto _cse29 = _cse1 * _interp6;
        const auto _cse30 = _interp2 * _interp4;
        const auto _cse31 = _cse29 + _cse30;
        const auto _cse32 = powr<6>(k); // clang-format off
_acc += fma(1.5,
_cse12 * (-0.5 * _cse13 + _cse20 + _cse21 + _cse24) *
powr<-2>(0.5 * _cse13 + _cse20 + _cse21 + _cse24) * powr<-1>(_cse31) * powr<-1>(_cse8) *
_interp25 * _interp26 *
(_interp27 * _interp28 + _interp13 * _interp29 +
_interp28 * (-50. * _interp13 + 50. * _interp30)) *
powr<-1>(_interp4 * _interp7 + _cse8 * _interp8) * _interp9 *
powr<-1>(_cse1 + _cse2 + _cse3 + 2. * cosl1p2 * l1 * p) *
(powr<2>(_cse1) * (-1. + _cse7) +
_cse1 * _cse3 *
(-1.833333333333333 + 0.833333333333333 * _cse4 + 0.6666666666666666 * powr<2>(_cse4) +
0.6666666666666666 * _cse5 * cosl1p2 - 1.333333333333333 * cosl1p1 * cosl1p2) +
powr<2>(_cse3) * (-0.5 + _cse7 + 1. * cosl1p1 * cosl1p2) +
(_cse6 - 2.5 * cosl1p1 - cosl1p2 + 0.6666666666666666 * _cse4 * cosl1p2) * powr<3>(l1) *
p +
(_cse6 - 2. * cosl1p1 - 1.666666666666667 * cosl1p2 +
2.333333333333333 * _cse4 * cosl1p2) *
l1 * powr<3>(p)),
fma(0.05555555555555555,
powr<-2>(_cse31) * _interp11 *
powr<-1>(_cse23 * powr<2>(_interp12) +
2. * sqrt(_cse26) * _interp12 * _interp13 * _interp14 +
_cse26 * powr<2>(_interp14) + powr<2>(_interp15)) *
_interp19 *
powr<-1>(_cse23 * powr<2>(_interp21) +
2. * sqrt(_cse28) * _interp13 * _interp21 * _interp22 +
_cse28 * powr<2>(_interp22) + powr<2>(_interp23)) *
_interp24 *
(_interp3 * _interp4 + _interp2 * (_interp1 + 50. * _cse32 * powr<-1>(1. + _cse32) *
(-_interp4 + _interp5))) *
((5. - 2. * _cse4) * _interp15 * _interp23 +
(-sqrt(powr<-1>(_cse26)) * _interp12 * _interp13 - _interp14) *
(-sqrt(powr<-1>(_cse28)) * _interp13 * _interp21 - _interp22) *
(_cse1 * (-1. + 4. * _cse4) -
_cse3 * (-1. + 2. * _cse4 + 4. * cosl1p1 * cosl1p2 + 4. * powr<2>(cosl1p2)) +
(cosl1p1 + 2. * cosl1p2) * l1 * p)),
0.));
        // clang-format on
      }
      { // subkernel 3
        const auto _interp33 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
        const auto _interp34 =
          ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
        const auto _interp35 = RFdot(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
        const auto _interp36 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
        const auto _interp37 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _interp38 = Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _cse1 = powr<2>(l1);
        const auto _cse2 = cosl1p1 + cosl1p2;
        const auto _cse3 = 2. * _cse2 * l1 * p;
        const auto _cse4 = powr<2>(p);
        const auto _cse5 = _cse1 + _cse3 + _cse4;
        const auto _cse6 = 2. * cosl1p2 * l1 * p;
        const auto _cse7 = _cse1 + _cse4 + _cse6;
        const auto _cse8 = powr<2>(_interp13);
        const auto _cse9 = sqrt(_cse5);
        const auto _cse10 = powr<2>(cosl1p1);
        const auto _cse11 = 2. * _cse10;
        const auto _cse12 = powr<2>(_interp31);
        const auto _cse13 = powr<2>(_interp32);
        const auto _cse14 = _cse13 * _cse5;
        const auto _cse15 = 2. * _cse9 * _interp13 * _interp28 * _interp32;
        const auto _cse16 = powr<2>(_interp28);
        const auto _cse17 = _cse16 * _cse8;
        const auto _cse18 = sqrt(_cse7);
        const auto _cse19 = sqrt(powr<-1>(_cse7));
        const auto _cse20 = _cse12 + _cse14 + _cse15 + _cse17;
        const auto _cse21 = powr<2>(_interp37);
        const auto _cse22 = powr<2>(_interp38);
        const auto _cse23 = _cse22 * _cse7;
        const auto _cse24 = 2. * _cse18 * _interp13 * _interp36 * _interp38;
        const auto _cse25 = powr<2>(_interp36);
        const auto _cse26 = _cse25 * _cse8;
        const auto _cse27 = _cse21 + _cse23 + _cse24 + _cse26;
        const auto _cse28 = -5. + _cse11;
        const auto _cse29 = _cse18 * _interp38;
        const auto _cse30 = _interp13 * _interp36;
        const auto _cse31 = _cse29 + _cse30;
        const auto _cse32 = sqrt(powr<-1>(_cse5));
        const auto _cse33 = 4. * _cse10; // clang-format off
_acc += fma(0.0833333333333333,
powr<-2>(_cse20) * powr<-1>(_cse27) * _cse32 * _interp25 *
(-_interp27 * _interp28 - _interp13 * _interp29 +
_interp28 * (50. * _interp13 - 50. * _interp30)) *
_interp33 * _interp34 * powr<-1>(_interp2 * _interp4 + _cse1 * _interp6) *
(1.333333333333333 * _cse28 * _cse9 * _interp31 *
(_interp13 * _interp28 + _cse9 * _interp32) * _interp37 -
0.6666666666666666 * (-_cse12 + _cse14 + _cse15 + _cse17) * _cse19 * _cse31 *
(_cse1 * (-1. + _cse33) -
_cse4 * (-1. + _cse11 + 4. * cosl1p1 * cosl1p2 + 4. * powr<2>(cosl1p2)) -
(cosl1p1 + 2. * cosl1p2) * l1 * p)),
fma(0.0833333333333333,
_cse19 *powr<-1>(_cse20) * powr<-2>(_cse27) * _interp25 * _interp33 * _interp34 *
(-_interp13 * _interp35 -
(-50. * _interp13 + _interp27 + 50. * _interp30) * _interp36) *
powr<-1>(_interp2 * _interp4 + _cse1 * _interp6) *
(1.333333333333333 * _cse18 * _cse28 * _cse31 * _interp31 * _interp37 -
0.6666666666666666 * (-_cse21 + _cse23 + _cse24 + _cse26) * _cse32 *
(_interp13 * _interp28 + _cse9 * _interp32) *
(_cse1 * (-1. + _cse33) -
_cse4 * (-1. + _cse11 + 4. * cosl1p1 * cosl1p2 + 4. * powr<2>(cosl1p2)) -
(cosl1p1 + 2. * cosl1p2) * l1 * p)),
0.));
        // clang-format on
      }
      return _acc;
    }

    // clang-format off
static KOKKOS_FORCEINLINE_FUNCTION auto constant(const double& p, const double& k, const double& Nf, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
{
using namespace DiFfRG;
      // clang-format on
      using namespace DiFfRG::compute;
      return 0.;
    }

  private:
    static KOKKOS_FORCEINLINE_FUNCTION auto RB(const auto &k2, const auto &p2) { return Regulator::RB(k2, p2); }

    static KOKKOS_FORCEINLINE_FUNCTION auto RF(const auto &k2, const auto &p2) { return Regulator::RF(k2, p2); }

    static KOKKOS_FORCEINLINE_FUNCTION auto RBdot(const auto &k2, const auto &p2) { return Regulator::RBdot(k2, p2); }

    static KOKKOS_FORCEINLINE_FUNCTION auto RFdot(const auto &k2, const auto &p2) { return Regulator::RFdot(k2, p2); }

    static KOKKOS_FORCEINLINE_FUNCTION auto dq2RB(const auto &k2, const auto &p2) { return Regulator::dq2RB(k2, p2); }

    static KOKKOS_FORCEINLINE_FUNCTION auto dq2RF(const auto &k2, const auto &p2) { return Regulator::dq2RF(k2, p2); }
  };
} // namespace DiFfRG
using DiFfRG::ZAqbq1_kernel;