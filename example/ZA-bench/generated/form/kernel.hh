#pragma once

#include "DiFfRG/physics/interpolation.hh"
#include "DiFfRG/physics/physics.hh"

namespace DiFfRG
{
  template <typename _Regulator> class ZA_kernel
  {
  public:
    using Regulator = _Regulator;

    static KOKKOS_FORCEINLINE_FUNCTION auto
    kernel(const double &l1, const double &cos1, const double &p, const double &k, const double &Nf,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZA3,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZAcbc,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZA4,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZAqbq1,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZAqbq4,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZAqbq7,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZAAqbq1,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZAAqbq2,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &dtZc,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &Zc,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &dtZA,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZA,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &dtZq,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &Zq,
           const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &Mq)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      const auto _repl1 = Zq(k);
      const auto _repl2 = RF(powr<2>(k), powr<2>(l1));
      const auto _repl3 = Zq(sqrt(powr<2>(l1) + (-2.) * ((cos1) * ((l1) * (p))) + powr<2>(p)));
      const auto _repl4 = RF(powr<2>(k), powr<2>(l1) + (-2.) * ((cos1) * ((l1) * (p))) + powr<2>(p));
      const auto _repl5 = Zq(l1);
      const auto _repl6 =
          ZAqbq7((0.816496580927726) * (sqrt(powr<2>(l1) + (-1.) * ((cos1) * ((l1) * (p))) + powr<2>(p))));
      const auto _repl7 =
          ZAqbq4((0.816496580927726) * (sqrt(powr<2>(l1) + (-1.) * ((cos1) * ((l1) * (p))) + powr<2>(p))));
      const auto _repl8 = Mq(l1);
      const auto _repl9 =
          ZAqbq1((0.816496580927726) * (sqrt(powr<2>(l1) + (-1.) * ((cos1) * ((l1) * (p))) + powr<2>(p))));
      const auto _repl10 = Mq(sqrt(powr<2>(l1) + (-2.) * ((cos1) * ((l1) * (p))) + powr<2>(p)));
      const auto _repl11 = RB(powr<2>(k), powr<2>(l1));
      const auto _repl12 = ZA(pow(1. + powr<6>(k), 0.16666666666666666667));
      const auto _repl13 = Zc(k);
      const auto _repl14 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _repl15 = RB(powr<2>(k), powr<2>(l1) + (-2.) * ((cos1) * ((l1) * (p))) + powr<2>(p));
      const auto _repl16 = Zq((1.02) * (k));
      const auto _repl17 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _repl18 = dtZq(k);
      const auto _repl19 = ZA(l1);
      const auto _repl20 = ZA((1.02) * (pow(1. + powr<6>(k), 0.16666666666666666667)));
      const auto _repl21 = dtZA(pow(1. + powr<6>(k), 0.16666666666666666667));
      const auto _repl22 = powr<2>(l1);
      const auto _repl23 = powr<2>(p);
      const auto _repl24 = powr<2>(k);
      const auto _repl25 = powr<6>(k);
      const auto _repl26 = powr<2>((_repl1) * (_repl2) + (_repl5) * (l1));
      const auto _repl27 = powr<2>(_repl8);
      const auto _repl28 = powr<2>(_repl6);
      const auto _repl29 = powr<2>(cos1);
      const auto _repl30 = powr<2>(_repl7);
      const auto _repl31 = powr<-2>(p);
      const auto _repl32 = powr<2>(_repl9);
      const auto _repl33 = powr<3>(l1);
      const auto _repl34 = powr<-2>((_repl11) * (_repl12) + (_repl19) * (_repl22));
      const auto _repl35 = powr<3>(p);
      const auto _repl36 = powr<2>(_repl1);
      const auto _repl37 = powr<4>(cos1);
      return (-1.333333333333333) *
                 ((-1.) * ((_repl1) * (_repl17)) + ((50.) * (_repl1) + (-50.) * (_repl16)) * (_repl2) +
                  (-1.) * ((_repl18) * (_repl2))) *
                 ((powr<-2>(_repl26 + _repl27)) *
                  ((-3.) *
                       ((_repl1) * (_repl4) + (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                       ((_repl23) * ((_repl26) * ((_repl28) * (_repl33)))) +
                   (2.) * ((-1.) * (_repl26) + _repl27) *
                       ((_repl23) *
                        ((_repl1) * (_repl4) + (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                        ((_repl28) * ((_repl33) * (_repl37)))) +
                   (_repl22) * ((-1.) * (_repl26) + _repl27) *
                       ((_repl28) *
                        ((_repl1) * (_repl4) + (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                        ((_repl35) * (powr<3>(cos1)))) +
                   (_repl26) *
                       ((_repl1) * (_repl4) + (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                       ((_repl32) * (l1)) +
                   (cos1) *
                       ((_repl22) *
                            ((_repl1) * (_repl4) +
                             (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                            ((_repl23) * ((_repl26) * (_repl28))) +
                        (-1.) * ((_repl22) * ((_repl27) * (_repl28)) + (-3.) * ((_repl26) * (_repl30))) *
                            (((_repl1) * (_repl4) +
                              (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                             (_repl23)) +
                        (12.) * ((_repl1) * (_repl2) + (_repl5) * (l1)) *
                            ((_repl7) *
                             ((_repl1) * (_repl4) +
                              (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                             ((_repl8) * ((_repl9) * (l1)))) +
                        (3.) *
                            ((2.) * ((_repl10) * ((_repl7) * ((_repl9) * (sqrt(_repl22 + _repl23 +
                                                                               (-2.) * ((cos1) * ((l1) * (p)))))))) +
                             ((-1.) * ((_repl23) * (_repl30)) + _repl32) *
                                 ((_repl1) * (_repl4) +
                                  (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p))))))) *
                            (_repl27) +
                        (-3.) *
                            ((2.) * ((_repl10) *
                                     ((_repl7) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p))))))) +
                             ((_repl1) * (_repl4) +
                              (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                                 (_repl9)) *
                            ((_repl26) * (_repl9))) *
                       (p) +
                   (_repl23) *
                       ((3.) *
                            ((_repl1) * (_repl4) +
                             (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                            ((_repl22) * ((_repl27) * (_repl28))) +
                        (2.) * ((_repl1) * (_repl2) + (_repl5) * (l1)) *
                            ((_repl6) *
                             ((_repl10) * ((_repl6) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) +
                              (-4.) *
                                  ((_repl1) * (_repl4) +
                                   (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                                  (_repl7)) *
                             ((_repl8) * (l1))) +
                        (-1.) *
                            ((4.) * ((_repl10) * ((_repl6) * ((_repl7) * (sqrt(_repl22 + _repl23 +
                                                                               (-2.) * ((cos1) * ((l1) * (p)))))))) +
                             (-1.) * (_repl30 + (-4.) * ((_repl6) * (_repl9))) *
                                 ((_repl1) * (_repl4) +
                                  (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p))))))) *
                            (_repl26)) *
                       (l1) +
                   (((6.) * ((_repl1) * (_repl2) + (_repl5) * (l1)) *
                         ((_repl10) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) +
                     (-1.) *
                         ((_repl1) * (_repl4) +
                          (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                         ((_repl8) * (l1))) *
                        (_repl32) +
                    ((-6.) * ((_repl1) * (_repl2) + (_repl5) * (l1)) *
                         (((_repl10) * ((_repl7) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) +
                           (2.) *
                               ((_repl1) * (_repl4) +
                                (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                               (_repl9)) *
                          (_repl7)) +
                     (_repl8) *
                         ((4.) * ((_repl10) * ((_repl6) * ((_repl7) * (sqrt(_repl22 + _repl23 +
                                                                            (-2.) * ((cos1) * ((l1) * (p)))))))) +
                          (-1.) * (_repl30 + (-4.) * ((_repl6) * (_repl9))) *
                              ((_repl1) * (_repl4) +
                               (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p))))))) *
                         (l1)) *
                        (_repl23)) *
                       (_repl8) +
                   (_repl29) *
                       ((5.) *
                            ((_repl1) * (_repl4) +
                             (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                            ((_repl22) * ((_repl23) * ((_repl26) * (_repl28)))) +
                        (2.) *
                            ((_repl1) * (_repl4) +
                             (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                            ((_repl26) * (_repl32)) +
                        (-2.) *
                            ((2.) * ((_repl10) *
                                     ((_repl23) * ((_repl6) * ((_repl7) * (sqrt(_repl22 + _repl23 +
                                                                                (-2.) * ((cos1) * ((l1) * (p))))))))) +
                             ((-2.) * ((_repl23) * (_repl30)) + _repl32 + (2.) * ((_repl23) * ((_repl6) * (_repl9)))) *
                                 ((_repl1) * (_repl4) +
                                  (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p))))))) *
                            (_repl27) +
                        (-1.) *
                            ((5.) *
                                 ((_repl1) * (_repl4) +
                                  (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                                 ((_repl22) * ((_repl27) * (_repl28))) +
                             (2.) * ((_repl1) * (_repl2) + (_repl5) * (l1)) *
                                 ((_repl6) *
                                  ((_repl10) *
                                       ((_repl6) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) +
                                   (-4.) *
                                       ((_repl1) * (_repl4) +
                                        (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                                       (_repl7)) *
                                  ((_repl8) * (l1))) +
                             (-4.) *
                                 ((_repl10) * ((_repl6) * ((_repl7) * (sqrt(_repl22 + _repl23 +
                                                                            (-2.) * ((cos1) * ((l1) * (p))))))) +
                                  (-1.) * (_repl30 + (-1.) * ((_repl6) * (_repl9))) *
                                      ((_repl1) * (_repl4) +
                                       (_repl3) * (sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p))))))) *
                                 (_repl26)) *
                            (_repl23)) *
                       (l1)) *
                  ((_repl31) *
                   ((Nf) *
                    ((sqrt(powr<-1>(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))) *
                     (powr<-1>(powr<2>(_repl10) + (_repl36) * (powr<2>(_repl4)) +
                               (2.) * ((_repl1) * ((_repl3) * ((_repl4) * (sqrt(_repl22 + _repl23 +
                                                                                (-2.) * ((cos1) * ((l1) * (p)))))))) +
                               (_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))) * (powr<2>(_repl3)))))))) +
             (-12.) *
                 ((_repl11) * ((-50.) * (_repl12) + (50.) * (_repl20)) * (_repl25) +
                  (_repl12) * (1. + (1.) * (_repl25)) * (_repl14) + (_repl11) * (1. + (1.) * (_repl25)) * (_repl21)) *
                 ((powr<-1>(1. + _repl25)) *
                  ((_repl22) *
                       (-2.666666666666667 + (2.333333333333333) * (_repl29) + (0.3333333333333333) * (_repl37)) *
                       (_repl23) +
                   (_repl35) * (2. + (-2.) * (_repl29)) * ((cos1) * (l1)) + (-1. + (1.) * (_repl29)) * (powr<4>(l1)) +
                   (_repl33) * (2. + (-2.) * (_repl29)) * ((cos1) * (p)) + (-1. + (1.) * (_repl29)) * (powr<4>(p))) *
                  ((_repl31) *
                   ((_repl34) * ((powr<-1>(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p))))) *
                                 ((powr<-1>((_repl12) * (_repl15) +
                                            (_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))) *
                                                (ZA(sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))))))) *
                                  (powr<2>(ZA3((0.816496580927726) *
                                               (sqrt(_repl22 + _repl23 + (-1.) * ((cos1) * ((l1) * (p))))))))))))) +
             (-1.) *
                 ((_repl12) * (_repl14) +
                  (_repl21 + (50.) * ((-1.) * (_repl12) + _repl20) * ((_repl25) * (powr<-1>(1. + _repl25)))) *
                      (_repl11)) *
                 ((_repl31) * (7. + (-1.) * (_repl29)) *
                  ((_repl34) * (ZA4((0.7071067811865475) * (sqrt(_repl22 + _repl23)))))) +
             (-8.) *
                 ((_repl1) * (_repl17) + ((-50.) * (_repl1) + (50.) * (_repl16)) * (_repl2) + (_repl18) * (_repl2)) *
                 ((_repl31) * ((_repl1) * (_repl2) + (_repl5) * (l1)) *
                  ((_repl8) * ((powr<-2>(_repl27 + (powr<2>(_repl2)) * (_repl36) + (_repl22) * (powr<2>(_repl5)) +
                                         (2.) * ((_repl1) * ((_repl2) * ((_repl5) * (l1)))))) *
                               ((Nf) * (ZAAqbq1((0.7071067811865475) * (sqrt(_repl22 + _repl23)))))))) +
             (2.) * (-1. + _repl29) *
                 ((_repl22) *
                  ((_repl13) * (_repl14) + (_repl11) * (dtZc(k)) +
                   ((-50.) * (_repl13) + (50.) * (Zc((1.02) * (k)))) * (_repl11)) *
                  ((_repl31) * ((powr<2>(ZAcbc((0.816496580927726) *
                                               (sqrt(_repl22 + _repl23 + (-1.) * ((cos1) * ((l1) * (p)))))))) *
                                ((powr<-2>((_repl11) * (_repl13) + (_repl22) * (Zc(l1)))) *
                                 (powr<-1>((_repl13) * (_repl15) +
                                           (_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p)))) *
                                               (Zc(sqrt(_repl22 + _repl23 + (-2.) * ((cos1) * ((l1) * (p))))))))))));
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto
    constant(const double &p, const double &k, const double &Nf,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZA3,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZAcbc,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZA4,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZAqbq1,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZAqbq4,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZAqbq7,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZAAqbq1,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZAAqbq2,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &dtZc,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &Zc,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &dtZA,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &ZA,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &dtZq,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &Zq,
             const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory> &Mq)
    {
      using namespace DiFfRG;
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
using DiFfRG::ZA_kernel;