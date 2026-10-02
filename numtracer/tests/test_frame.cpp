/// @file test_frame.cpp
/// @brief The user-facing layer: @ref numtracer::Frame (named symbols and momenta, automatic projector
///        atoms, checked evaluation), typed labels, @ref numtracer::SUN, and @ref numtracer::interpret.
///
/// Each check is a closed form a physicist can verify by hand, evaluated at several random points.
/// The last group lowers a contraction to a straight-line program and runs it, so the lowering is
/// graded against the polynomial it came from.
#include "numtracer/numtracer.hpp"

#include <cmath>
#include <cstdio>
#include <functional>
#include <random>
#include <stdexcept>
#include <string>

namespace nt = numtracer;
using nt::Cx;

namespace
{
  int fails = 0;
  void check(bool ok, const std::string &what)
  {
    std::printf("  %-62s %s\n", what.c_str(), ok ? "ok" : "FAIL");
    if (!ok) ++fails;
  }
  bool close(double a, double b) { return std::fabs(a - b) <= 1e-11 * (1.0 + std::fabs(b)); }

  /// Does @p f throw (anything)? Returns the message for a substring check.
  std::string throws(const std::function<void()> &f)
  {
    try {
      f();
    } catch (const std::exception &e) {
      return e.what();
    }
    return "";
  }
} // namespace

int main()
{
  std::mt19937 rng(20261002);
  std::uniform_real_distribution<double> U(0.2, 1.7), C(-0.95, 0.95);

  std::printf("== Frame: traces against closed forms ==\n");
  {
    // tr(p/ q/) = 4 p.q and tr(gamma^mu p/ gamma_mu q/) = -2 tr(p/ q/), with fully symbolic momenta.
    nt::Frame F;
    nt::Symbol s[8] = {F.symbol("p0"), F.symbol("p1"), F.symbol("p2"), F.symbol("p3"),
                       F.symbol("q0"), F.symbol("q1"), F.symbol("q2"), F.symbol("q3")};
    auto p = F.momentum(s[0], s[1], s[2], s[3]);
    auto q = F.momentum(s[4], s[5], s[6], s[7]);
    auto [mu, nu] = F.indices<2>();
    const nt::Poly trPQ = F.trace({nt::slash(p), nt::slash(q)});
    const nt::Poly trG = F.trace({nt::gamma(mu), nt::slash(p), nt::gamma(nu), nt::slash(q)}, nt::metric(mu, nu));
    bool ok1 = true, ok2 = true;
    for (int t = 0; t < 5; ++t) {
      double x[8];
      for (double &v : x) v = U(rng);
      const nt::Point pt = F.at(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]);
      const double pq = x[0] * x[4] + x[1] * x[5] + x[2] * x[6] + x[3] * x[7];
      ok1 = ok1 && close(F.eval(trPQ, pt).re, 4 * pq);
      ok2 = ok2 && close(F.eval(trG, pt).re, -8 * pq);
    }
    check(ok1, "tr(p/ q/) = 4 p.q");
    check(ok2, "tr(gamma^mu p/ gamma_mu q/) = -8 p.q");
  }
  {
    // p.P_T(l).p = p^2 (1 - cos^2) in a one-angle frame; the projector's 1/l^2 is assigned automatically.
    nt::Frame F;
    auto P = F.symbol("p"), L = F.symbol("l");
    auto [Cs, Sn] = F.angle("theta");
    auto p = F.momentum(P, 0, 0, 0);
    auto l = F.momentum(L * Cs, L * Sn, 0, 0);
    auto [mu, nu] = F.indices<2>();
    const nt::Poly v = F.contract(nt::vec(mu, p) * nt::projT(mu, nu, l) * nt::vec(nu, p));
    bool ok = true;
    for (int t = 0; t < 5; ++t) {
      const double pv = U(rng), lv = U(rng), c = C(rng);
      ok = ok && close(F.eval(v, F.at({{P, pv}, {L, lv}, {Cs, c}})).re, pv * pv * (1 - c * c));
    }
    check(ok, "p.P_T(l).p = p^2 (1 - cos^2)");
    check(F.denominators().size() == 1, "one projector momentum -> one registered 1/k^2");
    // the quark self-energy numerator tr[p/ g^mu q/ g^nu] P_T(l)_{mu nu}, q = l - p (tutorial step 5)
    auto q = l - p;
    const nt::Poly T = F.trace({nt::slash(p), nt::gamma(mu), nt::slash(q), nt::gamma(nu)}, nt::projT(mu, nu, l));
    bool okT = true;
    for (int t = 0; t < 5; ++t) {
      const double pv = U(rng), lv = U(rng), c = C(rng);
      const double closed = 4.0 * pv * (-3.0 * c * lv + pv + 2.0 * c * c * pv);
      okT = okT && close(F.eval(T, F.at(pv, lv, c)).re, closed);
    }
    check(okT, "tr[p/ g^mu q/ g^nu] P_T(l) = 4p(-3 c l + p + 2 c^2 p)");
    check(F.denominators().size() == 1, "projT on l again reuses its atom");

    // ---- lowering: the emitted program, run in-process, equals the polynomial ----
    nt::GlobalEnv g;
    const nt::GenProg prog = nt::to_genprog(T, g);
    bool okL = true;
    for (int t = 0; t < 5; ++t) {
      const nt::Point pt = F.at(U(rng), U(rng), C(rng));
      const std::vector<double> f = F.fill_values(g, pt);
      okL = okL && close(nt::interpret(prog, f.data()).re, F.eval(T, pt).re);
    }
    check(okL, "interpret(to_genprog(T)) = eval(T)");
  }

  std::printf("== SU(N): typed labels ==\n");
  {
    nt::SUN su3(3);
    auto [a, b, c] = su3.adjoint<3>();
    auto [i, j] = su3.fundamental<2>();
    check(nt::approx(su3.value(su3.T(a, i, j) * su3.T(a, j, i)), Cx{4, 0}), "tr(T^a T^a) = 4");
    check(nt::approx(su3.value(su3.f(a, b, c) * su3.f(a, b, c)), Cx{24, 0}), "f^abc f^abc = 24");
    // two SU(3) groups in one net: their labels never contract with each other
    nt::SUN fl3(3);
    auto [x] = fl3.adjoint<1>();
    const Cx both = nt::sun_value(su3.T(a, i, j) * su3.T(a, j, i) * fl3.delta(x, x));
    check(nt::approx(both, Cx{32, 0}), "SU(3) (x) SU(3): 4 * delta^xx = 32");
    check(!throws([&] { (void)su3.T(x, i, j); }).empty(), "label of another group -> throws");
    check(!throws([&] { (void)su3.value({fl3.delta(x, x)}); }).empty(), "SUN::value on a foreign factor -> throws");
  }

  std::printf("== guards ==\n");
  {
    nt::Frame F;
    auto P = F.symbol("p");
    auto p = F.momentum(P, 0, 0, 0);
    (void)p;
    check(throws([&] { (void)F.symbol("late"); }).find("declare every symbol") != std::string::npos,
          "symbol after the frame froze -> throws");
    check(!throws([&] { (void)F.at(1.0, 2.0); }).empty(), "at() with too many values -> throws");
    auto [mu] = F.indices<1>();
    check(throws([&] { (void)F.trace({nt::gamma(mu), nt::slash(p)}); }).find("is OPEN") != std::string::npos,
          "open Lorentz index -> throws");
    const nt::Momentum ghost{{{1.0, 7}}};
    check(throws([&] { (void)F.trace({nt::slash(ghost), nt::slash(p)}); }).find("not in the frame") != std::string::npos,
          "momentum not in the frame -> throws");
    nt::Frame G;
    auto Q = G.symbol("q");
    const nt::Poly other = Q;
    check(!throws([&] { (void)F.eval(other, G.at(1.0)); }).empty(), "point of another frame -> throws");
  }

  {
    nt::Frame F, G;
    auto X = G.symbol("x");
    auto Y = F.symbol("y");
    (void)Y;
    check(!throws([&] { (void)F.momentum(X, 0, 0, 0); }).empty(), "momentum from another frame's symbol -> throws");
    nt::Frame H;
    check(!throws([&] { (void)H.symbol("not a name"); }).empty(), "symbol name that is not a C++ identifier -> throws");
    check(throws([&] { (void)H.eval(nt::Poly{}, H.at()); }).empty(), "a default (zero) Poly evaluates to 0");
  }
  {
    nt::SUN su3(3);
    auto [a, b] = su3.adjoint<2>();
    check(!throws([&] { (void)su3.diag(a, b, {0, 1}); }).empty(), "diag with the wrong component count -> throws");
    check(throws([&] { (void)su3.value(su3.diag(a, b, std::vector<int>(8, 0)) * su3.delta(b, a)); })
                  .find("sun_value_dressed") != std::string::npos,
          "SUN::value on a diag factor -> points to sun_value_dressed");
    check(!throws([] { nt::SUN bad(0); }).empty(), "SU(0) -> throws");
  }

  std::printf(fails ? "FAILED (%d)\n" : "ALL TESTS PASSED\n", fails);
  return fails ? 1 : 0;
}
