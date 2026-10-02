// step-05 — Hand-coding a full diagram, end to end.
//
// The capstone: the Dirac-trace numerator of the quark self-energy,
//   T = tr[ p/ gamma^mu q/ gamma^nu ] P_{mu nu}(l),   q = l - p,
// contracted, checked against its closed form, then LOWERED to straight-line C++ — and the lowered
// program is run too, so all three agree at one point:
//   (1) the contraction        F.trace(...)                    -> a polynomial T
//   (2) the closed form        4 p (-3 c l + p + 2 c^2 p)      (c = cos theta)
//   (3) the lowered program    to_genprog(T), run by interpret -> the arithmetic a kernel performs
#include <numtracer.hpp> // the whole NumTracer API

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

namespace nt = numtracer;

int main() {
  // @snip begin: frame
  // One-angle frame: p along axis 0, l at angle theta to it, q = l - p.
  nt::Frame F;
  auto P = F.symbol("p"), L = F.symbol("l");
  auto [C, S] = F.angle("cos");
  nt::Momentum p = F.momentum(P, 0, 0, 0);
  nt::Momentum l = F.momentum(L * C, L * S, 0, 0);
  nt::Momentum q = l - p;
  auto [mu, nu] = F.indices<2>();
  // @snip end: frame

  // @snip begin: tokens
  // (1) The closed chain p/ gamma^mu q/ gamma^nu, its free indices mu, nu tied by P_T(l).
  nt::Poly T = F.trace({nt::slash(p), nt::gamma(mu), nt::slash(q), nt::gamma(nu)}, nt::projT(mu, nu, l));
  // @snip end: tokens

  // @snip begin: closed
  // (2) The closed form at one point.
  const double pv = 1.3, lv = 0.86, c = 0.58;
  const nt::Point pt = F.at(pv, lv, c);
  const double numeric = F.eval(T, pt).re;
  const double closed = 4.0 * pv * (-3.0 * c * lv + pv + 2.0 * c * c * pv);
  // @snip end: closed

  // @snip begin: lower
  // (3) Lower T to a straight-line program: Horner factoring + common-subexpression elimination over
  // an array f[] of "env" values (the symbols and the 1/l^2 the polynomial needs). The GlobalEnv `g`
  // records which f[i] holds what. Print the program as C++, together with the fill() that computes
  // f[] from p, l, cos ...
  nt::GlobalEnv g;
  nt::GenProg prog = nt::to_genprog(T, g);
  F.emit_fill(std::cout, g);        // static inline void fill(double* f, double p, double l, double cos)
  nt::emit_cpp(std::cout, prog, "T"); // static inline double T(const double* f)

  // ... and RUN the same program in-process, on the f[] values at our point. This is the guarantee
  // that lowering preserved the value: a Horner or CSE bug would show up here and nowhere else.
  const std::vector<double> f = F.fill_values(g, pt);
  const double lowered = nt::interpret(prog, f.data()).re;
  // @snip end: lower

  std::printf("\nnumeric  T = %.12f\n", numeric);
  std::printf("closed   T = %.12f   (4 p (-3 c l + p + 2 c^2 p))\n", closed);
  std::printf("lowered  T = %.12f   (the emitted program, interpreted)\n", lowered);

  // The lowered program fuses multiply-adds, so it agrees to rounding, not bit for bit.
  const bool ok = std::fabs(numeric - closed) < 1e-10 && std::fabs(lowered - closed) < 1e-10;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
