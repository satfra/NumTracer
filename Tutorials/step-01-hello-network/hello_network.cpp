// step-01 — hello, tensor network.
//
// The smallest possible use of the engine, with no physics vocabulary yet: contract two index
// networks over 4-vectors and read off the number each collapses to. Everything here is pure C++
// against the NumTracer library — no Mathematica, no external tools.
//
// The one rule the engine runs on: two factors are summed over an index exactly when they carry
// the same index label (Einstein convention). A network with no free labels left is a scalar.
//
// We do two contractions and check each against its closed form:
//   (1) a . b        = a_mu delta_{mu nu} b_nu       (a metric tying two vectors)
//   (2) a . P(k) . a = a^2 - (a.k)^2 / k^2           (a transverse projector on k)
#include <numtracer.hpp> // the whole NumTracer API

#include <cmath>
#include <cstdio>

namespace nt = numtracer;

// The three 4-vectors, as plain numbers. We use them both to build the network and to compute the
// closed-form check, so the two can never silently drift apart.
static constexpr double a[4] = {1.0, 0.5, -0.3, 0.2};
static constexpr double b[4] = {0.8, -0.4, 1.2, 0.1};
static constexpr double k[4] = {0.5, 0.7, 0.0, 0.0};

int main() {
  // @snip begin: frame
  // A Frame says what the momenta are. Here every component is a plain number, so the frame has no
  // symbols at all, and every contraction will collapse to a constant.
  nt::Frame F;
  nt::Momentum va = F.momentum(a[0], a[1], a[2], a[3]);
  nt::Momentum vb = F.momentum(b[0], b[1], b[2], b[3]);
  nt::Momentum vk = F.momentum(k[0], k[1], k[2], k[3]);

  // Two fresh, distinct Lorentz index labels. Reusing a label is how you say "sum these together".
  auto [mu, nu] = F.indices<2>();
  // @snip end: frame

  // @snip begin: dot
  // (1) a . b: the vector a on index mu, the metric delta_{mu nu}, the vector b on index nu.
  // `*` multiplies tensors; the shared labels mu and nu are summed when the network is contracted.
  nt::Poly dot = F.contract(nt::vec(mu, va) * nt::metric(mu, nu) * nt::vec(nu, vb));
  const double ab = F.eval(dot, F.at()).re; // no symbols, so the point is empty
  // @snip end: dot

  double abClosed = 0;
  for (int c = 0; c < 4; ++c) abClosed += a[c] * b[c];

  // @snip begin: proj
  // (2) a . P(k) . a with the transverse projector P(k)_{mu nu} = delta_{mu nu} - k_mu k_nu / k^2.
  // The frame knows k, so it supplies the 1/k^2 itself.
  nt::Poly proj = F.contract(nt::vec(mu, va) * nt::projT(mu, nu, vk) * nt::vec(nu, va));
  const double apa = F.eval(proj, F.at()).re;
  // @snip end: proj

  double k2 = 0, a2 = 0, ak = 0;
  for (int c = 0; c < 4; ++c) k2 += k[c] * k[c], a2 += a[c] * a[c], ak += a[c] * k[c];
  const double apaClosed = a2 - ak * ak / k2;

  std::printf("a . b        = %g   (expect %g)\n", ab, abClosed);
  std::printf("a . P(k) . a = %g   (expect a^2 - (a.k)^2/k^2 = %g)\n", apa, apaClosed);

  const bool ok = std::fabs(ab - abClosed) < 1e-12 && std::fabs(apa - apaClosed) < 1e-12;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
