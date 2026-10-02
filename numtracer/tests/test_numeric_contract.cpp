// Part B (task #24): correctness gate for the numeric (matrix-product) contraction backend.
//
// Validates the engine in et/numeric/* against independent numeric ground truth:
//   A) all-slash Dirac chains tr(p̸₁…p̸_K), K=2..10, vs the chiral γ-trace oracle chiral_gamma_trace<K>;
//   B) a free-leg chain closed by a metric, tr(γ^μ p̸ γ_μ q̸) = Σ_μ tr(γ^μ p̸ γ^μ q̸), vs a direct
//      numeric 4×4 contraction;
//   C) a transverse projector P(k) contracting two free legs — the monomial-k² case (unit-direction
//      loop → the 1/k² atom CANCELS, no surviving atom) and the non-monomial-k² case (shifted line →
//      the atom SURVIVES) — value-checked vs a direct numeric contraction.
//
// Symbolic momentum components are Poly variables; the engine result is evaluated at random points
// and compared to the numeric truth to ≤ 1e-12. Build via the test CMake (adds -I include).
#include "oracle/dense_trace.hpp"          // chiral_gamma_trace<K> (test-only γ-trace oracle)
#include "numtracer/dirac/dirac_data.hpp"  // kGamma, kGamma5, kC
#include "engine_test_util.hpp"
using namespace numtracer; // the typed builders (LorentzIndex, Momentum, …) the tests spell unqualified
#include "numtracer/numeric/frame.hpp"
#include "numtracer/numeric/numeric_contract.hpp"

#include <cmath>
#include <array>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

using numtracer::Cx;
using numtracer::dirac::kGamma;
namespace nm = numtracer::numeric;
namespace network = numtracer::network;

// ---- tiny numeric 4×4 complex matrix ground truth -------------------------------------------
struct NM {
  Cx a[4][4]{};
};
static NM nGamma(int mu)
{
  NM M;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      M.a[i][j] = kGamma[mu][i][j];
  return M;
}
static NM nGamma5()
{
  NM M;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      M.a[i][j] = numtracer::dirac::kGamma5[i][j];
  return M;
}
static NM nCmat()
{
  NM M;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      M.a[i][j] = numtracer::dirac::kC[i][j];
  return M;
}
static NM nSlash(const double p[4])
{
  NM M;
  for (int mu = 0; mu < 4; ++mu)
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        const Cx g = kGamma[mu][i][j];
        M.a[i][j] = M.a[i][j] + Cx{g.re * p[mu], g.im * p[mu]};
      }
  return M;
}
static NM nMul(const NM &A, const NM &B)
{
  NM C;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      Cx s{0, 0};
      for (int k = 0; k < 4; ++k)
        s = s + A.a[i][k] * B.a[k][j];
      C.a[i][j] = s;
    }
  return C;
}
static Cx nTrace(const NM &A)
{
  Cx s{0, 0};
  for (int i = 0; i < 4; ++i)
    s = s + A.a[i][i];
  return s;
}
// bare commutator [X,Y] = X·Y − Y·X of two explicit 4×4 matrices — the reference for the comm fold
// (the engine token carries NO i/2; the σ^{μν}=(i/2)[γ^μ,γ^ν] normalization lives in the scalar).
static NM nCommM(const NM &X, const NM &Y)
{
  NM xy = nMul(X, Y), yx = nMul(Y, X), C;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      C.a[i][j] = Cx{xy.a[i][j].re - yx.a[i][j].re, xy.a[i][j].im - yx.a[i][j].im};
  return C;
}

static double cdiff(Cx a, Cx b) { return std::abs(a.re - b.re) + std::abs(a.im - b.im); }

int main()
{
  std::mt19937 rng(2024);
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  int fails = 0;

  // ---- A) all-slash chains vs chiral_gamma_trace<K> ----
  std::printf("== A: all-slash chains tr(p1..pK) vs chiral_gamma_trace ==\n");
  for (int K : {2, 3, 4, 5, 6, 7, 8, 9, 10}) {
    const int nsym = 4 * K;
    nm::Frame env(ntest::names(nsym));
    std::vector<std::array<nm::Poly, 4>> comp(static_cast<std::size_t>(K));
    for (int m = 0; m < K; ++m)
      for (int mu = 0; mu < 4; ++mu)
        comp[static_cast<std::size_t>(m)][static_cast<std::size_t>(mu)] = env.var(4 * m + mu);
    network::DiracChain chain;
    for (int m = 0; m < K; ++m)
      chain.push_back(network::slash(Momentum{{{1.0, m}}}));
    nm::Poly tr = ntest::contract(env, chain, /*lorentz*/ {}, comp, /*atomDen*/ {});
    // random point
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    Cx sym = nm::ndetail::eval(tr, x, {});
    // ground truth
    double P[12][4];
    for (int m = 0; m < K; ++m)
      for (int mu = 0; mu < 4; ++mu)
        P[m][mu] = x[static_cast<std::size_t>(4 * m + mu)];
    auto toCx = [](std::complex<double> z) { return Cx{z.real(), z.imag()}; };
    Cx num{0, 0};
    switch (K) {
    case 2:
      num = toCx(numtracer::dirac::chiral_gamma_trace<2>(P));
      break;
    case 3:
      num = toCx(numtracer::dirac::chiral_gamma_trace<3>(P));
      break;
    case 4:
      num = toCx(numtracer::dirac::chiral_gamma_trace<4>(P));
      break;
    case 5:
      num = toCx(numtracer::dirac::chiral_gamma_trace<5>(P));
      break;
    case 6:
      num = toCx(numtracer::dirac::chiral_gamma_trace<6>(P));
      break;
    case 7:
      num = toCx(numtracer::dirac::chiral_gamma_trace<7>(P));
      break;
    case 8:
      num = toCx(numtracer::dirac::chiral_gamma_trace<8>(P));
      break;
    case 9:
      num = toCx(numtracer::dirac::chiral_gamma_trace<9>(P));
      break;
    case 10:
      num = toCx(numtracer::dirac::chiral_gamma_trace<10>(P));
      break;
    }
    const double err = cdiff(sym, num);
    std::printf("  K=%2d  monomials=%-4d |sym-num|=%.2e  %s\n", K, tr.size(), err, err < 1e-12 ? "ok" : "FAIL");
    if (!(err < 1e-12)) ++fails;
  }

  // ---- B) free legs closed by a metric: tr(γ^μ p̸ γ_μ q̸) ----
  std::printf("\n== B: tr(g^mu p g_mu q) via metric closure ==\n");
  {
    const int nsym = 8; // p:0..3, q:4..7
    nm::Frame env(ntest::names(nsym));
    std::vector<std::array<nm::Poly, 4>> comp(2);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][static_cast<std::size_t>(mu)] = env.var(mu);
      comp[1][static_cast<std::size_t>(mu)] = env.var(4 + mu);
    }
    network::DiracChain chain = {network::gamma(LorentzIndex{100}), network::slash(Momentum{{{1.0, 0}}}), network::gamma(LorentzIndex{101}), network::slash(Momentum{{{1.0, 1}}})};
    LorentzNet lor = {LorentzTerm{Cx{1, 0}, {ntest::fmet(100, 101)}}};
    nm::Poly tr = ntest::contract(env, chain, lor, comp, {});
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    Cx sym = nm::ndetail::eval(tr, x, {});
    double p[4], q[4];
    for (int mu = 0; mu < 4; ++mu) {
      p[mu] = x[static_cast<std::size_t>(mu)];
      q[mu] = x[static_cast<std::size_t>(4 + mu)];
    }
    Cx num{0, 0};
    for (int mu = 0; mu < 4; ++mu)
      num = num + nTrace(nMul(nMul(nMul(nGamma(mu), nSlash(p)), nGamma(mu)), nSlash(q)));
    const double err = cdiff(sym, num);
    std::printf("  monomials=%d |sym-num|=%.2e  %s\n", tr.size(), err, err < 1e-12 ? "ok" : "FAIL");
    if (!(err < 1e-12)) ++fails;
  }

  // ---- C) projector contraction: monomial-k² (cancels) vs non-monomial-k² (survives) ----
  std::printf("\n== C: P(k) contracting tr(g^mu p g^nu q) ==\n");
  auto runProj = [&](bool monomial) {
    const int nsym = 10; // p:0..3, q:4..7, l1:8, l2:9
    nm::Frame env(ntest::names(nsym));
    std::vector<std::array<nm::Poly, 4>> comp(3);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][static_cast<std::size_t>(mu)] = env.var(mu);
      comp[1][static_cast<std::size_t>(mu)] = env.var(4 + mu);
    }
    comp[2][0] = env.var(8);                // k_0 = l1
    if (!monomial) comp[2][1] = env.var(9); // k_1 = l2 (shifted-like → non-monomial k²)
    // atomDen[0] = k² = Σ comp[2][μ]²
    nm::Poly k2 = env.zero();
    for (int mu = 0; mu < 4; ++mu)
      k2 = k2 + comp[2][static_cast<std::size_t>(mu)] * comp[2][static_cast<std::size_t>(mu)];
    std::vector<nm::Poly> atomDen = {k2};
    network::DiracChain chain = {network::gamma(LorentzIndex{100}), network::slash(Momentum{{{1.0, 0}}}), network::gamma(LorentzIndex{101}), network::slash(Momentum{{{1.0, 1}}})};
    LorentzNet lor = {LorentzTerm{Cx{1, 0}, {ntest::fprojT(100, 101, {{1.0, 2}}, 0)}}};
    nm::Poly tr = ntest::contract(env, chain, lor, comp, atomDen);
    // does any monomial still carry the atom?
    bool hasAtom = false;
    for (const auto &[m, c] : tr.terms)
      if (!m.atoms.empty()) {
        hasAtom = true;
        break;
      }
    // evaluate
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    double k2v = 0;
    {
      double kc[4] = {x[8], monomial ? 0.0 : x[9], 0, 0};
      for (int mu = 0; mu < 4; ++mu)
        k2v += kc[mu] * kc[mu];
    }
    std::vector<double> atomVal = {1.0 / k2v};
    Cx sym = nm::ndetail::eval(tr, x, atomVal);
    // ground truth
    double p[4], q[4], kc[4] = {x[8], monomial ? 0.0 : x[9], 0, 0};
    for (int mu = 0; mu < 4; ++mu) {
      p[mu] = x[static_cast<std::size_t>(mu)];
      q[mu] = x[static_cast<std::size_t>(4 + mu)];
    }
    Cx num{0, 0};
    for (int mu = 0; mu < 4; ++mu)
      for (int nu = 0; nu < 4; ++nu) {
        const double proj = (mu == nu ? 1.0 : 0.0) - kc[mu] * kc[nu] / k2v;
        num = num + nTrace(nMul(nMul(nMul(nGamma(mu), nSlash(p)), nGamma(nu)), nSlash(q))) * Cx{proj, 0};
      }
    const double err = cdiff(sym, num);
    const bool atomOk = monomial ? !hasAtom : hasAtom;
    std::printf("  %-13s monomials=%-3d hasAtom=%d (expect %d) |sym-num|=%.2e  %s\n",
                monomial ? "monomial-k2" : "nonmono-k2", tr.size(), hasAtom, monomial ? 0 : 1, err,
                (err < 1e-12 && atomOk) ? "ok" : "FAIL");
    if (!(err < 1e-12 && atomOk)) ++fails;
  };
  runProj(true);
  runProj(false);

  // ---- D) the network::LorentzNet adapter path (what the generator uses): proj/met builders ----
  std::printf("\n== D: numeric_value_netval over network::LorentzNet (generator path) ==\n");
  {
    const int nsym = 10;
    nm::Frame env(ntest::names(nsym));
    std::vector<std::array<nm::Poly, 4>> comp(3);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][static_cast<std::size_t>(mu)] = env.var(mu);
      comp[1][static_cast<std::size_t>(mu)] = env.var(4 + mu);
    }
    comp[2][0] = env.var(8);
    comp[2][1] = env.var(9); // non-monomial k² → atom survives
    network::DiracChain chain = {network::gamma(LorentzIndex{100}), network::slash(Momentum{{{1.0, 0}}}), network::gamma(LorentzIndex{101}), network::slash(Momentum{{{1.0, 1}}})};
    // Lorentz net via the SAME inv builders the generator emits: P(k)_{100,101}, k = vec id 2, 1/k² env id 7.
    network::LorentzNet lor = ntest::iprojT(100, 101, 2, 7);
    std::vector<nm::Poly> atomDen = ntest::collect_atom_denoms(env, {lor}, comp);
    nm::Poly tr = ntest::contract(env, chain, lor, comp, atomDen);
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    double k2v = x[8] * x[8] + x[9] * x[9];
    std::vector<double> atomVal(8, 0.0);
    atomVal[7] = 1.0 / k2v;
    Cx sym = nm::ndetail::eval(tr, x, atomVal);
    double p[4], q[4], kc[4] = {x[8], x[9], 0, 0};
    for (int mu = 0; mu < 4; ++mu) {
      p[mu] = x[static_cast<std::size_t>(mu)];
      q[mu] = x[static_cast<std::size_t>(4 + mu)];
    }
    Cx num{0, 0};
    for (int mu = 0; mu < 4; ++mu)
      for (int nu = 0; nu < 4; ++nu) {
        const double proj = (mu == nu ? 1.0 : 0.0) - kc[mu] * kc[nu] / k2v;
        num = num + nTrace(nMul(nMul(nMul(nGamma(mu), nSlash(p)), nGamma(nu)), nSlash(q))) * Cx{proj, 0};
      }
    const double err = cdiff(sym, num);
    std::printf("  monomials=%d |sym-num|=%.2e  %s\n", tr.size(), err, err < 1e-12 ? "ok" : "FAIL");
    if (!(err < 1e-12)) ++fails;
  }

  // ---- E) 3-free-leg chain vs the proven dirac_value expansion (guards the variable-elimination
  //         transpose bug: a rank>=2 ASYMMETRIC intermediate factor — only a >=3-free-leg Dirac
  //         tensor produces one; pure-Lorentz nets use symmetric δ/projectors and never trip it) ----
  std::printf("\n== E: 3-free-leg chain vs dirac_value reference (asymmetric intermediate) ==\n");
  {
    const int nsym = 12; // p0:0-3, p1:4-7, k:8 (loop dir), m:9..11 spare
    nm::Frame env(ntest::names(nsym));
    std::vector<std::array<nm::Poly, 4>> comp(3);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][static_cast<std::size_t>(mu)] = env.var(mu);
      comp[1][static_cast<std::size_t>(mu)] = env.var(4 + mu);
    }
    comp[2][0] = env.var(8);
    comp[2][1] = env.var(9); // non-monomial k
    network::DiracChain chain = {network::gamma(LorentzIndex{50}),         network::slash(Momentum{{{1.0, 0}}}), network::gamma(LorentzIndex{51}),
                           network::slash(Momentum{{{1.0, 1}}}), network::gamma(LorentzIndex{52}),         network::slash(Momentum{{{1.0, 0}}})};
    // net contracting legs 50,51 (one projector) and 52 (projector to an aux leg closed by a vector).
    network::LorentzNet net = mul(ntest::iprojT(50, 51, 2, 90), ntest::iprojT(52, 60, 2, 91), ntest::ivec(60, 0));
    std::vector<nm::Poly> atomDen = ntest::collect_atom_denoms(env, {net}, comp);
    nm::Poly mine = ntest::contract(env, chain, net, comp, atomDen);
    network::LorentzNet full = mul(network::dirac_value(chain, 900000), net);
    nm::Poly ref = ntest::contract(env, network::DiracChain{}, full, comp, atomDen);
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    std::vector<double> av(92, 0.0);
    for (int a : {90, 91}) {
      double dv = nm::ndetail::eval(atomDen[static_cast<std::size_t>(a)], x, {}).re;
      av[static_cast<std::size_t>(a)] = dv != 0 ? 1.0 / dv : 0.0;
    }
    double err = cdiff(nm::ndetail::eval(mine, x, av), nm::ndetail::eval(ref, x, av));
    std::printf("  monomials mine=%d ref=%d  |mine-ref|=%.2e  %s\n", mine.size(), ref.size(), err,
                err < 1e-10 ? "ok" : "FAIL");
    if (!(err < 1e-10)) ++fails;
  }

  // ---- F) γ5 (axial) chains — block-diagonal γ5 in the chiral fold vs the direct 4×4 truth.
  //         Covers a leading γ5 with 4 slashes (the nonzero ε structure), a mid-chain γ5 that must
  //         trace to 0, and γ5 with free legs closed by a metric (the pion/axial vertex shape). ----
  std::printf("\n== F: gamma5 chains vs direct 4x4 (nGamma5) ==\n");
  // Generic "trace this closed chain and compare with an explicit 4x4 product at random momenta".
  // Not gamma5-specific — section H reuses it for the charge-conjugation token.
  auto closedChain = [&](const char *tag, const network::DiracChain &chain, int K, auto buildTruth) {
    const int nsym = 4 * K;
    nm::Frame env(ntest::names(nsym));
    std::vector<std::array<nm::Poly, 4>> comp(static_cast<std::size_t>(K));
    for (int m = 0; m < K; ++m)
      for (int mu = 0; mu < 4; ++mu)
        comp[static_cast<std::size_t>(m)][static_cast<std::size_t>(mu)] = env.var(4 * m + mu);
    nm::Poly tr = ntest::contract(env, chain, {}, comp, {});
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    Cx sym = nm::ndetail::eval(tr, x, {});
    double P[8][4];
    for (int m = 0; m < K; ++m)
      for (int mu = 0; mu < 4; ++mu)
        P[m][mu] = x[static_cast<std::size_t>(4 * m + mu)];
    Cx num = buildTruth(P);
    const double err = cdiff(sym, num);
    std::printf("  %-22s monomials=%-3d |sym|=%.2e |sym-num|=%.2e  %s\n", tag, tr.size(),
                std::abs(sym.re) + std::abs(sym.im), err, err < 1e-12 ? "ok" : "FAIL");
    if (!(err < 1e-12)) ++fails;
  };
  // F1: tr(γ5 p̸ q̸ r̸ s̸) — leading γ5, the nonzero axial (ε) structure.
  closedChain(
      "g5 p q r s",
      {network::gamma5(), network::slash(Momentum{{{1.0, 0}}}), network::slash(Momentum{{{1.0, 1}}}), network::slash(Momentum{{{1.0, 2}}}), network::slash(Momentum{{{1.0, 3}}})},
      4, [&](double P[8][4]) {
        return nTrace(nMul(nMul(nMul(nMul(nGamma5(), nSlash(P[0])), nSlash(P[1])), nSlash(P[2])), nSlash(P[3])));
      });
  // F2: tr(p̸ γ5 q̸) — mid-chain γ5, two slashes → traces to 0.
  closedChain("p g5 q (=0)", {network::slash(Momentum{{{1.0, 0}}}), network::gamma5(), network::slash(Momentum{{{1.0, 1}}})}, 2,
           [&](double P[8][4]) { return nTrace(nMul(nMul(nSlash(P[0]), nGamma5()), nSlash(P[1]))); });
  // F3: free legs + γ5 closed by a metric, with enough slashes to stay NONZERO:
  //     tr(γ5 γ^μ p̸ γ_μ q̸ r̸ s̸) = Σ_μ tr(γ5 γ^μ p̸ γ^μ q̸ r̸ s̸)  (= −2 tr(γ5 p̸ q̸ r̸ s̸) ≠ 0).
  {
    const int nsym = 16; // p:0-3 q:4-7 r:8-11 s:12-15
    nm::Frame env(ntest::names(nsym));
    std::vector<std::array<nm::Poly, 4>> comp(4);
    for (int m = 0; m < 4; ++m)
      for (int mu = 0; mu < 4; ++mu)
        comp[static_cast<std::size_t>(m)][static_cast<std::size_t>(mu)] = env.var(4 * m + mu);
    network::DiracChain chain = {network::gamma5(),
                           network::gamma(LorentzIndex{100}),
                           network::slash(Momentum{{{1.0, 0}}}),
                           network::gamma(LorentzIndex{101}),
                           network::slash(Momentum{{{1.0, 1}}}),
                           network::slash(Momentum{{{1.0, 2}}}),
                           network::slash(Momentum{{{1.0, 3}}})};
    LorentzNet lor = {LorentzTerm{Cx{1, 0}, {ntest::fmet(100, 101)}}};
    nm::Poly tr = ntest::contract(env, chain, lor, comp, {});
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    Cx sym = nm::ndetail::eval(tr, x, {});
    double P[4][4];
    for (int m = 0; m < 4; ++m)
      for (int mu = 0; mu < 4; ++mu)
        P[m][mu] = x[static_cast<std::size_t>(4 * m + mu)];
    Cx num{0, 0};
    for (int mu = 0; mu < 4; ++mu)
      num =
          num + nTrace(nMul(nMul(nMul(nMul(nMul(nMul(nGamma5(), nGamma(mu)), nSlash(P[0])), nGamma(mu)), nSlash(P[1])),
                                 nSlash(P[2])),
                            nSlash(P[3])));
    const double err = cdiff(sym, num);
    std::printf("  %-22s monomials=%-3d |sym|=%.2e |sym-num|=%.2e  %s\n", "g5 g^mu p g_mu q r s", tr.size(),
                std::abs(sym.re) + std::abs(sym.im), err, err < 1e-12 ? "ok" : "FAIL");
    if (!(err < 1e-12)) ++fails;
  }

  // ---- G) the BARE commutator token `comm` vs an explicit `[X,Y]` 4×4, in all three leg modes the
  //         front-end emits: both-free (`comm`, σ^{μν} with two open legs), both-slash (`comm`,
  //         the struct-7 external projector [B̸,Q̸]), and free-slash (`comm`, a loop vertex
  //         σ^{μν}k_ν). Each is ONE chain token (no commutator split). Reference: tr(comm · p̸ q̸). ----
  std::printf("\n== G: commutator token comm vs explicit [X,Y] 4x4 (free/free, slash/slash, free/slash) ==\n");
  {
    const int nsym = 12; // p:0-3, q:4-7, r:8-11  (p,q close the chain; r is leg-B's slash momentum)
    nm::Frame env(ntest::names(nsym));
    std::vector<std::array<nm::Poly, 4>> comp(3);
    for (int mu = 0; mu < 4; ++mu)
      for (int m = 0; m < 3; ++m)
        comp[static_cast<std::size_t>(m)][static_cast<std::size_t>(mu)] = env.var(4 * m + mu);
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    double p[4], q[4], r[4];
    for (int mu = 0; mu < 4; ++mu) {
      p[mu] = x[static_cast<std::size_t>(mu)];
      q[mu] = x[static_cast<std::size_t>(4 + mu)];
      r[mu] = x[static_cast<std::size_t>(8 + mu)];
    }

    // G1) both-free [γ^a, γ^b]: rank-2 tensor T^{ab} = tr([γ^a,γ^b] p̸ q̸), 16 entries.
    {
      network::DiracChain chain = {network::comm(LorentzIndex{100}, LorentzIndex{101}), network::slash(Momentum{{{1.0, 0}}}), network::slash(Momentum{{{1.0, 1}}})};
      auto T = nm::numeric_dirac(nsym, chain, comp);
      double err = 0, maxabs = 0;
      for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
          Cx sym = nm::ndetail::eval(T.entries[static_cast<std::size_t>(a * 4 + b)], x, {});
          Cx num = nTrace(nMul(nMul(nCommM(nGamma(a), nGamma(b)), nSlash(p)), nSlash(q)));
          maxabs = std::max(maxabs, std::abs(num.re) + std::abs(num.im));
          err = std::max(err, cdiff(sym, num));
        }
      std::printf("  G1 free/free  tr([g^a,g^b] p q): 16 entries |max|=%.2e max|sym-num|=%.2e %s\n", maxabs, err,
                  err < 1e-12 ? "ok" : "FAIL");
      if (!(err < 1e-12)) ++fails;
    }
    // G2) both-slash [r̸, p̸]: scalar T = tr([r̸,p̸] p̸ q̸) (no open legs).
    {
      network::DiracChain chain = {network::comm(Momentum{{{1.0, 2}}}, Momentum{{{1.0, 0}}}), network::slash(Momentum{{{1.0, 0}}}), network::slash(Momentum{{{1.0, 1}}})};
      auto T = nm::numeric_dirac(nsym, chain, comp);
      Cx sym = nm::ndetail::eval(T.entries[0], x, {});
      Cx num = nTrace(nMul(nMul(nCommM(nSlash(r), nSlash(p)), nSlash(p)), nSlash(q)));
      double err = cdiff(sym, num);
      std::printf("  G2 slash/slash tr([r,p] p q): |num|=%.2e |sym-num|=%.2e %s\n", std::abs(num.re) + std::abs(num.im),
                  err, err < 1e-12 ? "ok" : "FAIL");
      if (!(err < 1e-12)) ++fails;
    }
    // G3) free-slash [γ^a, r̸]: rank-1 tensor T^{a} = tr([γ^a,r̸] p̸ q̸), 4 entries.
    {
      network::DiracChain chain = {network::comm(LorentzIndex{100}, Momentum{{{1.0, 2}}}), network::slash(Momentum{{{1.0, 0}}}), network::slash(Momentum{{{1.0, 1}}})};
      auto T = nm::numeric_dirac(nsym, chain, comp);
      double err = 0, maxabs = 0;
      for (int a = 0; a < 4; ++a) {
        Cx sym = nm::ndetail::eval(T.entries[static_cast<std::size_t>(a)], x, {});
        Cx num = nTrace(nMul(nMul(nCommM(nGamma(a), nSlash(r)), nSlash(p)), nSlash(q)));
        maxabs = std::max(maxabs, std::abs(num.re) + std::abs(num.im));
        err = std::max(err, cdiff(sym, num));
      }
      std::printf("  G3 free/slash tr([g^a,r] p q): 4 entries |max|=%.2e max|sym-num|=%.2e %s\n", maxabs, err,
                  err < 1e-12 ? "ok" : "FAIL");
      if (!(err < 1e-12)) ++fails;
    }
  }

  // ---- H) the charge-conjugation matrix C, and the word-reversal identity it licenses ----------
  // C is a plain constexpr table (dirac/dirac_data.hpp), so nothing else would notice if an entry
  // were wrong -- this block is its safety net, exactly as section A is the gamma table's.
  //
  // The last two checks are the load-bearing ones for the engine:
  //   * tr(Mn..M1) == tr(M1..Mn) for an even gamma word. This is WHY orderDiracFacs may walk a
  //     closed spinor loop in either direction (a uniformly reversed traversal is routine -- the
  //     walk seeds on a symmetric spinor delta whose label order is arbitrary). It follows from
  //     C gamma^T C^-1 = -gamma: a reversal costs (-1)^(#gamma), and odd words already trace to 0.
  //   * conjugating a TRANSPOSED word by X = C.gamma5 reverses it with NO per-gamma sign and no
  //     momentum flip, because X gamma_mu^T X^-1 = +gamma_mu. Bare C would cost a sign per gamma
  //     (equivalently p -> -p on the reversed segment); X is what a diquark vertex supplies anyway.
  std::printf("\n== H: charge conjugation C = gamma^2 gamma^4 ==\n");
  {
    int fH = 0;
    auto nId = [] {
      NM I;
      for (int i = 0; i < 4; ++i)
        I.a[i][i] = Cx{1, 0};
      return I;
    };
    auto nT = [](const NM &A) {
      NM R;
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          R.a[i][j] = A.a[j][i];
      return R;
    };
    auto nNeg = [](const NM &A) {
      NM R;
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          R.a[i][j] = Cx{-A.a[i][j].re, -A.a[i][j].im};
      return R;
    };
    auto nEq = [](const NM &A, const NM &B) {
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
          if (std::abs(A.a[i][j].re - B.a[i][j].re) > 1e-12) return false;
          if (std::abs(A.a[i][j].im - B.a[i][j].im) > 1e-12) return false;
        }
      return true;
    };
    auto say = [&](const char *what, bool ok) {
      std::printf("  %-46s %s\n", what, ok ? "ok" : "FAIL");
      if (!ok) ++fH;
    };
    const NM C = nCmat(), g5 = nGamma5();
    // C^T = C^-1 = -C^dagger = -C  (C is real here, so the dagger check is the transpose check)
    say("C^T == -C", nEq(nT(C), nNeg(C)));
    say("C * (-C) == 1  (i.e. C^-1 == -C)", nEq(nMul(C, nNeg(C)), nId()));
    const NM Cinv = nNeg(C);
    // C gamma_mu^T C^-1 = -gamma_mu  and  C gamma5^T C^-1 = +gamma5
    bool okg = true;
    for (int mu = 0; mu < 4; ++mu)
      okg = okg && nEq(nMul(nMul(C, nT(nGamma(mu))), Cinv), nNeg(nGamma(mu)));
    say("C gamma_mu^T C^-1 == -gamma_mu (all mu)", okg);
    say("C gamma5^T C^-1 == +gamma5", nEq(nMul(nMul(C, nT(g5)), Cinv), g5));
    // C is block-DIAGONAL in the Weyl split: like gamma5 it must not flip the antidiagonal parity
    // that numeric_dirac uses to zero odd traces, and it costs no full 2x2 multiply.
    bool blockDiag = true;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        if ((i < 2) != (j < 2) && (C.a[i][j].re != 0 || C.a[i][j].im != 0)) blockDiag = false;
    say("C is block-diagonal in the Weyl split", blockDiag);
    // C.gamma5 antisymmetric -> the scalar diquark bilinear q^T C gamma5 q survives once the
    // antisymmetric colour and flavour factors are included.
    const NM X = nMul(C, g5);
    say("(C gamma5)^T == -(C gamma5)", nEq(nT(X), nNeg(X)));
    // Reversing a closed EVEN gamma word leaves the trace alone: tr(gd gc gb ga) == tr(ga gb gc gd).
    bool okrev = true;
    for (int a = 0; a < 4 && okrev; ++a)
      for (int b = 0; b < 4 && okrev; ++b)
        for (int c = 0; c < 4 && okrev; ++c)
          for (int d = 0; d < 4 && okrev; ++d) {
            const Cx fwd = nTrace(nMul(nMul(nGamma(a), nGamma(b)), nMul(nGamma(c), nGamma(d))));
            const Cx rev = nTrace(nMul(nMul(nGamma(d), nGamma(c)), nMul(nGamma(b), nGamma(a))));
            if (std::abs(fwd.re - rev.re) > 1e-12 || std::abs(fwd.im - rev.im) > 1e-12) okrev = false;
          }
    say("tr(reversed 4-gamma word) == tr(word)", okrev);
    // X W^T X^-1 == reverse(W), exhaustively over all 3-letter words. X^-1 = gamma5^-1 C^-1 = -gamma5 C.
    const NM Xinv = nNeg(nMul(g5, C));
    say("X * X^-1 == 1  (X = C gamma5)", nEq(nMul(X, Xinv), nId()));
    bool okw = true;
    for (int a = 0; a < 4 && okw; ++a)
      for (int b = 0; b < 4 && okw; ++b)
        for (int c = 0; c < 4 && okw; ++c) {
          const NM W = nMul(nMul(nGamma(a), nGamma(b)), nGamma(c));
          const NM R = nMul(nMul(nGamma(c), nGamma(b)), nGamma(a));
          if (!nEq(nMul(nMul(X, nT(W)), Xinv), R)) okw = false;
        }
    say("X W^T X^-1 == reverse(W), all 3-words", okw);
    // kC must BE gamma^2 gamma^4 — the table is typed out by hand, so pin it against the gamma
    // table rather than trusting the two to agree. (mu is 0-based: gamma^2 is mu=1, gamma^4 is mu=3.)
    say("kC == gamma^2 gamma^4", nEq(C, nMul(nGamma(1), nGamma(3))));
    fails += fH;
  }

  // ---- I) the DFac::C ENGINE token, traced by numeric_dirac, vs the same explicit 4x4 product. C is
  //         block-diagonal like gamma5, so it must not flip the antidiagonal parity — I4 is the
  //         check that it does not (an odd gamma/slash count still traces to 0 with a C present). ---
  std::printf("\n== I: DFac::C engine token vs direct 4x4 ==\n");
  // I1: tr(C p q) — leading C.
  closedChain("C p q", {network::chargeC(), network::slash(Momentum{{{1.0, 0}}}), network::slash(Momentum{{{1.0, 1}}})}, 2,
              [&](double P[8][4]) { return nTrace(nMul(nMul(nCmat(), nSlash(P[0])), nSlash(P[1]))); });
  // I2: tr(p C q C r s) — two mid-chain C's, so the running product visits both Weyl blocks.
  closedChain("p C q C r s",
              {network::slash(Momentum{{{1.0, 0}}}), network::chargeC(), network::slash(Momentum{{{1.0, 1}}}), network::chargeC(),
               network::slash(Momentum{{{1.0, 2}}}), network::slash(Momentum{{{1.0, 3}}})},
              4, [&](double P[8][4]) {
                return nTrace(nMul(nMul(nMul(nMul(nMul(nSlash(P[0]), nCmat()), nSlash(P[1])), nCmat()),
                                        nSlash(P[2])),
                                   nSlash(P[3])));
              });
  // I3: tr(C g5 p q) — the X = C.gamma5 combination a diquark vertex actually supplies; both
  //     block-diagonal tokens adjacent, exercising the gamma5 sign flip on top of the C multiply.
  closedChain("C g5 p q",
              {network::chargeC(), network::gamma5(), network::slash(Momentum{{{1.0, 0}}}), network::slash(Momentum{{{1.0, 1}}})}, 2,
              [&](double P[8][4]) {
                return nTrace(nMul(nMul(nMul(nCmat(), nGamma5()), nSlash(P[0])), nSlash(P[1])));
              });
  // I4: tr(C p) — ONE antidiagonal factor, so the parity rule must still zero it. If C were
  //     miscounted as antidiagonal this would come back nonzero (and every even chain would vanish).
  closedChain("C p (=0)", {network::chargeC(), network::slash(Momentum{{{1.0, 0}}})}, 1,
              [&](double P[8][4]) { return nTrace(nMul(nCmat(), nSlash(P[0]))); });

  // ---- J) TRANSPOSED tokens: the engine vs an independent INDEX-NETWORK contractor ---------------
  //
  // This is the correctness gate for `DFac::transposed`, and it is deliberately a property test
  // rather than a derivation. A closed spinor loop is a cycle in a degree-2 index network, so
  //     sum_{l0..l(n-1)} M1[l0,l1] M2[l1,l2] .. Mn[l(n-1),l0]
  // is a trace of a matrix product ONLY where each factor's declared (din,dout) order matches the
  // traversal; where it does not, the network wants M^T. The oracle below evaluates that sum
  // DIRECTLY -- brute force over all 4^n label assignments, no walk, no orientation model, no gamma
  // algebra -- so it shares nothing with the engine path it grades.
  //
  // Nets are randomised over token kind (gamma, gamma5, C, generic slash, sigma) AND over which of
  // each factor's two slots is declared first, which is exactly the "scrambled index order" that
  // makes the fast path dangerous. Measured before this feature existed: the no-transpose rule is
  // wrong on ~36% of such nets.
  std::printf("\n== J: transposed tokens vs a direct index-network contraction ==\n");
  {
    int fJ = 0, tested = 0, nonzero = 0;
    std::mt19937 jr(20260907);
    auto pick4 = [&] { return static_cast<int>(jr() % 4); };
    for (int trial = 0; trial < 3000; ++trial) {
      const int n = 2 + static_cast<int>(jr() % 4); // 2..5 factors, so 4^n <= 1024 assignments
      // Build a cycle over labels 0..n-1: factor k joins label k to label k+1 (mod n). For each
      // factor, randomly decide which of those two labels is its DECLARED din.
      std::vector<NM> mats(static_cast<std::size_t>(n));
      std::vector<std::array<int, 2>> labs(static_cast<std::size_t>(n));
      network::DiracChain chain;
      bool anyTransposed = false;
      for (int k = 0; k < n; ++k) {
        const int a = k, b = (k + 1) % n;
        const bool flip = (jr() % 2) == 0; // declared order reversed w.r.t. the traversal a -> b
        labs[static_cast<std::size_t>(k)] = flip ? std::array<int, 2>{b, a} : std::array<int, 2>{a, b};
        network::DFac tok;
        const int kind = static_cast<int>(jr() % 5);
        if (kind == 0) {
          // A FREE gamma leg would leave an open Lorentz index; use the unit slash of momentum mu
          // instead, whose components comp[mu] are the unit vector e_mu, so it IS gamma^mu.
          const int mu = pick4();
          mats[static_cast<std::size_t>(k)] = nGamma(mu);
          tok = network::slash(Momentum{{{1.0, mu}}});
        } else if (kind == 1) {
          mats[static_cast<std::size_t>(k)] = nGamma5();
          tok = network::gamma5();
        } else if (kind == 2) {
          mats[static_cast<std::size_t>(k)] = nCmat();
          tok = network::chargeC();
        } else if (kind == 3) {
          double comp[4];
          for (int mu = 0; mu < 4; ++mu)
            comp[mu] = static_cast<double>(static_cast<int>(jr() % 9) - 4);
          mats[static_cast<std::size_t>(k)] = nSlash(comp);
          std::vector<std::pair<double, int>> vlc;
          for (int mu = 0; mu < 4; ++mu)
            if (comp[mu] != 0.0) vlc.push_back({comp[mu], mu});
          if (vlc.empty()) vlc.push_back({1.0, 0});
          // rebuild the truth from the vlc actually emitted, so the two sides cannot drift
          double c2[4] = {0, 0, 0, 0};
          for (const auto &pr : vlc)
            c2[pr.second] += pr.first;
          mats[static_cast<std::size_t>(k)] = nSlash(c2);
          tok = network::slash(Momentum{vlc});
        } else {
          const int ma = pick4(), mb = pick4();
          double ca[4] = {0, 0, 0, 0}, cb[4] = {0, 0, 0, 0};
          ca[ma] = 1.0;
          cb[mb] = 1.0;
          mats[static_cast<std::size_t>(k)] = nCommM(nSlash(ca), nSlash(cb));
          tok = network::comm(Momentum{{{1.0, ma}}}, Momentum{{{1.0, mb}}});
        }
        // The traversal below always walks a -> b, so the factor is transposed exactly when its
        // declared order was flipped.
        if (flip) {
          tok = network::transposed(tok);
          anyTransposed = true;
        }
        chain.push_back(tok);
      }
      if (!anyTransposed) continue; // nothing to grade
      // --- oracle: sum over every label assignment, straight from the declared index order -------
      Cx truth{0, 0};
      int total = 1;
      for (int k = 0; k < n; ++k)
        total *= 4;
      for (int asg = 0; asg < total; ++asg) {
        std::vector<int> idx(static_cast<std::size_t>(n));
        int t = asg;
        for (int k = 0; k < n; ++k) {
          idx[static_cast<std::size_t>(k)] = t % 4;
          t /= 4;
        }
        Cx pr{1, 0};
        for (int k = 0; k < n; ++k) {
          const auto &lb = labs[static_cast<std::size_t>(k)];
          pr = pr * mats[static_cast<std::size_t>(k)]
                       .a[idx[static_cast<std::size_t>(lb[0])]][idx[static_cast<std::size_t>(lb[1])]];
          if (pr.re == 0 && pr.im == 0) break;
        }
        truth = truth + pr;
      }
      // --- engine: the same chain in traversal order, with the transpose flags ------------------
      const int nsym = 4;
      nm::Frame env(ntest::names(nsym));
      std::vector<std::array<nm::Poly, 4>> comp(4);
      for (int m = 0; m < 4; ++m)
        for (int mu = 0; mu < 4; ++mu)
          comp[static_cast<std::size_t>(m)][static_cast<std::size_t>(mu)] =
              (m == mu) ? nm::PolyFactory::constant(nsym, Cx{1, 0}) : nm::PolyFactory::zero(nsym);
      nm::Poly tr = ntest::contract(env, chain, {}, comp, {});
      const std::vector<double> x(static_cast<std::size_t>(nsym), 0.0);
      const Cx got = nm::ndetail::eval(tr, x, {});
      ++tested;
      if (std::abs(truth.re) + std::abs(truth.im) > 1e-9) ++nonzero;
      if (cdiff(got, truth) > 1e-9) {
        if (fJ < 3)
          std::printf("  FAIL n=%d  engine=(%.6g,%.6g) oracle=(%.6g,%.6g)\n", n, got.re, got.im,
                      truth.re, truth.im);
        ++fJ;
      }
    }
    std::printf("  %d scrambled nets with >=1 transposed token (%d non-zero): %d mismatch(es) %s\n",
                tested, nonzero, fJ, fJ == 0 ? "ok" : "FAIL");
    if (nonzero < 200) {
      std::printf("  FAIL: only %d non-zero nets -- the generator is not exercising the engine\n", nonzero);
      ++fJ;
    }
    fails += fJ;
  }

  // ---- K) a TRANSPOSED SLOT: dtrslot(k) vs the same chain written out by hand -------------------
  //
  // A slot is a sum of chains, so (A1..An)^T = An^T..A1^T -- transposing one means reversing the
  // option's token list and transposing each token, which happens at the dress_enumerate splice
  // (numeric_contract.hpp) rather than in the front end, because the slot caches key on the slot
  // expression. Grade that splice against the explicitly reversed-and-transposed chain, which goes
  // through the ordinary (non-dressed) path and so shares none of the splice logic.
  std::printf("\n== K: transposed slot splice vs an explicit reversed chain ==\n");
  {
    const int nsym = 16; // four momenta, 0-3
    nm::Frame env(ntest::names(nsym));
    std::vector<std::array<nm::Poly, 4>> comp(4);
    for (int m = 0; m < 4; ++m)
      for (int mu = 0; mu < 4; ++mu)
        comp[static_cast<std::size_t>(m)][static_cast<std::size_t>(mu)] = env.var(4 * m + mu);
    // A slot holding a TWO-token chain, spliced into a loop with two more slashes. Four slashes in
    // total, so the trace is non-zero (a gamma5 with only two would vanish and prove nothing), and
    // reversing the slot genuinely reorders the product.
    nm::DSlot slot{nm::DSlotOpt{Cx{1, 0}, {}, {network::slash(Momentum{{{1.0, 0}}}), network::slash(Momentum{{{1.0, 1}}})}, {}}};
    const std::vector<nm::DChainTok> chainT = {nm::dtrslot(0), nm::dtfix(network::slash(Momentum{{{1.0, 2}}})),
                                               nm::dtfix(network::slash(Momentum{{{1.0, 3}}}))};
    const nm::Poly viaSlot = ntest::contract_structural(env, chainT, {slot}, {}, comp, {});
    // by hand: (p0 p1)^T = p1^T p0^T -- reversed, every token transposed.
    const network::DiracChain explicitChain = {network::transposed(network::slash(Momentum{{{1.0, 1}}})),
                                             network::transposed(network::slash(Momentum{{{1.0, 0}}})),
                                             network::slash(Momentum{{{1.0, 2}}}), network::slash(Momentum{{{1.0, 3}}})};
    const nm::Poly viaExplicit = ntest::contract(env, explicitChain, {}, comp, {});
    // the UNtransposed slot, as a control that the two spellings are not trivially equal
    const std::vector<nm::DChainTok> chainU = {nm::dtslot(0), nm::dtfix(network::slash(Momentum{{{1.0, 2}}})),
                                               nm::dtfix(network::slash(Momentum{{{1.0, 3}}}))};
    const nm::Poly viaPlain = ntest::contract_structural(env, chainU, {slot}, {}, comp, {});
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    const Cx a = nm::ndetail::eval(viaSlot, x, {}), b = nm::ndetail::eval(viaExplicit, x, {}), c = nm::ndetail::eval(viaPlain, x, {});
    const double err = cdiff(a, b);
    std::printf("  dtrslot vs explicit reversed chain: |a|=%.3e |a-b|=%.3e %s\n",
                std::abs(a.re) + std::abs(a.im), err, err < 1e-12 ? "ok" : "FAIL");
    if (!(err < 1e-12)) ++fails;
    if (std::abs(a.re) + std::abs(a.im) < 1e-12) {
      std::printf("  FAIL: the transposed-slot trace is zero -- the check proves nothing\n");
      ++fails;
    }
    std::printf("  control: untransposed slot differs from transposed: |a-c|=%.3e %s\n", cdiff(a, c),
                cdiff(a, c) > 1e-12 ? "ok" : "FAIL");
    if (!(cdiff(a, c) > 1e-12)) ++fails;
  }

  std::printf(fails == 0 ? "\nALL TESTS PASSED\n" : "\nTESTS FAILED\n");
  return fails ? 1 : 0;
}
