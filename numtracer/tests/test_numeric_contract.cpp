// Part B (task #24): correctness gate for the numeric (matrix-product) contraction backend.
//
// Validates the engine in et/numeric/* against independent numeric ground truth:
//   A) all-slash Dirac chains tr(p̸₁…p̸_K), K=2..10, vs the dense oracle chiral_gamma_trace<K>;
//   B) a free-leg chain closed by a metric, tr(γ^μ p̸ γ_μ q̸) = Σ_μ tr(γ^μ p̸ γ^μ q̸), vs a direct
//      numeric 4×4 contraction;
//   C) a transverse projector P(k) contracting two free legs — the monomial-k² case (unit-direction
//      loop → the 1/k² atom CANCELS, no surviving atom) and the non-monomial-k² case (shifted line →
//      the atom SURVIVES) — value-checked vs a direct numeric contraction.
//
// Symbolic momentum components are MPoly variables; the engine result is evaluated at random points
// and compared to the numeric truth to ≤ 1e-12. Build via the test CMake (adds -I include).
#include "numtracer/dirac/dense_trace.hpp" // chiral_gamma_trace<K>
#include "numtracer/dirac/dirac_data.hpp"  // kGamma
#include "numtracer/numeric/numeric_contract.hpp"

#include <array>
#include <complex>
#include <cstdio>
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
// bare commutator [X,Y] = X·Y − Y·X of two explicit 4×4 matrices — the reference for the dcomm fold
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
    std::vector<std::array<nm::MPoly, 4>> comp(static_cast<std::size_t>(K));
    for (int m = 0; m < K; ++m)
      for (int mu = 0; mu < 4; ++mu)
        comp[static_cast<std::size_t>(m)][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 4 * m + mu);
    network::DiracNet chain;
    for (int m = 0; m < K; ++m)
      chain.push_back(network::dslash({{1.0, m}}));
    nm::MPoly tr = nm::numeric_value(nsym, chain, /*lorentz*/ {}, comp, /*atomDen*/ {});
    // random point
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    Cx sym = nm::eval(tr, x, {});
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
    std::printf("  K=%2d  monomials=%-4zu |sym-num|=%.2e  %s\n", K, tr.size(), err, err < 1e-12 ? "ok" : "FAIL");
    if (!(err < 1e-12)) ++fails;
  }

  // ---- B) free legs closed by a metric: tr(γ^μ p̸ γ_μ q̸) ----
  std::printf("\n== B: tr(g^mu p g_mu q) via metric closure ==\n");
  {
    const int nsym = 8; // p:0..3, q:4..7
    std::vector<std::array<nm::MPoly, 4>> comp(2);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, mu);
      comp[1][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 4 + mu);
    }
    network::DiracNet chain = {network::dgamma(100), network::dslash({{1.0, 0}}), network::dgamma(101), network::dslash({{1.0, 1}})};
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 101)}}};
    nm::MPoly tr = nm::numeric_value(nsym, chain, lor, comp, {});
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    Cx sym = nm::eval(tr, x, {});
    double p[4], q[4];
    for (int mu = 0; mu < 4; ++mu) {
      p[mu] = x[static_cast<std::size_t>(mu)];
      q[mu] = x[static_cast<std::size_t>(4 + mu)];
    }
    Cx num{0, 0};
    for (int mu = 0; mu < 4; ++mu)
      num = num + nTrace(nMul(nMul(nMul(nGamma(mu), nSlash(p)), nGamma(mu)), nSlash(q)));
    const double err = cdiff(sym, num);
    std::printf("  monomials=%zu |sym-num|=%.2e  %s\n", tr.size(), err, err < 1e-12 ? "ok" : "FAIL");
    if (!(err < 1e-12)) ++fails;
  }

  // ---- C) projector contraction: monomial-k² (cancels) vs non-monomial-k² (survives) ----
  std::printf("\n== C: P(k) contracting tr(g^mu p g^nu q) ==\n");
  auto runProj = [&](bool monomial) {
    const int nsym = 10; // p:0..3, q:4..7, l1:8, l2:9
    std::vector<std::array<nm::MPoly, 4>> comp(3);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, mu);
      comp[1][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 4 + mu);
    }
    comp[2][0] = nm::MPoly::var(nsym, 8);                // k_0 = l1
    if (!monomial) comp[2][1] = nm::MPoly::var(nsym, 9); // k_1 = l2 (shifted-like → non-monomial k²)
    // atomDen[0] = k² = Σ comp[2][μ]²
    nm::MPoly k2(nsym);
    for (int mu = 0; mu < 4; ++mu)
      k2 = k2 + comp[2][static_cast<std::size_t>(mu)] * comp[2][static_cast<std::size_t>(mu)];
    std::vector<nm::MPoly> atomDen = {k2};
    network::DiracNet chain = {network::dgamma(100), network::dslash({{1.0, 0}}), network::dgamma(101), network::dslash({{1.0, 1}})};
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nprojT(100, 101, {{1.0, 2}}, 0)}}};
    nm::MPoly tr = nm::numeric_value(nsym, chain, lor, comp, atomDen);
    // does any monomial still carry the atom?
    bool hasAtom = false;
    for (const auto &[m, c] : tr.t)
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
    Cx sym = nm::eval(tr, x, atomVal);
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
    std::printf("  %-13s monomials=%-3zu hasAtom=%d (expect %d) |sym-num|=%.2e  %s\n",
                monomial ? "monomial-k2" : "nonmono-k2", tr.size(), hasAtom, monomial ? 0 : 1, err,
                (err < 1e-12 && atomOk) ? "ok" : "FAIL");
    if (!(err < 1e-12 && atomOk)) ++fails;
  };
  runProj(true);
  runProj(false);

  // ---- D) the network::NetVal adapter path (what the generator uses): proj/met builders ----
  std::printf("\n== D: numeric_value_netval over network::NetVal (generator path) ==\n");
  {
    const int nsym = 10;
    std::vector<std::array<nm::MPoly, 4>> comp(3);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, mu);
      comp[1][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 4 + mu);
    }
    comp[2][0] = nm::MPoly::var(nsym, 8);
    comp[2][1] = nm::MPoly::var(nsym, 9); // non-monomial k² → atom survives
    network::DiracNet chain = {network::dgamma(100), network::dslash({{1.0, 0}}), network::dgamma(101), network::dslash({{1.0, 1}})};
    // Lorentz net via the SAME inv builders the generator emits: P(k)_{100,101}, k = vec id 2, 1/k² env id 7.
    network::NetVal lor = network::projT(100, 101, 2, 7);
    std::vector<nm::MPoly> atomDen = nm::collect_atom_denoms(nsym, {lor}, comp);
    nm::MPoly tr = nm::numeric_value_netval(nsym, chain, lor, comp, atomDen);
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    double k2v = x[8] * x[8] + x[9] * x[9];
    std::vector<double> atomVal(8, 0.0);
    atomVal[7] = 1.0 / k2v;
    Cx sym = nm::eval(tr, x, atomVal);
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
    std::printf("  monomials=%zu |sym-num|=%.2e  %s\n", tr.size(), err, err < 1e-12 ? "ok" : "FAIL");
    if (!(err < 1e-12)) ++fails;
  }

  // ---- E) 3-free-leg chain vs the proven dirac_value expansion (guards the variable-elimination
  //         transpose bug: a rank>=2 ASYMMETRIC intermediate factor — only a >=3-free-leg Dirac
  //         tensor produces one; pure-Lorentz nets use symmetric δ/projectors and never trip it) ----
  std::printf("\n== E: 3-free-leg chain vs dirac_value reference (asymmetric intermediate) ==\n");
  {
    const int nsym = 12; // p0:0-3, p1:4-7, k:8 (loop dir), m:9..11 spare
    std::vector<std::array<nm::MPoly, 4>> comp(3);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, mu);
      comp[1][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 4 + mu);
    }
    comp[2][0] = nm::MPoly::var(nsym, 8);
    comp[2][1] = nm::MPoly::var(nsym, 9); // non-monomial k
    network::DiracNet chain = {network::dgamma(50),         network::dslash({{1.0, 0}}), network::dgamma(51),
                           network::dslash({{1.0, 1}}), network::dgamma(52),         network::dslash({{1.0, 0}})};
    // net contracting legs 50,51 (one projector) and 52 (projector to an aux leg closed by a vector).
    network::NetVal net = network::contract(network::projT(50, 51, 2, 90), network::projT(52, 60, 2, 91), network::vec(60, 0));
    std::vector<nm::MPoly> atomDen = nm::collect_atom_denoms(nsym, {net}, comp);
    nm::MPoly mine = nm::numeric_value_netval(nsym, chain, net, comp, atomDen);
    network::NetVal full = network::contract(network::dirac_value(chain, 900000), net);
    nm::MPoly ref = nm::numeric_value_netval(nsym, network::DiracNet{}, full, comp, atomDen);
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    std::vector<double> av(92, 0.0);
    for (int a : {90, 91}) {
      double dv = nm::eval(atomDen[static_cast<std::size_t>(a)], x, {}).re;
      av[static_cast<std::size_t>(a)] = dv != 0 ? 1.0 / dv : 0.0;
    }
    double err = cdiff(nm::eval(mine, x, av), nm::eval(ref, x, av));
    std::printf("  monomials mine=%zu ref=%zu  |mine-ref|=%.2e  %s\n", mine.size(), ref.size(), err,
                err < 1e-10 ? "ok" : "FAIL");
    if (!(err < 1e-10)) ++fails;
  }

  // ---- F) γ5 (axial) chains — block-diagonal γ5 in the chiral fold vs the direct 4×4 truth.
  //         Covers a leading γ5 with 4 slashes (the nonzero ε structure), a mid-chain γ5 that must
  //         trace to 0, and γ5 with free legs closed by a metric (the pion/axial vertex shape). ----
  std::printf("\n== F: gamma5 chains vs direct 4x4 (nGamma5) ==\n");
  auto g5closed = [&](const char *tag, const network::DiracNet &chain, int K, auto buildTruth) {
    const int nsym = 4 * K;
    std::vector<std::array<nm::MPoly, 4>> comp(static_cast<std::size_t>(K));
    for (int m = 0; m < K; ++m)
      for (int mu = 0; mu < 4; ++mu)
        comp[static_cast<std::size_t>(m)][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 4 * m + mu);
    nm::MPoly tr = nm::numeric_value(nsym, chain, {}, comp, {});
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    Cx sym = nm::eval(tr, x, {});
    double P[8][4];
    for (int m = 0; m < K; ++m)
      for (int mu = 0; mu < 4; ++mu)
        P[m][mu] = x[static_cast<std::size_t>(4 * m + mu)];
    Cx num = buildTruth(P);
    const double err = cdiff(sym, num);
    std::printf("  %-22s monomials=%-3zu |sym|=%.2e |sym-num|=%.2e  %s\n", tag, tr.size(),
                std::abs(sym.re) + std::abs(sym.im), err, err < 1e-12 ? "ok" : "FAIL");
    if (!(err < 1e-12)) ++fails;
  };
  // F1: tr(γ5 p̸ q̸ r̸ s̸) — leading γ5, the nonzero axial (ε) structure.
  g5closed(
      "g5 p q r s",
      {network::dg5(), network::dslash({{1.0, 0}}), network::dslash({{1.0, 1}}), network::dslash({{1.0, 2}}), network::dslash({{1.0, 3}})},
      4, [&](double P[8][4]) {
        return nTrace(nMul(nMul(nMul(nMul(nGamma5(), nSlash(P[0])), nSlash(P[1])), nSlash(P[2])), nSlash(P[3])));
      });
  // F2: tr(p̸ γ5 q̸) — mid-chain γ5, two slashes → traces to 0.
  g5closed("p g5 q (=0)", {network::dslash({{1.0, 0}}), network::dg5(), network::dslash({{1.0, 1}})}, 2,
           [&](double P[8][4]) { return nTrace(nMul(nMul(nSlash(P[0]), nGamma5()), nSlash(P[1]))); });
  // F3: free legs + γ5 closed by a metric, with enough slashes to stay NONZERO:
  //     tr(γ5 γ^μ p̸ γ_μ q̸ r̸ s̸) = Σ_μ tr(γ5 γ^μ p̸ γ^μ q̸ r̸ s̸)  (= −2 tr(γ5 p̸ q̸ r̸ s̸) ≠ 0).
  {
    const int nsym = 16; // p:0-3 q:4-7 r:8-11 s:12-15
    std::vector<std::array<nm::MPoly, 4>> comp(4);
    for (int m = 0; m < 4; ++m)
      for (int mu = 0; mu < 4; ++mu)
        comp[static_cast<std::size_t>(m)][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 4 * m + mu);
    network::DiracNet chain = {network::dg5(),
                           network::dgamma(100),
                           network::dslash({{1.0, 0}}),
                           network::dgamma(101),
                           network::dslash({{1.0, 1}}),
                           network::dslash({{1.0, 2}}),
                           network::dslash({{1.0, 3}})};
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 101)}}};
    nm::MPoly tr = nm::numeric_value(nsym, chain, lor, comp, {});
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    Cx sym = nm::eval(tr, x, {});
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
    std::printf("  %-22s monomials=%-3zu |sym|=%.2e |sym-num|=%.2e  %s\n", "g5 g^mu p g_mu q r s", tr.size(),
                std::abs(sym.re) + std::abs(sym.im), err, err < 1e-12 ? "ok" : "FAIL");
    if (!(err < 1e-12)) ++fails;
  }

  // ---- G) the BARE commutator token `dcomm` vs an explicit `[X,Y]` 4×4, in all three leg modes the
  //         front-end emits: both-free (`dcomm`, σ^{μν} with two open legs), both-slash (`dcomm_ss`,
  //         the struct-7 external projector [B̸,Q̸]), and free-slash (`dcomm_fs`, a loop vertex
  //         σ^{μν}k_ν). Each is ONE chain token (no commutator split). Reference: tr(comm · p̸ q̸). ----
  std::printf("\n== G: commutator token dcomm vs explicit [X,Y] 4x4 (free/free, slash/slash, free/slash) ==\n");
  {
    const int nsym = 12; // p:0-3, q:4-7, r:8-11  (p,q close the chain; r is leg-B's slash momentum)
    std::vector<std::array<nm::MPoly, 4>> comp(3);
    for (int mu = 0; mu < 4; ++mu)
      for (int m = 0; m < 3; ++m)
        comp[static_cast<std::size_t>(m)][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 4 * m + mu);
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
      network::DiracNet chain = {network::dcomm(100, 101), network::dslash({{1.0, 0}}), network::dslash({{1.0, 1}})};
      auto T = nm::numeric_dirac(nsym, chain, comp);
      double err = 0, maxabs = 0;
      for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
          Cx sym = nm::eval(T.v[static_cast<std::size_t>(a * 4 + b)], x, {});
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
      network::DiracNet chain = {network::dcomm_ss({{1.0, 2}}, {{1.0, 0}}), network::dslash({{1.0, 0}}), network::dslash({{1.0, 1}})};
      auto T = nm::numeric_dirac(nsym, chain, comp);
      Cx sym = nm::eval(T.v[0], x, {});
      Cx num = nTrace(nMul(nMul(nCommM(nSlash(r), nSlash(p)), nSlash(p)), nSlash(q)));
      double err = cdiff(sym, num);
      std::printf("  G2 slash/slash tr([r,p] p q): |num|=%.2e |sym-num|=%.2e %s\n", std::abs(num.re) + std::abs(num.im),
                  err, err < 1e-12 ? "ok" : "FAIL");
      if (!(err < 1e-12)) ++fails;
    }
    // G3) free-slash [γ^a, r̸]: rank-1 tensor T^{a} = tr([γ^a,r̸] p̸ q̸), 4 entries.
    {
      network::DiracNet chain = {network::dcomm_fs(100, {{1.0, 2}}), network::dslash({{1.0, 0}}), network::dslash({{1.0, 1}})};
      auto T = nm::numeric_dirac(nsym, chain, comp);
      double err = 0, maxabs = 0;
      for (int a = 0; a < 4; ++a) {
        Cx sym = nm::eval(T.v[static_cast<std::size_t>(a)], x, {});
        Cx num = nTrace(nMul(nMul(nCommM(nGamma(a), nSlash(r)), nSlash(p)), nSlash(q)));
        maxabs = std::max(maxabs, std::abs(num.re) + std::abs(num.im));
        err = std::max(err, cdiff(sym, num));
      }
      std::printf("  G3 free/slash tr([g^a,r] p q): 4 entries |max|=%.2e max|sym-num|=%.2e %s\n", maxabs, err,
                  err < 1e-12 ? "ok" : "FAIL");
      if (!(err < 1e-12)) ++fails;
    }
  }

  std::printf(fails == 0 ? "\nALL TESTS PASSED\n" : "\nTESTS FAILED\n");
  return fails ? 1 : 0;
}
