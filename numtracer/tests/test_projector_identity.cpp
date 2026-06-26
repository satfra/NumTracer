// Correctness gate for the longitudinal projector in the numeric contraction backend.
//
// The transverse projector P_T(k)_{μν} = δ_{μν} − k_μ k_ν/k² and the new longitudinal projector
// P_L(k)_{μν} = k_μ k_ν/k² are built as numeric NElem factors (nprojT / nprojL) and contracted with
// the pure-Lorentz path (empty DiracNet, every index closed by a test vector or a metric so the net
// collapses to a scalar). Each projector identity is checked two ways at ~200 random kinematic points:
//
//   1. engine-vs-oracle: the engine scalar equals an explicit double contraction of the 4×4 oracle
//      matrices PT[i][j] = (i==j) − k_i k_j/k², PL[i][j] = k_i k_j/k².
//   2. engine-only identities (independent of the oracle):
//        completeness     a·P_T·b + a·P_L·b = a·b
//        traces           tr P_T = 3,  tr P_L = 1
//        idempotency      a·P_T·P_T·b = a·P_T·b,   a·P_L·P_L·b = a·P_L·b
//        orthogonality    a·P_T·P_L·b = 0
//        transversality   k·P_T·b = 0
//        longitudinality  k·P_L·b = k·b
//        symmetry         a·P·b = b·P·a
//   The finite-T electric P_E = P_T − P_M and magnetic P_M (spatial-transverse, component 0 = heat
//   bath) projectors (nprojE / nprojM, spatial atom |k⃗|²) add: engine-vs-oracle for P_E/P_M, the
//   decomposition P_E+P_M=P_T, traces (tr P_E=1, tr P_M=2), idempotency, orthogonality P_E·P_M=0,
//   and temporal nullity u·P_M=0 with u=(1,0,0,0).
//
// k, a, b carry symbolic MPoly components; each net is contracted once and evaluated at many points.
// Build via the test CMake (adds -I include). Prints ALL TESTS PASSED / exits non-zero on failure.
#include "numtracer/numeric/numeric_contract.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using numtracer::Cx;
namespace nm = numtracer::numeric;
namespace network = numtracer::network;

// 12 symbols: k = vars 0..3 (vid 0), a = vars 4..7 (vid 1), b = vars 8..11 (vid 2).
constexpr int nsym = 12;
constexpr int kVid = 0, aVid = 1, bVid = 2;

// net builders (single product term, coefficient 1). The projector's momentum is k (vid 0),
// its full 1/k² atom is id 0; the finite-T spatial 1/|k⃗|² atom is id 1.
static nm::NElem vecLeg(int id, int vid) { return nm::nvec(id, {{1.0, vid}}); }
static nm::NElem projT(int a, int b) { return nm::nprojT(a, b, {{1.0, kVid}}, 0); }
static nm::NElem projL(int a, int b) { return nm::nprojL(a, b, {{1.0, kVid}}, 0); }
static nm::NElem projE(int a, int b) { return nm::nprojE(a, b, {{1.0, kVid}}, 0, 1); }
static nm::NElem projM(int a, int b) { return nm::nprojM(a, b, {{1.0, kVid}}, 1); }
static nm::NNet net(std::initializer_list<nm::NElem> e)
{
  return nm::NNet{nm::NTerm{Cx{1, 0}, std::vector<nm::NElem>(e)}};
}

int main()
{
  // a 4th momentum u = (1,0,0,0) (the heat-bath / temporal unit vector) to probe P_M's temporal rows.
  constexpr int uVid = 3;
  std::vector<std::array<nm::MPoly, 4>> comp(4);
  for (int mu = 0; mu < 4; ++mu) {
    comp[kVid][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 0 + mu);
    comp[aVid][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 4 + mu);
    comp[bVid][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 8 + mu);
    comp[uVid][static_cast<std::size_t>(mu)] = nm::MPoly::constant(nsym, Cx{mu == 0 ? 1.0 : 0.0, 0});
  }
  // full 1/k² atom (id 0): atomDen[0] = k² = Σ_μ comp[k][μ]²;
  // spatial 1/|k⃗|² atom (id 1): atomDen[1] = |k⃗|² = Σ_{μ=1..3} comp[k][μ]² (component 0 = temporal).
  nm::MPoly k2(nsym), ks2(nsym);
  for (int mu = 0; mu < 4; ++mu)
    k2 = k2 + comp[kVid][static_cast<std::size_t>(mu)] * comp[kVid][static_cast<std::size_t>(mu)];
  for (int mu = 1; mu < 4; ++mu)
    ks2 = ks2 + comp[kVid][static_cast<std::size_t>(mu)] * comp[kVid][static_cast<std::size_t>(mu)];
  const std::vector<nm::MPoly> atomDen = {k2, ks2};

  // contract each identity's net ONCE (symbolic in k/a/b); evaluate at many points below.
  // index convention: 10/11 are the outer (closed) legs, 12 the inner contracted leg.
  const nm::MPoly aPTb = nm::numeric_value(nsym, {}, net({vecLeg(10, aVid), projT(10, 11), vecLeg(11, bVid)}), comp, atomDen);
  const nm::MPoly aPLb = nm::numeric_value(nsym, {}, net({vecLeg(10, aVid), projL(10, 11), vecLeg(11, bVid)}), comp, atomDen);
  const nm::MPoly adotb = nm::numeric_value(nsym, {}, net({vecLeg(10, aVid), vecLeg(10, bVid)}), comp, atomDen);
  const nm::MPoly trPT = nm::numeric_value(nsym, {}, net({projT(10, 11), nm::nmet(10, 11)}), comp, atomDen);
  const nm::MPoly trPL = nm::numeric_value(nsym, {}, net({projL(10, 11), nm::nmet(10, 11)}), comp, atomDen);
  const nm::MPoly aPTPTb =
      nm::numeric_value(nsym, {}, net({vecLeg(10, aVid), projT(10, 12), projT(12, 11), vecLeg(11, bVid)}), comp, atomDen);
  const nm::MPoly aPLPLb =
      nm::numeric_value(nsym, {}, net({vecLeg(10, aVid), projL(10, 12), projL(12, 11), vecLeg(11, bVid)}), comp, atomDen);
  const nm::MPoly aPTPLb =
      nm::numeric_value(nsym, {}, net({vecLeg(10, aVid), projT(10, 12), projL(12, 11), vecLeg(11, bVid)}), comp, atomDen);
  const nm::MPoly kPTb = nm::numeric_value(nsym, {}, net({vecLeg(10, kVid), projT(10, 11), vecLeg(11, bVid)}), comp, atomDen);
  const nm::MPoly kPLb = nm::numeric_value(nsym, {}, net({vecLeg(10, kVid), projL(10, 11), vecLeg(11, bVid)}), comp, atomDen);
  const nm::MPoly kdotb = nm::numeric_value(nsym, {}, net({vecLeg(10, kVid), vecLeg(10, bVid)}), comp, atomDen);
  const nm::MPoly bPTa = nm::numeric_value(nsym, {}, net({vecLeg(10, bVid), projT(10, 11), vecLeg(11, aVid)}), comp, atomDen);

  // finite-T electric / magnetic projectors
  const nm::MPoly aPEb = nm::numeric_value(nsym, {}, net({vecLeg(10, aVid), projE(10, 11), vecLeg(11, bVid)}), comp, atomDen);
  const nm::MPoly aPMb = nm::numeric_value(nsym, {}, net({vecLeg(10, aVid), projM(10, 11), vecLeg(11, bVid)}), comp, atomDen);
  const nm::MPoly trPE = nm::numeric_value(nsym, {}, net({projE(10, 11), nm::nmet(10, 11)}), comp, atomDen);
  const nm::MPoly trPM = nm::numeric_value(nsym, {}, net({projM(10, 11), nm::nmet(10, 11)}), comp, atomDen);
  const nm::MPoly aPEPEb =
      nm::numeric_value(nsym, {}, net({vecLeg(10, aVid), projE(10, 12), projE(12, 11), vecLeg(11, bVid)}), comp, atomDen);
  const nm::MPoly aPMPMb =
      nm::numeric_value(nsym, {}, net({vecLeg(10, aVid), projM(10, 12), projM(12, 11), vecLeg(11, bVid)}), comp, atomDen);
  const nm::MPoly aPEPMb =
      nm::numeric_value(nsym, {}, net({vecLeg(10, aVid), projE(10, 12), projM(12, 11), vecLeg(11, bVid)}), comp, atomDen);
  const nm::MPoly uPMb = nm::numeric_value(nsym, {}, net({vecLeg(10, uVid), projM(10, 11), vecLeg(11, bVid)}), comp, atomDen);

  std::mt19937 rng(2024);
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  constexpr double tol = 1e-12;
  int fails = 0, npts = 200;

  // worst-case error per named check, accumulated over all points.
  struct Chk {
    const char *name;
    double worst = 0;
  };
  Chk checks[] = {
      {"engine==oracle a.PT.b"}, {"engine==oracle a.PL.b"}, {"engine==oracle tr PT"},
      {"engine==oracle tr PL"},  {"completeness PT+PL=d"},  {"trace PT=3"},
      {"trace PL=1"},            {"idempotent PT.PT=PT"},   {"idempotent PL.PL=PL"},
      {"orthogonal PT.PL=0"},    {"transverse k.PT=0"},     {"longitudinal k.PL=k.b"},
      {"symmetry a.P.b=b.P.a"},
      // finite-T electric/magnetic
      {"engine==oracle a.PE.b"}, {"engine==oracle a.PM.b"}, {"engine==oracle tr PE"},
      {"engine==oracle tr PM"},  {"decomposition PE+PM=PT"}, {"trace PE=1"},
      {"trace PM=2"},            {"idempotent PE.PE=PE"},   {"idempotent PM.PM=PM"},
      {"orthogonal PE.PM=0"},    {"temporal u.PM=0"},
  };
  auto rec = [&](int i, double err) {
    if (err > checks[static_cast<std::size_t>(i)].worst) checks[static_cast<std::size_t>(i)].worst = err;
  };
  auto re = [](Cx z) { return z.re; }; // all identity scalars are real

  for (int it = 0; it < npts; ++it) {
    std::vector<double> x(nsym);
    for (double &v : x)
      v = U(rng);
    double kk[4], aa[4], bb[4], k2v = 0, ks2v = 0;
    for (int mu = 0; mu < 4; ++mu) {
      kk[mu] = x[static_cast<std::size_t>(0 + mu)];
      aa[mu] = x[static_cast<std::size_t>(4 + mu)];
      bb[mu] = x[static_cast<std::size_t>(8 + mu)];
      k2v += kk[mu] * kk[mu];
      if (mu > 0) ks2v += kk[mu] * kk[mu]; // |k⃗|² (spatial atom)
    }
    const std::vector<double> atomVal = {1.0 / k2v, 1.0 / ks2v};
    auto E = [&](const nm::MPoly &p) { return re(nm::eval(p, x, atomVal)); };

    // oracle 4×4 projector matrices and double contractions (ks2v = |k⃗|² computed above)
    double PT[4][4], PL[4][4], PM[4][4], PE[4][4];
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        PL[i][j] = kk[i] * kk[j] / k2v;
        PT[i][j] = (i == j ? 1.0 : 0.0) - PL[i][j];
        PM[i][j] = (i == j && i > 0 ? 1.0 : 0.0) - ((i > 0 && j > 0) ? kk[i] * kk[j] / ks2v : 0.0);
        PE[i][j] = PT[i][j] - PM[i][j];
      }
    auto vMv = [](const double u[4], const double M[4][4], const double w[4]) {
      double s = 0;
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          s += u[i] * M[i][j] * w[j];
      return s;
    };
    double trPTo = 0, trPLo = 0, adotbo = 0;
    for (int i = 0; i < 4; ++i) {
      trPTo += PT[i][i];
      trPLo += PL[i][i];
      adotbo += aa[i] * bb[i];
    }

    // 1) engine vs oracle
    rec(0, std::abs(E(aPTb) - vMv(aa, PT, bb)));
    rec(1, std::abs(E(aPLb) - vMv(aa, PL, bb)));
    rec(2, std::abs(E(trPT) - trPTo));
    rec(3, std::abs(E(trPL) - trPLo));

    // 2) engine-only identities
    rec(4, std::abs(E(aPTb) + E(aPLb) - E(adotb)));               // completeness
    rec(5, std::abs(E(trPT) - 3.0));                              // tr P_T = 3
    rec(6, std::abs(E(trPL) - 1.0));                              // tr P_L = 1
    rec(7, std::abs(E(aPTPTb) - E(aPTb)));                        // P_T idempotent
    rec(8, std::abs(E(aPLPLb) - E(aPLb)));                        // P_L idempotent
    rec(9, std::abs(E(aPTPLb)));                                  // P_T · P_L = 0
    rec(10, std::abs(E(kPTb)));                                   // k transverse to P_T
    rec(11, std::abs(E(kPLb) - E(kdotb)));                        // k longitudinal: k·P_L·b = k·b
    rec(12, std::abs(E(aPTb) - E(bPTa)));                         // symmetry

    // finite-T electric/magnetic
    rec(13, std::abs(E(aPEb) - vMv(aa, PE, bb)));                 // engine vs oracle a.PE.b
    rec(14, std::abs(E(aPMb) - vMv(aa, PM, bb)));                 // engine vs oracle a.PM.b
    double trPEo = 0, trPMo = 0;
    for (int i = 0; i < 4; ++i) {
      trPEo += PE[i][i];
      trPMo += PM[i][i];
    }
    rec(15, std::abs(E(trPE) - trPEo));                           // engine vs oracle tr PE
    rec(16, std::abs(E(trPM) - trPMo));                           // engine vs oracle tr PM
    rec(17, std::abs(E(aPEb) + E(aPMb) - E(aPTb)));               // decomposition P_E+P_M=P_T
    rec(18, std::abs(E(trPE) - 1.0));                             // tr P_E = 1
    rec(19, std::abs(E(trPM) - 2.0));                             // tr P_M = 2
    rec(20, std::abs(E(aPEPEb) - E(aPEb)));                       // P_E idempotent
    rec(21, std::abs(E(aPMPMb) - E(aPMb)));                       // P_M idempotent
    rec(22, std::abs(E(aPEPMb)));                                 // P_E · P_M = 0
    rec(23, std::abs(E(uPMb)));                                   // temporal: u·P_M·b = 0
    (void)adotbo;
  }

  std::printf("== projector identities (%d random points) ==\n", npts);
  for (const Chk &c : checks) {
    const bool ok = c.worst < tol;
    std::printf("  %-26s worst=%.2e  %s\n", c.name, c.worst, ok ? "ok" : "FAIL");
    if (!ok) ++fails;
  }

  std::printf(fails ? "\nTESTS FAILED (%d)\n" : "\nALL TESTS PASSED\n", fails);
  return fails ? 1 : 0;
}
