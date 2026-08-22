// Correctness gate for the dressing-collection layer (symbolic dressing collection).
//
// Validates:
//   A) DPoly arithmetic — add / operator+ / operator* collect over dressing monomials, and eval()
//      equals the manual sum.
//   B) Dressing-free reduction — a dressed chain with NO slots produces a DPoly with a single empty
//      dressing monomial whose MPoly is bit-equal to numeric_value(...), and to_genprog(DPoly) emits
//      the SAME program (instruction stream + root) as to_genprog(MPoly). This guarantees no-regression
//      for every existing (un-dressed) flow.
//   C) Dressed contraction — a quark line with dressed numerators S(p)=Mq·δ + Z(p)·p̸ kept EAGER and
//      collected into one DPoly equals the explicit distributed sum Σ_combos (∏dressings)·trace, both
//      evaluated at random points (≤1e-10), with the trace computed ONCE instead of per combination.
#include "numtracer/numeric/env.hpp"
#include "numtracer/numeric/numeric_contract.hpp"
#include "numtracer/numeric/trace_fold.hpp" // fold_net_dressed (lever (b) assembly)

#include <cstdio>
#include <random>
#include <vector>

using numtracer::Cx;
namespace nm = numtracer::numeric;
namespace network = numtracer::network;

static double cdiff(Cx a, Cx b) { return std::abs(a.re - b.re) + std::abs(a.im - b.im); }

int main()
{
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  int fails = 0;

  // ---- A) DPoly arithmetic + eval ----
  std::printf("== A: DPoly arithmetic ==\n");
  {
    const int nsym = 2;
    nm::LorentzEnv env(nsym);
    // a = d0 * (x0)  +  {} * (3)
    nm::DPoly a = env.dzero();
    a.add(nm::DMono{0}, env.var(0));
    a.add(nm::DMono{}, env.constant(Cx{3, 0}));
    // b = d0 * (x1) + d1 * (2)
    nm::DPoly b = env.dzero();
    b.add(nm::DMono{0}, env.var(1));
    b.add(nm::DMono{1}, env.constant(Cx{2, 0}));
    nm::DPoly sum = a + b;   // d0*(x0+x1) + d1*2 + {}*3
    nm::DPoly prod = a * b;  // collects d0*d0, d0*d1, {}*d0, {}*d1
    std::vector<double> x = {0.4, -0.7};
    std::vector<double> atomVal;            // none
    std::vector<double> drVal = {1.5, -2.0}; // d0=1.5, d1=-2.0
    Cx esum = nm::eval(sum, x, atomVal, drVal);
    Cx eprod = nm::eval(prod, x, atomVal, drVal);
    // manual
    double av = 1.5 * x[0] + 3.0;
    double bv = 1.5 * x[1] + (-2.0) * 2.0;
    double msum = (1.5 * (x[0] + x[1])) + (-2.0) * 2.0 + 3.0;
    double mprod = av * bv;
    double ea = cdiff(esum, Cx{msum, 0}) + cdiff(eprod, Cx{mprod, 0});
    std::printf("  sum terms=%d prod terms=%d  |err|=%.2e  %s\n", sum.size(), prod.size(), ea,
                ea < 1e-12 ? "ok" : "FAIL");
    if (!(ea < 1e-12)) ++fails;
  }

  // ---- B) dressing-free reduction == MPoly path (no-regression guarantee) ----
  std::printf("\n== B: dressing-free DPoly == MPoly path ==\n");
  {
    const int nsym = 8; // p:0..3, q:4..7
    nm::LorentzEnv env(nsym);
    std::vector<std::array<nm::MPoly, 4>> comp(2);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][mu] = env.var(mu);
      comp[1][mu] = env.var(4 + mu);
    }
    // tr(γ^μ p̸ γ_μ q̸) closed by a metric
    network::DiracNet chain = {network::dgamma(100), network::dslash({{1.0, 0}}), network::dgamma(101),
                               network::dslash({{1.0, 1}})};
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 101)}}};
    nm::MPoly mp = env.numeric_value(chain, lor, comp, {});
    // same chain expressed as a slot-free dressed chain
    std::vector<nm::DChainTok> dchain;
    for (const auto &f : chain)
      dchain.push_back(nm::dtfix(f));
    nm::DPoly dp = env.numeric_value_dressed(dchain, /*slots*/ {}, lor, comp, {});
    bool oneTerm = (dp.size() == 1) && dp.terms[0].first.empty();
    bool mpEqual = oneTerm && (dp.terms[0].second.terms.size() == mp.terms.size());
    if (mpEqual)
      for (std::size_t i = 0; i < mp.terms.size(); ++i) {
        const auto &dterm = dp.terms[0].second.terms[i];
        const auto &mterm = mp.terms[i];
        if (!(dterm.first == mterm.first) || dterm.second.re != mterm.second.re ||
            dterm.second.im != mterm.second.im) {
          mpEqual = false;
          break;
        }
      }
    // to_genprog must emit the SAME program
    network::GlobalEnv g1, g2;
    network::GenProg pm = nm::to_genprog(mp, g1);
    network::GenProg pd = nm::to_genprog(dp, g2);
    bool progEqual = (pm.ins.size() == pd.ins.size()) && (pm.root == pd.root) && (pm.rootIm == pd.rootIm);
    std::printf("  oneTerm=%d mpEqual=%d progEqual=%d (ins %zu vs %zu)  %s\n", oneTerm, mpEqual, progEqual,
                pm.ins.size(), pd.ins.size(), (oneTerm && mpEqual && progEqual) ? "ok" : "FAIL");
    if (!(oneTerm && mpEqual && progEqual)) ++fails;
  }

  // ---- C) dressed quark line: S(p)=Mq·δ + Z(p)·p̸ collected vs distributed ----
  std::printf("\n== C: dressed numerators collected vs distributed ==\n");
  {
    // Chain: γ^μ · S(p) · γ^ν · S(q), closed by metric δ_{μν}. Each S is a dressed slot:
    //   option 0: identity δ with coeff Mq (dressing atom 0 = "Mq")
    //   option 1: slash with momentum, coeff 1, dressing atom 1 (for p) / 2 (for q) = "Z(p)"/"Z(q)"
    const int nsym = 8; // p:0..3, q:4..7
    nm::LorentzEnv env(nsym);
    std::vector<std::array<nm::MPoly, 4>> comp(2);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][mu] = env.var(mu);
      comp[1][mu] = env.var(4 + mu);
    }
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 101)}}};
    nm::DSlot sP = {nm::DSlotOpt{Cx{1, 0}, {0}, {}, {}}, nm::DSlotOpt{Cx{1, 0}, {1}, {network::dslash({{1.0, 0}})}, {}}};
    nm::DSlot sQ = {nm::DSlotOpt{Cx{1, 0}, {0}, {}, {}}, nm::DSlotOpt{Cx{1, 0}, {2}, {network::dslash({{1.0, 1}})}, {}}};
    std::vector<nm::DChainTok> dchain = {nm::dtfix(network::dgamma(100)), nm::dtslot(0),
                                         nm::dtfix(network::dgamma(101)), nm::dtslot(1)};
    nm::DPoly dp = env.numeric_value_dressed(dchain, {sP, sQ}, lor, comp, {});

    // distributed reference: enumerate the 2×2 structure choices explicitly
    auto refTrace = [&](int cp, int cq, const std::vector<double> &x) {
      network::DiracNet c = {network::dgamma(100)};
      if (cp == 1) c.push_back(network::dslash({{1.0, 0}}));
      c.push_back(network::dgamma(101));
      if (cq == 1) c.push_back(network::dslash({{1.0, 1}}));
      nm::MPoly t = env.numeric_value(c, lor, comp, {});
      return nm::eval(t, x, {});
    };

    int worst = 0;
    double maxerr = 0.0;
    for (int it = 0; it < 5000; ++it) {
      std::vector<double> x(nsym);
      for (double &v : x)
        v = U(rng);
      std::vector<double> drVal = {U(rng), U(rng), U(rng)}; // Mq, Z(p), Z(q)
      Cx collected = nm::eval(dp, x, {}, drVal);
      // distributed: Σ (dressing product) · trace(choice)
      Cx dist{0, 0};
      for (int cp = 0; cp < 2; ++cp)
        for (int cq = 0; cq < 2; ++cq) {
          double w = (cp == 0 ? drVal[0] : drVal[1]) * (cq == 0 ? drVal[0] : drVal[2]);
          Cx tr = refTrace(cp, cq, x);
          dist = dist + Cx{tr.re * w, tr.im * w};
        }
      double e = cdiff(collected, dist);
      maxerr = std::max(maxerr, e);
      if (e >= 1e-10) ++worst;
    }
    std::printf("  dp terms=%d  worst |collected-distributed|=%.2e  (%d/5000 fail)  %s\n", dp.size(), maxerr,
                worst, worst == 0 ? "ok" : "FAIL");
    if (worst != 0) ++fails;
  }

  // ---- D) σ (dcomm) vertex + dressed slot collected vs distributed (the 1/4/7 case) ----
  std::printf("\n== D: sigma (dcomm) + dressed slot collected vs distributed ==\n");
  {
    // tr( [γ^100,γ^101] · S(p) · γ^102 · γ^103 ), free legs closed by metrics 100-102, 101-103.
    const int nsym = 4; // p:0..3
    nm::LorentzEnv env(nsym);
    std::vector<std::array<nm::MPoly, 4>> comp(1);
    for (int mu = 0; mu < 4; ++mu)
      comp[0][mu] = env.var(mu);
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 102), nm::nmet(101, 103)}}};
    nm::DSlot sP = {nm::DSlotOpt{Cx{1, 0}, {0}, {}, {}}, nm::DSlotOpt{Cx{1, 0}, {1}, {network::dslash({{1.0, 0}})}, {}}};
    std::vector<nm::DChainTok> dchain = {nm::dtfix(network::dcomm(100, 101)), nm::dtslot(0),
                                         nm::dtfix(network::dgamma(102)), nm::dtfix(network::dgamma(103))};
    nm::DPoly dp = env.numeric_value_dressed(dchain, {sP}, lor, comp, {});
    auto refTrace = [&](int cp, const std::vector<double> &x) {
      network::DiracNet c = {network::dcomm(100, 101)};
      if (cp == 1) c.push_back(network::dslash({{1.0, 0}}));
      c.push_back(network::dgamma(102));
      c.push_back(network::dgamma(103));
      return nm::eval(env.numeric_value(c, lor, comp, {}), x, {});
    };
    int worst = 0;
    double maxerr = 0.0;
    for (int it = 0; it < 5000; ++it) {
      std::vector<double> x(nsym);
      for (double &v : x)
        v = U(rng);
      std::vector<double> drVal = {U(rng), U(rng)};
      Cx collected = nm::eval(dp, x, {}, drVal);
      Cx dist{0, 0};
      for (int cp = 0; cp < 2; ++cp) {
        double w = (cp == 0 ? drVal[0] : drVal[1]);
        Cx tr = refTrace(cp, x);
        dist = dist + Cx{tr.re * w, tr.im * w};
      }
      double e = cdiff(collected, dist);
      maxerr = std::max(maxerr, e);
      if (e >= 1e-10) ++worst;
    }
    std::printf("  dp terms=%d worst |collected-distributed|=%.2e (%d/5000 fail)  %s\n", dp.size(), maxerr, worst,
                worst == 0 ? "ok" : "FAIL");
    if (worst != 0) ++fails;
  }

  // ---- E) 147-like chain: multi-gamma + σ(dcomm) + TWO dressed slots, collected vs distributed ----
  std::printf("\n== E: multi-gamma + sigma + 2 dressed slots (147-like) ==\n");
  {
    // tr( γ^100 γ^101 · S(p) · [γ^102,γ^103] · S(q) ), free legs closed by metrics 100-102, 101-103.
    const int nsym = 8; // p:0..3, q:4..7
    nm::LorentzEnv env(nsym);
    std::vector<std::array<nm::MPoly, 4>> comp(2);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][mu] = env.var(mu);
      comp[1][mu] = env.var(4 + mu);
    }
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 102), nm::nmet(101, 103)}}};
    nm::DSlot sP = {nm::DSlotOpt{Cx{1, 0}, {0}, {}, {}}, nm::DSlotOpt{Cx{1, 0}, {1}, {network::dslash({{1.0, 0}})}, {}}};
    nm::DSlot sQ = {nm::DSlotOpt{Cx{1, 0}, {0}, {}, {}}, nm::DSlotOpt{Cx{1, 0}, {2}, {network::dslash({{1.0, 1}})}, {}}};
    std::vector<nm::DChainTok> dchain = {nm::dtfix(network::dgamma(100)), nm::dtfix(network::dgamma(101)),
                                         nm::dtslot(0), nm::dtfix(network::dcomm(102, 103)), nm::dtslot(1)};
    nm::DPoly dp = env.numeric_value_dressed(dchain, {sP, sQ}, lor, comp, {});
    auto refTrace = [&](int cp, int cq, const std::vector<double> &x) {
      network::DiracNet c = {network::dgamma(100), network::dgamma(101)};
      if (cp == 1) c.push_back(network::dslash({{1.0, 0}}));
      c.push_back(network::dcomm(102, 103));
      if (cq == 1) c.push_back(network::dslash({{1.0, 1}}));
      return nm::eval(env.numeric_value(c, lor, comp, {}), x, {});
    };
    int worst = 0;
    double maxerr = 0.0;
    for (int it = 0; it < 5000; ++it) {
      std::vector<double> x(nsym);
      for (double &v : x)
        v = U(rng);
      std::vector<double> drVal = {U(rng), U(rng), U(rng)};
      Cx collected = nm::eval(dp, x, {}, drVal);
      Cx dist{0, 0};
      for (int cp = 0; cp < 2; ++cp)
        for (int cq = 0; cq < 2; ++cq) {
          double w = (cp == 0 ? drVal[0] : drVal[1]) * (cq == 0 ? drVal[0] : drVal[2]);
          Cx tr = refTrace(cp, cq, x);
          dist = dist + Cx{tr.re * w, tr.im * w};
        }
      double e = cdiff(collected, dist);
      maxerr = std::max(maxerr, e);
      if (e >= 1e-10) ++worst;
    }
    std::printf("  dp terms=%d worst |collected-distributed|=%.2e (%d/5000 fail)  %s\n", dp.size(), maxerr, worst,
                worst == 0 ? "ok" : "FAIL");
    if (worst != 0) ++fails;
  }

  // ---- F) σ with a SLASHED leg (dcomm_fs) + dressed slots — the 147 loop-σ case ----
  std::printf("\n== F: sigma with slashed leg (dcomm_fs) + dressed slots ==\n");
  {
    // tr( γ^100 · S(p) · [γ^101, r̸] · S(q) ), free legs closed by metric 100-101.  r = comp[2].
    const int nsym = 12; // p:0..3, q:4..7, r:8..11
    nm::LorentzEnv env(nsym);
    std::vector<std::array<nm::MPoly, 4>> comp(3);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][mu] = env.var(mu);
      comp[1][mu] = env.var(4 + mu);
      comp[2][mu] = env.var(8 + mu);
    }
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 101)}}};
    nm::DSlot sP = {nm::DSlotOpt{Cx{1, 0}, {0}, {}, {}}, nm::DSlotOpt{Cx{1, 0}, {1}, {network::dslash({{1.0, 0}})}, {}}};
    nm::DSlot sQ = {nm::DSlotOpt{Cx{1, 0}, {0}, {}, {}}, nm::DSlotOpt{Cx{1, 0}, {2}, {network::dslash({{1.0, 1}})}, {}}};
    std::vector<nm::DChainTok> dchain = {nm::dtfix(network::dgamma(100)), nm::dtslot(0),
                                         nm::dtfix(network::dcomm_fs(101, {{1.0, 2}})), nm::dtslot(1)};
    nm::DPoly dp = env.numeric_value_dressed(dchain, {sP, sQ}, lor, comp, {});
    auto refTrace = [&](int cp, int cq, const std::vector<double> &x) {
      network::DiracNet c = {network::dgamma(100)};
      if (cp == 1) c.push_back(network::dslash({{1.0, 0}}));
      c.push_back(network::dcomm_fs(101, {{1.0, 2}}));
      if (cq == 1) c.push_back(network::dslash({{1.0, 1}}));
      return nm::eval(env.numeric_value(c, lor, comp, {}), x, {});
    };
    int worst = 0;
    double maxerr = 0.0;
    for (int it = 0; it < 5000; ++it) {
      std::vector<double> x(nsym);
      for (double &v : x)
        v = U(rng);
      std::vector<double> drVal = {U(rng), U(rng), U(rng)};
      Cx collected = nm::eval(dp, x, {}, drVal);
      Cx dist{0, 0};
      for (int cp = 0; cp < 2; ++cp)
        for (int cq = 0; cq < 2; ++cq) {
          double w = (cp == 0 ? drVal[0] : drVal[1]) * (cq == 0 ? drVal[0] : drVal[2]);
          Cx tr = refTrace(cp, cq, x);
          dist = dist + Cx{tr.re * w, tr.im * w};
        }
      double e = cdiff(collected, dist);
      maxerr = std::max(maxerr, e);
      if (e >= 1e-10) ++worst;
    }
    std::printf("  dp terms=%d worst |collected-distributed|=%.2e (%d/5000 fail)  %s\n", dp.size(), maxerr, worst,
                worst == 0 ? "ok" : "FAIL");
    if (worst != 0) ++fails;
  }

  // ---- G) OPEN gluon-leg vertex slot {γ^μ, σ^{μν}p̸_ν} collected vs distributed (Stage 4, CP1) ----
  // The quark-gluon vertex is a SUM of structures sharing ONE open gluon axis μ. Unlike a propagator
  // numerator (every option internally contracted), here every option carries the SAME open μ, closed
  // by the surrounding net — so the collected DPoly must equal the distributed structure sum.
  std::printf("\n== G: open gluon-leg vertex slot {gamma^mu, sigma^{mu nu}} collected vs distributed ==\n");
  {
    const int nsym = 4; // p:0..3
    nm::LorentzEnv env(nsym);
    std::vector<std::array<nm::MPoly, 4>> comp(1);
    for (int mu = 0; mu < 4; ++mu)
      comp[0][mu] = env.var(mu);
    // net closes the shared gluon leg μ=200 against the fixed γ^101 (playing the external gluon line).
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(200, 101)}}};
    // vertex slot: opt0 = γ^200 (T1, dressing atom 0); opt1 = σ^{200,ν}p̸_ν (T7, dressing atom 1).
    nm::DSlotOpt o0;
    o0.coeff = Cx{1, 0}; o0.dress = {0}; o0.toks = {network::dgamma(200)};                  // T1: γ^μ
    nm::DSlotOpt o1;
    o1.coeff = Cx{1, 0}; o1.dress = {1}; o1.toks = {network::dcomm_fs(200, {{1.0, 0}})};      // T7: σ^{μν}p̸_ν
    nm::DSlot sV = {o0, o1};
    std::vector<nm::DChainTok> dchain = {nm::dtslot(0), nm::dtfix(network::dgamma(101))};
    nm::DPoly dp = env.numeric_value_dressed(dchain, {sV}, lor, comp, {});
    auto refTrace = [&](int cp, const std::vector<double> &x) {
      network::DiracNet c;
      if (cp == 0) c.push_back(network::dgamma(200));            // γ^μ
      else c.push_back(network::dcomm_fs(200, {{1.0, 0}}));      // σ^{μν}p̸_ν
      c.push_back(network::dgamma(101));
      return nm::eval(env.numeric_value(c, lor, comp, {}), x, {});
    };
    int worst = 0;
    double maxerr = 0.0;
    for (int it = 0; it < 5000; ++it) {
      std::vector<double> x(nsym);
      for (double &v : x)
        v = U(rng);
      std::vector<double> drVal = {U(rng), U(rng)};
      Cx collected = nm::eval(dp, x, {}, drVal);
      Cx dist{0, 0};
      for (int cp = 0; cp < 2; ++cp) {
        double w = (cp == 0 ? drVal[0] : drVal[1]);
        Cx tr = refTrace(cp, x);
        dist = dist + Cx{tr.re * w, tr.im * w};
      }
      double e = cdiff(collected, dist);
      maxerr = std::max(maxerr, e);
      if (e >= 1e-10) ++worst;
    }
    std::printf("  dp terms=%d worst |collected-distributed|=%.2e (%d/5000 fail)  %s\n", dp.size(), maxerr, worst,
                worst == 0 ? "ok" : "FAIL");
    if (worst != 0) ++fails;
  }

  // ---- H) VecMu open-leg vector p^μ (T4) collected vs distributed (Stage 4, CP2) ----
  // The full quark-gluon vertex sum Σ c_a T_a shares ONE open gluon axis μ: T1 = γ^μ, T4 = p^μ·(δ|slash)
  // — the open index rides a Lorentz VECTOR routed into the net — and T7 = σ^{μν}p̸_ν. Collect all four
  // families in ONE slot and check it equals the explicit distributed sum. Two surrounding chains of
  // OPPOSITE Dirac parity are used so every family is nonzero in at least one (the wrong-parity ones
  // trace to zero on BOTH sides — the collection must reproduce that too).
  std::printf("\n== H: VecMu open-leg vector p^mu (T4) collected vs distributed ==\n");
  {
    const int nsym = 20; // 5 momenta (p:0..3, q:4..7, r:8..11, a:12..15, b:16..19)
    nm::LorentzEnv env(nsym);
    std::vector<std::array<nm::MPoly, 4>> comp(5);
    for (int v = 0; v < 5; ++v)
      for (int mu = 0; mu < 4; ++mu)
        comp[v][mu] = env.var(4 * v + mu);
    const std::vector<std::pair<double, int>> pMom = {{1.0, 0}}; // T4 vector p^μ (comp 0)
    const std::vector<std::pair<double, int>> qMom = {{1.0, 1}}; // slash / σ momentum   (comp 1)
    const std::vector<std::pair<double, int>> rMom = {{1.0, 2}}; // external gluon-leg partner (comp 2)
    const std::vector<std::pair<double, int>> aMom = {{1.0, 3}}; // external partner of γ^101 (comp 3)
    const std::vector<std::pair<double, int>> bMom = {{1.0, 4}}; // external partner of γ^102 (comp 4)

    // a Lorentz-net Vector factor p^μ on leg μ (an open leg routed into the net).
    auto netVec = [](int mu, std::vector<std::pair<double, int>> vlc) {
      return network::Elem{network::Elem::Vector, mu, -1, -1, -1, std::move(vlc)};
    };
    // one slot carrying all four vertex structures, sharing the open gluon axis μ = 200.
    nm::DSlotOpt oT1;  oT1.coeff = Cx{1, 0};  oT1.dress = {0}; oT1.toks = {network::dgamma(200)};        // T1: γ^μ
    nm::DSlotOpt oT4d; oT4d.coeff = Cx{1, 0}; oT4d.dress = {1}; oT4d.netFacs = {netVec(200, pMom)};       // T4: p^μ·δ (open leg on a vector)
    nm::DSlotOpt oT4s; oT4s.coeff = Cx{1, 0}; oT4s.dress = {2}; oT4s.toks = {network::dslash(qMom)};
    oT4s.netFacs = {netVec(200, pMom)};                                                                   // T4: p^μ·(γ·q̸)
    nm::DSlotOpt oT7;  oT7.coeff = Cx{1, 0};  oT7.dress = {3}; oT7.toks = {network::dcomm_fs(200, qMom)};  // T7: σ^{μν}p̸_ν
    nm::DSlot sV = {oT1, oT4d, oT4s, oT7};

    // build the concrete chain + net for ONE structure choice (distributed reference). `base` carries the
    // external gluon-leg partner r^200; T4 appends its own p^200 vector so μ closes as (r·p).
    auto optChainNet = [&](int cp, const network::DiracNet &pre, const network::DiracNet &post,
                           const nm::NNet &base) {
      network::DiracNet c = pre;
      nm::NNet net = base;
      if (cp == 0) c.push_back(network::dgamma(200));                     // T1: γ^μ (Dirac free leg 200)
      else if (cp == 1) for (auto &t : net) t.e.push_back(nm::nvec(200, pMom)); // T4-δ: p^200 into the net
      else if (cp == 2) { c.push_back(network::dslash(qMom)); for (auto &t : net) t.e.push_back(nm::nvec(200, pMom)); } // T4-slash
      else c.push_back(network::dcomm_fs(200, qMom));                     // T7: σ^{μν}p̸_ν (Dirac free leg 200)
      for (const network::DFac &d : post) c.push_back(d);
      return std::make_pair(c, net);
    };

    struct Cfg { const char *name; network::DiracNet pre, post; nm::NNet base; };
    std::vector<Cfg> cfgs = {
        // EVEN surrounding chain [γ101, ·, γ102] ⇒ T4-δ, T7 nonzero (T1, T4-slash trace to 0). The outer
        // γ's close against DISTINCT external vectors a^101, b^102 (not a shared metric) so σ does not
        // collapse via γ^αγ_α = 4 · tr(σ) = 0 — this keeps T7 genuinely nonzero.
        {"even [g101,slot,g102]", {network::dgamma(101)}, {network::dgamma(102)},
         {nm::NTerm{Cx{1, 0}, {nm::nvec(101, aMom), nm::nvec(102, bMom), nm::nvec(200, rMom)}}}},
        // ODD surrounding chain [γ101, ·] ⇒ T1, T4-slash nonzero (T4-δ, T7 trace to 0)
        {"odd  [g101,slot]", {network::dgamma(101)}, {},
         {nm::NTerm{Cx{1, 0}, {nm::nvec(101, qMom), nm::nvec(200, rMom)}}}},
    };

    for (const Cfg &cfg : cfgs) {
      std::vector<nm::DChainTok> dchain = {nm::dtfix(cfg.pre[0]), nm::dtslot(0)};
      for (const network::DFac &d : cfg.post) dchain.push_back(nm::dtfix(d));
      nm::DPoly dp = env.numeric_value_dressed(dchain, {sV}, cfg.base, comp, {});
      int worst = 0;
      double maxerr = 0.0;
      for (int it = 0; it < 5000; ++it) {
        std::vector<double> x(nsym);
        for (double &v : x)
          v = U(rng);
        std::vector<double> drVal = {U(rng), U(rng), U(rng), U(rng)};
        Cx collected = nm::eval(dp, x, {}, drVal);
        Cx dist{0, 0};
        for (int cp = 0; cp < 4; ++cp) {
          auto cn = optChainNet(cp, cfg.pre, cfg.post, cfg.base);
          Cx tr = nm::eval(env.numeric_value(cn.first, cn.second, comp, {}), x, {});
          dist = dist + Cx{tr.re * drVal[cp], tr.im * drVal[cp]};
        }
        double e = cdiff(collected, dist);
        maxerr = std::max(maxerr, e);
        if (e >= 1e-10) ++worst;
      }
      std::printf("  %-22s dp terms=%d worst |collected-distributed|=%.2e (%d/5000 fail)  %s\n", cfg.name,
                  dp.size(), maxerr, worst, worst == 0 ? "ok" : "FAIL");
      if (worst != 0) ++fails;
    }
  }

  // ---- I) MULTI-TOKEN option: open leg is NOT the only token (real T4/T7 shape) (Stage 4 R0) ----
  // The real quark-gluon vertex options are Dirac-token CHAINS, not single tokens: T4 = p̸₁·γ^μ (a slash
  // THEN the open γ), and the open leg can sit anywhere in the chain. Collect a slot whose options are
  // such chains and check it equals the distributed sum. Two fixed surroundings of opposite parity make
  // the single-token (T1) and multi-token (T4) options each nonzero in one.
  std::printf("\n== I: multi-token slot options {gamma^mu, slash.gamma^mu, slash.gamma^mu.slash} (R0) ==\n");
  {
    const int nsym = 20; // 5 momenta
    nm::LorentzEnv env(nsym);
    std::vector<std::array<nm::MPoly, 4>> comp(5);
    for (int v = 0; v < 5; ++v)
      for (int mu = 0; mu < 4; ++mu)
        comp[v][mu] = env.var(4 * v + mu);
    const std::vector<std::pair<double, int>> f1 = {{1.0, 0}}, f2 = {{1.0, 1}}, p1 = {{1.0, 2}},
                                              rr = {{1.0, 3}};
    // slot: T1 = γ^μ (1 token); T4 = p̸₁ γ^μ (2 tokens); T4' = p̸₁ γ^μ p̸₂ (3 tokens, open leg in the middle).
    nm::DSlotOpt oT1;  oT1.dress = {0}; oT1.toks = {network::dgamma(200)};
    nm::DSlotOpt oT4;  oT4.dress = {1}; oT4.toks = {network::dslash(p1), network::dgamma(200)};
    nm::DSlotOpt oT4b; oT4b.dress = {2}; oT4b.toks = {network::dslash(p1), network::dgamma(200), network::dslash(f2)};
    nm::DSlot sV = {oT1, oT4, oT4b};
    nm::NNet base = {nm::NTerm{Cx{1, 0}, {nm::nvec(200, rr)}}}; // close μ=200 against external momentum r
    struct Cfg { const char *name; network::DiracNet pre; };
    std::vector<Cfg> cfgs = {
        {"even [slash.f1, slash.f2, slot]", {network::dslash(f1), network::dslash(f2)}},
        {"odd  [slash.f1, slot]", {network::dslash(f1)}}};
    for (const Cfg &cfg : cfgs) {
      std::vector<nm::DChainTok> dchain;
      for (const network::DFac &d : cfg.pre) dchain.push_back(nm::dtfix(d));
      dchain.push_back(nm::dtslot(0));
      nm::DPoly dp = env.numeric_value_dressed(dchain, {sV}, base, comp, {});
      int worst = 0; double maxerr = 0.0;
      for (int it = 0; it < 5000; ++it) {
        std::vector<double> x(nsym);
        for (double &v : x) v = U(rng);
        std::vector<double> drVal = {U(rng), U(rng), U(rng)};
        Cx collected = nm::eval(dp, x, {}, drVal);
        Cx dist{0, 0};
        for (int cp = 0; cp < 3; ++cp) {
          network::DiracNet c = cfg.pre;
          for (const network::DFac &d : sV[cp].toks) c.push_back(d);
          Cx tr = nm::eval(env.numeric_value(c, base, comp, {}), x, {});
          dist = dist + Cx{tr.re * drVal[cp], tr.im * drVal[cp]};
        }
        double e = cdiff(collected, dist);
        maxerr = std::max(maxerr, e);
        if (e >= 1e-10) ++worst;
      }
      std::printf("  %-33s dp terms=%d worst=%.2e (%d/5000)  %s\n", cfg.name, dp.size(), maxerr, worst,
                  worst == 0 ? "ok" : "FAIL");
      if (worst != 0) ++fails;
    }
  }

  // ---- J) TWO open legs per option (k=2, an AAqbq-like vertex) (Stage 4 R1) ----
  // A two-gluon quark vertex exposes TWO open Lorentz legs {μ,ν}. Options may carry them on the Dirac
  // side (two open γ's), on the net (a metric g^{μν}·δ), or mixed (γ^μ + a net vector on ν). All are
  // collected in one slot and closed by a 2-leg net; collected must equal the distributed sum.
  std::printf("\n== J: two-open-leg slot options {gamma^mu gamma^nu, g^{mu nu}.delta, gamma^mu . vec^nu} (R1) ==\n");
  {
    const int nsym = 20;
    nm::LorentzEnv env(nsym);
    std::vector<std::array<nm::MPoly, 4>> comp(5);
    for (int v = 0; v < 5; ++v)
      for (int mu = 0; mu < 4; ++mu)
        comp[v][mu] = env.var(4 * v + mu);
    const std::vector<std::pair<double, int>> f1 = {{1.0, 0}}, f2 = {{1.0, 1}}, s = {{1.0, 2}},
                                              q = {{1.0, 3}};
    auto netVec = [](int mu, std::vector<std::pair<double, int>> vlc) {
      return network::Elem{network::Elem::Vector, mu, -1, -1, -1, std::move(vlc)};
    };
    auto netMet = [](int a, int b) { return network::Elem{network::Elem::Metric, a, b, -1, -1, {}}; };
    // local Elem→NElem for the distributed reference (public builders; elem_to_nelem is body-only).
    auto toNElem = [](const network::Elem &e) {
      return e.kind == network::Elem::Metric ? nm::nmet(e.a, e.b) : nm::nvec(e.a, e.vlc);
    };
    const int MU = 200, NU = 201;
    // options (all k=2, legs {MU,NU}): two open γ's; a net metric with δ spinor; γ^MU with a net vec on NU.
    nm::DSlotOpt oA; oA.dress = {0}; oA.toks = {network::dgamma(MU), network::dgamma(NU)};
    nm::DSlotOpt oB; oB.dress = {1}; oB.netFacs = {netMet(MU, NU)};                 // g^{μν}·δ (both legs on the net)
    nm::DSlotOpt oC; oC.dress = {2}; oC.toks = {network::dslash(s), network::dgamma(MU)};
    oC.netFacs = {netVec(NU, q)};                                                  // p̸·γ^μ · p^ν (mixed: γ leg + net-vec leg)
    nm::DSlot sV = {oA, oB, oC};
    // fixed surrounding: two slashes (even, non-collapsing) so no option zeroes by parity or loses tr(1).
    network::DiracNet pre = {network::dslash(f1), network::dslash(f2)};
    // 2-leg net: contract the two gluon legs together (g_{μν}) — a closed two-gluon sub-net.
    nm::NNet base = {nm::NTerm{Cx{1, 0}, {nm::nmet(MU, NU)}}};
    std::vector<nm::DChainTok> dchain = {nm::dtfix(pre[0]), nm::dtfix(pre[1]), nm::dtslot(0)};
    nm::DPoly dp = env.numeric_value_dressed(dchain, {sV}, base, comp, {});
    int worst = 0; double maxerr = 0.0; std::vector<double> mag(3, 0.0);
    for (int it = 0; it < 5000; ++it) {
      std::vector<double> x(nsym);
      for (double &v : x) v = U(rng);
      std::vector<double> drVal = {U(rng), U(rng), U(rng)};
      Cx collected = nm::eval(dp, x, {}, drVal);
      Cx dist{0, 0};
      for (int cp = 0; cp < 3; ++cp) {
        network::DiracNet c = pre;
        for (const network::DFac &d : sV[cp].toks) c.push_back(d);
        nm::NNet net = base;
        for (const network::Elem &e : sV[cp].netFacs) net[0].e.push_back(toNElem(e));
        Cx tr = nm::eval(env.numeric_value(c, net, comp, {}), x, {});
        mag[cp] = std::max(mag[cp], std::abs(tr.re) + std::abs(tr.im));
        dist = dist + Cx{tr.re * drVal[cp], tr.im * drVal[cp]};
      }
      double e = cdiff(collected, dist);
      maxerr = std::max(maxerr, e);
      if (e >= 1e-10) ++worst;
    }
    std::printf("  dp terms=%d worst=%.2e (%d/5000)  per-option |tr| max = {%.2f, %.2f, %.2f}  %s\n", dp.size(),
                maxerr, worst, mag[0], mag[1], mag[2], worst == 0 ? "ok" : "FAIL");
    if (worst != 0) ++fails;
  }

  // ---- J2) an option's netFacs against an EMPTY surrounding net (with_slot_facs regression) ----
  // Every other netFacs case above hands the option a net that already has a term, so the splice
  // (`append facs to every term of lor`) always had something to append to. When the Lorentz REST of
  // a sub-term is trivial the front end emits a net with ZERO terms, and appending into it is a
  // no-op — the option's factors used to be dropped silently, after which the legs they would have
  // closed abort in close_free_legs (or, on an even multiplicity, self-contract to a wrong number).
  //
  // The shape is the one that bit ZAAqbq1: the FIXED chain carries γ^μ … γ^ν, and the slot is the
  // δ_{μν}δ_Dirac + [γ_μ,γ_ν] vertex branch — so the commutator option closes μ,ν on the Dirac side
  // (paired, legal with no net at all) while the metric option can only close them through netFacs.
  // A pion/σ exchange is what makes the surrounding net empty in the first place, hence the γ5 pair.
  std::printf("\n== J2: netFacs against an EMPTY Lorentz net (metric must survive the splice) ==\n");
  {
    const int nsym = 8; // f1:0..3, f2:4..7
    nm::LorentzEnv env(nsym);
    std::vector<std::array<nm::MPoly, 4>> comp(2);
    for (int v = 0; v < 2; ++v)
      for (int mu = 0; mu < 4; ++mu)
        comp[v][mu] = env.var(4 * v + mu);
    auto met = [](int a, int b) { return network::Elem{network::Elem::Metric, a, b, -1, -1, {}}; };
    const std::vector<std::pair<double, int>> f1 = {{1.0, 0}}, f2 = {{1.0, 1}};
    const int MU = 300, NU = 301;
    // γ^MU S(f1) γ5 [slot] γ5 S(f2) γ^NU — even γ count under BOTH options, so neither branch is
    // killed by parity and both contribute.
    std::vector<nm::DChainTok> chain = {nm::dtfix(network::dgamma(MU)), nm::dtfix(network::dslash(f1)),
                                        nm::dtfix(network::dg5()),     nm::dtslot(0),
                                        nm::dtfix(network::dg5()),     nm::dtfix(network::dslash(f2)),
                                        nm::dtfix(network::dgamma(NU))};
    nm::DSlotOpt oComm; oComm.dress = {0}; oComm.toks = {network::dgamma(MU), network::dgamma(NU)};
    nm::DSlotOpt oMet;  oMet.dress = {1};  oMet.netFacs = {met(MU, NU)};   // the branch that needs the splice
    nm::DSlot sV = {oComm, oMet};
    // the whole point: NO surrounding Lorentz structure.
    const nm::NNet emptyNNet = {};
    const network::NetVal emptyNetVal = {};
    nm::DPoly dpN = env.numeric_value_dressed(chain, {sV}, emptyNNet, comp, {});
    nm::DPoly dpV = env.numeric_value_dressed_netval(chain, {sV}, emptyNetVal, comp, {});
    // distributed reference: the commutator option keeps the empty net, the metric option's factor is
    // written into a one-term net by hand — that IS the distributed diagram, built with public builders.
    // Splice IN PLACE while walking the chain: the slot sits between the two γ5, and a γ chain does not
    // commute, so appending the option's tokens at the end would be a different trace.
    network::DiracNet cComm, cMet;
    for (const nm::DChainTok &t : chain) {
      if (!t.isSlot) { cComm.push_back(t.fac); cMet.push_back(t.fac); continue; }
      for (const network::DFac &d : oComm.toks) cComm.push_back(d);
    }
    const nm::NNet refMetNet = {nm::NTerm{Cx{1, 0}, {nm::nmet(MU, NU)}}};
    nm::MPoly trComm = env.numeric_value(cComm, emptyNNet, comp, {});
    nm::MPoly trMet = env.numeric_value(cMet, refMetNet, comp, {});
    int worst = 0; double maxerr = 0.0; double magComm = 0.0, magMet = 0.0;
    for (int it = 0; it < 5000; ++it) {
      std::vector<double> x(nsym);
      for (double &v : x) v = U(rng);
      std::vector<double> drVal = {U(rng), U(rng)};
      Cx a = nm::eval(trComm, x, {}), b = nm::eval(trMet, x, {});
      magComm = std::max(magComm, std::abs(a.re) + std::abs(a.im));
      magMet = std::max(magMet, std::abs(b.re) + std::abs(b.im));
      Cx dist{a.re * drVal[0] + b.re * drVal[1], a.im * drVal[0] + b.im * drVal[1]};
      double e = std::max(cdiff(nm::eval(dpN, x, {}, drVal), dist), cdiff(nm::eval(dpV, x, {}, drVal), dist));
      maxerr = std::max(maxerr, e);
      if (e >= 1e-10) ++worst;
    }
    // A zero metric channel would make the comparison vacuous — that is exactly the dropped-factor
    // symptom this case exists to catch, so assert the channel is alive before trusting the match.
    const bool alive = magMet > 1e-8 && magComm > 1e-8;
    std::printf("  dp terms=%d/%d worst=%.2e (%d/5000)  |tr| max = {comm %.2f, metric %.2f}  %s\n",
                dpN.size(), dpV.size(), maxerr, worst, magComm, magMet,
                (worst == 0 && alive) ? "ok" : "FAIL");
    if (worst != 0 || !alive) ++fails;
  }

  // ---- K) LEVER (b): structural MPoly traces + carried dressing assemble to numeric_value_dressed ----
  // The generator no longer keys its trace table on the dressing: it strips each option's dressing
  // (coeff→1, dress→{}) into a per-sub-term scalar (`sc`, the numeric part) and monomial (`dmono`, the
  // dressing-atom ids), contracts each STRUCTURAL combination once into a plain MPoly
  // (numeric_value_dressed_netval_mp), and folds the table back into a DPoly per net (fold_net_dressed:
  // Σ_j sc[j]·T[k_j] ⊗ dmono[j]). This must reproduce the collected DPoly EXACTLY — a half-carried
  // dressing would silently drop a channel. Validate the whole decomposition against the reference
  // numeric_value_dressed_netval on the 147-like chain (multi-γ + σ + two dressed slots).
  std::printf("\n== K: lever (b) structural-trace + carried-dressing assembly == numeric_value_dressed ==\n");
  {
    const int nsym = 8; // p:0..3, q:4..7
    nm::LorentzEnv env(nsym);
    std::vector<std::array<nm::MPoly, 4>> comp(2);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][mu] = env.var(mu);
      comp[1][mu] = env.var(4 + mu);
    }
    // tr( γ^100 · S(p) · γ^101 · S(q) ), closed by metric δ_{100,101}. Both the δδ and the slash-slash
    // combos survive (even parity), so ≥2 distinct dressing channels are assembled and compared.
    // NetVal (the generator's Lorentz representation) rather than NNet, since _mp is the netval path.
    auto met = [](int a, int b) { return network::Elem{network::Elem::Metric, a, b, -1, -1, {}}; };
    network::NetVal lor = {network::PTerm{Cx{1, 0}, {met(100, 101)}}};
    // Two dressed slots, each with a NON-trivial complex option coeff so the (coeff·dress) split is
    // exercised (σ's −i lands in a slot coeff in the real flow); atoms 0/1 (p), 0/2 (q) as in case C.
    nm::DSlot sP = {nm::DSlotOpt{Cx{2, 0}, {0}, {}, {}},
                    nm::DSlotOpt{Cx{0, -1}, {1}, {network::dslash({{1.0, 0}})}, {}}};
    nm::DSlot sQ = {nm::DSlotOpt{Cx{1, 0}, {0}, {}, {}},
                    nm::DSlotOpt{Cx{3, 0}, {2}, {network::dslash({{1.0, 1}})}, {}}};
    std::vector<nm::DSlot> slots = {sP, sQ};
    std::vector<nm::DChainTok> chain = {nm::dtfix(network::dgamma(100)), nm::dtslot(0),
                                        nm::dtfix(network::dgamma(101)), nm::dtslot(1)};

    // reference: the collected DPoly.
    nm::DPoly ref = env.numeric_value_dressed_netval(chain, slots, lor, comp, {});

    // lever (b): enumerate the structural combinations, contract each to a PLAIN MPoly, and carry the
    // dressing (coeff → sc, atoms → dmono) exactly as Codegen.m does. Then fold_net_dressed reassembles.
    std::vector<nm::MPoly> T;
    std::vector<int> sidx;
    std::vector<Cx> sc;
    std::vector<nm::DMono> dmono;
    const int nSlots = static_cast<int>(slots.size());
    std::vector<int> ch(nSlots, 0);
    bool done = false;
    while (!done) {
      // structural slots (dressing stripped), and this combo's (dressCx, dressIds).
      std::vector<nm::DSlot> ss(nSlots);
      Cx dressCx{1, 0};
      nm::DMono dressIds;
      for (int s = 0; s < nSlots; ++s) {
        const nm::DSlotOpt &o = slots[s][ch[s]];
        dressCx = Cx{dressCx.re * o.coeff.re - dressCx.im * o.coeff.im,
                     dressCx.re * o.coeff.im + dressCx.im * o.coeff.re};
        dressIds.insert(dressIds.end(), o.dress.begin(), o.dress.end());
        ss[s] = {nm::DSlotOpt{Cx{1, 0}, {}, o.toks, o.netFacs}}; // one structural option per slot
      }
      nm::MPoly tr = env.numeric_value_dressed_netval_mp(chain, ss, lor, comp, {});
      sidx.push_back(static_cast<int>(T.size()));
      T.push_back(std::move(tr));
      sc.push_back(dressCx);
      dmono.push_back(nm::dmono_sorted(std::move(dressIds)));
      int k = 0;
      for (; k < nSlots; ++k) {
        if (++ch[k] < static_cast<int>(slots[k].size())) break;
        ch[k] = 0;
      }
      if (k == nSlots) done = true;
    }
    // nCache == T.size() ⇒ every trace resident, the recompute lambda is never called.
    auto never = [&](int) { return env.zero(); };
    nm::DPoly asmb =
        nm::fold_net_dressed(nsym, sidx, sc, dmono, T, static_cast<long>(T.size()), never);

    // structural comparison: same dressing channels, coefficients equal (bit-exact where the only extra
    // factor is a power of two — the tr(1)=4 collapse; a general dressCx keeps them within round-off).
    bool sameKeys = (ref.size() == asmb.size());
    double coeffErr = 0.0;
    if (sameKeys)
      for (std::size_t i = 0; i < ref.terms.size(); ++i) {
        if (ref.terms[i].first != asmb.terms[i].first) { sameKeys = false; break; }
        const auto &ra = ref.terms[i].second, &ab = asmb.terms[i].second;
        if (ra.terms.size() != ab.terms.size()) { sameKeys = false; break; }
        for (std::size_t m = 0; m < ra.terms.size(); ++m) {
          if (!(ra.terms[m].first == ab.terms[m].first)) { sameKeys = false; break; }
          coeffErr = std::max(coeffErr, cdiff(ra.terms[m].second, ab.terms[m].second));
        }
      }
    // value comparison: eval both at random points (the hard gate — never trust a byte/structural check
    // alone, per the codegen-test-gap note).
    int worst = 0;
    double maxerr = 0.0;
    for (int it = 0; it < 5000; ++it) {
      std::vector<double> x(nsym);
      for (double &v : x)
        v = U(rng);
      std::vector<double> drVal = {U(rng), U(rng), U(rng)};
      double e = cdiff(nm::eval(ref, x, {}, drVal), nm::eval(asmb, x, {}, drVal));
      maxerr = std::max(maxerr, e);
      if (e >= 1e-10) ++worst;
    }
    bool ok = sameKeys && worst == 0;
    std::printf("  traces=%zu ref terms=%d asm terms=%d sameKeys=%d coeffErr=%.2e  value worst=%.2e (%d/5000)  %s\n",
                T.size(), ref.size(), asmb.size(), sameKeys, coeffErr, maxerr, worst, ok ? "ok" : "FAIL");
    if (!ok) ++fails;
  }

  // ---- L) MonoExp HEAP-FALLBACK path: nsym > kInlineSyms (24) and degree > kMaxExp (31) ----
  // MonoExp packs exponents into two 64-bit words for nsym <= 24 / degree <= 31 and spills to a heap
  // list outside that range. Every committed flow has nsym <= 6, so the spill path had NO coverage at
  // all — yet it is what makes nsym/degree unbounded. Exercise both triggers and check the two
  // properties the packed representation must preserve: the lexicographic monomial ORDER (which is
  // what keeps emitted kernels byte-identical) and the arithmetic itself.
  std::printf("\n== L: MonoExp heap fallback (nsym=30 > 24 inline, and degree > 31) ==\n");
  {
    const int nsym = 30; // past kInlineSyms, so high symbols live on the heap
    nm::LorentzEnv env(nsym);
    std::mt19937 rg(20260723);
    std::uniform_real_distribution<double> UU(-2.0, 2.0);

    // (i) products spanning the inline/heap boundary: (x0 + x25)*(x3 + x29) must expand to the four
    //     cross terms and evaluate as the factored form does.
    nm::MPoly a = env.var(0) + env.var(25);
    nm::MPoly b = env.var(3) + env.var(29);
    nm::MPoly ab = a * b;

    // (ii) degree past kMaxExp: x7^40 — 40 > 31 forces the spill on a LOW symbol index too.
    nm::MPoly hi = env.var(7);
    for (int i = 1; i < 40; ++i)
      hi = hi * env.var(7);

    int worst = 0;
    double maxerr = 0.0;
    for (int it = 0; it < 2000; ++it) {
      std::vector<double> x(nsym);
      for (double &v : x)
        v = UU(rg);
      const Cx got = nm::eval(ab, x, {});
      const double want = (x[0] + x[25]) * (x[3] + x[29]);
      maxerr = std::max(maxerr, std::fabs(got.re - want));
      if (std::fabs(got.re - want) >= 1e-10 || std::fabs(got.im) >= 1e-12) ++worst;
    }
    // x7^40 at x7 = 1 and -1: exact, and the sign proves the exponent survived the spill intact.
    std::vector<double> one(nsym, 0.0), neg(nsym, 0.0);
    one[7] = 1.0;
    neg[7] = -1.0;
    const bool powOk = std::fabs(nm::eval(hi, one, {}).re - 1.0) < 1e-12 &&
                       std::fabs(nm::eval(hi, neg, {}).re - 1.0) < 1e-12 && hi.size() == 1;

    // (iii) ORDER: x_i must sort ascending in i across the inline→heap boundary, exactly as the old
    //       element-wise exponent-vector comparison did. Build Σ x_i and read the term order back.
    nm::MPoly sum = env.var(0);
    for (int i = 1; i < nsym; ++i)
      sum = sum + env.var(i);
    bool orderOk = sum.size() == nsym;
    for (int i = 0; orderOk && i < nsym; ++i)
      orderOk = sum.terms[static_cast<std::size_t>(i)].first.e[nsym - 1 - i] == 1; // descending exponent order

    const bool ok = worst == 0 && powOk && orderOk;
    std::printf("  ab terms=%d worst=%.2e (%d/2000)  x7^40 exact=%d  cross-boundary order=%d  %s\n", ab.size(),
                maxerr, worst, powOk, orderOk, ok ? "ok" : "FAIL");
    if (!ok) ++fails;
  }

  std::printf("\n%s\n", fails == 0 ? "ALL TESTS PASSED" : "TESTS FAILED");
  return fails == 0 ? 0 : 1;
}
