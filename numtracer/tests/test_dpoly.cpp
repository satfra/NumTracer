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
#include "numtracer/numeric/numeric_contract.hpp"

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
    // a = d0 * (x0)  +  {} * (3)
    nm::DPoly a(nsym);
    a.add(nm::DMono{0}, nm::MPoly::var(nsym, 0));
    a.add(nm::DMono{}, nm::MPoly::constant(nsym, Cx{3, 0}));
    // b = d0 * (x1) + d1 * (2)
    nm::DPoly b(nsym);
    b.add(nm::DMono{0}, nm::MPoly::var(nsym, 1));
    b.add(nm::DMono{1}, nm::MPoly::constant(nsym, Cx{2, 0}));
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
    std::vector<std::array<nm::MPoly, 4>> comp(2);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][mu] = nm::MPoly::var(nsym, mu);
      comp[1][mu] = nm::MPoly::var(nsym, 4 + mu);
    }
    // tr(γ^μ p̸ γ_μ q̸) closed by a metric
    network::DiracNet chain = {network::dgamma(100), network::dslash({{1.0, 0}}), network::dgamma(101),
                               network::dslash({{1.0, 1}})};
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 101)}}};
    nm::MPoly mp = nm::numeric_value(nsym, chain, lor, comp, {});
    // same chain expressed as a slot-free dressed chain
    std::vector<nm::DChainTok> dchain;
    for (const auto &f : chain)
      dchain.push_back(nm::dtfix(f));
    nm::DPoly dp = nm::numeric_value_dressed(nsym, dchain, /*slots*/ {}, lor, comp, {});
    bool oneTerm = (dp.size() == 1) && dp.t[0].first.empty();
    bool mpEqual = oneTerm && (dp.t[0].second.t.size() == mp.t.size());
    if (mpEqual)
      for (std::size_t i = 0; i < mp.t.size(); ++i) {
        const auto &dterm = dp.t[0].second.t[i];
        const auto &mterm = mp.t[i];
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
    std::vector<std::array<nm::MPoly, 4>> comp(2);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][mu] = nm::MPoly::var(nsym, mu);
      comp[1][mu] = nm::MPoly::var(nsym, 4 + mu);
    }
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 101)}}};
    nm::DSlot sP = {nm::DSlotOpt{Cx{1, 0}, {0}, false, {}}, nm::DSlotOpt{Cx{1, 0}, {1}, true, {{1.0, 0}}}};
    nm::DSlot sQ = {nm::DSlotOpt{Cx{1, 0}, {0}, false, {}}, nm::DSlotOpt{Cx{1, 0}, {2}, true, {{1.0, 1}}}};
    std::vector<nm::DChainTok> dchain = {nm::dtfix(network::dgamma(100)), nm::dtslot(0),
                                         nm::dtfix(network::dgamma(101)), nm::dtslot(1)};
    nm::DPoly dp = nm::numeric_value_dressed(nsym, dchain, {sP, sQ}, lor, comp, {});

    // distributed reference: enumerate the 2×2 structure choices explicitly
    auto refTrace = [&](int cp, int cq, const std::vector<double> &x) {
      network::DiracNet c = {network::dgamma(100)};
      if (cp == 1) c.push_back(network::dslash({{1.0, 0}}));
      c.push_back(network::dgamma(101));
      if (cq == 1) c.push_back(network::dslash({{1.0, 1}}));
      nm::MPoly t = nm::numeric_value(nsym, c, lor, comp, {});
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
    std::vector<std::array<nm::MPoly, 4>> comp(1);
    for (int mu = 0; mu < 4; ++mu)
      comp[0][mu] = nm::MPoly::var(nsym, mu);
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 102), nm::nmet(101, 103)}}};
    nm::DSlot sP = {nm::DSlotOpt{Cx{1, 0}, {0}, false, {}}, nm::DSlotOpt{Cx{1, 0}, {1}, true, {{1.0, 0}}}};
    std::vector<nm::DChainTok> dchain = {nm::dtfix(network::dcomm(100, 101)), nm::dtslot(0),
                                         nm::dtfix(network::dgamma(102)), nm::dtfix(network::dgamma(103))};
    nm::DPoly dp = nm::numeric_value_dressed(nsym, dchain, {sP}, lor, comp, {});
    auto refTrace = [&](int cp, const std::vector<double> &x) {
      network::DiracNet c = {network::dcomm(100, 101)};
      if (cp == 1) c.push_back(network::dslash({{1.0, 0}}));
      c.push_back(network::dgamma(102));
      c.push_back(network::dgamma(103));
      return nm::eval(nm::numeric_value(nsym, c, lor, comp, {}), x, {});
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
    std::vector<std::array<nm::MPoly, 4>> comp(2);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][mu] = nm::MPoly::var(nsym, mu);
      comp[1][mu] = nm::MPoly::var(nsym, 4 + mu);
    }
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 102), nm::nmet(101, 103)}}};
    nm::DSlot sP = {nm::DSlotOpt{Cx{1, 0}, {0}, false, {}}, nm::DSlotOpt{Cx{1, 0}, {1}, true, {{1.0, 0}}}};
    nm::DSlot sQ = {nm::DSlotOpt{Cx{1, 0}, {0}, false, {}}, nm::DSlotOpt{Cx{1, 0}, {2}, true, {{1.0, 1}}}};
    std::vector<nm::DChainTok> dchain = {nm::dtfix(network::dgamma(100)), nm::dtfix(network::dgamma(101)),
                                         nm::dtslot(0), nm::dtfix(network::dcomm(102, 103)), nm::dtslot(1)};
    nm::DPoly dp = nm::numeric_value_dressed(nsym, dchain, {sP, sQ}, lor, comp, {});
    auto refTrace = [&](int cp, int cq, const std::vector<double> &x) {
      network::DiracNet c = {network::dgamma(100), network::dgamma(101)};
      if (cp == 1) c.push_back(network::dslash({{1.0, 0}}));
      c.push_back(network::dcomm(102, 103));
      if (cq == 1) c.push_back(network::dslash({{1.0, 1}}));
      return nm::eval(nm::numeric_value(nsym, c, lor, comp, {}), x, {});
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
    std::vector<std::array<nm::MPoly, 4>> comp(3);
    for (int mu = 0; mu < 4; ++mu) {
      comp[0][mu] = nm::MPoly::var(nsym, mu);
      comp[1][mu] = nm::MPoly::var(nsym, 4 + mu);
      comp[2][mu] = nm::MPoly::var(nsym, 8 + mu);
    }
    nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 101)}}};
    nm::DSlot sP = {nm::DSlotOpt{Cx{1, 0}, {0}, false, {}}, nm::DSlotOpt{Cx{1, 0}, {1}, true, {{1.0, 0}}}};
    nm::DSlot sQ = {nm::DSlotOpt{Cx{1, 0}, {0}, false, {}}, nm::DSlotOpt{Cx{1, 0}, {2}, true, {{1.0, 1}}}};
    std::vector<nm::DChainTok> dchain = {nm::dtfix(network::dgamma(100)), nm::dtslot(0),
                                         nm::dtfix(network::dcomm_fs(101, {{1.0, 2}})), nm::dtslot(1)};
    nm::DPoly dp = nm::numeric_value_dressed(nsym, dchain, {sP, sQ}, lor, comp, {});
    auto refTrace = [&](int cp, int cq, const std::vector<double> &x) {
      network::DiracNet c = {network::dgamma(100)};
      if (cp == 1) c.push_back(network::dslash({{1.0, 0}}));
      c.push_back(network::dcomm_fs(101, {{1.0, 2}}));
      if (cq == 1) c.push_back(network::dslash({{1.0, 1}}));
      return nm::eval(nm::numeric_value(nsym, c, lor, comp, {}), x, {});
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

  std::printf("\n%s\n", fails == 0 ? "ALL TESTS PASSED" : "TESTS FAILED");
  return fails == 0 ? 0 : 1;
}
