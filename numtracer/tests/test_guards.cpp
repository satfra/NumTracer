// NumTracer — guards against silent-corruption traps on code paths no current flow exercises.
//
// These cover latent bugs the end-to-end FORM/dense oracles cannot catch because no committed flow
// triggers them:
//   (A1) network::add collapses a same-index vector sum into one compound-vector leaf whose weights
//        are REAL (Elem::vlc is std::pair<double,int>). A complex coefficient on such a term would
//        have its imaginary part silently dropped; the collapse must refuse it loudly instead.
//   (B3) an EMPTY colour network is the identity factor (1), not a UnionFind(-1)/degenerate path.
//   (C1) dirac_value's gamma-parity count must treat a Comm (sigma) as TWO gammas and a LoopSep as
//        NONE — and, since its Wick recursion can trace neither, must refuse both loudly. The old
//        `!= Gamma5` count called each ONE gamma, which inverted the verdict in both directions:
//        {sigma,gamma,gamma} (4 gammas, nonzero trace) was reported as an odd chain and returned a
//        structural ZERO, while {sigma,gamma} (3 gammas, vanishing trace) passed the even gate into
//        trace_rec, whose pair_factor read the Comm's empty vlc and silently collapsed it. This
//        function is the cross-validation ORACLE in test_numeric_contract.cpp, so a wrong verdict
//        here would "confirm" a wrong engine result.
//   (D1) core/axplan.hpp's kMaxAxisRank=16 scratch arrays were entirely unbounds-checked, though the header
//        itself says an overflow "silently corrupts the contraction". Checking the OPERAND ranks is
//        not enough: the result rank is nFreeA+nFreeB, so two in-range rank-10 operands sharing no
//        axis give RR=20>kMaxAxisRank and overrun rid/rdim. This is the dense validation ORACLE's planner, so
//        a silent corruption here confirms a wrong engine result.
//   (D2) MPoly::atom narrowed the atom id to int16 with no check; aid>=32768 wrapped to a different
//        (or negative) id, so the term carried somebody else's denominator into the cancellation.
//   (D3) collect_atom_denoms overwrote atomDen[id] unconditionally, so two projectors sharing an id
//        with different momenta silently gave the second one's k² to both.
//   (D4) snap_coeff passed NaN/Inf through. RBuilder::ieq compares constants with k==k, so a NaN
//        never dedups (unbounded SSA slots) and gets printed as a literal `nan` into a committed
//        kernel header.
//   (E1) an OPEN (once-occurring) index made the engine return a plausible wrong number in silence.
//        A finished network is a scalar, but neither contraction core checked that the network it was
//        handed actually closes:
//          * Lorentz: `eliminate` sums EVERY live id over 0..3, so a lone index was contracted against
//            (1,1,1,1). `p·P(l)` with ν left open returned Σ_ν (p·P)_ν = 0.2459… — no diagnostic. Only
//            a free DIRAC leg with NO Lorentz net was caught (close_free_legs); with any Lorentz net
//            present the leg joined the elimination and was summed like the rest.
//          * SU(N) adjoint: dense-summed the same way — `f^{abc}f^{abd}` with c,d open returned 24
//            instead of the rank-2 `N δ^{cd}`.
//          * SU(N) δ legs: worse — union-find identifies a δ's two labels, the class is then "touched
//            only by δ" and counts as a CLOSED LOOP, so a lone δ^{ab} returned N²−1 = 8 and a lone
//            δ^{ij} returned N = 3. Only the open fundamental chain THROUGH A GENERATOR was caught.
//        Both cores now count occurrences up front and refuse count == 1. Counts ≥ 3 are deliberately
//        left alone (the front-end's NumTrace::badlabel owns those), so no valid net changes value —
//        which is what the anchors below pin down.
#include "numtracer/core/axplan.hpp"
#include "numtracer/network/dirac.hpp"
#include "numtracer/network/network.hpp"
#include "numtracer/network/sun_net.hpp"
#include "numtracer/numeric/env.hpp"
#include "numtracer/numeric/numeric_contract.hpp"
#include "numtracer/sun/sun_data.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

using namespace numtracer;

int main() {
  int fail = 0;
  auto ok = [&](const char *name, bool pass) {
    std::printf("  %s: %s\n", pass ? "ok  " : "FAIL", name);
    if (!pass) ++fail;
  };

  // ---- A1: real-coefficient vector sums collapse; complex ones throw ------------------------
  {
    using namespace numtracer::network;
    // two momentum legs on the same Lorentz index, real weights → one compound-vector leaf.
    NetVal s = add(scale(Cx{2, 0}, vec(0, 5)), scale(Cx{-1, 0}, vec(0, 7)));
    ok("real vecsum collapses to one term", s.size() == 1 && s[0].e.size() == 1);
    ok("compound vector carries both momenta", !s.empty() && s[0].e[0].vlc.size() == 2);

    bool threw = false;
    try {
      (void)add(scale(Cx{0, 1}, vec(0, 5)), scale(Cx{0, 1}, vec(0, 7)));
    } catch (const std::exception &) {
      threw = true;
    }
    ok("complex vecsum coefficient throws (no silent im-drop)", threw);
  }

  // ---- B3: the empty colour network is the identity factor ----------------------------------
  {
    using namespace numtracer::network;
    const Cx e = sun_value_cx(SUNNet{});
    ok("sun_value_cx(empty) == 1", e.re == 1.0 && e.im == 0.0);

    const SUNPoly p = sun_value_dressed(SUNNet{});
    const bool unit = (p.size() == 1 && p[0].dress.empty() && p[0].coeff.re == 1.0 && p[0].coeff.im == 0.0);
    ok("sun_value_dressed(empty) == single unit term", unit);
  }

  // ---- C1: dirac_value refuses tokens its Wick recursion cannot trace ------------------------
  {
    using namespace numtracer::network;
    auto refuses = [&](const DiracNet &ch) {
      try {
        (void)dirac_value(ch, 900000);
        return false;
      } catch (const std::exception &) {
        return true;
      }
    };
    // physically 4 gammas (nonzero trace) — the old count said "odd" and returned a structural zero
    ok("dirac_value refuses {sigma,gamma,gamma}", refuses({dcomm(0, 1), dgamma(2), dgamma(3)}));
    // physically 3 gammas (vanishing trace) — the old count said "even" and fell into trace_rec
    ok("dirac_value refuses {sigma,gamma}", refuses({dcomm(0, 1), dgamma(2)}));
    ok("dirac_value refuses a LoopSep chain", refuses({dgamma(0), dloopsep(), dgamma(1)}));

    // plain gamma/slash chains must be completely unaffected by the guard + recount
    ok("even gamma chain still traces", !dirac_value({dgamma(0), dgamma(1)}, 900000).empty());
    ok("odd gamma chain still vanishes", dirac_value({dgamma(0), dgamma(1), dgamma(2)}, 900000).empty());
    ok("slash pair still traces",
       !dirac_value({dslash({{1.0, 0}}), dslash({{1.0, 4}})}, 900000).empty());
  }

  auto throws = [](auto &&f) {
    try {
      f();
      return false;
    } catch (const std::exception &) {
      return true;
    }
  };

  // ---- D1: axplan kMaxAxisRank bounds ----------------------------------------------------------------
  // NOTE the pre-fix behaviour of the overflow cases is UB (an out-of-bounds write into a stack
  // array), not a wrong-but-defined value, so these were confirmed red by hand under
  // -fsanitize=address rather than by expecting a particular garbage result.
  {
    using namespace numtracer::core;
    // Rank-8 + rank-8 sharing 4 axes: RR = 8, comfortably in range. Must NOT throw, and must give
    // the same plan as before the guard — this is the "don't break the legitimate case" anchor.
    {
      std::array<int, 8> ida{}, adim{}, idb{}, bdim{};
      for (int i = 0; i < 8; ++i) {
        ida[i] = i;
        adim[i] = 2;
        idb[i] = (i < 4) ? i : 100 + i; // ids 0..3 shared, 4..7 free on each side
        bdim[i] = 2;
      }
      bool threw = throws([&] { (void)make_eplan(ida, adim, 8, idb, bdim, 8); });
      const EPlan p = make_eplan(ida, adim, 8, idb, bdim, 8);
      ok("axplan: rank-8 x rank-8 with 4 shared does not throw", !threw);
      ok("axplan: ...and still plans RR=8, nSh=4", p.RR == 8 && p.nSh == 4);
    }
    // Operand rank past kMaxAxisRank.
    {
      std::array<int, 17> id{}, dim{};
      for (int i = 0; i < 17; ++i) {
        id[i] = i;
        dim[i] = 2;
      }
      std::array<int, 1> id1{{999}}, dim1{{2}};
      ok("axplan: RA=17 throws", throws([&] { (void)make_eplan(id, dim, 17, id1, dim1, 1); }));
      ok("axplan: RB=17 throws", throws([&] { (void)make_eplan(id1, dim1, 1, id, dim, 17); }));
      ok("axplan: make_eaddplan R=17 throws", throws([&] { (void)make_eaddplan(id, dim, 17, id, dim); }));
    }
    // THE case operand-rank checks miss: two in-range rank-10 operands, NO shared axis, RR = 20.
    {
      std::array<int, 10> ida{}, adim{}, idb{}, bdim{};
      for (int i = 0; i < 10; ++i) {
        ida[i] = i;
        adim[i] = 2;
        idb[i] = 500 + i; // disjoint ids => pure outer product
        bdim[i] = 2;
      }
      ok("axplan: rank-10 (x) rank-10 disjoint => RR=20 throws",
         throws([&] { (void)make_eplan(ida, adim, 10, idb, bdim, 10); }));
    }
    // Shared identity with disagreeing extents was a comment, not a check.
    {
      std::array<int, 1> ida{{7}}, adim{{4}}, idb{{7}}, bdim{{3}};
      ok("axplan: shared axis with mismatched extent throws",
         throws([&] { (void)make_eplan(ida, adim, 1, idb, bdim, 1); }));
    }
  }

  // ---- D2 / D3 / D4: numeric-engine silent-corruption guards ---------------------------------
  {
    using namespace numtracer::numeric;
    LorentzEnv env(4);
    ok("MPoly::atom accepts the largest in-range id", !throws([&] { (void)env.atom(32767); }));
    ok("MPoly::atom rejects an id that would wrap int16", throws([&] { (void)env.atom(32768); }));
    ok("MPoly::atom rejects a negative id", throws([&] { (void)env.atom(-1); }));

    // D3: two ProjT on the SAME inv id. Same vid => the ordinary repeated-projector case, silent.
    // Different vid => the corruption, must throw.
    {
      using namespace numtracer::network;
      std::vector<std::array<MPoly, 4>> comp(2);
      for (int v = 0; v < 2; ++v)
        for (int mu = 0; mu < 4; ++mu) {
          std::vector<int> e(4, 0);
          e[mu] = 1;
          // vid 0 -> symbols as-is; vid 1 -> a different momentum (scaled), so k^2 differs.
          comp[static_cast<std::size_t>(v)][static_cast<std::size_t>(mu)] =
              env.mono(e, Cx{v == 0 ? 1.0 : 2.0, 0});
        }
      auto projNet = [](int vid, int inv, int a, int b) {
        Elem el{};
        el.kind = Elem::ProjT;
        el.a = a;
        el.b = b;
        el.vid = vid;
        el.inv = inv;
        el.invS = -1;
        PTerm t;
        t.coeff = Cx{1, 0};
        t.e.push_back(el);
        NetVal nv;
        nv.push_back(t);
        return nv;
      };
      std::vector<NetVal> sameVid{projNet(0, 3, 0, 1), projNet(0, 3, 1, 2)};
      std::vector<NetVal> diffVid{projNet(0, 3, 0, 1), projNet(1, 3, 1, 2)};
      ok("collect_atom_denoms: same inv id + same momentum is silent",
         !throws([&] { (void)collect_atom_denoms(4, sameVid, comp); }));
      ok("collect_atom_denoms: same inv id + DIFFERENT momentum throws",
         throws([&] { (void)collect_atom_denoms(4, diffVid, comp); }));
    }

    // D4: a non-finite coefficient must abort the lowering, not reach the emitted kernel.
    ok("snap_coeff rejects NaN", throws([] { (void)snap_coeff(std::numeric_limits<double>::quiet_NaN()); }));
    ok("snap_coeff rejects +Inf", throws([] { (void)snap_coeff(std::numeric_limits<double>::infinity()); }));
    ok("snap_coeff still passes finite values", std::abs(snap_coeff(1.5) - 1.5) < 1e-15);
    ok("snap_coeff still passes zero", snap_coeff(0.0) == 0.0);
  }

  // ---- E1a: an open LORENTZ index must throw, not be summed over 0..3 ------------------------
  {
    using namespace numtracer::numeric;
    using namespace numtracer::network;
    // The step-04 frame: p = (p0,0,0,0) along axis 0, l = (l0,l1,0,0) in the 0-1 plane.
    const int nsym = 3;
    LorentzEnv env(nsym);
    std::vector<std::array<MPoly, 4>> comp(2, {env.zero(), env.zero(), env.zero(), env.zero()});
    comp[0][0] = env.var(0);
    comp[1][0] = env.var(1);
    comp[1][1] = env.var(2);
    MPoly l2 = env.zero();
    for (int i = 0; i < 4; ++i)
      l2 = l2 + comp[1][i] * comp[1][i];
    const std::vector<MPoly> atomDen = {l2};
    const double Pm = 1.3, l0 = 0.5, l1 = 0.7;
    const std::vector<double> x = {Pm, l0, l1};
    const double l2v = l0 * l0 + l1 * l1;
    const std::vector<double> av = {1.0 / l2v};
    enum { mu, nu };

    // ANCHOR: the closed contraction p·P(l)·p = p²(1−cos²θ) is untouched by the guard.
    {
      NNet lor = {NTerm{Cx{1, 0}, {nvec(mu, {{1.0, 0}}), nprojT(mu, nu, {{1.0, 1}}, 0), nvec(nu, {{1.0, 0}})}}};
      const double got = eval(env.numeric_value(DiracNet{}, lor, comp, atomDen), x, av).re;
      const double cth = l0 / std::sqrt(l2v);
      ok("open-index guard: closed p.P(l).p still exact", std::fabs(got - Pm * Pm * (1 - cth * cth)) < 1e-12);
    }
    // THE trap: drop the second p and ν occurs once. Pre-guard this returned Σ_ν (p·P)_ν silently.
    {
      NNet lor = {NTerm{Cx{1, 0}, {nvec(mu, {{1.0, 0}}), nprojT(mu, nu, {{1.0, 1}}, 0)}}};
      ok("open Lorentz index throws (not summed over 0..3)",
         throws([&] { (void)env.numeric_value(DiracNet{}, lor, comp, atomDen); }));
    }
    // ANCHOR: a SELF-PAIRED free leg (γ^μ … γ^μ, both slots on one Dirac factor) is a legitimate
    // count-2 contraction that close_free_legs routes through contract_factors — it must still work.
    // tr(γ^μ p̸ γ_μ q̸) = −8 p·q in d = 4.
    {
      const DiracNet ch = {dgamma(mu), dslash({{1.0, 0}}), dgamma(mu), dslash({{1.0, 1}})};
      const double got = eval(env.numeric_value(ch, NNet{}, comp, atomDen), x, av).re;
      ok("open-index guard: self-paired gamma legs still contract", std::fabs(got - (-8.0 * Pm * l0)) < 1e-12);
    }
  }

  // ---- E1b: an open SU(N) leg must throw, in every sector -------------------------------------
  {
    using namespace numtracer::network;
    SUNEnv sun3(3); // Adim = 8
    enum { a, b, c, d, i, j, k };

    // ANCHORS: the closed nets keep their values.
    ok("open-leg guard: f^{abc}f^{abc} still 24", sun_value({sun3.f(a, b, c), sun3.f(a, b, c)}) == 24.0);
    ok("open-leg guard: tr(T^a T^a) still 4", std::fabs(sun_value({sun3.T(a, i, j), sun3.T(a, j, i)}) - 4.0) < 1e-12);
    // a label used TWICE WITHIN ONE factor is a legal closed loop (δ^{aa} = N²−1), not an open leg.
    ok("open-leg guard: delta^{aa} still 8", sun_value({sun3.deltaAdj(a, a)}) == 8.0);

    // THE traps. Pre-guard values in comments — every one of them silently wrong.
    ok("open adjoint delta leg throws (was 8)", throws([&] { (void)sun_value({sun3.deltaAdj(a, b)}); }));
    ok("open adjoint f legs throw (was 24)",
       throws([&] { (void)sun_value({sun3.f(a, b, c), sun3.f(a, b, d)}); }));
    ok("open adjoint generator legs throw (was 4)",
       throws([&] { (void)sun_value({sun3.T(a, i, j), sun3.T(b, j, i)}); }));
    ok("open fundamental delta leg throws (was 3)", throws([&] { (void)sun_value({sun3.deltaFund(i, j)}); }));
    // the one case that was already caught (extract_cycles) — now caught earlier, still caught.
    ok("open fundamental generator chain still throws",
       throws([&] { (void)sun_value({sun3.T(a, i, j), sun3.T(a, j, k)}); }));

    // the DRESSED fold is a second contraction core with the same union-find; guard it too.
    auto arr = [](int dim) {
      std::vector<int> v(dim);
      for (int n = 0; n < dim; ++n) v[n] = n;
      return v;
    };
    ok("dressed: closed diagAdj loop still evaluates",
       !throws([&] { (void)sun_value_dressed({sun3.diagAdj(a, b, arr(8)), sun3.deltaAdj(b, a)}); }));
    ok("dressed: open diagAdj leg throws",
       throws([&] { (void)sun_value_dressed({sun3.diagAdj(a, b, arr(8)), sun3.deltaAdj(b, c)}); }));
  }

  std::printf("\n%s (%d failure%s)\n", fail ? "TESTS FAILED" : "ALL TESTS PASSED", fail, fail == 1 ? "" : "s");
  return fail ? 1 : 0;
}
