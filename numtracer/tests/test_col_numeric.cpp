// NumTracer — the build-time numeric SU(N) colour/flavour net (network/sun_net.hpp).
//
// sun_net.hpp contracts a fully-contracted SU(N) network to a number at codegen time, so a net too
// large to instantiate as a dense tensor never reaches the compiler. This test validates the
// generalized net against known analytic group factors — covering adjoint f's, fundamental
// generator traces, fundamental deltas, the large-fundamental (OOM-escape) regime, two coexisting
// groups, an arbitrary-rank (untabulated) SU(5) group, and the backward-compatible adjoint path
// used by the committed generators.
#include "numtracer/network/sun_net.hpp"
#include "numtracer/sun/sun_data.hpp"

#include <cstdio>

using namespace numtracer;
using namespace numtracer::network;

int main() {
  int fail = 0;
  // One SUNEnv per distinct group rank; every SU(N) builder routes through its rank's env.
  SUNEnv sun2(2), sun3(3), sun5(5);
  auto rep = [&](const char *name, Cx got, Cx want) {
    const bool ok = approx(got, want);
    std::printf("  %s: %-40s = (%.6g, %.6g)\n", ok ? "ok  " : "FAIL", name, got.re, got.im);
    if (!ok) ++fail;
  };

  // ---- adjoint: numeric net == analytic ----------------------------------------------------
  std::printf("[adjoint nets — numeric vs analytic]\n");
  // f^{acd} f^{acd} = N(N^2-1): SU(3) -> 24, SU(2) -> 6.
  rep("SU(3) f.f closed (=24)", sun_value_cx({sun3.f(0, 1, 2), sun3.f(0, 1, 2)}), Cx{24, 0});
  rep("SU(2) f.f closed (=6)", sun_value_cx({sun2.f(0, 1, 2), sun2.f(0, 1, 2)}), Cx{6, 0});
  // f f f f = N^2(N^2-1): SU(2) -> 12.
  {
    SUNNet ffff = {sun2.f(0, 1, 2), sun2.f(1, 2, 3), sun2.f(3, 4, 5), sun2.f(4, 5, 0)};
    rep("SU(2) f f f f (=12)", sun_value_cx(ffff), Cx{12, 0});
  }

  // ---- fundamental generator traces -------------------------------------------------------
  std::printf("[fundamental nets — numeric vs analytic]\n");
  // f^{abc} tr(T^a T^b T^c) = i N(N^2-1)/4: SU(2) -> 1.5 i, SU(3) -> 6 i.
  {
    SUNNet fTTT2 = {sun2.f(0, 1, 2), sun2.T(0, 10, 11), sun2.T(1, 11, 12), sun2.T(2, 12, 10)};
    SUNNet fTTT3 = {sun3.f(0, 1, 2), sun3.T(0, 10, 11), sun3.T(1, 11, 12), sun3.T(2, 12, 10)};
    rep("SU(2) f T T T (=1.5 i)", sun_value_cx(fTTT2), Cx{0, 1.5});
    rep("SU(3) f T T T (=6 i)", sun_value_cx(fTTT3), Cx{0, 6});
  }
  // tr(T^a T^a) = C_F * N = (N^2-1)/2 [d^{ab} closes the adjoint sum]: SU(3) -> 4, SU(2) -> 1.5.
  {
    SUNNet trTT3 = {sun3.deltaAdj(0, 1), sun3.T(0, 10, 11), sun3.T(1, 11, 10)};
    SUNNet trTT2 = {sun2.deltaAdj(0, 1), sun2.T(0, 10, 11), sun2.T(1, 11, 10)};
    rep("SU(3) tr(T^a T^a) (=4)", sun_value_cx(trTT3), Cx{4, 0});
    rep("SU(2) tr(T^a T^a) (=1.5)", sun_value_cx(trTT2), Cx{1.5, 0});
  }
  // closed fundamental delta loop d^{ij} d^{ji} = N: SU(3) -> 3, SU(2) -> 2.
  rep("SU(3) d_fund loop (=3)", sun_value_cx({sun3.deltaFund(0, 1), sun3.deltaFund(1, 0)}), Cx{3, 0});
  rep("SU(2) d_fund loop (=2)", sun_value_cx({sun2.deltaFund(0, 1), sun2.deltaFund(1, 0)}), Cx{2, 0});

  // ---- large fundamental (OOM-escape regime): analytic + small-N et cross-check -----------
  std::printf("[large fundamental — tr(T^a T^b T^a T^b) = -(N^2-1)/(4N)]\n");
  {
    // adjoint labels 0,1 each appear twice (summed densely); fundamental chain 10..13 closes.
    SUNNet n3 = {sun3.T(0, 10, 11), sun3.T(1, 11, 12), sun3.T(0, 12, 13), sun3.T(1, 13, 10)};
    SUNNet n2 = {sun2.T(0, 10, 11), sun2.T(1, 11, 12), sun2.T(0, 12, 13), sun2.T(1, 13, 10)};
    rep("SU(3) (= -2/3)", sun_value_cx(n3), Cx{-2.0 / 3.0, 0});
    rep("SU(2) (= -0.375)", sun_value_cx(n2), Cx{-0.375, 0});
  }

  // ---- two groups in one net factorise -----------------------------------------------------
  std::printf("[two-group numeric net = product of per-group values]\n");
  {
    // SU(3) f.f (24) x SU(2) d_fund loop (2) = 48.
    SUNNet mix = {sun3.f(0, 1, 2), sun3.f(0, 1, 2), sun2.deltaFund(50, 51), sun2.deltaFund(51, 50)};
    rep("SU(3) f.f x SU(2) d_fund (=48)", sun_value_cx(mix), Cx{48, 0});
    const Cx c3 = sun_value_cx({sun3.f(0, 1, 2), sun3.f(0, 1, 2)});
    const Cx c2 = sun_value_cx({sun2.deltaFund(50, 51), sun2.deltaFund(51, 50)});
    rep("  == product of singles", sun_value_cx(mix), c3 * c2);
  }

  // ---- arbitrary rank: an untabulated SU(5) exercises the runtime generalized-Gell-Mann builder --
  std::printf("[arbitrary rank — untabulated SU(5) built on demand]\n");
  {
    // f^{abc} f^{abc} = N(N^2-1): SU(5) -> 120  (validates the runtime-built structure constants).
    SUNNet ff5 = {sun5.f(0, 1, 2), sun5.f(0, 1, 2)};
    rep("SU(5) f.f closed (=120)", sun_value_cx(ff5), Cx{120, 0});
    // tr(T^a T^a) = (N^2-1)/2: SU(5) -> 12  (validates the runtime-built generators + their traces).
    SUNNet trTT5 = {sun5.deltaAdj(0, 1), sun5.T(0, 10, 11), sun5.T(1, 11, 10)};
    rep("SU(5) tr(T^a T^a) (=12)", sun_value_cx(trTT5), Cx{12, 0});
    // closed fundamental delta loop d^{ij} d^{ji} = N: SU(5) -> 5.
    rep("SU(5) d_fund loop (=5)", sun_value_cx({sun5.deltaFund(0, 1), sun5.deltaFund(1, 0)}), Cx{5, 0});
  }

  // ---- backward compat: cf/cdelta (g defaulted to 3) ------
  std::printf("[backward-compat — committed adjoint nets via sun3.f()/sun3.deltaAdj()]\n");
  {
    // The 5 nonzero SU(3) colour nets emitted by gen_qcd_za_inv.cpp:60. These must keep their
    // values across the refactor (byte-identical kernel regeneration is the companion gate).
    const SUNNet za[5] = {
        {sun3.deltaAdj(0, 1), sun3.deltaAdj(2, 3), sun3.deltaAdj(4, 5), sun3.deltaAdj(6, 3), sun3.deltaAdj(6, 7), sun3.f(0, 5, 2), sun3.f(1, 7, 4)},
        {sun3.deltaAdj(0, 1), sun3.deltaAdj(2, 3), sun3.deltaAdj(4, 3), sun3.deltaAdj(4, 5), sun3.f(6, 1, 0), sun3.f(6, 5, 2)},
        {sun3.deltaAdj(0, 1), sun3.deltaAdj(2, 3), sun3.deltaAdj(4, 3), sun3.deltaAdj(4, 5), sun3.f(6, 0, 2), sun3.f(6, 1, 5)},
        {sun3.deltaAdj(0, 1), sun3.deltaAdj(2, 3), sun3.deltaAdj(4, 3), sun3.deltaAdj(4, 5), sun3.f(6, 0, 5), sun3.f(6, 1, 2)},
        {sun3.deltaAdj(0, 1), sun3.deltaAdj(2, 3), sun3.deltaAdj(4, 5), sun3.deltaAdj(6, 7), sun3.deltaAdj(7, 2), sun3.f(0, 3, 4), sun3.f(1, 5, 6)}};
    const double want[5] = {-24, 0, 24, 24, -24}; // SU(3) ZA colour factors (the refactor must preserve these)
    for (int i = 0; i < 5; ++i)
      rep("ZA colnet", Cx{sun_value_cx(za[i]).re, 0}, Cx{want[i], 0});
  }

  // ---- group-diagonal dressings: the fold returns a SUNPoly Σ_a c_a D_a, not a single number ----
  std::printf("[group-diagonal dressings — sun_value_dressed -> SUNPoly]\n");
  {
    // A diag factor carries comp2dr: component v -> a dressing-id (-1 = drop). Here every component
    // gets its own id (0..dim-1), so a leaf id names its component; D maps an id -> real value.
    auto arr = [](int dim) { std::vector<int> v(dim); for (int i = 0; i < dim; ++i) v[i] = i; return v; };
    // evaluate a SUNPoly at a dressing-id assignment D(id) -> real value.
    auto evalPoly = [](const SUNPoly &p, auto &&D) {
      Cx s{0, 0};
      for (const auto &t : p) {
        Cx c = t.coeff;
        for (int id : t.dress) c = c * Cx{D(id), 0.0};
        s = s + c;
      }
      return s;
    };
    auto ones = [](int) { return 1.0; };

    // (a) all D_a == 1 reproduces the undressed delta — the δ-reduction invariant (fund + adj).
    {
      SUNPoly fp = sun_value_dressed({sun3.diagFund(10, 11, arr(3)), sun3.deltaFund(11, 10)});
      rep("SU(3) diagFund loop, all D=1 (=3)", evalPoly(fp, ones), Cx{3, 0});
      SUNPoly ap = sun_value_dressed({sun3.diagAdj(0, 1, arr(8)), sun3.deltaAdj(1, 0)});
      rep("SU(3) diagAdj loop, all D=1 (=8)", evalPoly(ap, ones), Cx{8, 0});
    }

    // (b) fundamental δ-loop with a diagonal dressing == Σ_i D_i (weighted trace).
    {
      SUNPoly p = sun_value_dressed({sun3.diagFund(10, 11, arr(3)), sun3.deltaFund(11, 10)});
      // id i -> i+1  ->  1+2+3 = 6
      rep("SU(3) Σ_i D_i  (D_i=i+1, =6)", evalPoly(p, [](int i) { return i + 1.0; }), Cx{6, 0});
      // a closed loop of 4 diagonal dressings on ONE flavour line (the in-group u/d structure), all
      // four sharing the same per-component ids: Σ_i D_i^4 — here SU(2), D = {2,5} -> 2^4+5^4 = 641.
      SUNPoly q = sun_value_dressed({sun2.diagFund(10, 11, arr(2)), sun2.diagFund(11, 12, arr(2)),
                                     sun2.diagFund(12, 13, arr(2)), sun2.diagFund(13, 10, arr(2))});
      auto d2 = [](int i) { return i == 0 ? 2.0 : 5.0; };
      rep("SU(2) Σ_i Π_k D_i  (=641)", evalPoly(q, d2), Cx{641, 0});
      rep("  same, all D=1 (=2=Nf)", evalPoly(q, ones), Cx{2, 0});
    }

    // (c) adjoint: Σ_c Z_c f^{acd}f^{bcd}δ^{ab}|diag.  Undressed = N(N²−1)=24; per-component c_a = N = 3.
    {
      SUNPoly p = sun_value_dressed({sun3.f(0, 1, 2), sun3.f(3, 1, 2), sun3.diagAdj(0, 3, arr(8))});
      rep("SU(3) Σ_a Z_a f f, all Z=1 (=24)", evalPoly(p, ones), Cx{24, 0});
      // id a -> a+1 -> 3 * Σ_{a=0}^{7}(a+1) = 3*36 = 108
      rep("SU(3) Σ_a Z_a f f  (Z_a=a+1, =108)", evalPoly(p, [](int a) { return a + 1.0; }), Cx{108, 0});
      bool allN = !p.empty();
      for (const auto &t : p) allN = allN && approx(t.coeff, Cx{3, 0});
      std::printf("  %s: %-40s\n", allN ? "ok  " : "FAIL", "each per-colour coeff == N (=3)");
      if (!allN) ++fail;
    }

    // (d) two groups still factorise: SU(3) f.f (24) × SU(2) diagFund loop (Σ_i D_i).
    {
      SUNPoly p = sun_value_dressed({sun3.f(0, 1, 2), sun3.f(0, 1, 2),
                                     sun2.diagFund(50, 51, arr(2)), sun2.deltaFund(51, 50)});
      rep("SU(3)ff × SU(2)diag, all D=1 (=48)", evalPoly(p, ones), Cx{48, 0});
      rep("  D_i=i+1: 24*(1+2)=72", evalPoly(p, [](int i) { return i + 1.0; }), Cx{72, 0});
    }

    // (e2) DROP: a component with comp2dr = -1 contributes nothing. Fund SU(3) loop dressing only
    //      components 0 and 2 (id -1 on component 1) folds to D_0 + D_2, not D_0+D_1+D_2.
    {
      SUNPoly p = sun_value_dressed({sun3.diagFund(10, 11, {0, -1, 2}), sun3.deltaFund(11, 10)});
      rep("SU(3) diagFund drop comp 1 (D_i=i+1, =4)", evalPoly(p, [](int i) { return i + 1.0; }), Cx{4, 0});
    }

    // (e) undressed nets routed through sun_value_dressed are a single term == sun_value_cx.
    {
      SUNNet n = {sun3.f(0, 1, 2), sun3.f(0, 1, 2)};
      SUNPoly p = sun_value_dressed(n);
      const bool single = (p.size() == 1 && p[0].dress.empty());
      rep("undressed via sun_value_dressed (=24)", single ? p[0].coeff : Cx{-1, 0}, sun_value_cx(n));
    }

    // ---- diagFund ON A GENERATOR LINE (the colour-background / quark-propagator case) ----------
    // A colour-diagonal quark propagator always sits BETWEEN two T^a vertices, so its δ lands on a
    // generator cycle. loop_poly_dressed expands the trace over that index's components instead of
    // matrix-multiplying through it. Both closed forms below follow from the Fierz identity
    //     Σ_a (T^a)_{AB}(T^a)_{BA} = ½(1 − δ_AB/N)      (no sum on A,B).
    {
      // (f) self-energy topology  tr(T^a D T^a) = C_F · Σ_B D_B.  SU(3): C_F = 4/3.
      SUNNet se = {sun3.deltaAdj(0, 1), sun3.T(0, 10, 11), sun3.diagFund(11, 12, arr(3)), sun3.T(1, 12, 10)};
      SUNPoly p = sun_value_dressed(se);
      // all D=1 must collapse back to the undressed tr(T^a T^a) = C_F N = 4 — the invariant that
      // makes a zero background reproduce the colour-blind flow exactly.
      rep("SU(3) tr(T^a D T^a), all D=1 (=4)", evalPoly(p, ones), Cx{4, 0});
      // D=(1,2,3): C_F · 6 = 8
      rep("SU(3) tr(T^a D T^a)  (D=1,2,3, =8)", evalPoly(p, [](int i) { return i + 1.0; }), Cx{8, 0});

      // (g) TWO dressed propagators on one fermion loop:
      //     tr(T^a D T^a D) = ½[ (Σ D)² − (1/N)·Σ D² ].
      // The second term is a genuine same-colour correlation — it is exactly what an independent
      // per-propagator average would drop, and why the gluon-polarisation quark loop needs the
      // resolved trace rather than a colour average.
      SUNNet lp = {sun3.deltaAdj(0, 1), sun3.T(0, 10, 11), sun3.diagFund(11, 12, arr(3)),
                   sun3.T(1, 12, 13), sun3.diagFund(13, 10, arr(3))};
      SUNPoly q = sun_value_dressed(lp);
      rep("SU(3) tr(T^a D T^a D), all D=1 (=4)", evalPoly(q, ones), Cx{4, 0});
      // D=(1,2,3): ½[36 − 14/3] = 47/3
      rep("SU(3) tr(T^a D T^a D) (D=1,2,3, =47/3)", evalPoly(q, [](int i) { return i + 1.0; }),
          Cx{47.0 / 3.0, 0});

      // (g2) COMPONENT-SENSITIVE check. Everything above sums symmetrically over components, so it
      // passes even if the per-component values are PERMUTED — tr(T^a T^a) = 1/2 for every a. This
      // one reads the generator diagonal directly: tr(T^a D) with D = diag(1,2,3).
      // NumTracer builds GENERALIZED Gell-Mann generators with the diagonal (Cartan) ones LAST, so
      // 0-based gen[6] is the standard Gell-Mann T^3 = diag(1/2,-1/2,0)  -> tr(T^3 D) = -1/2
      // and    gen[7] is             Gell-Mann T^8 = diag(1,1,-2)/(2√3) -> tr(T^8 D) = -√3/2.
      // gen[0..5] are off-diagonal and give 0 against any diagonal D. FunKitAdapter's
      // ntCartanComponent relies on exactly this ordering, so pin it down here.
      {
        auto keepId = [](int c, int id) { std::vector<int> v(8, -1); v[c] = id; return v; };
        auto D = [](int id) { return id == 9 ? 1.0 : id + 1.0; }; // id 9 = the unit (pin) dressing
        auto trAD = [&](int comp) {
          return evalPoly(sun_value_dressed({sun3.T(0, 10, 11), sun3.diagFund(11, 10, arr(3)),
                                             sun3.diagAdj(0, 0, keepId(comp, 9))}), D);
        };
        for (int c = 0; c < 6; ++c) rep("off-diagonal gen vs diagonal D (=0)", trAD(c), Cx{0, 0});
        rep("gen[6] is Gell-Mann T^3: tr(T^3 D)", trAD(6), Cx{-0.5, 0});
        rep("gen[7] is Gell-Mann T^8: tr(T^8 D)", trAD(7), Cx{-0.8660254037844386, 0});
      }

      // (h) a diagFund that KEEPS ONE component is a colour projector P_c on the quark line.
      //     tr(T^a P_c T^a) = Σ_A ½(1 − δ_Ac/N) = ½(N − 1/N) = C_F, for every choice of c.
      for (int c = 0; c < 3; ++c) {
        std::vector<int> keep(3, -1);
        keep[c] = 0; // component c carries dressing id 0, the others are dropped
        SUNPoly pc = sun_value_dressed({sun3.deltaAdj(0, 1), sun3.T(0, 10, 11),
                                        sun3.diagFund(11, 12, keep), sun3.T(1, 12, 10)});
        rep("SU(3) tr(T^a P_c T^a) = C_F (=4/3)", evalPoly(pc, ones), Cx{4.0 / 3.0, 0});
      }
    }
  }

  std::printf("\n%s (%d failure%s)\n", fail ? "TESTS FAILED" : "ALL TESTS PASSED", fail, fail == 1 ? "" : "s");
  return fail ? 1 : 0;
}
