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
  auto rep = [&](const char *name, Cx got, Cx want) {
    const bool ok = approx(got, want);
    std::printf("  %s: %-40s = (%.6g, %.6g)\n", ok ? "ok  " : "FAIL", name, got.re, got.im);
    if (!ok) ++fail;
  };

  // ---- adjoint: numeric net == analytic ----------------------------------------------------
  std::printf("[adjoint nets — numeric vs analytic]\n");
  // f^{acd} f^{acd} = N(N^2-1): SU(3) -> 24, SU(2) -> 6.
  rep("SU(3) f.f closed (=24)", sun_value_cx({SUN::f(3, 0, 1, 2), SUN::f(3, 0, 1, 2)}), Cx{24, 0});
  rep("SU(2) f.f closed (=6)", sun_value_cx({SUN::f(2, 0, 1, 2), SUN::f(2, 0, 1, 2)}), Cx{6, 0});
  // f f f f = N^2(N^2-1): SU(2) -> 12.
  {
    SUNNet ffff = {SUN::f(2, 0, 1, 2), SUN::f(2, 1, 2, 3), SUN::f(2, 3, 4, 5), SUN::f(2, 4, 5, 0)};
    rep("SU(2) f f f f (=12)", sun_value_cx(ffff), Cx{12, 0});
  }

  // ---- fundamental generator traces -------------------------------------------------------
  std::printf("[fundamental nets — numeric vs analytic]\n");
  // f^{abc} tr(T^a T^b T^c) = i N(N^2-1)/4: SU(2) -> 1.5 i, SU(3) -> 6 i.
  {
    SUNNet fTTT2 = {SUN::f(2, 0, 1, 2), SUN::T(2, 0, 10, 11), SUN::T(2, 1, 11, 12), SUN::T(2, 2, 12, 10)};
    SUNNet fTTT3 = {SUN::f(3, 0, 1, 2), SUN::T(3, 0, 10, 11), SUN::T(3, 1, 11, 12), SUN::T(3, 2, 12, 10)};
    rep("SU(2) f T T T (=1.5 i)", sun_value_cx(fTTT2), Cx{0, 1.5});
    rep("SU(3) f T T T (=6 i)", sun_value_cx(fTTT3), Cx{0, 6});
  }
  // tr(T^a T^a) = C_F * N = (N^2-1)/2 [d^{ab} closes the adjoint sum]: SU(3) -> 4, SU(2) -> 1.5.
  {
    SUNNet trTT3 = {SUN::deltaAdj(3, 0, 1), SUN::T(3, 0, 10, 11), SUN::T(3, 1, 11, 10)};
    SUNNet trTT2 = {SUN::deltaAdj(2, 0, 1), SUN::T(2, 0, 10, 11), SUN::T(2, 1, 11, 10)};
    rep("SU(3) tr(T^a T^a) (=4)", sun_value_cx(trTT3), Cx{4, 0});
    rep("SU(2) tr(T^a T^a) (=1.5)", sun_value_cx(trTT2), Cx{1.5, 0});
  }
  // closed fundamental delta loop d^{ij} d^{ji} = N: SU(3) -> 3, SU(2) -> 2.
  rep("SU(3) d_fund loop (=3)", sun_value_cx({SUN::deltaFund(3, 0, 1), SUN::deltaFund(3, 1, 0)}), Cx{3, 0});
  rep("SU(2) d_fund loop (=2)", sun_value_cx({SUN::deltaFund(2, 0, 1), SUN::deltaFund(2, 1, 0)}), Cx{2, 0});

  // ---- large fundamental (OOM-escape regime): analytic + small-N et cross-check -----------
  std::printf("[large fundamental — tr(T^a T^b T^a T^b) = -(N^2-1)/(4N)]\n");
  {
    // adjoint labels 0,1 each appear twice (summed densely); fundamental chain 10..13 closes.
    SUNNet n3 = {SUN::T(3, 0, 10, 11), SUN::T(3, 1, 11, 12), SUN::T(3, 0, 12, 13), SUN::T(3, 1, 13, 10)};
    SUNNet n2 = {SUN::T(2, 0, 10, 11), SUN::T(2, 1, 11, 12), SUN::T(2, 0, 12, 13), SUN::T(2, 1, 13, 10)};
    rep("SU(3) (= -2/3)", sun_value_cx(n3), Cx{-2.0 / 3.0, 0});
    rep("SU(2) (= -0.375)", sun_value_cx(n2), Cx{-0.375, 0});
  }

  // ---- two groups in one net factorise -----------------------------------------------------
  std::printf("[two-group numeric net = product of per-group values]\n");
  {
    // SU(3) f.f (24) x SU(2) d_fund loop (2) = 48.
    SUNNet mix = {SUN::f(3, 0, 1, 2), SUN::f(3, 0, 1, 2), SUN::deltaFund(2, 50, 51), SUN::deltaFund(2, 51, 50)};
    rep("SU(3) f.f x SU(2) d_fund (=48)", sun_value_cx(mix), Cx{48, 0});
    const Cx c3 = sun_value_cx({SUN::f(3, 0, 1, 2), SUN::f(3, 0, 1, 2)});
    const Cx c2 = sun_value_cx({SUN::deltaFund(2, 50, 51), SUN::deltaFund(2, 51, 50)});
    rep("  == product of singles", sun_value_cx(mix), c3 * c2);
  }

  // ---- arbitrary rank: an untabulated SU(5) exercises the runtime generalized-Gell-Mann builder --
  std::printf("[arbitrary rank — untabulated SU(5) built on demand]\n");
  {
    // f^{abc} f^{abc} = N(N^2-1): SU(5) -> 120  (validates the runtime-built structure constants).
    SUNNet ff5 = {SUN::f(5, 0, 1, 2), SUN::f(5, 0, 1, 2)};
    rep("SU(5) f.f closed (=120)", sun_value_cx(ff5), Cx{120, 0});
    // tr(T^a T^a) = (N^2-1)/2: SU(5) -> 12  (validates the runtime-built generators + their traces).
    SUNNet trTT5 = {SUN::deltaAdj(5, 0, 1), SUN::T(5, 0, 10, 11), SUN::T(5, 1, 11, 10)};
    rep("SU(5) tr(T^a T^a) (=12)", sun_value_cx(trTT5), Cx{12, 0});
    // closed fundamental delta loop d^{ij} d^{ji} = N: SU(5) -> 5.
    rep("SU(5) d_fund loop (=5)", sun_value_cx({SUN::deltaFund(5, 0, 1), SUN::deltaFund(5, 1, 0)}), Cx{5, 0});
  }

  // ---- backward compat: cf/cdelta (g defaulted to 3) + the double sun_value() shim ------
  std::printf("[backward-compat — committed adjoint nets via SUN::f(3,)/SUN::deltaAdj(3,)]\n");
  {
    // The 5 nonzero SU(3) colour nets emitted by gen_qcd_za_inv.cpp:60. These must keep their
    // values across the refactor (byte-identical kernel regeneration is the companion gate).
    const SUNNet za[5] = {
        {SUN::deltaAdj(3,0, 1), SUN::deltaAdj(3,2, 3), SUN::deltaAdj(3,4, 5), SUN::deltaAdj(3,6, 3), SUN::deltaAdj(3,6, 7), SUN::f(3,0, 5, 2), SUN::f(3,1, 7, 4)},
        {SUN::deltaAdj(3,0, 1), SUN::deltaAdj(3,2, 3), SUN::deltaAdj(3,4, 3), SUN::deltaAdj(3,4, 5), SUN::f(3,6, 1, 0), SUN::f(3,6, 5, 2)},
        {SUN::deltaAdj(3,0, 1), SUN::deltaAdj(3,2, 3), SUN::deltaAdj(3,4, 3), SUN::deltaAdj(3,4, 5), SUN::f(3,6, 0, 2), SUN::f(3,6, 1, 5)},
        {SUN::deltaAdj(3,0, 1), SUN::deltaAdj(3,2, 3), SUN::deltaAdj(3,4, 3), SUN::deltaAdj(3,4, 5), SUN::f(3,6, 0, 5), SUN::f(3,6, 1, 2)},
        {SUN::deltaAdj(3,0, 1), SUN::deltaAdj(3,2, 3), SUN::deltaAdj(3,4, 5), SUN::deltaAdj(3,6, 7), SUN::deltaAdj(3,7, 2), SUN::f(3,0, 3, 4), SUN::f(3,1, 5, 6)}};
    const double want[5] = {-24, 0, 24, 24, -24}; // SU(3) ZA colour factors (the refactor must preserve these)
    for (int i = 0; i < 5; ++i)
      rep("ZA colnet", Cx{sun_value(za[i]), 0}, Cx{want[i], 0}); // double shim, as the generators use it
  }

  // ---- group-diagonal dressings: the fold returns a SUNPoly Σ_a c_a D_a, not a single number ----
  std::printf("[group-diagonal dressings — sun_value_dressed -> SUNPoly]\n");
  {
    // evaluate a SUNPoly at a per-component assignment D(dr,component) -> real value.
    auto evalPoly = [](const SUNPoly &p, auto &&D) {
      Cx s{0, 0};
      for (const auto &t : p) {
        Cx c = t.coeff;
        for (const auto &dc : t.dress) c = c * Cx{D(dc.first, dc.second), 0.0};
        s = s + c;
      }
      return s;
    };
    auto ones = [](int, int) { return 1.0; };

    // (a) all D_a == 1 reproduces the undressed delta — the δ-reduction invariant (fund + adj).
    {
      SUNPoly fp = sun_value_dressed({SUN::diagFund(3, 10, 11, 0), SUN::deltaFund(3, 11, 10)});
      rep("SU(3) diagFund loop, all D=1 (=3)", evalPoly(fp, ones), Cx{3, 0});
      SUNPoly ap = sun_value_dressed({SUN::diagAdj(3, 0, 1, 0), SUN::deltaAdj(3, 1, 0)});
      rep("SU(3) diagAdj loop, all D=1 (=8)", evalPoly(ap, ones), Cx{8, 0});
    }

    // (b) fundamental δ-loop with a diagonal dressing == Σ_i D_i (weighted trace).
    {
      SUNPoly p = sun_value_dressed({SUN::diagFund(3, 10, 11, 0), SUN::deltaFund(3, 11, 10)});
      // D_i = i+1  ->  1+2+3 = 6
      rep("SU(3) Σ_i D_i  (D_i=i+1, =6)", evalPoly(p, [](int, int i) { return i + 1.0; }), Cx{6, 0});
      // a closed loop of 4 diagonal dressings on ONE flavour line (the in-group u/d structure):
      // Σ_i D0_i D1_i D2_i D3_i — here SU(2), all four dressings equal, D_i = {2,5} -> 2^4+5^4 = 641.
      SUNPoly q = sun_value_dressed({SUN::diagFund(2, 10, 11, 0), SUN::diagFund(2, 11, 12, 1),
                                     SUN::diagFund(2, 12, 13, 2), SUN::diagFund(2, 13, 10, 3)});
      auto d2 = [](int, int i) { return i == 0 ? 2.0 : 5.0; };
      rep("SU(2) Σ_i Π_k D_i  (=641)", evalPoly(q, d2), Cx{641, 0});
      rep("  same, all D=1 (=2=Nf)", evalPoly(q, ones), Cx{2, 0});
    }

    // (c) adjoint: Σ_c Z_c f^{acd}f^{bcd}δ^{ab}|diag.  Undressed = N(N²−1)=24; per-component c_a = N = 3.
    {
      SUNPoly p = sun_value_dressed({SUN::f(3, 0, 1, 2), SUN::f(3, 3, 1, 2), SUN::diagAdj(3, 0, 3, 0)});
      rep("SU(3) Σ_a Z_a f f, all Z=1 (=24)", evalPoly(p, ones), Cx{24, 0});
      // Z_a = a+1 -> 3 * Σ_{a=0}^{7}(a+1) = 3*36 = 108
      rep("SU(3) Σ_a Z_a f f  (Z_a=a+1, =108)", evalPoly(p, [](int, int a) { return a + 1.0; }), Cx{108, 0});
      bool allN = !p.empty();
      for (const auto &t : p) allN = allN && approx(t.coeff, Cx{3, 0});
      std::printf("  %s: %-40s\n", allN ? "ok  " : "FAIL", "each per-colour coeff == N (=3)");
      if (!allN) ++fail;
    }

    // (d) two groups still factorise: SU(3) f.f (24) × SU(2) diagFund loop (Σ_i D_i).
    {
      SUNPoly p = sun_value_dressed({SUN::f(3, 0, 1, 2), SUN::f(3, 0, 1, 2),
                                     SUN::diagFund(2, 50, 51, 7), SUN::deltaFund(2, 51, 50)});
      rep("SU(3)ff × SU(2)diag, all D=1 (=48)", evalPoly(p, ones), Cx{48, 0});
      rep("  D_i=i+1: 24*(1+2)=72", evalPoly(p, [](int, int i) { return i + 1.0; }), Cx{72, 0});
    }

    // (e) undressed nets routed through sun_value_dressed are a single term == sun_value_cx.
    {
      SUNNet n = {SUN::f(3, 0, 1, 2), SUN::f(3, 0, 1, 2)};
      SUNPoly p = sun_value_dressed(n);
      const bool single = (p.size() == 1 && p[0].dress.empty());
      rep("undressed via sun_value_dressed (=24)", single ? p[0].coeff : Cx{-1, 0}, sun_value_cx(n));
    }
  }

  std::printf("\n%s (%d failure%s)\n", fail ? "TESTS FAILED" : "ALL TESTS PASSED", fail, fail == 1 ? "" : "s");
  return fail ? 1 : 0;
}
