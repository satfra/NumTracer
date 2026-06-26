// ZA-bench — full-QCD (AAqbq) gluon self-energy ZA: one NumTracer variant vs FormTracer.
//
// This is the flow whose NumTracer "Generate" output blew up (375 KB + 212 KB, 273 trace
// functions) versus the FormTracer analogue (~19 KB, one collected kernel).
//
// The TU is parameterized so run_bench.sh compiles it once per NumTracer variant
// (nt_before = pre-rework per-diagram path, nt = reworked global-collection path):
//   -DNT_HDR='"nt/kernel.hh"'  -DNT_CLASS=ZA_nt_kernel
// Each variant ships its own `za_qcd*` support namespace, so including FormTracer (which has no
// support namespace) + exactly ONE NumTracer variant never collides.
//
// It reports NumTracer-vs-FormTracer agreement (informational: the QCD_NumTracer derivation is a
// work in progress, so a residual mismatch here is pre-existing, not introduced by the rework)
// and ns/eval, and dumps NumTracer integrand values on a deterministic grid to argv[1]. The
// rework's faithfulness gate is `diff` of the nt vs nt_before dumps (must agree to ~1e-12): the
// global-collection rework is an output-preserving algebraic regrouping.
//
// Nf=2 is baked into the NumTracer kernel (SetNf[2] -> no Nf parameter); the FormTracer oracle
// keeps Nf a runtime double, so we call it with Nf=2.0 to match. std::real is defensive.
#include "form/kernel.hh"
#include NT_HDR
#include "za_bench_shim.hpp"

#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <random>
#include <vector>

int main(int argc, char **argv)
{
  using Form = DiFfRG::ZA_kernel<DiFfRG::ShimRegulator>;
  using NT = DiFfRG::NT_CLASS<DiFfRG::ShimRegulator>;
  DressingSetZA d;
  const double Nf = 2.0; // matches SetNf[2] baked into the NumTracer kernel

  auto callForm = [&](double l1, double c1, double p, double k) {
    return Form::kernel(l1, c1, p, k, Nf, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4, d.ZAqbq7,
                        d.ZAAqbq1, d.ZAAqbq2, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq);
  };
  auto callNT = [&](double l1, double c1, double p, double k) {
    return std::real(NT::kernel(l1, c1, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4, d.ZAqbq7,
                                d.ZAAqbq1, d.ZAAqbq2, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
  const int Np = 200000;
  std::vector<double> L(Np), C1(Np), P(Np), K(Np);
  for (int i = 0; i < Np; ++i) { L[i] = U(rng); C1[i] = Uc(rng); P[i] = U(rng); K[i] = U(rng); }

  // NumTracer-vs-FormTracer (informational): pointwise + cos1-integrated
  double pw = 0;
  for (int i = 0; i < Np; ++i) {
    double r = callForm(L[i], C1[i], P[i], K[i]);
    pw = std::max(pw, std::fabs(callNT(L[i], C1[i], P[i], K[i]) - r) / (1e-300 + std::fabs(r)));
  }
  auto ang1 = [&](auto kfn, double l1, double p, double k) {
    const int M = 321; double s = 0, h = 2.0 / (M - 1);
    for (int i = 0; i < M; ++i) {
      double c1 = -1.0 + i * h, w = (i == 0 || i == M - 1) ? 0.5 : 1.0;
      s += w * kfn(l1, c1, p, k);
    }
    return s * h;
  };
  double ie = 0;
  for (double l1 : {0.6, 1.3, 2.4})
    for (double p : {0.7, 1.8})
      for (double k : {0.5, 1.4})
        ie = std::max(ie, std::fabs(ang1(callNT, l1, p, k) - ang1(callForm, l1, p, k)) /
                              (1e-300 + std::fabs(ang1(callForm, l1, p, k))));

  std::printf("ZA (full-QCD AAqbq) [%s] vs FormTracer:\n", NT_NAME);
  std::printf("  NT-vs-Form pointwise max rel err   = %.3e\n", pw);
  std::printf("  NT-vs-Form cos1-integrated max err = %.3e   (informational; derivation WIP)\n", ie);

  // deterministic value dump for the rework faithfulness gate (nt vs nt_before)
  if (argc > 1) {
    if (FILE *fp = std::fopen(argv[1], "w")) {
      for (double l1 = 0.2; l1 < 3.0; l1 += 0.37)
        for (double c1 = -0.9; c1 < 0.95; c1 += 0.31)
          for (double p = 0.3; p < 2.6; p += 0.5)
            for (double k = 0.4; k < 2.4; k += 0.7)
              std::fprintf(fp, "%.17g\n", callNT(l1, c1, p, k));
      std::fclose(fp);
      std::printf("  wrote NT value dump -> %s\n", argv[1]);
    }
  }

  auto bench = [&](auto fn, int n) {
    volatile double sink = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < n; ++i) sink += fn(L[i], C1[i], P[i], K[i]);
    auto t1 = std::chrono::steady_clock::now();
    (void)sink;
    return std::chrono::duration<double, std::nano>(t1 - t0).count() / n;
  };
  const double tf = bench(callForm, Np);
  const double tn = bench(callNT, Np);
  std::printf("\ntiming over %d random points:\n  Form %12.2f ns/eval\n  NT   %12.2f ns/eval  (%.2fx Form)\n",
              Np, tf, tn, tn / tf);
  std::printf("RESULT variant=%s ns_form=%.2f ns_nt=%.2f ratio=%.3f nt_vs_form_relerr=%.3e\n",
              NT_NAME, tf, tn, tn / tf, ie);
  return 0;
}
