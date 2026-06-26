# NumTracer — performance summary

Generated numeric kernels vs their FormTracer (FORM) reference, on this machine
(`-O3 -march=native`, `bench_aqbq147` over 200k random points). `num/FORM` is the ratio of
kernel evaluation time; correctness is the max relative error against FORM.

## Symbolic dressing collection (default-on, no gate)

Keeping a dressed propagator-numerator sum eager as one collected `DPoly` trace (instead of
distributing into `2^D` diagrams) is now a **pure win** — faster than the previous GlobalCollect
path on both runtime and generation, so it stays default with no gate.

Fundamental-colour quark–gluon vertex `ZAqbq{1,4,7}_147` (struct 1 / 4 / 7):

| metric                | GlobalCollect (old) | collection, unfixed | **collection (current)** |
|-----------------------|---------------------|---------------------|--------------------------|
| runtime, num/FORM     | 1.42 / 1.45 / 1.10  | 1.74 / 1.72 / 1.35  | **1.35 / 1.37 / 1.05**   |
| generation (struct 1) | 27 s                | 45 s                | **27 s**                 |
| rel-err vs FORM       | —                   | —                   | 2.5e-7 / 5.1e-9 / 3.5e-8 |

What made it a win (both fixes in `mathematica/Codegen.m`, source-only — non-dressed kernels
regenerate byte-identical):

- **Runtime:** the imaginary-vanishing probe (`numericImagProbeRealQ`) is now dressed-aware. It
  previously failed to compile against a dressed `fill()` (missing the `dr_<id>` dressing-atom
  args), so every dressed kernel was conservatively kept complex and computed-then-discarded its
  imaginary half. It now stubs and passes the dressing atoms, confirms `Im≈0`, and emits the
  lossless **RePart** double kernel (re/im split; the dead imaginary halves are DCE'd by the
  compiler) — exactly what the distributed path already got.
- **Generation:** the dressed generator's main TU compiles at `-O1` (env `NT_GEN_MAIN_OPT`) — it
  runs once (~1 s) so `-O2`'s extra cost on the `DPoly`/`dch`/`dsl` templates was wasted (28.7 s
  vs 17.5 s, identical run). The `lc/dc` net-CSE also now runs for the dressed path.

`to_genprog(DPoly)` prunes round-off **per dressing channel** (each kinematic poly against its own
max), not against one global tolerance — a channel reweighted small at generation time can dominate
at runtime once a regulator suppresses the others, so a global tolerance could silently delete its
genuine terms.

## Other flows

All committed numeric kernels match their FORM / dense / equivalence oracle — the full ctest
suite passes, including the per-flavour / per-component / disconnected flows
(`flow_flavour_split`, `flow_flavour_ingroup`, `flow_gluon_condensate`, `discdirac_num`,
`flow_sigl_dirac`). For indicative ratios of the larger gluon kernels (quark box `ZA4` ≈ 1.15×,
three-gluon `ZA3` 1/4/7 within a small factor of FORM), see the per-flow `compare_*` drivers and
`tests/refshim/`. The generation pipeline (the disposable build-time generator program) is
described in the top-level `CLAUDE.md` and `numtracer/README.md`.

## Reproducing

```bash
# regenerate a flow (needs wolframscript + FunKit, and FORM for the YM/QCD tensor bases):
wolfram -script tests/gen/gen_qcd_aqbq147_numeric.wls
# runtime + rel-err vs FORM (not a ctest — compile by hand):
g++ -O3 -march=native -std=c++20 -Iinclude -Itests/gen -Itests/refshim \
    tests/refshim/bench_aqbq147.cpp -o /tmp/b && /tmp/b
# the full numeric gate:
ctest --test-dir build --output-on-failure
```
