# Trace tables — production design (draft, 2026-11-11)

Background and measurements: `NUMTRACER_TRACE_TABLES_PROBE.md`. The prototype (post-processing tool plus DiFfRG
overlay) gave 1.40x on YangMills, and its tables are accurate to float32 precision at every S0/k. This document is the
plan to make it a NumTracer and DiFfRG feature.

## Principle

The kinematics of a trace live on two kinds of variables:
- **fixed per run:** quadrature nodes in units of the box (q̂ = q/k, angles) and the external coordinates that are not
  a momentum scale (S1, SPhi);
- **changing per call:** the RG scale k (through q = k·q̂) and the dimensionful external coordinate S0.

After clearing the propagator denominators, every trace is a homogeneous polynomial in the per-call variables, with
coefficients that depend only on the fixed ones:

    T = N · ∏_j inv_j^{m_j},    N = Σ_{n=n0..n1} c_n(q̂, angles, shape) · k^n · S0^{D−n}

The table holds c_n at every fixed node. The kernel evaluates Horner in r = min(k,S0)/max(k,S0) and multiplies by
s^d · ∏_j (s^{w_j} inv_j)^{m_j}, with s = max(k, S0).

## 1. NumTracer generator (C++): exact coefficients

The generator already has everything this needs:
- each trace exactly as `acc`, a `Poly` over the frame variables and the inv atoms, in the
  `fold_groups_streaming` callback just before `to_genprog`;
- each denominator Q_j as a `Poly` over the frame variables (`frame.denominators()`).

New step, per tabulated trace (in place of, or next to, `to_genprog(acc)`):
1. **m_j** = the max power of inv_j in `acc`.
2. **N** = Σ_terms c · (frame monomial) · ∏_j Q_j^{m_j − p_j}. This is exact Poly arithmetic, and identically vanishing
   coefficients come out as exact zeros. (It is the probing-free version of what the prototype did numerically.)
3. Collect N by powers of the radial variable (and of the scale variable). Homogeneity fixes the exponent of S0, so
   only n0..n1 survive. Each c_n is a `Poly` over the fixed variables, and is checked to be homogeneous; a violation is
   a hard error.
4. Lower each c_n with the existing `to_genprog`/`best_into` into a table environment that contains only the fixed
   slots. Emit:
   - `tab_fill(double *out, <fixed args>)`: device-callable, fp64, writes all entries of a node;
   - `rtrN(row, stride, f, s, r, kle)`: the reconstruction;
   - `constexpr` metadata: `tab_size`, the per-trace entries and inline flop count, and the balance decision.

Dressed traces (`DressedPoly` = Σ_α D^α · Poly_α) follow the same steps per α, and the reconstruction sums
D^α · Horner_α. Section 4 of the probe shows this pays only for some flows. The balance switch decides per (trace, α).

## 2. Wolfram front end

- **Option on `MakeNTKernel` / `MakeNTKernelDiFfRG`:** a new entry in `Options[mkGenerateKernel]` and `ntGenOptions`,
  forwarded like `"Hoist"`. The name and shape are open (Q1).
- **Variable roles:** the radial variable is known from `"IntegrationVariables"` and the integrator type. The scale
  variable(s) are inferred by checking which subset of `"CoordinateArguments"` makes every frame vector homogeneous of
  degree 1. An explicit override is possible (Q2).
- **Emission of `kernel.hh`:**
  - `has_tab`, `tab_size` (read from `kernels.hh`) and `tab_ext_mask`;
  - `tab()`, which calls `ns::tab_fill` with q̂ = l1/k;
  - `kernel_tab()` / `post_tab()`: the existing bodies with `trN(fenv)` replaced by `rtrN(...)`. The prototype already
    does this textually; here it comes from `ntAssembleIntegrand`.

## 3. DiFfRG

From the prototype patch (`gpubench/tools/tracetab/diffrg_tab.patch`), cleaned up:
- **Kernel concept `kernel_has_tab`:** `has_tab`, `tab_size`, `tab_ext_mask`, and `tab()` as a device function (Q3).
- **`TabTable`:**
  - keyed on the coordinates and the non-radial box geometry, never on k or the interpolator epoch;
  - layout shape-fastest, which is required on GPUs (Kokkos iterates the external index fastest);
  - **filled on the device** in one `parallel_for` over (shape, loop node). That is about 2% of one RHS evaluation,
    versus 49 s on the host in the prototype.
- **`map()` branch** before the hoisted and plain paths; `map_points()` and the other integrators (1 and 3 angles,
  finite T) follow.
- **Transform wrappers** forward the interface and apply the measure.

## 4. Scope of a first version

| in | later |
|---|---|
| vacuum integrators with one radial scale k and one dimensionful external coordinate | finite T (Matsubara axis plus T as a second scale: bivariate tables) |
| dressing-free traces | dressed traces (per α), once measured worthwhile |
| `map()` | `map_points()`, `get()` |

## 5. Verification and performance budget

- **Off by default:** emitted kernels are byte-identical with the option off (`regen_check`, `git status -- tests/gen`).
- **Exactness:** a new codegen test reconstructs every tabulated trace from its coefficients at random nodes and
  compares with the interpreted inline program (`codegen/interpret.hpp`) in fp64. Tolerance about 1e-12 of the
  integrand scale.
- **Precision:** the double-literal reference of the probe (ZA3, Zcbcbcc) has to be reproduced: table ≤ 2e-6 at every
  S0/k.
- **Speed (YangMills `full_AAcbc_sigma_split`, 4070):**
  - ≥ 1.35x per RHS (the prototype reached 1.40x);
  - the fill must cost < 5% of one RHS;
  - remeasure on A100/H100 before making the option a default anywhere.
- **Compile time:** the coefficient code is about the size of the traces. With tables on, the inline `trN` can be
  dropped (Q4), so device code should not grow.

## 6. Decisions (2026-11-11)

| question | decision |
|---|---|
| Q1 option | `"TraceTables" -> False \| True \| Automatic \| <\|"Balance" -> B, "Scales" -> {...}\|>`. The default stays `False` until measured on A100/H100. |
| Q2 scales | The radial variable comes from `"IntegrationVariables"`. Dimensionful external coordinates are inferred from frame homogeneity; `"Scales"` overrides; an ambiguous or failed inference is a hard error. |
| Q3 balance | **Generation time only.** The generator decides per trace (inline flop vs B · bytes per node) and emits only that path for `map()`. Retuning needs a regeneration. |
| Q4 DiFfRG | A new branch of `/mnt/data/Documents/Uni/Code/DiFfRG`, in a worktree so the uncommitted work stays untouched, installed to a separate prefix until merged. |

**Consequence of Q3.**
- `kernel()` / `post()` stay emitted with inline traces, because `get()`, `map_points()` and the other integrators
  still need them in v1.
- With `has_tab`, DiFfRG's `map()` takes the table path at compile time, with no run-time fallback. So the inline traces
  are not instantiated in the map TUs, and device compile cost does not grow.
- The prototype's `DIFFRG_NO_TAB` A/B switch therefore becomes "regenerate with `"TraceTables" -> False`".

## 7. Implementation stages

Each stage closes with its gate.

1. **Generator: exact coefficients.**
   - Coefficient extraction in C++ (`codegen/` next to `lower_poly.hpp`): clear denominators, collect by powers,
     check homogeneity, lower each c_n.
   - Emit `tab_fill` + `rtrN` + metadata behind a generator flag.
   - Gate: a unit test against the interpreted inline program at random nodes (fp64, ≤ 1e-12); flag off is
     byte-identical.
2. **Wolfram:**
   - the `"TraceTables"` option;
   - scale inference;
   - generator flag and metadata hand-off;
   - `kernel.hh` emission (`has_tab`, `tab()`, `kernel_tab()`, `post_tab()`) from `ntAssembleIntegrand`.
   - Gate: `ctest -L codegen` green with real generation; `tests/gen` byte-identical with the option off.
3. **DiFfRG branch:**
   - `kernel_has_tab`, `TabTable` (device fill, shape-fastest), `map()` branch, 4D 1/2/3-angle Transform forwarding;
   - a DiFfRG unit test with a toy tabulated kernel against its inline twin.
4. **End to end, YangMills `full_AAcbc_sigma_split` in a scratch copy:**
   - Speed: ≥ 1.35x per RHS on the 4070, fill < 5% of one RHS.
   - Physics: snapshots to t = 1.25 agree with the inline run to ≤ 1e-5 above S0 ≈ 1, and the ZA3/Zcbcbcc values at
     small S0 match the double-literal reference.
   - Then A100/H100 on bwHPC.
5. **Later:** `map_points()`, finite T, dressed traces per α.

## 8. Stage 1 — done (2026-11-12)

**Design change found while implementing: the table has no radial axis, and k never appears.**
- Every trace is homogeneous in the loop momentum l and the external scale S jointly.
- So the cleared numerator is expanded directly as N = Σ_n a_n · l^n · S^(D−n), and the a_n depend only on the
  dimensionless symbols (angles, shape).
- **Table axes:** the non-radial box variables × the non-scale external coordinates. For YangMills that is
  36 × 60 = 2160 nodes, about 0.4 MB for ZA3, instead of 51 840 nodes.
- **Reconstruction:** Horner in r = min(l, S)/max(l, S), as validated.
- It works for any radial grid, and the name of the RG scale is irrelevant.
- For DiFfRG this means a mask on the box variables too, not only on the external axes.

**Code:**
- **`numtracer/include/numtracer/codegen/trace_table.hpp`** (compiled into the library):
  - `table_form(trace, frame.denominators(), {radial, external})`: exact, no fitting. It groups the terms by their
    atom multiset, multiplies by cached Q_j powers, prunes round-off with the lowering's rule, checks homogeneity
    (throws if the scales are wrong), and collects by the power of l. Returns nullopt for complex traces or
    denominators.
  - `eval_table_trace`: host reference.
  - `lower_table_fill` / `emit_table_fill`: one fp64 program for all entries, with shared CSE. A float kernel also
    gets a double `powr` overload.
  - `intern_table_slots` / `emit_table_trace`: the reconstruction `rtrN(row, stride, f)` in the kernel precision.
- **`numtracer/tests/test_trace_table.cpp`** (ctest `trace_table`):
  - real projector traces (one and two `P_T`) and a single-scale case;
  - reconstruction vs `Frame::eval` ≤ 5e-15, interpreted fill vs coefficients ≤ 5e-15;
  - wrong or missing scales throw;
  - emission shape.
- **End to end:** hand-patched copies of the ZA3 and Zcbcbcc generator mains call the new API.
  - The kernel part of the output is byte-identical to the committed kernels. ZA3 differs only in `noinline`, which
    comes from the `NT_GEN_DEVICE` environment the build driver sets.
  - Against the double-precision reference: table ≤ 1.5e-6 at every S0/k, i.e. float32 precision. Inline float32
    reaches 9.3e-3 (ZA3) and 7.9e4 (Zcbcbcc) on the same cells.
- **Fast ctest subset:** frame, lowering_norm, generator_phases, poly_cancel, guards, trace_table all pass.

**Not in Stage 1:**
- dressed traces (`DressedPoly`);
- complex traces;
- the generator-main emission from Wolfram (Stage 2); the hand patch shows the four insertion points:
  `table_form` in the fold callback, `intern_table_slots` before the fill is printed, then `emit_table_fill` +
  `emit_table_trace` + `tab_size` at the end.

## 9. Stage 2 — done (2026-11-12)

**Option.** `"TraceTables"` on `MakeNTKernel` (normalised by `ntTraceTableOption`) and on `MakeNTKernelDiFfRG`, where
it is forwarded in the single-flow path; the kernel-group path does not have it yet.

| value | behaviour |
|---|---|
| `False` (default) | byte-identical to before |
| `Automatic` | tabulates where the inline trace costs more than `$ntTraceTableBalance` = 1 flop per table byte; otherwise falls back with a log line |
| `True` | tabulates every trace it can; aborts when the flow cannot be tabulated |
| `<|"Balance" -> B, "Scales" -> {...}|>` | explicit balance and/or scales |

**Scale inference (`CodegenTables.m`).**
- The radial variable is the first kernel argument.
- The external scale is the one `"CoordinateArguments"` entry that is a frame symbol of its own and makes every
  momentum component exactly degree 1, checked on the exact component polynomials (`ncomp["compExpr"]`, new).
- Ambiguous or failed inference is reported; `"Scales"` overrides.
- v1 refuses dressed, complex and finite-T flows.

**Generator (`CodegenGenerator.m`).**
- Phase B also calls `table_form` per group.
- The balance compares the lowered trace's instruction count with B × entry bytes.
- The slots are interned before the fill is printed.
- After the traces it emits `tab_fill` (fp64, the full fill argument list; scales unused), one `ttrN(row, stride, f)`
  per group (a reconstruction, or a forwarder to `trN`), and `tab_size`.

**`kernel.hh` (`ntTableKernelText` / `ntTraceTableTraits`).**
- Traits: `has_tab` (= `tab_size > 0`), `tab_size`, `tab_box_mask` and `tab_ext_mask`, false exactly on the scales.
- `tab(double *out, std::size_t stride, <kernel params>)`.
- `kernel_tab` / `post_tab`: copies of `kernel` / `post` with `(row, stride)` inserted after the leading parameters
  and every `ns::trN(fenv)` turned into `ns::ttrN(row_, stride_, fenv)`.

**Checks.**
- `tests/test_trace_tables.wls`: 24 checks (option, inference cases, traits, kernel text). Standalone, like the other
  `.wls` unit tests.
- **End to end:** ZA3 and Zcbcbcc of `full_AAcbc_sigma_split`, regenerated from the app's own `YangMills.wl` with
  `"TraceTables" -> True`.
  - The trace part of `kernels.hh` is byte-identical to the same regeneration with tables off.
  - With tables off, `kernel.hh` is byte-identical to the committed one for both flows, and so is Zcbcbcc's
    `kernels.hh`.
  - ZA3's `kernels.hh` differs from the committed one with tables off too, so that is drift of the committed kernels
    against today's engine, not this change.
  - Table accuracy against the double-precision reference: ≤ 1.5e-6 at every S0/k, as in Stage 1.
  - The app's map and get TUs of both flows compile against the current DiFfRG, which does not use the tables yet.
- **Note:** the app's `numtrace` build target compiles generators against the *installed* NumTracer, so table-enabled
  flows need `cmake --install` of this branch before they build there.

## 10. Stage 3 — done (2026-11-12)

**DiFfRG branch `trace-tables`** (worktree `/mnt/data/Documents/Uni/Code/DiFfRG_tables`, commit 6c411727, installed to
`~/.local/share/DiFfRG_tables`, whose `bundled` is a symlink to the shared bundle):
- `physics/integration/tab.hh`: the `kernel_has_tab` concept and `TabTable`.
- `QuadratureIntegrator::map()`:
  - builds the table on the device: fp64 `tab()`, cast to ctype, one launch per grid;
  - keys it on the coordinates and the masked grids only;
  - stores it shape-fastest;
  - calls `post_tab` / `kernel_tab`.
- `set_trace_tables(false)` or `DIFFRG_NO_TRACE_TABLES=1` turned it off (removed in Stage 4: no runtime switch).
- The vacuum Transform wrappers forward the interface.
- Test `tests/physics/integration/tabulated_vacuum.cc`: 602 assertions. `hoisted_vacuum`, `quadrature_integrator` and
  `map_points` still pass.
- **Trap:** indexing the kernel's `static constexpr std::array` mask with a runtime index in device code odr-uses it,
  and nvcc rejects it ("undefined in device code"). Copy it into a captured `device::array` first.

**NumTracer follow-up.**
- The tadpole frames carry unit directions e = p/|p| as frame vectors. Scale inference now accepts every frame vector
  homogeneous of degree 1 (momentum) or 0 (direction), and prefers an external scale over none.
- A trace that is nevertheless not homogeneous in the inferred scales now stays inline under `Automatic`, with a
  generator note; it is still an error under `True`. Correctness never depends on the inference: `table_form`
  re-checks exactly.

**End to end, YangMills `full_AAcbc_sigma_split`, RTX 4070, float kernels.**
- Every flow regenerated with `"TraceTables" -> Automatic`. Every `map()` flow is tabulated; the three lam* flows use
  `map_points` and are not.
- The app is built against both branch installs (`-DDiFfRG_DIR`, `-DNumTracer_DIR`).
- A/B in one binary (`DIFFRG_NO_TRACE_TABLES=1`), nsys, alternating:

| | off | tables | speedup |
|---|--:|--:|--:|
| per RHS (wall clock, incl. table fill) | 575 ms | 382 ms | **1.50x** |
| ZA4tadpoleP1 | 136.7 | 64.4 | 2.12x |
| ZA3 | 77.2 | 36.0 | 2.15x |
| hA | 107.5 | 70.4 | 1.53x |
| Zcbcbcc | 60.3 | 40.3 | 1.49x |
| hc / ZAcbc | 59.4 / 29.5 | 45.3 / 22.4 | 1.31x / 1.32x |
| ZAAcbcT | 64.5 | 63.3 | 1.02x |

- ZA4tadpoleP1 beats its stubbed-traces "ceiling" (1.70x): that ablation still ran `fill()` over the full trace
  environment, while the table path needs only the scales and denominators.
- **Physics, t = 0 … 1.25:** ZA, Zc, ZA4, ZAcbc, ZA4tadpole, ZAAcbcT and m2A agree to ≤ 3e-6, hA/hc to 1e-4. ZA3
  and Zcbcbcc differ only at the lowest S0 nodes, where the double-precision reference showed the float32 inline traces
  are wrong.
- **Full-run time to t = 1.25:** 176 s → 110 s.

**Not done.** `map_points()`, the finite-T integrators, the kernel-group path (`MakeNTKernelGroupDiFfRG`), dressed and
complex traces (done in Stage 4), and measurements on A100/H100.

## 11. Stage 4 — dressed and complex traces (2026-11-12)

**Generator (`codegen/trace_table.hpp`).**
- A trace is a `DressedPoly` Σ_α D^α P_α. `tabulate()` splits every P_α into its real and imaginary part (with
  the inline lowering's noise rule) and takes the table form of each: a `TabulatedTrace` is a list of
  `TablePiece {dressing monomial, imag, TableTrace}`. A real-consumed trace (`realOnly`) drops its imaginary parts.
- The reconstruction `ttrN` returns the inline trace's type (`nt_complex_t` for a complex one) and sums
  D^α · piece, with the dressings read from the kernel env as before. Complex dressings or denominators stay inline.
- `select_table_traces`: keep a table iff inline instructions > reconstruction instructions + B · bytes · entries
  (B ≤ 0 keeps all). The reconstruction term is new; it leaves tiny traces inline (YM: Gsig*, Usig*, h2, ZA, Zc).
- `emit_trace_tables` prints the whole table half: one fp64 `tab_fillN` per tabulated trace (noinline above the
  usual threshold) called from `tab_fill`; traces whose inline body duplicates an earlier one forward to its `ttr`
  and share its entries.
- **Overflow fix (affects Stage 2/3 kernels too).** The old reconstruction formed `powr<d>(s)`, and for a
  single-power piece `l^(d-(D-n)) S^(D-n)`, on their own. In float these overflow at s ~ 1e3, d ≥ 13 (Nf2,
  Λ = 200, has d up to 14), and at small l for negative exponents. Now every piece is built from bounded factors
  (Horner in r, r powers, (s^w inv)^m, the dressings) and `tab_scale<d>(h, s)` multiplies s in last, one factor at
  a time, so a piece overflows only if its value does.

**DiFfRG: no runtime switches.** Whether a flow uses tables (and hoisting) is decided when it is generated
(`"TraceTables"`, `"Hoist"`). `set_trace_tables()`, `DIFFRG_NO_TRACE_TABLES` and `set_hoist()` are gone, and
`map()` / `get()` instantiate one functor, fixed at compile time: tables, else hoisted, else plain. A tabulated
flow's `map()` therefore no longer compiles its inline traces or the other three kernel bodies. `get()` and
`map_points()` still evaluate `kernel()` / `post()`. An A/B is two generations (or, for a harness, the same kernel
class with `has_tab` hidden). To patch on merge: qcd-codes `with_mesons/eval_batch.cc` calls `set_hoist`.

**Wolfram.** `ntTraceTableSpec` no longer refuses dressed or complex flows (only a complex frame and finite T);
`tab_fill` takes the kinematic fill arguments only; `tab()` is emitted after every `kernel()`, i.e. once per
verdict branch of a complex kernel.

**End to end: Nf2 `with_mesons_fullAqbq_FP32` ZA4tadpole** (dressed: 18 dressing slots, 192 of 223 traces complex,
verdict RePart). Regenerated with `"TraceTables" -> Automatic`:
- 221 of 223 traces tabulated (the other 2 duplicate earlier ones), 5877 entries per node; table 36 angle × 40
  shape nodes, 34 MB in float. `kernels.hh` 31 → 51 MB (the fp64 fill).
- Host check, every `ttrN` against `trN` over 400 random nodes: ≤ 2.0e-6 for S0/l1 ∈ [0.1, 10], ≤ 6e-5 out to
  S0/l1 = 1e±3; no non-finite value (612 before the overflow fix).
- GPU, RTX 4070, `map()` on the production tadpole grid with synthetic smooth inputs
  (`gpubench`-style harness, scratch): the generated kernel vs the same kernel class with `has_tab` hidden, i.e.
  the inline path of a tables-off generation; two binaries, alternating:

| k | inline | tables | speedup | max \|tab − inline\| / max \|inline\| |
|---|--:|--:|--:|--:|
| 200 (= Λ) | 681 ms | 85 ms | 8.0x | 3.4e-7 |
| 0.7 | 681 | 86 | 7.9x | 4.2e-6 |

  (20 and 0.05 gave the same speedup and ≤ 4.2e-6 in an earlier one-binary build.) The first `map()` builds the
  table: +50 ms. The committed production kernel takes 684 ms.
- **Compile time of the flow's map TU** (4070 build, sm_89): inline only 20.3 min / 11.5 GB peak; tables only
  7.6 min / 6.4 GB. With the runtime switches (all four paths in one binary) it was 35 min / 12.7 GB.
- Not measured: the flow's share of the full Nf2 RHS, and physics with the full app.


## 12. A100 / H100 (bwHPC, 2026-11-12)

Same sources as the 4070 runs, built for sm_80;90, two binaries per case (tables vs inline), alternating runs on
one GPU per job (jobs 8583061 H100, 8583062 A100 80GB PCIe). Runs repeated twice; spread < 1%.

**Nf2 ZA4tadpole `map()`** (harness of §11, median of 9):

| GPU | k | inline | tables | speedup |
|---|---|--:|--:|--:|
| H100 | 200 | 278 ms | 43.2 ms | 6.4x |
| H100 | 0.7 | 278 | 40.5 | 6.9x |
| A100 | 200 | 397 | 80.8 | 4.9x |
| A100 | 0.7 | 397 | 76.8 | 5.2x |
| 4070 | 200 / 0.7 | 681 | 85 | 8.0x |

First `map()` (table build) +40 ms. Max |tab − inline| / max |inline|: 3.4e-7 (k = 200), 4.2e-6 (k = 0.7), as on
the 4070.

**YangMills `full_AAcbc_sigma_split`**, flow to t = 1.25 (k = 286.5; the run then stops with "Diverging result",
exit code 1, because it has not reached the IR — expected for a truncated flow), average time per RHS:

| GPU | inline | tables | speedup |
|---|--:|--:|--:|
| H100 | 204 ms | 130 ms | 1.57x |
| A100 | 293 | 200 | 1.46x |
| 4070 | 575 | 382 | 1.50x |

Results at t = 1.25: m2A −489.93677 vs −489.93723 (9e-7 relative), Zc(0) 0.96155681 vs 0.96155629 (5e-7);
the tabulated run needs 288 instead of 306 RHS evaluations (different adaptive step sequence).

Both flows use float kernels; the tables win on all three GPUs. The gain is smallest on the A100, but far from break-even.
