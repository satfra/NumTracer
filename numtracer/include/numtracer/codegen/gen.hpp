/// @file gen.hpp
/// @brief Build-time generator support for the numeric kernels: lower several diagrams
///        into a **shared** fundamental-symbol environment, and emit each as a straight-line C++
///        function. Runtime-only (uses `<ostream>`) — this is the codegen side, not the kernel.
///
/// Contracting a large trace at compile time is RAM/time-prohibitive in GCC (it retains every
/// intermediate), whereas the same contraction runs in a fraction of the memory and time at
/// runtime. So the contraction is a codegen step: a small C++ program builds the network (the
/// contraction stays in C++), contracts each diagram numerically → lowers it (see
/// @ref numtracer::numeric::to_genprog), and **prints** a committed header of straight-line
/// `trN(const double* f)` functions — the same flat form as a FORM reference kernel. The kernel
/// computes the shared `f[]` once per call and invokes them.
#pragma once

#include "numtracer/core/export.hpp"   // NUMTRACER_FUNC / NUMTRACER_DEFINE_BODIES (compiled vs header-only)
#include "numtracer/core/envvar.hpp"   // env_flag / env_int — the single truth test for NT_* switches
#include "numtracer/codegen/lower.hpp"
#include "numtracer/network/network.hpp" // NetVal / Elem / GenProg

#include <climits>
#include <cstdint> // SIZE_MAX (the NT_GEN_NOINLINE_MIN "off" sentinel)
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iomanip>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

namespace numtracer::network
{

  /// @brief A shared fundamental-symbol environment: assigns one global env id (`f[]` index) per
  ///        distinct symbol (a scalar product `sp(a,b)` or an inverse `inv(id)`), across all the
  ///        diagrams of a kernel — so the ~few fundamental symbols are computed once per call.
  struct GlobalEnv {
    std::vector<std::tuple<int, int, int>> syms; ///< env id i → (kind 0=sp/1=inv, a, b); inv is (1,id,0).
    std::vector<int> bucket;                     ///< open-addressed index into @ref syms (or -1).
    std::size_t mask = 0;                        ///< bucket.size()-1 (power of two); 0 while empty.

    static std::uint64_t hsym(int k, int a, int b)
    {
      return hash_combine(
          hash_combine(hash_combine(splitmix64_finalise(2), static_cast<std::uint64_t>(static_cast<unsigned>(k))),
                       static_cast<std::uint64_t>(static_cast<unsigned>(a))),
          static_cast<std::uint64_t>(static_cast<unsigned>(b)));
    }
    void rehash(std::size_t cap)
    {
      bucket.assign(cap, -1);
      mask = cap - 1;
      for (std::size_t s = 0; s < syms.size(); ++s) {
        auto [k, a, b] = syms[s];
        std::size_t p = hsym(k, a, b) & mask;
        while (bucket[p] != -1)
          p = (p + 1) & mask;
        bucket[p] = static_cast<int>(s);
      }
    }
    /// First-seen lookup of symbol `(k,a,b)`; appends on miss so env ids stay in first-seen order.
    int intern(int k, int a, int b)
    {
      const std::uint64_t h = hsym(k, a, b);
      if (mask) {
        std::size_t p = h & mask;
        while (bucket[p] != -1) {
          if (syms[bucket[p]] == std::tuple<int, int, int>{k, a, b}) return bucket[p];
          p = (p + 1) & mask;
        }
      }
      if ((syms.size() + 1) * 10 >= (mask + 1) * 7) rehash(mask == 0 ? 16 : (mask + 1) * 2);
      const int s = static_cast<int>(syms.size());
      syms.push_back({k, a, b});
      std::size_t p = h & mask;
      while (bucket[p] != -1)
        p = (p + 1) & mask;
      bucket[p] = s;
      return s;
    }
    int inv_id(int v) { return intern(1, v, 0); }
    /// A raw USER-SYMBOL leaf (kind 3): a kernel argument (a momentum component / angle) the numeric
    /// backend interns directly. It fills its `f[]` slot from the argument verbatim (see
    /// @ref FillFormulas::var) and rides the polynomial as a monomial variable — the inv backend never
    /// emits kind 3, so its env layout / fill are byte-identical.
    int var_id(int v) { return intern(3, v, 0); }
    /// A DRESSING symbol (kind 2): an opaque runtime leaf (a propagator dressing / regulator call)
    /// the kernel evaluates once and stores in `f[]`. It enters the polynomial only as a monomial
    /// factor, so it lowers exactly like an `inv` symbol — the difference is purely how the kernel
    /// fills its slot (a dressing C++ expression vs `1/q²`); see @ref FillFormulas::dress.
    int dr_id(int id) { return intern(2, id, 0); }
  };

  /// @brief Sentinel for @ref GenProg::rootIm: the program has no imaginary part (it is purely real).
  inline constexpr int kRealProgram = INT_MIN;

  /// @brief A lowered diagram program over the shared env (its `RVAR` ids are global `f[]` indices).
  ///
  /// `rootIm == kRealProgram` marks a purely REAL program (`root` is the result slot; emitted as
  /// `double`). Otherwise the program is COMPLEX: the real and imaginary parts share one instruction
  /// stream `ins`, `root` is the real-part slot and `rootIm` the imaginary-part slot (either `-1` if
  /// that part is structurally zero); emitted as `std::complex<double>`. A complex `tr_i` arises when a
  /// complex colour factor (an imaginary non-abelian-vertex T-trace) was folded into the polynomial —
  /// the kernel multiplies it by the (also complex) runtime dressing coefficient and the consumer takes
  /// Re, so the imaginary part must survive to the kernel rather than being dropped at this lowering.
  struct GenProg {
    std::vector<RInstr> ins;
    int root = -1;
    int rootIm = kRealProgram;
  };

  /// @brief Several diagram programs lowered into ONE shared instruction stream (`CrossTraceCSE`).
  ///
  /// Emitted as `void trace_all(const double* f, <T>* t)` rather than N independent `trN()`, so a
  /// subexpression reached by more than one trace is computed once. This is the only way to exploit
  /// CROSS-DIAGRAM monomial duplication: FormTracer sums every diagram into one polynomial before
  /// expanding, so identical monomials collect; NumTracer keeps one polynomial per trace and cannot.
  /// A compiler can never recover it either — collecting a monomial out of N traces is a
  /// floating-point reassociation across function boundaries (measured: force-inlining all 108 traces
  /// of ZAqbq1_147 recovered only -3.5% instructions).
  ///
  /// `root[i]` / `rootIm[i]` are the same pair as @ref GenProg, per trace; `rootIm[i] ==
  /// kRealProgram` marks trace `i` real. The emitted array element type is `std::complex<double>` if
  /// ANY trace is complex (uniform, so the kernel's `tarr[i]` reads keep the type `trN(fenv)`
  /// returned), otherwise `double`.
  /// `offset` is the global trace index of `root[0]`: a fused program may cover a CONTIGUOUS SLICE of
  /// the traces rather than all of them (see `NT_GEN_CC_CHUNK`). Fusing everything into one function
  /// maximises sharing but produces a single enormous basic block that spills; chunking trades a little
  /// sharing for register pressure. `offset` lets each chunk store into the caller's `t[]` directly.
  struct FusedProg {
    std::vector<RInstr> ins;
    std::vector<int> root;
    std::vector<int> rootIm;
    int offset = 0;
  };

  /// @brief Lower a polynomial's monomials into the shared env via CSE (`RBuilder`) + Horner.
  ///
  /// The input is first sorted into a **canonical order** so the emitted program (env-id layout, Horner
  /// pivot choices, hence the operation count) depends only on the polynomial as a *set*, not on the
  /// order the upstream reduction happened to produce it. This makes the kernel reproducible across the
  /// serial and (term-reordering) sharded-parallel rebase paths — the same op count either way.
  // Lowering internals (Horner ordering sweep) — ODR-used only via to_genprog (the library TU).
#if NUMTRACER_DEFINE_BODIES
  namespace gdetail
  {
    /// The deterministic monomial orderings the greedy (order-sensitive) Horner is tried on; the caller
    /// keeps the factorization with the fewest emitted ops — reproducible and never worse than any single
    /// order. Only the first @p maxOrders variants are BUILT (each is a full deep copy of the monomial
    /// set, and best_into sweeps 3 or 1 for mid/large polynomials — building all 8 regardless was pure
    /// waste). The variant SEQUENCE is frozen: reordering it would change which ordering wins a tie on
    /// op count, and with it the emitted kernel bytes.
    inline std::vector<std::vector<LMono>> make_orderings(const std::vector<LMono> &monos,
                                                          std::size_t maxOrders = 8)
    {
      auto compareAscending = [](const LMono &x, const LMono &y) { return x.vp < y.vp; };
      auto degree = [](const LMono &m) {
        int d = 0;
        for (auto [id, e] : m.vp)
          d += e;
        return d;
      };
      std::vector<std::vector<LMono>> orders;
      auto want = [&] { return orders.size() < maxOrders; };
      if (want()) orders.push_back(monos); // as-built (canonical)
      if (want()) {
        auto v = monos;
        std::sort(v.begin(), v.end(), compareAscending);
        orders.push_back(std::move(v));
      }
      if (want()) {
        auto v = monos;
        std::sort(v.begin(), v.end(), compareAscending);
        std::reverse(v.begin(), v.end());
        orders.push_back(std::move(v));
      }
      if (want()) {
        auto v = monos;
        std::sort(v.begin(), v.end(), [](const LMono &x, const LMono &y) {
          if (x.vp.size() != y.vp.size()) return x.vp.size() > y.vp.size(); // most-factors first
          return x.vp < y.vp;
        });
        orders.push_back(std::move(v));
      }
      if (want()) {
        auto v = monos;
        std::sort(v.begin(), v.end(), [](const LMono &x, const LMono &y) {
          if (x.vp.size() != y.vp.size()) return x.vp.size() < y.vp.size(); // fewest-factors first
          return x.vp < y.vp;
        });
        orders.push_back(std::move(v));
      }
      if (want()) {
        auto v = monos;
        std::sort(v.begin(), v.end(), [&](const LMono &x, const LMono &y) {
          const int dx = degree(x), dy = degree(y);
          if (dx != dy) return dx > dy; // highest total degree first
          return x.vp < y.vp;
        });
        orders.push_back(std::move(v));
      }
      if (want()) {
        auto v = monos;
        std::sort(v.begin(), v.end(), [&](const LMono &x, const LMono &y) {
          const int dx = degree(x), dy = degree(y);
          if (dx != dy) return dx < dy; // lowest total degree first
          return x.vp < y.vp;
        });
        orders.push_back(std::move(v));
      }
      if (want()) {
        auto v = monos;
        std::sort(v.begin(), v.end(), [](const LMono &x, const LMono &y) { return x.vp > y.vp; });
        orders.push_back(std::move(v));
      } // reverse-lex
      return orders;
    }
  } // namespace gdetail

  namespace gdetail
  {
    /// Pick the cheapest Horner ordering of @p monos (costed on scratch builders), then replay it into
    /// @p builder (so several parts — e.g. a trace's real and imaginary halves — share one CSE stream). Returns
    /// the result slot in @p builder. Takes the monomials by value so the no-sweep path (numOrderings <= 1,
    /// i.e. every >2000-monomial trace) hands them straight to horner without a deep copy.
    inline int best_into(std::vector<LMono> monos, rdetail::RBuilder &builder)
    {
      // The greedy Horner is order-sensitive, so we cost several deterministic orderings and keep the
      // cheapest. But each trial is a FULL horner pass: for the dense 1/4/7 traces (tens of thousands of
      // monomials) the 8-way sweep dominates GENERATION while changing the op count only ~1% (the op count
      // tracks the monomial count, not the pivot order — see the noise-prune notes). So scale the sweep
      // down with the polynomial size; NT_GEN_HORNER_ORDERS=<n> overrides.
      std::size_t numOrderings = 8;
      if (const long v = env_int("NT_GEN_HORNER_ORDERS", 0); v > 0)
        numOrderings = static_cast<std::size_t>(v);
      else if (monos.size() > 2000)
        numOrderings = 1;
      else if (monos.size() > 500)
        numOrderings = 3;
      if (numOrderings <= 1) return horner(builder, std::move(monos)); // canonical (as-built) order only — no sweep
      auto orders = make_orderings(monos, numOrderings);
      if (numOrderings > orders.size()) numOrderings = orders.size();
      std::size_t bestIdx = 0, bestOps = 0;
      bool have = false;
      for (std::size_t i = 0; i < numOrderings; ++i) {
        rdetail::RBuilder scratch; // cost this ordering on a throwaway builder
        horner(scratch, orders[i]);
        if (!have || scratch.ins.size() < bestOps) {
          bestOps = scratch.ins.size();
          bestIdx = i;
          have = true;
        }
      }
      return horner(builder, std::move(orders[bestIdx]));
    }
  } // namespace gdetail
#endif // NUMTRACER_DEFINE_BODIES

  // Code-printing entry points, ODR-used by the generator TU: declared always, defined once
  // (library TU / header-only build). See core/export.hpp.
  NUMTRACER_FUNC void emit_cpp(std::ostream &out, const GenProg &p, const std::string &name,
                               const std::string &decor = "static inline");
  NUMTRACER_FUNC void emit_cpp_fused(std::ostream &out, const std::vector<FusedProg> &ps,
                                     const std::string &name,
                                     const std::string &decor = "static inline");
  NUMTRACER_FUNC void emit_env_layout(std::ostream &out, const GlobalEnv &g);

#if NUMTRACER_DEFINE_BODIES
  namespace edetail
  {
    /// Per-program statement-emission plan, shared by the single-output (@ref emit_cpp) and fused
    /// (@ref emit_cpp_fused) writers so their statement forms cannot drift. Two statement-level
    /// rewrites, both pure emission (the SSA program itself is untouched, slot numbering included):
    ///
    ///  - fmaFold: a MUL whose ONLY use is one ADD is emitted inside that consumer as
    ///    `fma(a,b,c)` instead of its own `const double` line. gcc/nvcc contract these pairs anyway
    ///    (-ffp-contract=fast / -fmad=true), so the value is (a) the line count — the pair collapses
    ///    to one statement — and (b) making contraction GUARANTEED where the default doesn't reach:
    ///    clang's -ffp-contract=on stops at statement boundaries, and the pair spans two statements.
    ///    Ported from FunKit COEN's fmaRestructure (the consumer kernels already emit fma()).
    ///  - cInline: an RCONST used exactly once is emitted as a (parenthesised) literal at its use
    ///    site instead of occupying its own declaration line.
    ///
    /// Both rewrites are unconditional. They used to sit behind NT_GEN_NO_FMA / NT_GEN_NO_CONST_INLINE
    /// A/B controls; both measured as wins, both hatches were referenced nowhere (no test, no
    /// fixture, no doc), and both were spelled `getenv(x) != nullptr` — under which `NT_GEN_NO_FMA=0`
    /// turned fma folding OFF. Reinstating an A/B here means deleting the branch, not setting a
    /// variable.
    struct EmitPlan
    {
      std::vector<int> use;       ///< total references: instruction operands + result roots
      std::vector<int> consumer;  ///< single consuming instruction, or -2 (root / >1 consumers)
      std::vector<char> fmaFold;  ///< MUL folded into its single ADD consumer
      std::vector<char> cInline;  ///< RCONST inlined at its single use site
    };

    inline EmitPlan make_plan(const std::vector<RInstr> &ins, const int *roots, std::size_t nroots)
    {
      const int n = static_cast<int>(ins.size());
      EmitPlan pl;
      pl.use.assign(ins.size(), 0);
      pl.consumer.assign(ins.size(), -2);
      pl.fmaFold.assign(ins.size(), 0);
      pl.cInline.assign(ins.size(), 0);
      auto touch = [&](int r, int c) {
        if (r < 0 || r >= n) return; // kRealProgram / -1 roots
        pl.consumer[r] = (++pl.use[r] == 1) ? c : -2;
      };
      for (int i = 0; i < n; ++i) {
        const RInstr &in = ins[i];
        if (in.op == RADD || in.op == RMUL) {
          touch(in.a, i);
          touch(in.b, i);
        } else if (in.op == RNEG)
          touch(in.a, i);
      }
      for (std::size_t r = 0; r < nroots; ++r) touch(roots[r], -2);
      // At most ONE folded MUL per consumer (an fma has one addend); prefer operand `a`.
      auto foldable = [&](int r, int c) {
        return r >= 0 && r < n && ins[r].op == RMUL && pl.use[r] == 1 && pl.consumer[r] == c;
      };
      for (int i = 0; i < n; ++i) {
        if (ins[i].op != RADD) continue;
        if (foldable(ins[i].a, i))
          pl.fmaFold[ins[i].a] = 1;
        else if (foldable(ins[i].b, i))
          pl.fmaFold[ins[i].b] = 1;
      }
      for (int i = 0; i < n; ++i)
        if (ins[i].op == RCONST && pl.use[i] == 1) pl.cInline[i] = 1;
      return pl;
    }

    /// Print one operand: an inlined single-use constant as a parenthesised literal (parentheses are
    /// load-bearing: `s5--46.7` would lex as a decrement), anything else as its slot name.
    inline void emit_operand(std::ostream &out, const std::vector<RInstr> &ins, const EmitPlan &pl, int r)
    {
      if (r < 0) {
        out << "0.0";
        return;
      }
      if (pl.cInline[static_cast<std::size_t>(r)])
        out << "(" << ins[static_cast<std::size_t>(r)].value << ")";
      else
        out << "s" << r;
    }

    /// Emit instruction @p i as a full `const double s<i> = <rhs>;` statement — or nothing, when the
    /// plan folded it into its consumer. The opcode set lives here so the two writers cannot drift.
    inline void emit_stmt(std::ostream &out, const std::vector<RInstr> &ins, std::size_t i, const EmitPlan &pl)
    {
      if (pl.fmaFold[i] || pl.cInline[i]) return;
      const RInstr &in = ins[i];
      auto opnd = [&](int r) { emit_operand(out, ins, pl, r); };
      // Emit `fma(a, b, c)` for a folded MUL (`mul` = a*b) and its addend. Only the `+` pattern
      // exists: the SSA has no subtract opcode — a difference is RADD against an RNEG — so there is
      // no `a*b - c` or `c - a*b` form to spell.
      auto fma3 = [&](int mul, int addend) {
        out << "fma(";
        opnd(ins[static_cast<std::size_t>(mul)].a);
        out << ", ";
        opnd(ins[static_cast<std::size_t>(mul)].b);
        out << ", ";
        opnd(addend);
        out << ")";
      };
      out << (pl.use[i] ? "  const double s" : "  [[maybe_unused]] const double s") << i << " = ";
      switch (in.op) {
      case RCONST:
        out << in.value;
        break;
      case RVAR:
        out << "f[" << in.a << "]";
        break;
      case RADD:
        if (in.a >= 0 && pl.fmaFold[static_cast<std::size_t>(in.a)])
          fma3(in.a, in.b);
        else if (in.b >= 0 && pl.fmaFold[static_cast<std::size_t>(in.b)])
          fma3(in.b, in.a);
        else {
          opnd(in.a);
          out << "+";
          opnd(in.b);
        }
        break;
      case RMUL:
        opnd(in.a);
        out << "*";
        opnd(in.b);
        break;
      default: // RNEG
        out << "-";
        opnd(in.a);
        break;
      }
      out << ";\n";
    }

    /// @brief Decide the effective decorator for one emitted trace/chunk function of @p nInstr SSA
    ///        instructions, out-of-lining it when that is faster (see @ref emit_cpp).
    ///
    /// SIZE-GATED, DEVICE ONLY. Measured on sm_89 (RTX 4070), GPU runtime, sweep over 12 flows: whether a
    /// trace function should be `inline` tracks its OWN instruction count, not the whole-kernel size. A
    /// function inlined into the kernel adds its SSA temporaries to the register pool; past ~500
    /// instructions that spills (up to 26 KB/thread on the dense 4-point flows) and the kernel goes
    /// memory-bound, so out-of-lining is runtime-FASTER — ZA4 (655 instr/fn) 0.73x, ZAAqbq1 (772) 0.70x —
    /// as well as ~3x cheaper to compile. BELOW the threshold inlining wins decisively via cross-trace CSE
    /// and no call overhead: ZAqbq1_147's 108 small (125-instr) functions run 1.66x FASTER inlined, even
    /// though that kernel is large. So the gate is per-function, and it never picks the losing side on the
    /// swept flows (>=500: noinline wins or ties; <500: inline wins).
    ///
    /// Only for DEVICE code: the host has no 255-register cliff, and its all-inline emission stays
    /// byte-identical.
    ///
    /// HOW DEVICE-NESS IS DECIDED — and why it is no longer a string sniff. This used to test the
    /// decorator for a literal `__device__` and call that "the CONTRACT". It was silently broken:
    /// `ntKokkosDecor` (Codegen.m) rewrites `__host__ __device__ inline` to `KOKKOS_INLINE_FUNCTION`
    /// BEFORE emission, and DiFfRG_compat.m hands the Kokkos spelling directly for a GPU target — so
    /// no production decorator has ever contained `__device__` and the gate NEVER fired, on any real
    /// flow, with no diagnostic (measured 2026-08-08: QCD_Nf2/no_mesons ZA4 has 75 trace functions,
    /// 7 of them over 500 lines, `tr0` at 2281, and zero `noinline`). The measured cost of missing
    /// it is the 0.70-0.73x above.
    ///
    /// Sniffing for `KOKKOS_` instead would be wrong in the other direction: those macros expand to
    /// plain `inline` on a host-only Kokkos build, where all-inline is correct. So the target is now
    /// stated EXPLICITLY by the caller via `NT_GEN_DEVICE=1` (Codegen.m sets it from the same
    /// condition that picks the decorator). The `__device__` test is kept only as a back-compat
    /// fallback for callers still passing the raw CUDA spelling.
    ///
    /// That fix closed the ONLINE path only, and every production flow is generated OFFLINE — the
    /// generator runs from a `cmake -P` step (NumTracerNumtraceRun.cmake) that inherits nothing of
    /// the emitting Wolfram kernel's environment, so the shell prefix never reached it and the gate
    /// stayed dead regardless. Closed 2026-08-11 by carrying a `"device"` field in numtrace.json —
    /// the same trip the NT_GEN_MAXW thread caps already make, for the same reason. A manifest
    /// without the field reads as false, so flows keep their all-inline emission until re-emitted.
    ///
    /// Overrides:
    ///   NT_GEN_NOINLINE_TRACES : force out-of-line for EVERY function (host+device) — the compile-cost
    ///                            lever for the 100k+-line kernels, and the A/B control.
    ///   NT_GEN_NOINLINE_MIN=N  : per-function threshold, default 500. Out-of-lines when `nInstr > N`,
    ///                            so N=0 out-of-lines every device function that has any instruction at
    ///                            all (a 0-instruction function stays inline — degenerate, and it has
    ///                            nothing to spill). `off` (or any negative value) disables the gate
    ///                            entirely — the escape hatch back to all-inline emission, which is
    ///                            NOT 0. Read ONCE per process and cached: emission must be
    ///                            consistent across every function in a run, so changing the variable
    ///                            mid-process deliberately has no effect.
    ///
    /// WHY THE THRESHOLD IS STILL 500, and what would move it. The 500 came from a 12-flow sm_89
    /// sweep whose companion figures ("ZA4 0.73x") were later re-measured at −4.8% — wrong by ~5x —
    /// so the sweep is not a sound basis for the number. What is solid is a static SASS sweep on nf2
    /// ZA4: on sm_90, ungated leaves 11,636 B of spill at 18% occupancy, min=300 leaves none, min=200
    /// reaches 25% and min=100 reaches 32%, at +2.6%/+5%/+7.6% ops respectively. That argues for
    /// 200-300 on the datacenter parts. It is NOT changed here because the one flow claimed to LOSE
    /// from out-of-lining (ZAqbq1_147, 108 functions of ~125 instructions) sits far below 500 and is
    /// therefore untouched at this threshold — lowering it is what would put that claim in play, and
    /// that claim has never been reproduced. Re-measure it before moving the default.

    /// Is this generation targeting device code? Authoritative signal is `NT_GEN_DEVICE` (set by
    /// Codegen.m from the same condition that chooses the decorator); the raw-CUDA spelling is
    /// honoured as a fallback. Read once per process: emission must be consistent across every
    /// function in a run.
    inline bool device_target(const std::string &decor)
    {
      static const bool envDevice = env_flag("NT_GEN_DEVICE");
      return envDevice || decor.find("__device__") != std::string::npos;
    }

    inline std::string eff_decor(const std::string &decor, std::size_t nInstr = 0)
    {
      std::string effDecor = decor;
      bool noinline = env_flag("NT_GEN_NOINLINE_TRACES");
      if (!noinline && device_target(decor)) {
        static const std::size_t noinlineMinInstr = [] {
          const char *e = std::getenv("NT_GEN_NOINLINE_MIN");
          if (e == nullptr || *e == '\0') return static_cast<std::size_t>(500); // empty means unset,
                                                                               // NOT the very
                                                                               // aggressive 0 below
          // "off" / any negative value disables the gate outright (all-inline emission, the
          // pre-2026-08-11 behaviour). It needs its own spelling because the natural guess, 0,
          // means the OPPOSITE here: the test is `nInstr > N`, so 0 out-of-lines everything.
          // Anything unparsable reads as "off" too: silently treating a typo as 0 would
          // out-of-line every device function in the kernel.
          if (std::strcmp(e, "off") == 0) return SIZE_MAX;
          const long v = env_int("NT_GEN_NOINLINE_MIN", -1);
          return v < 0 ? SIZE_MAX : static_cast<std::size_t>(v);
        }();
        noinline = nInstr > noinlineMinInstr;
      }
      if (noinline) {
        // KOKKOS_FORCEINLINE_FUNCTION expands to `__device__ __host__ __forceinline__`, so appending
        // `__attribute__((noinline))` to it emits a self-contradiction the compiler is free to
        // resolve either way. KOKKOS_INLINE_FUNCTION is only a plain `inline` (Kokkos_Macros.hpp:
        // KOKKOS_IMPL_INLINE_FUNCTION = inline) and would not strictly contradict it — but `inline`
        // still biases the inliner, and mixing the two spellings reads as an accident. Swap either
        // for KOKKOS_FUNCTION — same host/device qualification, no inline hint — and only then
        // attach the attribute.
        //
        // Worth knowing what this is actually overriding: nvcc is NOT force-inlining everything
        // today. On SP_EM ZA4 (55 traces, 435k SSA) it out-of-lines 26 of them on its own; the gate
        // takes that to 51. So the choice is between nvcc's heuristic and the size rule, not between
        // "one giant function" and "many small ones".
        for (const char *kokkosInline : {"KOKKOS_FORCEINLINE_FUNCTION", "KOKKOS_INLINE_FUNCTION"}) {
          const std::size_t at = effDecor.find(kokkosInline);
          if (at != std::string::npos) {
            effDecor.replace(at, std::strlen(kokkosInline), "KOKKOS_FUNCTION");
            break;
          }
        }
        const std::string kw = " inline";
        if (effDecor.size() >= kw.size() && effDecor.compare(effDecor.size() - kw.size(), kw.size(), kw) == 0)
          effDecor.replace(effDecor.size() - kw.size(), kw.size(), " __attribute__((noinline))");
        else
          effDecor += " __attribute__((noinline))";
      }
      return effDecor;
    }
  } // namespace edetail

  /// @brief Print a lowered program as a straight-line C++ function `name(const double* f)`.
  ///        `decor` is the function decorator/prefix (e.g. `"static KOKKOS_INLINE_FUNCTION"`
  ///        for a device-callable kernel); the default keeps the emitted bytes unchanged.
  NUMTRACER_FUNC void emit_cpp(std::ostream &out, const GenProg &p, const std::string &name,
                               const std::string &decor)
  {
    // Per-function inline decision (see @ref edetail::eff_decor): on the DEVICE target a trace function
    // above ~500 SSA instructions is out-of-lined (register isolation — measurably faster AND cheaper to
    // compile), below that it stays inline (cross-trace CSE wins). Host emission is unchanged (byte-identical).
    // NT_GEN_NOINLINE_TRACES forces out-of-line everywhere; `__attribute__((noinline))` is honoured by both
    // g++ and nvcc; fill()/powr stay inline.
    const std::string effDecor = edetail::eff_decor(decor, p.ins.size());
    // Statement plan: liveness ([[maybe_unused]] tagging of dead slots), single-use-constant
    // inlining, and MUL->ADD/SUB fma folding — shared with the fused writer via edetail.
    const int roots[2] = {p.root, p.rootIm}; // kRealProgram (INT_MIN) and -1 are filtered inside
    const edetail::EmitPlan pl = edetail::make_plan(p.ins, roots, 2);
    auto opnd = [&](int r) { edetail::emit_operand(out, p.ins, pl, r); };

    // COMPLEX trace (folded imaginary colour): real + imaginary halves share the instruction stream;
    // return std::complex<double>{re, im}. The kernel multiplies it by the (complex) dressing
    // coefficient and the consumer takes std::real — so the imaginary part reaches the kernel.
    if (p.rootIm != kRealProgram) {
      out << effDecor << " nt_complex_t " << name << "([[maybe_unused]] const double *f) {\n";
      out << std::setprecision(17);
      for (std::size_t i = 0; i < p.ins.size(); ++i) edetail::emit_stmt(out, p.ins, i, pl);
      out << "  return nt_complex_t{";
      opnd(p.root);
      out << ", ";
      opnd(p.rootIm);
      out << "};\n}\n";
      return;
    }
    out << effDecor << " double " << name << "([[maybe_unused]] const double *f) {\n";
    if (p.root < 0) {
      out << "  return 0.0;\n}\n";
      return;
    }
    out << std::setprecision(17);
    for (std::size_t i = 0; i < p.ins.size(); ++i) edetail::emit_stmt(out, p.ins, i, pl);
    out << "  return ";
    opnd(p.root);
    out << ";\n}\n";
  }

  /// @brief Print a FUSED multi-output program as `void name(const double* f, <T>* t)`.
  ///
  /// One instruction stream, N stores. The element type is `std::complex<double>` if any trace is
  /// complex and `double` otherwise — uniform across the array, because the kernel indexes it as
  /// `tarr[i]` wherever it used to call `tr<i>(fenv)` and those calls had a single return type per
  /// trace. A real trace in a complex array stores `{re, 0.0}`; the zero is a compile-time constant
  /// the consumer's arithmetic folds away.
  ///
  /// Liveness is seeded from EVERY root (a slot feeding only trace 7's result is live even though it
  /// feeds no other instruction) — the single-output writer's two `markUse` calls become a loop.
  namespace edetail
  {
    /// One fused chunk, written as `void <name>(const double* f, <T>* t)`. @p anyComplex and @p elemT
    /// are the WHOLE array's verdict, not this chunk's — chunks share the caller's `t[]`.
    inline void emit_fused_one(std::ostream &out, const FusedProg &p, const std::string &name,
                               const std::string &decor, bool anyComplex, const char *elemT)
    {
    const std::string effDecor = edetail::eff_decor(decor, p.ins.size());
    const std::size_t n = p.root.size();

    // Statement plan seeded from EVERY root (a slot feeding only trace 7's result is live even though
    // it feeds no other instruction); kRealProgram (INT_MIN) and -1 are filtered inside.
    std::vector<int> roots;
    roots.reserve(2 * n);
    roots.insert(roots.end(), p.root.begin(), p.root.end());
    roots.insert(roots.end(), p.rootIm.begin(), p.rootIm.end());
    const edetail::EmitPlan pl = edetail::make_plan(p.ins, roots.data(), roots.size());

    auto opnd = [&](int r) { edetail::emit_operand(out, p.ins, pl, r); };
    auto emitStore = [&](std::size_t i) {
      out << "  t[" << (static_cast<std::size_t>(p.offset) + i) << "] = ";
      if (!anyComplex)
        opnd(p.root[i]);
      else if (p.rootIm[i] == kRealProgram) {
        out << "nt_complex_t{";
        opnd(p.root[i]);
        out << ", 0.0}";
      } else {
        out << "nt_complex_t{";
        opnd(p.root[i]);
        out << ", ";
        opnd(p.rootIm[i]);
        out << "}";
      }
      out << ";\n";
    };

    // All N stores go at the END, after the whole instruction stream.
    //
    // REFUTED ALTERNATIVE (measured 2026-07-19, do not retry): sinking each `t[i] = …` to just after
    // its last input slot, so results don't stay live to the end. The register-pressure argument for it
    // is intuitive and WRONG here — on ZAqbq1_147 Mq-in it cut instructions 40.70e9 -> 39.83e9 but RAISED
    // cycles 17.90e9 -> 18.25e9, a consistent 2% runtime LOSS over repeated runs (3018 -> 3080 ns/eval).
    // Interleaving stores into the arithmetic constrains the scheduler more than the shortened live
    // ranges buy back.
    out << effDecor << " void " << name << "([[maybe_unused]] const double *f, " << elemT << " *t) {\n";
    out << std::setprecision(17);
    for (std::size_t i = 0; i < p.ins.size(); ++i) edetail::emit_stmt(out, p.ins, i, pl);
    for (std::size_t i = 0; i < n; ++i)
      emitStore(i);
    out << "}\n";
    }
  } // namespace edetail

  NUMTRACER_FUNC void emit_cpp_fused(std::ostream &out, const std::vector<FusedProg> &ps,
                                     const std::string &name, const std::string &decor)
  {
    // The element type and its published typedef are decided ONCE over all chunks: they share the
    // caller's `t[]`, so a per-chunk verdict could disagree and would redefine `<name>_t`.
    bool anyComplex = false;
    for (const FusedProg &q : ps)
      for (std::size_t i = 0; i < q.root.size(); ++i)
        if (q.rootIm[i] != kRealProgram) anyComplex = true;
    const char *elemT = anyComplex ? "nt_complex_t" : "double";
    out << "using " << name << "_t = " << elemT << ";\n";

    for (std::size_t ci = 0; ci < ps.size(); ++ci)
      edetail::emit_fused_one(out, ps[ci], ps.size() == 1 ? name : name + "_c" + std::to_string(ci),
                              decor, anyComplex, elemT);
    if (ps.size() == 1) return;
    // Wrapper, so the kernel's call site is identical whether or not the traces were chunked.
    out << edetail::eff_decor(decor) << " void " << name << "([[maybe_unused]] const double *f, "
        << elemT << " *t) {\n";
    for (std::size_t ci = 0; ci < ps.size(); ++ci)
      out << "  " << name << "_c" << ci << "(f, t);\n";
    out << "}\n";
  }

  /// @brief Print the shared env layout as a comment (which `f[i]` is which symbol) so the codegen /
  ///        kernel knows the formula to fill each slot with.
  NUMTRACER_FUNC void emit_env_layout(std::ostream &out, const GlobalEnv &g)
  {
    out << "// fundamental-symbol env layout (fill f[i] per call):\n";
    for (std::size_t i = 0; i < g.syms.size(); ++i) {
      auto [kind, a, b] = g.syms[i];
      if (kind == 0)
        out << "//   f[" << i << "] = sp(" << a << "," << b << ")\n";
      else if (kind == 1)
        out << "//   f[" << i << "] = inv(" << a << ")\n";
      else if (kind == 2)
        out << "//   f[" << i << "] = dress(" << a << ")\n";
      else
        out << "//   f[" << i << "] = var(" << a << ")\n";
    }
  }
#endif // NUMTRACER_DEFINE_BODIES

  /// @brief Formula providers for the fundamental symbols: given the symbol the reduction assigned to
  ///        an `f[]` slot, return the C++ expression that computes it from the kernel scalars.
  ///        `sp(a,b)` is a scalar product of fundamental momenta `a·b`; `inv(id)` is `1/q_id²`.
  struct FillFormulas {
    std::function<std::string(int /*a*/, int /*b*/)> sp; ///< C++ for the scalar product `sp(a,b)`.
    std::function<std::string(int /*invId*/)> inv;       ///< C++ for the inverse `inv(invId)`.
    std::function<std::string(int /*drId*/)> dress;      ///< C++ for the dressing/regulator call `dress(drId)`.
    std::function<std::string(int /*varId*/)> var;       ///< C++ for a raw user symbol (numeric backend, kind 3).
  };

  NUMTRACER_FUNC void emit_fill(std::ostream &out, const GlobalEnv &g, const std::string &name,
                                const std::string &argSig, const FillFormulas &fm,
                                const std::string &decor = "static inline");

#if NUMTRACER_DEFINE_BODIES
  /// @brief Print a `fill(double* f, <args>)` that fills every shared `f[]` slot from the kernel
  ///        scalars, using the supplied formula providers. Emitting this **from the generator**
  ///        (which alone knows the reduced layout) decouples the kernel from the `f[]` ordering: the
  ///        kernel just calls `fill(...)` then the `trN(f)`. `argSig` is the parameter list (e.g.
  ///        `"double l1, double cos1, double cos2, double p"`).
  NUMTRACER_FUNC void emit_fill(std::ostream &out, const GlobalEnv &g, const std::string &name,
                                const std::string &argSig, const FillFormulas &fm, const std::string &decor)
  {
    out << std::setprecision(17);
    out << decor << " void " << name << "(double *f, " << argSig << ") {\n";
    for (std::size_t i = 0; i < g.syms.size(); ++i) {
      auto [kind, a, b] = g.syms[i];
      out << "  f[" << i << "] = "
          << (kind == 0   ? fm.sp(a, b)
              : kind == 1 ? fm.inv(a)
              : kind == 2 ? fm.dress(a)
                          : fm.var(a))
          << ";\n";
    }
    out << "}\n";
  }
#endif // NUMTRACER_DEFINE_BODIES

} // namespace numtracer::network
