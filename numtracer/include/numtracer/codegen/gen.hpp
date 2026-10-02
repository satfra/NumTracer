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
#include "numtracer/core/hash.hpp"     // splitmix64_finalise / hash_combine (GlobalEnv)
#include "numtracer/core/intern_table.hpp"
#include "numtracer/codegen/lower.hpp"
#include "numtracer/codegen/precision.hpp" // float vs double emission
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

namespace numtracer::inline network
{

  /// @brief What an env slot holds. The integer values are hashed (@ref GlobalEnv::SymHash); keep them.
  enum class SymKind : int {
    sp = 0,    ///< scalar product `sp(a,b)` of fundamental momenta
    inv = 1,   ///< inverse `inv(id)` = `1/q_id²`
    dress = 2, ///< dressing / regulator call `dress(id)`
    var = 3,   ///< raw user symbol `var(id)` (numeric backend)
  };

  /// @brief A shared fundamental-symbol environment: assigns one global env id (`f[]` index) per
  ///        distinct symbol (a scalar product `sp(a,b)` or an inverse `inv(id)`), across all the
  ///        diagrams of a kernel — so the ~few fundamental symbols are computed once per call.
  struct GlobalEnv {
    using Sym = std::tuple<SymKind, int, int>; ///< (kind, a, b); every kind but sp has b = 0.
    std::vector<Sym> syms;                     ///< env id i → its symbol
    struct SymHash {
      std::uint64_t operator()(const Sym &s) const
      {
        const auto [k, a, b] = s;
        return hash_combine(
            hash_combine(hash_combine(splitmix64_finalise(2), static_cast<std::uint64_t>(static_cast<unsigned>(k))),
                         static_cast<std::uint64_t>(static_cast<unsigned>(a))),
            static_cast<std::uint64_t>(static_cast<unsigned>(b)));
      }
    };
    struct SymEq {
      bool operator()(const Sym &x, const Sym &y) const { return x == y; }
    };
    InternTable<Sym, SymHash, SymEq> index; ///< hash lookup into @ref syms

    /// First-seen lookup of symbol `(k,a,b)`; appends on miss so env ids stay in first-seen order.
    int intern(SymKind k, int a, int b) { return index.intern(syms, {k, a, b}); }
    int inv_id(int v) { return intern(SymKind::inv, v, 0); }
    /// A raw USER-SYMBOL leaf: a kernel argument (a momentum component / angle) the numeric backend
    /// interns directly. It fills its `f[]` slot from the argument verbatim (see @ref FillFormulas::var)
    /// and rides the polynomial as a monomial variable — the inv backend never emits it, so its env
    /// layout / fill are byte-identical.
    int var_id(int v) { return intern(SymKind::var, v, 0); }
    /// A DRESSING symbol: an opaque runtime leaf (a propagator dressing / regulator call) the kernel
    /// evaluates once and stores in `f[]`. It enters the polynomial only as a monomial factor, so it
    /// lowers exactly like an `inv` symbol — the difference is purely how the kernel fills its slot
    /// (a dressing C++ expression vs `1/q²`); see @ref FillFormulas::dress.
    int dr_id(int id) { return intern(SymKind::dress, id, 0); }
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

  // Lowering internals (the Horner ordering sweep), ODR-used only via to_genprog (the library TU).
  // The monomials arrive in MPoly's sorted order, so the emitted program (env-id layout, Horner pivots,
  // op count) depends only on the polynomial as a SET, not on the order the reduction produced it in.
#if NUMTRACER_DEFINE_BODIES
  namespace gdetail
  {
    /// The deterministic monomial orderings the greedy (order-sensitive) Horner is tried on; the caller
    /// keeps the factorization with the fewest emitted ops — reproducible and never worse than any single
    /// order. Only the first @p maxOrders variants are BUILT (each is a deep copy of the monomial set).
    /// The variant SEQUENCE is frozen: reordering it would change which ordering wins a tie on op
    /// count, and with it the emitted kernel bytes.
    inline std::vector<std::vector<LMono>> make_orderings(const std::vector<LMono> &monos,
                                                          std::size_t maxOrders = 8)
    {
      // Each ordering sorts monomial INDICES with the comparator on the monomials it refers to: the
      // same comparisons give the same permutation as sorting copies, without moving LMonos around.
      const std::size_t n = monos.size();
      std::vector<int> deg(n, 0);
      for (std::size_t i = 0; i < n; ++i)
        for (auto [id, e] : monos[i].vp)
          deg[i] += e;
      auto sortedBy = [&](auto less) {
        std::vector<std::uint32_t> idx(n);
        for (std::size_t i = 0; i < n; ++i) idx[i] = static_cast<std::uint32_t>(i);
        std::sort(idx.begin(), idx.end(), less);
        std::vector<LMono> v;
        v.reserve(n);
        for (std::uint32_t i : idx) v.push_back(monos[i]);
        return v;
      };
      auto vp = [&](std::uint32_t i) -> const auto & { return monos[i].vp; };
      std::vector<std::vector<LMono>> orders;
      auto want = [&] { return orders.size() < maxOrders; };
      if (want()) orders.push_back(monos); // as-built (canonical)
      if (want()) orders.push_back(sortedBy([&](std::uint32_t x, std::uint32_t y) { return vp(x) < vp(y); }));
      if (want()) { // the ascending order just built, reversed
        auto v = orders.back();
        std::reverse(v.begin(), v.end());
        orders.push_back(std::move(v));
      }
      if (want()) // most-factors first
        orders.push_back(sortedBy([&](std::uint32_t x, std::uint32_t y) {
          if (vp(x).size() != vp(y).size()) return vp(x).size() > vp(y).size();
          return vp(x) < vp(y);
        }));
      if (want()) // fewest-factors first
        orders.push_back(sortedBy([&](std::uint32_t x, std::uint32_t y) {
          if (vp(x).size() != vp(y).size()) return vp(x).size() < vp(y).size();
          return vp(x) < vp(y);
        }));
      if (want()) // highest total degree first
        orders.push_back(sortedBy([&](std::uint32_t x, std::uint32_t y) {
          if (deg[x] != deg[y]) return deg[x] > deg[y];
          return vp(x) < vp(y);
        }));
      if (want()) // lowest total degree first
        orders.push_back(sortedBy([&](std::uint32_t x, std::uint32_t y) {
          if (deg[x] != deg[y]) return deg[x] < deg[y];
          return vp(x) < vp(y);
        }));
      if (want()) // reverse-lex
        orders.push_back(sortedBy([&](std::uint32_t x, std::uint32_t y) { return vp(x) > vp(y); }));
      return orders;
    }
  } // namespace gdetail

  namespace gdetail
  {
    /// @brief Traces at or below this monomial count have BOTH lowerings costed; larger ones are
    ///        normalised unconditionally.
    ///
    /// Scalar normalisation is a large win on a trace with many repeated shapes and a small LOSS
    /// (~5-10%) on one with none; costing both arms keeps the lever self-guarding. The second Horner
    /// pass is spent only on small traces: big ones reliably win and are expensive to re-lower
    /// (measurements in docs/NUMTRACER_DESIGN_NOTES.md).
    inline constexpr std::size_t kNormGuardMax = 2000;

    /// Pick the cheapest Horner ordering of @p monos (costed on scratch builders), then replay it into
    /// @p builder (so several parts — e.g. a trace's real and imaginary halves — share one CSE stream). Returns
    /// the result slot in @p builder. Takes the monomials by value so the no-sweep path (numOrderings <= 1,
    /// i.e. every >2000-monomial trace) hands them straight to horner without a deep copy.
    ///
    /// The scalar of the normalised lowering (@ref NVal) is materialised HERE, via @ref scale_into, so
    /// callers keep dealing in plain slots. Doing it at this one seam rather than at each `lower_into`
    /// call site keeps every caller's root a plain slot. `scale_into` is likewise what keeps
    /// a trace whose polynomial is a pure CONSTANT correct: its shape slot is negative, and handing that
    /// to `rmul` — which reads a negative slot as structural zero — would silently emit `return 0.0;`.
    inline int best_into(std::vector<LMono> monos, rdetail::RBuilder &builder)
    {
      // The greedy Horner is order-sensitive, so we cost several deterministic orderings and keep the
      // cheapest. Each trial is a FULL horner pass, and on big polynomials the sweep dominates generation
      // while moving the op count only ~1%, so it shrinks with the polynomial size.
      std::size_t numOrderings = 8;
      if (monos.size() > 2000)
        numOrderings = 1;
      else if (monos.size() > 500)
        numOrderings = 3;

      std::vector<LMono> chosen;
      std::size_t sweptOps = 0; // the winning ordering's op count, normalised lowering; 0 = no sweep ran
      if (numOrderings <= 1)
        chosen = std::move(monos); // canonical (as-built) order only — no sweep, no deep copy
      else {
        auto orders = make_orderings(monos, numOrderings);
        if (numOrderings > orders.size()) numOrderings = orders.size();
        std::size_t bestIdx = 0, bestOps = 0;
        bool have = false;
        for (std::size_t i = 0; i < numOrderings; ++i) {
          rdetail::RBuilder scratch; // cost this ordering on a throwaway builder
          scale_into(scratch, horner(scratch, orders[i], true));
          if (!have || scratch.ins.size() < bestOps) {
            bestOps = scratch.ins.size();
            bestIdx = i;
            have = true;
          }
        }
        chosen = std::move(orders[bestIdx]);
        sweptOps = bestOps;
      }

      // Cost the two lowerings against each other and keep the smaller. Normalisation is a large win
      // on a trace with many repeated shapes and a small loss on one with none, and this is what makes
      // it self-guarding. The comparison is STRICT so a tie keeps the normalised form: flipping ties to
      // the plain lowering would change the emitted kernel on those traces for no gain at all.
      bool useNorm = true;
      if (chosen.size() <= kNormGuardMax) {
        // The sweep already costed `chosen` with the normalised lowering; reuse that count.
        std::size_t normOps = sweptOps;
        if (normOps == 0) {
          rdetail::RBuilder sn;
          scale_into(sn, horner(sn, chosen, true));
          normOps = sn.ins.size();
        }
        rdetail::RBuilder sp;
        scale_into(sp, horner(sp, chosen, false));
        if (sp.ins.size() < normOps) useNorm = false;
      }
      return scale_into(builder, horner(builder, std::move(chosen), useNorm));
    }
  } // namespace gdetail
#endif // NUMTRACER_DEFINE_BODIES

  // Code-printing entry points, ODR-used by the generator TU: declared always, defined once
  // (library TU / header-only build). See core/export.hpp.
  NUMTRACER_FUNC void emit_cpp(std::ostream &out, const GenProg &p, const std::string &name,
                               const std::string &decor = "static inline");
  NUMTRACER_FUNC void emit_env_layout(std::ostream &out, const GlobalEnv &g);

#if NUMTRACER_DEFINE_BODIES
  namespace edetail
  {
    /// Per-program statement-emission plan used by @ref emit_cpp. Two statement-level
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
    /// Both rewrites are unconditional (both measured as wins).
    struct EmitPlan
    {
      std::vector<int> use;       ///< total references: instruction operands + result roots
      std::vector<int> consumer;  ///< single consuming instruction, or -2 (root / >1 consumers)
      std::vector<char> fmaFold;  ///< MUL folded into its single ADD consumer
      std::vector<char> negFold;  ///< NEG folded into its single ADD consumer (spelled as a subtract)
      std::vector<char> cInline;  ///< RCONST inlined at its single use site
    };

    /// @brief Does @p r name an RNEG whose only consumer is the RADD @p c? Such a NEG is spelled as
    ///        the minus of a subtraction at @p c instead of taking a declaration line of its own.
    ///
    /// `make_plan` (which sets the flag) and `emit_stmt` (which prints it) MUST agree on every
    /// instance: a NEG marked folded but not actually consumed leaves a slot referenced and never
    /// emitted. Both therefore reproduce the same operand preference, which is why this predicate —
    /// and the order in which the two callers try `b` before `a` — is stated once, here.
    inline bool neg_foldable(const std::vector<RInstr> &ins, const EmitPlan &pl, int r, int c)
    {
      return r >= 0 && r < static_cast<int>(ins.size()) &&
             ins[static_cast<std::size_t>(r)].op == RNEG && pl.use[static_cast<std::size_t>(r)] == 1 &&
             pl.consumer[static_cast<std::size_t>(r)] == c;
    }

    inline EmitPlan make_plan(const std::vector<RInstr> &ins, const int *roots, std::size_t nroots)
    {
      const int n = static_cast<int>(ins.size());
      EmitPlan pl;
      pl.use.assign(ins.size(), 0);
      pl.consumer.assign(ins.size(), -2);
      pl.fmaFold.assign(ins.size(), 0);
      pl.negFold.assign(ins.size(), 0);
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
      // Then fold a single-use NEG addend into the same ADD, so `x + (-y)` is spelled `x-y` and
      // `fma(a,b,-y)` rather than costing a separate `const double s = -sy;` line. Runs AFTER the fma
      // pass because which operand is left over as the addend depends on that choice.
      for (int i = 0; i < n; ++i) {
        if (ins[i].op != RADD) continue;
        const int a = ins[i].a, b = ins[i].b;
        const bool fa = a >= 0 && pl.fmaFold[a], fb = b >= 0 && pl.fmaFold[b];
        if (fa || fb) {
          const int addend = fa ? b : a;
          if (neg_foldable(ins, pl, addend, i)) pl.negFold[addend] = 1;
        } else if (neg_foldable(ins, pl, b, i))
          pl.negFold[b] = 1;
        else if (neg_foldable(ins, pl, a, i))
          pl.negFold[a] = 1;
      }
      for (int i = 0; i < n; ++i)
        if (ins[i].op == RCONST && pl.use[i] == 1) pl.cInline[i] = 1;
      return pl;
    }

    /// Print a real constant in the emitted precision (see codegen/precision.hpp).
    inline void emit_real(std::ostream &out, double v)
    {
      if (codegen::emit_single())
        out << codegen::float_literal(v);
      else
        out << v;
    }

    inline const char *zero_literal() { return codegen::emit_single() ? "0.f" : "0.0"; }

    /// Print one operand: an inlined single-use constant as a parenthesised literal (parentheses are
    /// load-bearing: `s5--46.7` would lex as a decrement), anything else as its slot name.
    inline void emit_operand(std::ostream &out, const std::vector<RInstr> &ins, const EmitPlan &pl, int r)
    {
      if (r < 0) {
        out << zero_literal();
        return;
      }
      if (pl.cInline[static_cast<std::size_t>(r)]) {
        out << "(";
        emit_real(out, ins[static_cast<std::size_t>(r)].value);
        out << ")";
      } else
        out << "s" << r;
    }

    /// Emit instruction @p i as a full `const double s<i> = <rhs>;` statement — or nothing, when the
    /// plan folded it into its consumer. The opcode set lives here so the two writers cannot drift.
    inline void emit_stmt(std::ostream &out, const std::vector<RInstr> &ins, std::size_t i, const EmitPlan &pl)
    {
      if (pl.fmaFold[i] || pl.negFold[i] || pl.cInline[i]) return;
      const RInstr &in = ins[i];
      auto opnd = [&](int r) { emit_operand(out, ins, pl, r); };
      // Emit `fma(a, b, c)` for a folded MUL (`mul` = a*b) and its addend. The SSA has no subtract
      // opcode — a difference is an RADD against an RNEG — so when that RNEG is single-use the minus
      // is spelled here, as `fma(a, b, -c)`, instead of taking a line of its own.
      auto fma3 = [&](int mul, int addend) {
        const bool neg = addend >= 0 && pl.negFold[static_cast<std::size_t>(addend)];
        out << "fma(";
        opnd(ins[static_cast<std::size_t>(mul)].a);
        out << ", ";
        opnd(ins[static_cast<std::size_t>(mul)].b);
        out << ", ";
        if (neg) {
          out << "-";
          opnd(ins[static_cast<std::size_t>(addend)].a);
        } else
          opnd(addend);
        out << ")";
      };
      out << (pl.use[i] ? "  const " : "  [[maybe_unused]] const ") << codegen::emit_real_type() << " s" << i << " = ";
      switch (in.op) {
      case RCONST:
        emit_real(out, in.value);
        break;
      case RVAR:
        out << "f[" << in.a << "]";
        break;
      case RADD:
        if (in.a >= 0 && pl.fmaFold[static_cast<std::size_t>(in.a)])
          fma3(in.a, in.b);
        else if (in.b >= 0 && pl.fmaFold[static_cast<std::size_t>(in.b)])
          fma3(in.b, in.a);
        else if (in.b >= 0 && pl.negFold[static_cast<std::size_t>(in.b)]) {
          opnd(in.a);
          out << "-";
          opnd(ins[static_cast<std::size_t>(in.b)].a);
        } else if (in.a >= 0 && pl.negFold[static_cast<std::size_t>(in.a)]) {
          opnd(in.b);
          out << "-";
          opnd(ins[static_cast<std::size_t>(in.a)].a);
        } else {
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

    /// Is this generation targeting device code? The only signal is `NT_GEN_DEVICE`, set by
    /// CodegenBuild.m (online) or the numtrace manifest's "device" field (offline) from the same
    /// condition that chooses the decorator, so an explicit `"DeviceTarget" -> False` is honoured even
    /// with a `__device__` decorator. Read once per process: emission must be consistent across every
    /// function in a run.
    inline bool device_target()
    {
      static const bool envDevice = env_flag("NT_GEN_DEVICE");
      return envDevice;
    }

    /// @brief Decide the effective decorator for one emitted trace/chunk function of @p nInstr SSA
    ///        instructions, out-of-lining it when that is faster (see @ref emit_cpp).
    ///
    /// SIZE-GATED, DEVICE ONLY. A function inlined into the kernel adds its SSA temporaries to the
    /// register pool; past the threshold that spills and the kernel goes memory-bound, so
    /// out-of-lining is faster (and cheaper to compile). Below it inlining wins (cross-trace CSE, no
    /// call overhead). The decision is therefore per function, on its OWN instruction count. The host
    /// has no register cliff, so host emission stays all-inline and byte-identical.
    ///
    /// Device-ness comes ONLY from `NT_GEN_DEVICE` (see @ref device_target), never from sniffing the
    /// decorator: production decorators are Kokkos macros that never spell `__device__`, and the
    /// `KOKKOS_` macros expand to plain `inline` on a host-only build. A gate that silently never
    /// fires emits a valid all-inline kernel, so tests/test_noinline_gate.cpp asserts that it FIRES.
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
    /// The default 500 rests on a sweep later shown to be wrong; a static SASS sweep argues for
    /// 200-300 on sm_90. Re-measure before moving it (docs/NUMTRACER_DESIGN_NOTES.md).
    inline std::string eff_decor(const std::string &decor, std::size_t nInstr = 0)
    {
      std::string effDecor = decor;
      bool noinline = env_flag("NT_GEN_NOINLINE_TRACES");
      if (!noinline && device_target()) {
        static const std::size_t noinlineMinInstr = [] {
          const char *e = std::getenv("NT_GEN_NOINLINE_MIN");
          if (e == nullptr || *e == '\0') return static_cast<std::size_t>(500); // empty means unset,
                                                                               // NOT the very
                                                                               // aggressive 0 below
          // "off" / any negative value disables the gate outright (all-inline emission). It needs
          // its own spelling because the natural guess, 0, means the OPPOSITE here: the test is
          // `nInstr > N`, so 0 out-of-lines everything. Anything unparsable reads as "off" too:
          // silently treating a typo as 0 would out-of-line every device function in the kernel.
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
        // attach the attribute. (This overrides nvcc's own inlining heuristic, which already
        // out-of-lines some large functions.)
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
    // Per-function inline decision: see @ref edetail::eff_decor. `__attribute__((noinline))` is
    // honoured by both g++ and nvcc; fill()/powr stay inline.
    const std::string effDecor = edetail::eff_decor(decor, p.ins.size());
    // Statement plan: liveness ([[maybe_unused]] tagging of dead slots), single-use-constant
    // inlining, and MUL->ADD/SUB fma folding (see edetail::EmitPlan).
    const int roots[2] = {p.root, p.rootIm}; // kRealProgram (INT_MIN) and -1 are filtered inside
    const edetail::EmitPlan pl = edetail::make_plan(p.ins, roots, 2);
    auto opnd = [&](int r) { edetail::emit_operand(out, p.ins, pl, r); };

    // COMPLEX trace (folded imaginary colour): real + imaginary halves share the instruction stream;
    // return std::complex<double>{re, im}. The kernel multiplies it by the (complex) dressing
    // coefficient and the consumer takes std::real — so the imaginary part reaches the kernel.
    if (p.rootIm != kRealProgram) {
      out << effDecor << " nt_complex_t " << name << "([[maybe_unused]] const " << codegen::emit_real_type()
          << " *f) {\n";
      out << std::setprecision(17);
      for (std::size_t i = 0; i < p.ins.size(); ++i) edetail::emit_stmt(out, p.ins, i, pl);
      out << "  return nt_complex_t{";
      opnd(p.root);
      out << ", ";
      opnd(p.rootIm);
      out << "};\n}\n";
      return;
    }
    const char *realT = codegen::emit_real_type();
    out << effDecor << " " << realT << " " << name << "([[maybe_unused]] const " << realT << " *f) {\n";
    if (p.root < 0) {
      out << "  return " << edetail::zero_literal() << ";\n}\n";
      return;
    }
    out << std::setprecision(17);
    for (std::size_t i = 0; i < p.ins.size(); ++i) edetail::emit_stmt(out, p.ins, i, pl);
    out << "  return ";
    opnd(p.root);
    out << ";\n}\n";
  }

  /// @brief Print the shared env layout as a comment (which `f[i]` is which symbol) so the codegen /
  ///        kernel knows the formula to fill each slot with.
  NUMTRACER_FUNC void emit_env_layout(std::ostream &out, const GlobalEnv &g)
  {
    out << "// fundamental-symbol env layout (fill f[i] per call):\n";
    for (std::size_t i = 0; i < g.syms.size(); ++i) {
      auto [kind, a, b] = g.syms[i];
      switch (kind) {
      case SymKind::sp: out << "//   f[" << i << "] = sp(" << a << "," << b << ")\n"; break;
      case SymKind::inv: out << "//   f[" << i << "] = inv(" << a << ")\n"; break;
      case SymKind::dress: out << "//   f[" << i << "] = dress(" << a << ")\n"; break;
      case SymKind::var: out << "//   f[" << i << "] = var(" << a << ")\n"; break;
      }
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
    std::function<std::string(int /*varId*/)> var;       ///< C++ for a raw user symbol (numeric backend, SymKind::var).
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
    out << decor << " void " << name << "(" << codegen::emit_real_type() << " *f, " << argSig << ") {\n";
    for (std::size_t i = 0; i < g.syms.size(); ++i) {
      auto [kind, a, b] = g.syms[i];
      out << "  f[" << i << "] = "
          << (kind == SymKind::sp      ? fm.sp(a, b)
              : kind == SymKind::inv   ? fm.inv(a)
              : kind == SymKind::dress ? fm.dress(a)
                                       : fm.var(a))
          << ";\n";
    }
    out << "}\n";
  }
#endif // NUMTRACER_DEFINE_BODIES

} // namespace numtracer::network
