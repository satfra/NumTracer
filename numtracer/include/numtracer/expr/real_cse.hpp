/// @file real_cse.hpp
/// @brief Real-arithmetic lowering of the complex CSE DAG.
///
/// @ref numtracer::expr::eval_cse_unrolled evaluates the DAG in complex (@ref
/// numtracer::Cx) arithmetic: every multiply is 4 real multiplies + 2 adds, even
/// though the gamma entries are only `±1`/`±i`, the momentum/dressing leaves are
/// real, and a trace result is real. That complex overhead is the ~2x gap to the
/// FORM kernel (a real polynomial in scalar products).
///
/// This pass lowers the complex DAG to a **real SSA**: each complex slot becomes a
/// pair `(re, im)` of real slots, with a `-1` sentinel meaning "structurally
/// zero". Real leaves keep `im` zero, so a `(real)*(real)` product collapses to a
/// single real multiply and the imaginary machinery largely evaporates. The pass
/// then:
///   - value-numbers the real ops (CSE) with commutative canonicalisation, and
///   - dead-code-eliminates, keeping only real slots reachable from `Re(result)`
///     — the imaginary part of the (real) trace is never computed.
///
/// Result: straight-line real-`double` evaluation at roughly FORM's flop count.
/// Valid when every leaf is real (momenta, dressings), which is the case here.
#pragma once

#include "numtracer/core/config.hpp" // kUnrollMax / kChunkMax / kChunkBlock (lowering tiers)
#include "numtracer/expr/cse.hpp"

#include <array>
#include <vector>

namespace numtracer::expr
{

  /// @brief Opcode tag for a real-SSA instruction.
  enum ROp : int {
    RCONST, ///< Real constant (value in @ref RInstr::k).
    RVAR,   ///< Real variable (environment index in @ref RInstr::a).
    RADD,   ///< `a + b`.
    RSUB,   ///< `a - b`.
    RMUL,   ///< `a * b`.
    RNEG    ///< `-a`.
  };
  /// @brief One instruction of the lowered real SSA.
  struct RInstr {
    int op = RCONST; ///< Opcode (one of @ref ROp).
    int a = -1;      ///< First operand slot; for `RVAR` this is the environment id.
    int b = -1;      ///< Second operand slot (binary ops only).
    double k = 0;    ///< Constant value, for `RCONST`.
  };

  namespace rdetail
  {
    /// @brief Append an instruction, or reuse an identical existing one (value numbering).
    /// @param w The working instruction array.
    /// @param n In/out instruction count.
    /// @param e The instruction to find or add.
    /// @return The slot index of the (possibly newly added) instruction.
    constexpr int find_or_add(RInstr *w, int &n, RInstr e)
    {
      for (int j = 0; j < n; ++j)
        if (w[j].op == e.op && w[j].a == e.a && w[j].b == e.b && w[j].k == e.k) return j;
      w[n] = e;
      return n++;
    }
    /// @brief Emit a real constant, mapping `0` to the structural-zero sentinel.
    /// @param w The working instruction array.
    /// @param n In/out instruction count.
    /// @param k The constant value.
    /// @return The slot of the constant, or `-1` if `k == 0`.
    constexpr int rconst(RInstr *w, int &n, double k) { return k == 0.0 ? -1 : find_or_add(w, n, {RCONST, -1, -1, k}); }
    /// @brief Emit a real variable (environment lookup).
    /// @param w The working instruction array.
    /// @param n In/out instruction count.
    /// @param id The environment index.
    /// @return The slot of the variable.
    constexpr int rvar(RInstr *w, int &n, int id) { return find_or_add(w, n, {RVAR, id, -1, 0}); }
    /// @brief Emit a negation, folding `-(const)` and cancelling `-(-x)`.
    /// @param w The working instruction array.
    /// @param n In/out instruction count.
    /// @param y The operand slot (`-1` for structural zero).
    /// @return The slot of the result, or `-1` if `y` is zero.
    constexpr int rneg(RInstr *w, int &n, int y)
    {
      if (y < 0) return -1;
      if (w[y].op == RCONST) return find_or_add(w, n, {RCONST, -1, -1, -w[y].k}); // -(const)
      if (w[y].op == RNEG) return w[y].a;                                         // -(-x) = x
      return find_or_add(w, n, {RNEG, y, -1, 0});
    }
    /// @brief Emit a multiply, folding `×0`, `×1`, `×(-1)` and canonicalising operand order.
    /// @param w The working instruction array.
    /// @param n In/out instruction count.
    /// @param x First operand slot (`-1` for structural zero).
    /// @param y Second operand slot (`-1` for structural zero).
    /// @return The slot of the product, or `-1` if either operand is zero.
    constexpr int rmul(RInstr *w, int &n, int x, int y)
    {
      if (x < 0 || y < 0) return -1;
      // fold the gamma +-1 entries: x1 = identity, x(-1) = negate (kills trivial mults)
      if (w[x].op == RCONST) {
        if (w[x].k == 1.0) return y;
        if (w[x].k == -1.0) return rneg(w, n, y);
      }
      if (w[y].op == RCONST) {
        if (w[y].k == 1.0) return x;
        if (w[y].k == -1.0) return rneg(w, n, x);
      }
      return find_or_add(w, n, {RMUL, x < y ? x : y, x < y ? y : x, 0}); // commutative -> canonical order
    }
    /// @brief Emit an addition, dropping structural-zero operands and canonicalising order.
    /// @param w The working instruction array.
    /// @param n In/out instruction count.
    /// @param x First operand slot (`-1` for structural zero).
    /// @param y Second operand slot (`-1` for structural zero).
    /// @return The slot of the sum (or the surviving operand if one side is zero).
    constexpr int radd(RInstr *w, int &n, int x, int y)
    {
      if (x < 0) return y;
      if (y < 0) return x;
      return find_or_add(w, n, {RADD, x < y ? x : y, x < y ? y : x, 0});
    }
    /// @brief Emit a subtraction, handling structural-zero operands.
    /// @param w The working instruction array.
    /// @param n In/out instruction count.
    /// @param x Minuend slot (`-1` for structural zero).
    /// @param y Subtrahend slot (`-1` for structural zero).
    /// @return The slot of `x - y` (or `x`, or `-y`, when an operand is zero).
    constexpr int rsub(RInstr *w, int &n, int x, int y)
    {
      if (y < 0) return x;
      if (x < 0) return rneg(w, n, y);
      return find_or_add(w, n, {RSUB, x, y, 0});
    }

    // ---- hash-indexed buffer for the real SSA: no fixed capacity, O(1) amortised dedup -----
    // The pointer-based helpers above carry a fixed capacity (CAP) and are constexpr-friendly for
    // small fixed-size use (kept for expr/poly.hpp). This builder instead grows dynamically, so the
    // real lowering has no CAP *and* value-numbering stays O(w) (not O(w²)) as the real SSA grows.
    struct RBuilder {
      std::vector<RInstr> ins; ///< Real instructions, in emission order.
      std::vector<int> bucket; ///< Open-addressed index: slot, or -1 if empty.
      std::size_t mask = 0;    ///< `bucket.size()-1` (power of two); 0 while empty.

      static constexpr std::uint64_t ihash(const RInstr &e)
      {
        const std::uint64_t h = hcomb(hcomb(e.op, e.a), e.b);
        return hcomb(h, std::bit_cast<std::uint64_t>(e.k));
      }
      static constexpr bool ieq(const RInstr &a, const RInstr &b)
      {
        return a.op == b.op && a.a == b.a && a.b == b.b && a.k == b.k;
      }
      constexpr void rehash(std::size_t cap)
      {
        bucket.assign(cap, -1);
        mask = cap - 1;
        for (int s = 0; s < (int)ins.size(); ++s) {
          std::size_t p = ihash(ins[s]) & mask;
          while (bucket[p] != -1)
            p = (p + 1) & mask;
          bucket[p] = s;
        }
      }
      constexpr int find_or_add(RInstr e)
      {
        if (mask) {
          std::size_t p = ihash(e) & mask;
          while (bucket[p] != -1) {
            if (ieq(ins[bucket[p]], e)) return bucket[p];
            p = (p + 1) & mask;
          }
        }
        if ((ins.size() + 1) * 10 >= (mask + 1) * 7) rehash(mask == 0 ? 16 : (mask + 1) * 2);
        const int s = static_cast<int>(ins.size());
        ins.push_back(e);
        std::size_t p = ihash(e) & mask;
        while (bucket[p] != -1)
          p = (p + 1) & mask;
        bucket[p] = s;
        return s;
      }
    };
    constexpr int rconst(RBuilder &w, double k) { return k == 0.0 ? -1 : w.find_or_add({RCONST, -1, -1, k}); }
    constexpr int rvar(RBuilder &w, int id) { return w.find_or_add({RVAR, id, -1, 0}); }
    constexpr int rneg(RBuilder &w, int y)
    {
      if (y < 0) return -1;
      if (w.ins[y].op == RCONST) return w.find_or_add({RCONST, -1, -1, -w.ins[y].k});
      if (w.ins[y].op == RNEG) return w.ins[y].a;
      return w.find_or_add({RNEG, y, -1, 0});
    }
    constexpr int rmul(RBuilder &w, int x, int y)
    {
      if (x < 0 || y < 0) return -1;
      if (w.ins[x].op == RCONST) {
        if (w.ins[x].k == 1.0) return y;
        if (w.ins[x].k == -1.0) return rneg(w, y);
      }
      if (w.ins[y].op == RCONST) {
        if (w.ins[y].k == 1.0) return x;
        if (w.ins[y].k == -1.0) return rneg(w, x);
      }
      return w.find_or_add({RMUL, x < y ? x : y, x < y ? y : x, 0});
    }
    constexpr int radd(RBuilder &w, int x, int y)
    {
      if (x < 0) return y;
      if (y < 0) return x;
      return w.find_or_add({RADD, x < y ? x : y, x < y ? y : x, 0});
    }
    constexpr int rsub(RBuilder &w, int x, int y)
    {
      if (y < 0) return x;
      if (x < 0) return rneg(w, y);
      return w.find_or_add({RSUB, x, y, 0});
    }
  } // namespace rdetail

  /// @brief The lowered real SSA as a growable vector (intermediate of the two-step).
  struct RealVec {
    std::vector<RInstr> ins; ///< Real instructions, in evaluation order.
    int root = -1;           ///< Slot holding `Re(result)`; -1 if structurally zero.
  };

  /// @brief Lower the complex CSE DAG to a real SSA (split into re/im, CSE, then DCE).
  ///
  /// Each complex slot `(a+bi)` becomes a pair of real slots; a complex product uses
  /// `(a+bi)(c+di) = (ac-bd) + (ad+bc)i`. Only slots reachable from the real part of the
  /// result are retained, then compacted into topological order. A single forward pass
  /// growing into a `std::vector` — **no fixed capacity** (the old `12*CN+32` cap is gone).
  /// @tparam E The expression type to lower.
  /// @return The compacted real SSA (transient vector).
  template <class E> constexpr RealVec lower_real()
  {
    constexpr int cn = ssa_size<E>;
    rdetail::RBuilder w;
    std::vector<int> re(cn > 0 ? cn : 1, -1), im(cn > 0 ? cn : 1, -1);

    using namespace rdetail;
    for (int i = 0; i < cn; ++i) {
      const Node &nd = kSsa<E>.ins[i];
      switch (nd.op) {
      case OZERO:
        re[i] = -1;
        im[i] = -1;
        break;
      case OONE:
        re[i] = rconst(w, 1.0);
        im[i] = -1;
        break;
      case OLIT:
        re[i] = rconst(w, nd.lit.re);
        im[i] = rconst(w, nd.lit.im);
        break;
      case OVAR:
        re[i] = rvar(w, nd.varId);
        im[i] = -1;
        break; // leaves are real
      case OADD:
        re[i] = radd(w, re[nd.l], re[nd.r]);
        im[i] = radd(w, im[nd.l], im[nd.r]);
        break;
      case OMUL: { // (a+bi)(c+di) = (ac-bd) + (ad+bc)i
        const int ac = rmul(w, re[nd.l], re[nd.r]);
        const int bd = rmul(w, im[nd.l], im[nd.r]);
        const int ad = rmul(w, re[nd.l], im[nd.r]);
        const int bc = rmul(w, im[nd.l], re[nd.r]);
        re[i] = rsub(w, ac, bd);
        im[i] = radd(w, ad, bc);
      } break;
      }
    }
    const int rootRe = (cn > 0) ? re[kSsa<E>.root] : -1;
    const int workingSize = static_cast<int>(w.ins.size());

    // Dead-code elimination. Only the real part of the result is wanted, so the imaginary slots of
    // operations that don't feed it are dead. Mark reachability backward from Re(result): since every
    // operand index is strictly below its user (emission order is topological), one reverse sweep
    // marks each operand before — or when — its user is visited.
    std::vector<char> reach(workingSize, 0);
    if (rootRe >= 0) reach[rootRe] = 1;
    for (int i = workingSize - 1; i >= 0; --i)
      if (reach[i]) {
        const RInstr &in = w.ins[i];
        if (in.op == RADD || in.op == RSUB || in.op == RMUL) {
          reach[in.a] = 1;
          reach[in.b] = 1;
        } else if (in.op == RNEG) {
          reach[in.a] = 1;
        }
      }

    // compact reachable slots (topological order preserved), remap operands.
    RealVec out;
    std::vector<int> nidx(workingSize, -1);
    for (int i = 0; i < workingSize; ++i)
      if (reach[i]) {
        RInstr in = w.ins[i];
        if (in.op == RADD || in.op == RSUB || in.op == RMUL) {
          in.a = nidx[in.a];
          in.b = nidx[in.b];
        } else if (in.op == RNEG) {
          in.a = nidx[in.a];
        }
        nidx[i] = static_cast<int>(out.ins.size());
        out.ins.push_back(in);
      }
    out.root = (rootRe >= 0) ? nidx[rootRe] : -1;
    return out;
  }

  /// @brief Number of real instructions after lowering (phase 1 of the two-step).
  template <class E> consteval int real_count() { return static_cast<int>(lower_real<E>().ins.size()); }
  template <class E> inline constexpr int kRealN = real_count<E>();

  /// @brief The lowered real SSA for an expression: split, CSE'd, dead-code-eliminated.
  /// @tparam E The expression type this SSA was built from.
  template <class E> struct RealSsa {
    std::array<RInstr, (kRealN<E> > 0 ? kRealN<E> : 1)> ins{}; ///< Real instructions (length = @ref kRealN).
    int n = 0;                                                 ///< Number of real instructions kept.
    int root = -1;                                             ///< Slot holding `Re(result)`; -1 if structurally zero.
  };
  /// @brief Lower `E` and materialise it into the tight @ref RealSsa array (phase 2).
  template <class E> consteval RealSsa<E> build_real_ssa()
  {
    RealVec v = lower_real<E>();
    RealSsa<E> s;
    s.n = static_cast<int>(v.ins.size());
    s.root = v.root;
    for (int i = 0; i < s.n; ++i)
      s.ins[i] = v.ins[i];
    return s;
  }

  /// @brief The lowered real SSA for `E`, computed once per expression type.
  /// @tparam E The expression type.
  template <class E> inline constexpr auto kRealSsa = build_real_ssa<E>();
  /// @brief The number of real instructions in the lowered SSA of `E`.
  /// @tparam E The expression type.
  template <class E> inline constexpr int real_size = kRealSsa<E>.n;

  /// @brief Compute one slot of the unrolled real SSA (a single straight-line step).
  ///
  /// The slot's opcode is a compile-time constant (`if constexpr`), so the compiler
  /// emits a single real-`double` operation with no runtime dispatch.
  /// @tparam E The expression type.
  /// @tparam I The slot index to compute.
  /// @param renv The real environment array; `renv[id]` holds each `Var<id>`'s value.
  /// @param s The slot array, with operands for slot `I` already assigned.
  /// @return The value of slot `I`.
  template <class E, int I>
  constexpr double rstep(const double *renv, const std::array<double, (real_size<E> > 0 ? real_size<E> : 1)> &s)
  {
    constexpr RInstr in = kRealSsa<E>.ins[I];
    if constexpr (in.op == RCONST)
      return in.k;
    else if constexpr (in.op == RVAR)
      return renv[in.a];
    else if constexpr (in.op == RADD)
      return s[in.a] + s[in.b];
    else if constexpr (in.op == RSUB)
      return s[in.a] - s[in.b];
    else if constexpr (in.op == RMUL)
      return s[in.a] * s[in.b];
    else
      return -s[in.a]; // RNEG
  }
  /// @brief Fill every real slot in order via a fold over the slot index pack.
  /// @tparam E The expression type.
  /// @tparam I The slot indices (a `0..real_size<E>-1` pack).
  /// @param renv The real environment array.
  /// @return The value held in the root slot, or `0` if the result is structurally zero.
  template <class E, std::size_t... I> double real_unrolled_impl(const double *renv, std::index_sequence<I...>)
  {
    std::array<double, (real_size<E> > 0 ? real_size<E> : 1)> s{};
    ((s[I] = rstep<E, static_cast<int>(I)>(renv, s)), ...);
    return kRealSsa<E>.root < 0 ? 0.0 : s[kRealSsa<E>.root];
  }

  // The real lowering is a size-tiered ladder selected by the cutoffs in core/config.hpp:
  //   * up to `kUnrollMax`  — fully-unrolled straight-line code (@ref eval_real). One
  //     `rstep<E,I>` template per SSA slot: branch-free, register-resident, full superscalar
  //     throughput at run time. Compile cost grows **super-linearly** in the slot count
  //     (dominated by per-slot template instantiation over the expression type), so the cutoff
  //     trades run-time speed against compile time.
  //   * up to `kChunkMax`   — chunked straight-line code (@ref eval_real_chunked): a sequence of
  //     `noinline` blocks of `kChunkBlock` slots over a shared buffer. Much faster than the
  //     interpreter while keeping each compiled function small; compile cost is still
  //     super-linear, so this is not a bounded-compile guarantee.
  //   * above `kChunkMax`   — the flat interpreter below: one template instantiation total,
  //     constant compile cost, the backstop that keeps any expression compilable.
  // All three constants live in core/config.hpp and are overridable per-TU with `-D`.

  /// @brief Interpret the real SSA in a runtime loop over the (compile-time) instruction
  ///        array — one template instantiation total, independent of the slot count.
  ///
  /// The CSE/lowering is still done at compile time (`kRealSsa<E>` is a `constexpr`
  /// value); only the final walk is a loop. Slightly slower per eval than the unrolled
  /// path, but compilation stays flat regardless of expression size.
  /// @tparam E The expression type to evaluate.
  template <class E> double eval_real_interp(const double *renv)
  {
    constexpr int n = real_size<E>;
    std::array<double, (n > 0 ? n : 1)> s; // uninitialised: the SSA is topologically ordered,
                                           // so every slot is written before it is read
    for (int i = 0; i < n; ++i) {
      const RInstr in = kRealSsa<E>.ins[i];
      switch (in.op) {
      case RCONST:
        s[i] = in.k;
        break;
      case RVAR:
        s[i] = renv[in.a];
        break;
      case RADD:
        s[i] = s[in.a] + s[in.b];
        break;
      case RSUB:
        s[i] = s[in.a] - s[in.b];
        break;
      case RMUL:
        s[i] = s[in.a] * s[in.b];
        break;
      default:
        s[i] = -s[in.a];
        break; // RNEG
      }
    }
    return kRealSsa<E>.root < 0 ? 0.0 : s[kRealSsa<E>.root];
  }

  /// @brief Compute one chunked-tier slot from a raw slot buffer (operands read via `s[]`).
  ///
  /// Identical arithmetic to @ref rstep, but reads operands from a `const double*` buffer
  /// rather than a sized `std::array` — so one buffer is shared across the `noinline` blocks
  /// of @ref eval_real_chunked. The opcode and operand indices are compile-time constants.
  /// @tparam E The expression type.
  /// @tparam I The slot index to compute.
  template <class E, int I> constexpr double rstep_buf(const double *renv, const double *s)
  {
    constexpr RInstr in = kRealSsa<E>.ins[I];
    if constexpr (in.op == RCONST)
      return in.k;
    else if constexpr (in.op == RVAR)
      return renv[in.a];
    else if constexpr (in.op == RADD)
      return s[in.a] + s[in.b];
    else if constexpr (in.op == RSUB)
      return s[in.a] - s[in.b];
    else if constexpr (in.op == RMUL)
      return s[in.a] * s[in.b];
    else
      return -s[in.a]; // RNEG
  }
  /// @brief Write slot `I` into the buffer, guarded so out-of-range slots are never instantiated.
  ///
  /// The `if constexpr` guard is required: it stops a trailing partial block from instantiating
  /// `rstep_buf<E,I>` (hence the out-of-bounds constexpr `ins[I]`) for `I >= real_size<E>`.
  template <class E, int I> inline void rslot_buf(const double *renv, double *s)
  {
    if constexpr (I < real_size<E>) s[I] = rstep_buf<E, I>(renv, s);
  }
  /// @brief Fill one block of `kChunkBlock` consecutive slots via a fold over the index pack.
  template <class E, int Base, std::size_t... J>
  inline void rblock_impl(const double *renv, double *s, std::index_sequence<J...>)
  {
    (rslot_buf<E, Base + static_cast<int>(J)>(renv, s), ...);
  }
  /// @brief One `noinline` straight-line block — the unit that bounds the optimizer's basic block.
  ///
  /// `noinline` is essential: at `-O3` the blocks would otherwise be inlined back into one giant
  /// function, collapsing the chunked tier into the full-unroll tier (and its compile cost).
  template <class E, int Base> [[gnu::noinline]] void rblock(const double *renv, double *s)
  {
    rblock_impl<E, Base>(renv, s, std::make_index_sequence<kChunkBlock>{});
  }
  /// @brief Emit the blocks `[0,kChunkBlock), [kChunkBlock,2·kChunkBlock), …` covering all slots.
  ///
  /// `if constexpr (Base < real_size<E>)` instantiates exactly `ceil(real_size/kChunkBlock)`
  /// blocks and stops — no runtime block-count bookkeeping.
  template <class E, int Base = 0> void rchunk_rec(const double *renv, double *s)
  {
    if constexpr (Base < real_size<E>) {
      rblock<E, Base>(renv, s);
      rchunk_rec<E, Base + kChunkBlock>(renv, s);
    }
  }
  /// @brief Evaluate the real SSA as a sequence of `noinline` straight-line blocks (middle tier).
  ///
  /// Straight-line speed within each block (no per-op dispatch), but a bounded basic block per
  /// compiled function — ~9x faster than the interpreter for large diagrams (@ref kChunkMax).
  /// @tparam E The expression type to evaluate.
  template <class E> double eval_real_chunked(const double *renv)
  {
    constexpr int n = real_size<E>;
    std::array<double, (n > 0 ? n : 1)> s; // uninitialised; blocks write slots before reading them
    rchunk_rec<E>(renv, s.data());
    return kRealSsa<E>.root < 0 ? 0.0 : s[kRealSsa<E>.root];
  }

  /// @brief Evaluate `Re(expression)` in real arithmetic via a size-tiered ladder.
  ///
  /// Valid when all leaves are real (momenta, dressings) — the case throughout NumTracer.
  /// Dispatches at compile time on the real-SSA size:
  ///   - `<= kUnrollMax`: fully-unrolled straight-line code (fastest, ~FORM's flop count);
  ///   - `<= kChunkMax`:  chunked `noinline` straight-line blocks (~9x over the interpreter);
  ///   - otherwise:       the interpreter loop — the template-free, bounded-compile backstop
  ///                      so a pathologically large trace always compiles.
  /// @tparam E The expression type to evaluate.
  /// @param renv The real environment array; `renv[id]` holds the real value of each `Var<id>`.
  /// @return The real part of the expression's value.
  template <class E> double eval_real(const double *renv)
  {
    if constexpr (real_size<E> <= kUnrollMax)
      return real_unrolled_impl<E>(renv, std::make_index_sequence<real_size<E>>{});
    else if constexpr (real_size<E> <= kChunkMax)
      return eval_real_chunked<E>(renv);
    else
      return eval_real_interp<E>(renv);
  }

} // namespace numtracer::expr
