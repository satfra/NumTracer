/// @file cse.hpp
/// @brief Common-subexpression elimination: evaluate an expression as a DAG so
///        each distinct subexpression is computed exactly once.
///
/// The eager tree walk @ref numtracer::expr::evaluate recomputes a shared
/// subexpression once per occurrence. The fix is to evaluate the expression as a
/// *directed acyclic graph* (DAG): compute each distinct subexpression once and reuse it.
///
/// A distinct subexpression is, by construction, a distinct C++ type — so type identity
/// already *is* global value-numbering. We exploit that with **value-level hash-consing
/// inside a single `consteval`**, the technique compile-time engines like CTRE use:
///   1. a memoised structural hash @ref numtracer::expr::kHash (O(1)/type, like
///      @ref numtracer::expr::node_count) gives every distinct sub-type a key;
///   2. one `consteval` @ref numtracer::expr::walk builds a `std::vector<Node>`, emitting
///      each distinct sub-type once; a `find` by key **before** recursing prunes the
///      (exponentially larger) tree to O(DAG-edge) executions and shares nodes across the
///      two children of every node;
///   3. the result is materialised into a tight `std::array<Node, N>` via the C++20
///      two-step (`consteval` returns the count `N`, then fills) — **no fixed cap**.
///
/// Only the requested (root) type materialises a @ref numtracer::expr::Ssa; the ~O(DAG)
/// intermediate sub-types instantiate only the cheap @ref numtracer::expr::kHash and
/// @ref numtracer::expr::walk, never a per-type DAG array. This is what lets large traces
/// (the three-gluon triangle and beyond) build without the per-sub-type storage and
/// cumulative-merge cost of a value-numbered-array approach.
///
/// Steps 1–3 run once at compile time; only the slot loop runs per grid point.
/// @see numtracer::expr::eval_real for the real-arithmetic lowering of this DAG.
#pragma once

#include "numtracer/expr/expr.hpp"

#include <array>
#include <bit>
#include <cstdint>
#include <vector>

namespace numtracer::expr
{

  /// @brief Opcode tag for a flattened/SSA node.
  enum Op : int {
    OZERO, ///< Constant 0.
    OONE,  ///< Constant 1.
    OLIT,  ///< Complex literal (value in @ref Node::lit).
    OVAR,  ///< Runtime variable (environment index in @ref Node::varId).
    OADD,  ///< Addition of slots `l` and `r`.
    OMUL   ///< Multiplication of slots `l` and `r`.
  };

  /// @brief One instruction of the flattened expression / CSE DAG.
  struct Node {
    int op = OZERO; ///< Opcode (one of @ref Op).
    int l = -1;     ///< Left child slot; -1 if none.
    int r = -1;     ///< Right child slot; -1 if none.
    int varId = 0;  ///< Environment index, for `OVAR`.
    Cx lit{};       ///< Literal value, for `OLIT`.
  };

  // ---- per-type structural traits (op kind + payload + children) -------------

  /// @brief Maps an expression type to its opcode, payload, and (for binary nodes) children.
  /// @tparam E The expression type.
  template <class E> struct node_traits;
  template <> struct node_traits<Zero> {
    static constexpr int op = OZERO;
  };
  template <> struct node_traits<One> {
    static constexpr int op = OONE;
  };
  template <Cx C> struct node_traits<Lit<C>> {
    static constexpr int op = OLIT;
    static constexpr Cx lit = C; ///< The literal value.
  };
  template <int Id> struct node_traits<Var<Id>> {
    static constexpr int op = OVAR;
    static constexpr int var = Id; ///< The environment index.
  };
  template <class L, class R> struct node_traits<Add<L, R>> {
    static constexpr int op = OADD;
    using left = L;  ///< Left operand type.
    using right = R; ///< Right operand type.
  };
  template <class L, class R> struct node_traits<Mul<L, R>> {
    static constexpr int op = OMUL;
    using left = L;  ///< Left operand type.
    using right = R; ///< Right operand type.
  };

  // ---- memoised structural hash (O(1) per distinct sub-type, like node_count) -

  /// @brief splitmix64 finalizer — scrambles a 64-bit value.
  constexpr std::uint64_t mix64(std::uint64_t x)
  {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
  }
  /// @brief Combine two hashes into one (boost-style mix; order-sensitive).
  constexpr std::uint64_t hcomb(std::uint64_t a, std::uint64_t b)
  {
    return mix64(a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2)));
  }

  /// @brief Structural hash of an expression type — the hash-consing key.
  ///
  /// Bottom-up and memoised per distinct sub-type (the same variable-template shape as
  /// @ref numtracer::expr::node_count), so it costs O(1) per type. Encodes the opcode, the
  /// `Lit`/`Var` payload, and the children's hashes, so two types collide only by chance
  /// (~1e-13 at the scales here); a collision would merge two subexpressions and is caught
  /// loudly by the numeric oracle checks.
  /// @tparam E The expression type.
  template <class E> inline constexpr std::uint64_t kHash = mix64(0); // primary: unused for the node kinds below
  template <> inline constexpr std::uint64_t kHash<Zero> = mix64(1);
  template <> inline constexpr std::uint64_t kHash<One> = mix64(2);
  template <Cx C>
  inline constexpr std::uint64_t kHash<Lit<C>> =
      hcomb(hcomb(3, std::bit_cast<std::uint64_t>(C.re)), std::bit_cast<std::uint64_t>(C.im));
  template <int Id>
  inline constexpr std::uint64_t kHash<Var<Id>> = hcomb(4, static_cast<std::uint64_t>(static_cast<unsigned>(Id)));
  template <class L, class R> inline constexpr std::uint64_t kHash<Add<L, R>> = hcomb(hcomb(5, kHash<L>), kHash<R>);
  template <class L, class R> inline constexpr std::uint64_t kHash<Mul<L, R>> = hcomb(hcomb(6, kHash<L>), kHash<R>);

  // ---- value-level hash-consed DAG construction (one consteval, no caps) ------

  /// @brief Working state for the DAG build: the node list plus an open-addressed hash
  ///        index from structural key to slot, so dedup lookups are O(1) amortised (the
  ///        build is then O(DAG), not O(DAG²) — it scales to far larger expressions).
  struct Builder {
    std::vector<Node> ins;           ///< Emitted nodes, children before parents.
    std::vector<std::uint64_t> keys; ///< `keys[s]` is the @ref kHash of the type at slot `s`.
    std::vector<int> bucket;         ///< Open-addressed table: slot index, or -1 if empty.
    std::size_t mask = 0;            ///< `bucket.size()-1` (a power of two); 0 while empty.

    /// @brief Resize the hash table to `cap` (a power of two) and re-insert every key.
    constexpr void rehash(std::size_t cap)
    {
      bucket.assign(cap, -1);
      mask = cap - 1;
      for (int s = 0; s < (int)keys.size(); ++s) {
        std::size_t p = keys[s] & mask;
        while (bucket[p] != -1)
          p = (p + 1) & mask;
        bucket[p] = s;
      }
    }
    /// @brief The slot holding key `h`, or -1 if not present.
    constexpr int find(std::uint64_t h) const
    {
      if (mask == 0) return -1;
      std::size_t p = h & mask;
      while (bucket[p] != -1) {
        if (keys[bucket[p]] == h) return bucket[p];
        p = (p + 1) & mask;
      }
      return -1;
    }
    /// @brief Append node `nd` with key `h` and index it; returns its slot.
    constexpr int put(std::uint64_t h, Node nd)
    {
      if ((keys.size() + 1) * 10 >= (mask + 1) * 7) rehash(mask == 0 ? 16 : (mask + 1) * 2); // load < 0.7
      const int s = static_cast<int>(ins.size());
      keys.push_back(h);
      ins.push_back(nd);
      std::size_t p = h & mask;
      while (bucket[p] != -1)
        p = (p + 1) & mask;
      bucket[p] = s;
      return s;
    }
  };

  /// @brief Emit `E`'s sub-DAG into `b` (children first), returning E's slot.
  ///
  /// The `find` by @ref kHash runs **before** recursing, so an already-emitted sub-type
  /// returns its slot without re-descending — this prunes the type *tree* to O(DAG edges)
  /// executions and deduplicates nodes shared between a node's two children.
  /// @tparam E The expression type to emit.
  /// @param b The builder to append to.
  /// @return The slot index of `E`.
  template <class E> constexpr int walk(Builder &b)
  {
    constexpr std::uint64_t h = kHash<E>;
    if (const int s = b.find(h); s >= 0) return s; // already emitted — share it
    using NT = node_traits<E>;
    Node nd{};
    if constexpr (NT::op == OZERO)
      nd = Node{OZERO};
    else if constexpr (NT::op == OONE)
      nd = Node{OONE};
    else if constexpr (NT::op == OLIT)
      nd = Node{OLIT, -1, -1, 0, NT::lit};
    else if constexpr (NT::op == OVAR)
      nd = Node{OVAR, -1, -1, NT::var, {}};
    else {
      const int leftSlot = walk<typename NT::left>(b); // recurse only on first encounter
      const int rightSlot = walk<typename NT::right>(b);
      nd = Node{NT::op, leftSlot, rightSlot, 0, {}};
    }
    return b.put(h, nd);
  }

  /// @brief The number of distinct subexpressions (DAG size) of `E`.
  ///
  /// Phase 1 of the two-step: build into a transient `std::vector` and return the count,
  /// which then sizes the persistent array — so there is no hard capacity limit.
  /// @tparam E The expression type.
  template <class E> consteval int ssa_count()
  {
    Builder b;
    (void)walk<E>(b);
    return static_cast<int>(b.ins.size());
  }
  template <class E> inline constexpr int kSsaN = ssa_count<E>();

  /// @brief The SSA/DAG form of an expression: distinct subexpressions, deduplicated.
  /// @tparam E The expression type this SSA was built from. The DAG is in `ins[0..n)`,
  ///         in evaluation order, with the result in slot `root`.
  template <class E> struct Ssa {
    std::array<Node, (kSsaN<E> > 0 ? kSsaN<E> : 1)> ins{}; ///< DAG slots (length = @ref kSsaN).
    int n = 0;                                             ///< Number of distinct subexpressions.
    int root = 0;                                          ///< Slot holding the result.
  };

  /// @brief Lower an expression type to its deduplicated CSE DAG.
  ///
  /// Phase 2 of the two-step: re-run @ref walk into a transient `std::vector`, then copy
  /// into the tight array. The `std::vector` is destroyed before the `consteval` returns,
  /// so only the POD `std::array` escapes (respecting the C++20 transient-allocation rule).
  /// @tparam E The expression type.
  /// @return The @ref Ssa for `E`.
  template <class E> consteval Ssa<E> build_ssa()
  {
    Builder b;
    const int root = walk<E>(b);
    Ssa<E> s;
    s.n = static_cast<int>(b.ins.size());
    s.root = root;
    for (int i = 0; i < s.n; ++i)
      s.ins[i] = b.ins[i];
    return s;
  }
  /// @brief The CSE DAG for `E`, memoised — materialised once per requested type.
  template <class E> inline constexpr Ssa<E> kSsa = build_ssa<E>();
  /// @brief The number of distinct subexpressions (DAG size) of `E`.
  /// @tparam E The expression type.
  template <class E> inline constexpr int ssa_size = kSsaN<E>;

  /// @brief Evaluate the CSE DAG with a runtime interpreter loop.
  ///
  /// Each distinct subexpression is computed exactly once, in slot order, reading
  /// operands from already-filled slots.
  /// @tparam E The expression type to evaluate.
  /// @param env The runtime environment array; `env[id]` supplies each `Var<id>`.
  /// @return The complex value of the expression.
  template <class E> Cx eval_cse(const Cx *env)
  {
    constexpr const Ssa<E> &s = kSsa<E>;
    std::array<Cx, (ssa_size<E> > 0 ? ssa_size<E> : 1)> slots{};
    for (int i = 0; i < s.n; ++i) {
      const Node &in = s.ins[i];
      switch (in.op) {
      case OZERO:
        slots[i] = Cx{0, 0};
        break;
      case OONE:
        slots[i] = Cx{1, 0};
        break;
      case OLIT:
        slots[i] = in.lit;
        break;
      case OVAR:
        slots[i] = env[in.varId];
        break;
      case OADD:
        slots[i] = slots[in.l] + slots[in.r];
        break;
      case OMUL:
        slots[i] = slots[in.l] * slots[in.r];
        break;
      }
    }
    return slots[s.root];
  }

  /// @brief Compute one slot of the unrolled CSE DAG (a single straight-line step).
  ///
  /// The slot's opcode is a compile-time constant (`if constexpr`), so there is no runtime
  /// switch and the operand indices are constants — the compiler register-allocates the
  /// slots and optimizes like a hand-written straight-line kernel.
  /// @tparam E The expression type.
  /// @tparam I The slot index to compute.
  /// @param env The runtime environment array.
  /// @param slots The slot array, with operands for slot `I` already assigned.
  /// @return The value of slot `I`.
  template <class E, int I>
  constexpr Cx cse_step(const Cx *env, const std::array<Cx, (ssa_size<E> > 0 ? ssa_size<E> : 1)> &slots)
  {
    constexpr Node nd = kSsa<E>.ins[I];
    if constexpr (nd.op == OZERO)
      return Cx{0, 0};
    else if constexpr (nd.op == OONE)
      return Cx{1, 0};
    else if constexpr (nd.op == OLIT)
      return nd.lit;
    else if constexpr (nd.op == OVAR)
      return env[nd.varId];
    else if constexpr (nd.op == OADD)
      return slots[nd.l] + slots[nd.r];
    else
      return slots[nd.l] * slots[nd.r]; // OMUL
  }
  /// @brief Fill every slot in order via a fold over the slot index pack.
  /// @tparam E The expression type.
  /// @tparam I The slot indices (a `0..ssa_size<E>-1` pack).
  /// @param env The runtime environment array.
  /// @return The value held in the root slot.
  template <class E, std::size_t... I> Cx eval_unrolled_impl(const Cx *env, std::index_sequence<I...>)
  {
    std::array<Cx, (ssa_size<E> > 0 ? ssa_size<E> : 1)> slots{};
    ((slots[I] = cse_step<E, static_cast<int>(I)>(env, slots)), ...); // assign in slot order
    return slots[kSsa<E>.root];
  }
  /// @brief Evaluate the CSE DAG as straight-line code (the fast realization).
  ///
  /// Same DAG as @ref numtracer::expr::eval_cse, but fully unrolled with no interpreter overhead.
  /// @tparam E The expression type to evaluate.
  /// @param env The runtime environment array; `env[id]` supplies each `Var<id>`.
  /// @return The complex value of the expression.
  template <class E> Cx eval_cse_unrolled(const Cx *env)
  {
    return eval_unrolled_impl<E>(env, std::make_index_sequence<ssa_size<E>>{});
  }

} // namespace numtracer::expr
