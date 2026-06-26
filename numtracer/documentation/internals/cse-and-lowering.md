# CSE and real-arithmetic lowering

> Headers: `cse.hpp`, `real_cse.hpp` · Namespace: `numtracer::expr`

After the [numeric engine](numeric-engine.md) reduces a diagram, the result is a polynomial
built from expression types — a pruned sum-of-products of `Var` leaves. For a `Var`-free result you are done
(it already folded to a constant). For a result that carries runtime momenta and
dressings, you still need to turn that type into a *fast numeric kernel*. This page
covers the two passes that do it: common-subexpression elimination (**CSE**),
which removes redundant work, and real-arithmetic lowering, which drops complex
arithmetic the result does not need.

## Why the eager walk is not enough

`expr::evaluate` (from [the expression algebra](expression-algebra.md)) is a *tree*
walk. If a sub-expression appears in two branches of the tree, it is evaluated twice
— and a gamma trace's surviving polynomial shares many sub-expressions. The header
states the problem directly:

> The eager tree walk recomputes a shared subexpression once per occurrence. That is
> why a deferred expression can lose to a dense matmul chain, which collapses each
> intermediate to numbers (maximal sharing). The fix is to evaluate the expression
> as a *directed acyclic graph* (DAG): compute each distinct subexpression once and
> reuse the result.

So the goal is to evaluate the expression as a **DAG**, not a tree.

## CSE: building the DAG

A distinct sub-expression is, by construction, a distinct C++ *type* — so type identity
already *is* global value-numbering. We exploit that with **value-level hash-consing**
inside a single `consteval` (the technique compile-time engines like CTRE use). The
pipeline (`cse.hpp`) is:

1. **Hash** each type. A memoised structural hash `kHash<E>` (the same variable-template
   shape as `node_count`) gives every distinct sub-type an O(1) key.
2. **Walk** the type once. One `consteval` emits each distinct sub-expression exactly once
   into a `std::vector<Node>`, deduplicating by key through an open-addressed hash index.
   The result is *static-single-assignment* (**SSA**) form — a DAG in which every distinct
   sub-expression occupies exactly one slot.
3. **Materialise** the DAG into a tight `std::array<Node, N>` via the C++20 two-step (a
   `consteval` returns the size, then a second fills) — **no fixed-capacity cap**.
4. **Evaluate** with a flat loop over the slots, each computed once.

Only the requested (root) type materialises an `Ssa`; the O(DAG) intermediate sub-types
instantiate only the cheap `kHash`/`walk`, never a per-type array. That is what lets large
traces — the three-gluon triangle and beyond — build cheaply; see [the compile-cost
note](#compile-cost) at the end.

### Node and Op

A flattened/SSA instruction is a plain struct, tagged by an opcode:

```cpp
// cse.hpp
enum Op : int { OZERO, OONE, OLIT, OVAR, OADD, OMUL };
struct Node {
  int op = OZERO; // opcode
  int l = -1;     // left child SSA slot (-1 if none)
  int r = -1;     // right child SSA slot
  int varId = 0;  // env index, for OVAR
  Cx  lit{};      // literal value, for OLIT
};
```

### Hash each type

`kHash<E>` is a memoised structural hash — bottom-up and O(1) per distinct sub-type, the
same variable-template shape as `node_count`. It is the hash-consing key:

```cpp
// cse.hpp
template <class E> inline constexpr std::uint64_t kHash = mix64(0);            // primary (unused)
template <> inline constexpr std::uint64_t kHash<Zero> = mix64(1);
template <int Id> inline constexpr std::uint64_t kHash<Var<Id>> = hcomb(4, (std::uint64_t)Id);
template <class L, class R> inline constexpr std::uint64_t kHash<Add<L, R>> = hcomb(hcomb(5, kHash<L>), kHash<R>);
template <class L, class R> inline constexpr std::uint64_t kHash<Mul<L, R>> = hcomb(hcomb(6, kHash<L>), kHash<R>);
// (Lit hashes its bit-cast components.)
```

A collision would merge two distinct sub-expressions; at these scales that is ~1e-13 and
is caught loudly by the numeric oracle checks (the merge would corrupt the value).

### Walk the type once

`walk<E>` is a single `consteval` that emits each distinct sub-expression into a growing
`std::vector<Node>`. The crucial detail is the `find` by hash **before** recursing: an
already-emitted sub-type returns its slot without re-descending, which prunes the
(exponentially larger) type *tree* to O(DAG-edge) executions and shares nodes between a
node's two children. The hash index makes each lookup O(1), so the whole build is O(DAG):

```cpp
// cse.hpp — Builder holds `ins` (the nodes) + an open-addressed hash index
template <class E> constexpr int walk(Builder &b) {
  if (const int s = b.find(kHash<E>); s >= 0) return s;      // already emitted — share it
  using NT = node_traits<E>;                                 // op + (for binary) children
  Node nd{};
  if constexpr (NT::op == OVAR)      nd = Node{OVAR, -1, -1, NT::var, {}};
  else if constexpr (NT::op == OLIT) nd = Node{OLIT, -1, -1, 0, NT::lit};
  // … Zero / One similar …
  else {                                                     // Add / Mul
    const int li = walk<typename NT::left>(b);               // recurse only on first encounter
    const int ri = walk<typename NT::right>(b);
    nd = Node{NT::op, li, ri, 0, {}};
  }
  return b.put(kHash<E>, nd);                                // append + index
}
```

### Cap-free two-step storage

The `std::vector` is transient to the `consteval` (it never escapes, respecting the C++20
transient-allocation rule). The persistent store is sized to the *actual* DAG via the
two-step — phase 1 returns the count, phase 2 fills a tight `std::array`:

```cpp
template <class E> consteval int ssa_count() { Builder b; (void)walk<E>(b); return (int)b.ins.size(); }
template <class E> inline constexpr int kSsaN = ssa_count<E>();              // no fixed cap
template <class E> struct Ssa { std::array<Node, kSsaN<E>> ins{}; int n=0, root=0; };
template <class E> inline constexpr Ssa<E> kSsa = build_ssa<E>();           // phase 2: fill
template <class E> inline constexpr int ssa_size = kSsaN<E>;                // DAG size
```

The difference between `node_count<E>` (the tree size) and `ssa_size<E>` (the DAG size) is
exactly how much sharing CSE recovered.

### Evaluate the DAG — two realizations

The straightforward realization is an interpreter loop over the slots:

```cpp
template <class E> Cx eval_cse(const Cx *env) {
  constexpr const Ssa<E> &s = kSsa<E>;
  std::array<Cx, /*…*/> slots{};
  for (int i = 0; i < s.n; ++i) {
    const Node &in = s.ins[i];
    switch (in.op) {
      case OZERO: slots[i] = Cx{0,0}; break;
      case OONE:  slots[i] = Cx{1,0}; break;
      case OLIT:  slots[i] = in.lit;  break;
      case OVAR:  slots[i] = env[in.varId]; break;
      case OADD:  slots[i] = slots[in.l] + slots[in.r]; break;
      case OMUL:  slots[i] = slots[in.l] * slots[in.r]; break;
    }
  }
  return slots[s.root];
}
```

That still pays for a runtime `switch` per slot. The faster realization *unrolls*
the slot loop: each slot's opcode is a compile-time constant via `if constexpr`, so
the compiler emits one arithmetic operation per slot with constant operand indices,
register-allocates the slots, and optimizes like a hand-written straight-line
kernel:

```cpp
template <class E, int I> constexpr Cx cse_step(const Cx *env, const std::array<Cx, /*…*/> &slots) {
  constexpr Node nd = kSsa<E>.ins[I];
  if constexpr (nd.op == OZERO) return Cx{0,0};
  else if constexpr (nd.op == OONE) return Cx{1,0};
  else if constexpr (nd.op == OLIT) return nd.lit;
  else if constexpr (nd.op == OVAR) return env[nd.varId];
  else if constexpr (nd.op == OADD) return slots[nd.l] + slots[nd.r];
  else                              return slots[nd.l] * slots[nd.r]; // OMUL
}
// fill every slot in order, then return the root:
template <class E, std::size_t... I> Cx eval_unrolled_impl(const Cx *env, std::index_sequence<I...>) {
  std::array<Cx, /*…*/> slots{};
  ((slots[I] = cse_step<E, (int)I>(env, slots)), ...);   // assign in slot order
  return slots[kSsa<E>.root];
}
template <class E> Cx eval_cse_unrolled(const Cx *env) { /* … calls the above … */ }
```

This already computes each distinct sub-expression exactly once. But it is still
*complex* arithmetic.

## Real lowering: dropping the complex overhead

`eval_cse_unrolled` works in `Cx`: every complex multiply is 4 real multiplies + 2
adds. Yet for these traces, the gamma entries are only $\pm 1$ / $\pm i$, the leaves
(momenta, dressings) are *real*, and the final trace is *real*. All that complex
machinery is the ~2× gap to a FORM kernel (which is a real polynomial in scalar
products). `real_cse.hpp` closes it by lowering the complex DAG to a **real SSA**.

The idea: split every complex slot $(a + bi)$ into a *pair* of real slots
`(re, im)`, with a `-1` sentinel meaning "structurally zero". Real leaves keep `im`
absent, so a (real)·(real) product collapses to a single real multiply and almost
all the imaginary work evaporates. Then value-number the real ops and dead-code
eliminate everything not needed for the real part of the result.

### Real opcodes

```cpp
// real_cse.hpp
enum ROp : int { RCONST, RVAR, RADD, RSUB, RMUL, RNEG };
struct RInstr {
  int op = RCONST; // opcode
  int a = -1;      // first operand slot; for RVAR, the env id
  int b = -1;      // second operand slot
  double k = 0;    // constant value, for RCONST
};
```

Note `RSUB` and `RNEG` appear here but not in the complex DAG — they arise naturally
when a complex product is expanded.

### Splitting a complex multiply

The lowering walks the complex SSA and, for each node, produces `re[i]`/`im[i]` real
slots. Addition is componentwise; multiplication uses $(a+bi)(c+di) = (ac-bd) +
(ad+bc)i$:

```cpp
// real_cse.hpp — OMUL case of lower_real (w is the hash-indexed RBuilder)
case OMUL: { // (a+bi)(c+di) = (ac-bd) + (ad+bc)i
  const int ac = rmul(w, re[nd.l], re[nd.r]);
  const int bd = rmul(w, im[nd.l], im[nd.r]);
  const int ad = rmul(w, re[nd.l], im[nd.r]);
  const int bc = rmul(w, im[nd.l], re[nd.r]);
  re[i] = rsub(w, ac, bd);
  im[i] = radd(w, ad, bc);
} break;
```

The `rmul`/`radd`/`rsub`/`rneg`/`rconst`/`rvar` helpers do value-numbering *and*
algebraic folding as they emit, growing into the `RBuilder`'s `std::vector` (so the real
SSA is cap-free too, and its dedup is O(1) via the same open-addressed hashing as the
complex builder). They are where the savings concentrate:

* `rconst(0)` returns the `-1` structural-zero sentinel; any operation with a `-1`
  operand short-circuits (`rmul` returns `-1`, `radd` returns the other operand).
  This is how the imaginary halves of real leaves cost nothing.
* `rmul` folds `×1` to the operand and `×(-1)` to a negation — killing the trivial
  multiplies from the $\pm 1$ gamma entries — and canonicalises operand order so
  `a*b` and `b*a` dedup.
* `rneg` folds `-(const)` and cancels `-(-x)`.

```cpp
constexpr int rmul(RBuilder &w, int x, int y) {
  if (x < 0 || y < 0) return -1;                                  // ×0 -> structurally zero
  if (w.ins[x].op == RCONST) { if (w.ins[x].k == 1.0) return y; if (w.ins[x].k == -1.0) return rneg(w, y); }
  if (w.ins[y].op == RCONST) { if (w.ins[y].k == 1.0) return x; if (w.ins[y].k == -1.0) return rneg(w, x); }
  return w.find_or_add({RMUL, x < y ? x : y, x < y ? y : x, 0}); // commutative -> canonical
}
```

```{note}
CSE dedups structurally-identical subexpressions but never **collects like terms** or factors.
That is fine here because the [numeric engine](numeric-engine.md) already delivers a *collected*
polynomial in the frame's scalar symbols — the like-term collection happens during the numeric
contraction, before lowering, so the lowering only has to remove redundant work, not change
basis.
```

### Dead-code elimination

The trace result is real, so only the *real* part of the root matters. The pass
marks every slot reachable from `Re(result)` and discards the rest — the imaginary
part of the result is never computed — then compacts the survivors into topological
order, remapping operands:

```cpp
// real_cse.hpp — DCE: keep only slots reachable from Re(result)
std::vector<char> reach(wn, 0);
if (rootRe >= 0) reach[rootRe] = 1;
for (int i = wn - 1; i >= 0; --i)
  if (reach[i]) {
    const RInstr &in = w.ins[i];
    if (in.op == RADD || in.op == RSUB || in.op == RMUL) { reach[in.a] = 1; reach[in.b] = 1; }
    else if (in.op == RNEG) { reach[in.a] = 1; }
  }
```

### Evaluating the real kernel

The result is `RealSsa<E>`, memoised as `kRealSsa<E>` with size `real_size<E>`. Its
evaluator is the unrolled, straight-line, real-`double` kernel — the fast path —
taking a `const double *renv` environment of real leaf values:

```cpp
template <class E> double eval_real(const double *renv) { /* unrolled over real_size<E> slots */ }
```

This is "valid when every leaf is real (momenta, dressings)", which is the case
throughout NumTracer. The result is straight-line real arithmetic at roughly FORM's
floating-point operation count.

## The evaluation-path ladder

There are four ways to turn a surviving entry type into a number, fastest last:

| Path | Function | What it does | When to use |
|---|---|---|---|
| Eager tree | `expr::evaluate<E>` | recursive walk; recomputes shared sub-exprs | constants, tests |
| CSE interpreter | `expr::eval_cse<E>` | DAG, one slot loop with a runtime switch | debugging the DAG |
| CSE unrolled | `expr::eval_cse_unrolled<E>` | straight-line complex arithmetic, no switch | complex-valued results |
| Real lowered | `expr::eval_real<E>` | straight-line real arithmetic, DCE'd | the numerical hot loop |

Three node counts track how much the pipeline shrinks the work:

```cpp
constexpr int treeNodes = expr::node_count<T>; // pruned tree size  (what evaluate walks)
constexpr int cseNodes  = expr::ssa_size<T>;   // complex DAG size  (after CSE dedup)
constexpr int realNodes = expr::real_size<T>;  // real DAG size     (after lowering + DCE)
```

`treeNodes ≥ cseNodes ≥ realNodes`: CSE collapses repeated sub-expressions, and the real
lowering then both folds away the $\pm 1$/$\pm i$ triviality and deletes the unused imaginary
half. `eval_real` is the path the generated kernels use — straight-line real arithmetic at roughly FORM's
floating-point operation count.

(compile-cost)=
## Compile cost — why hash-consing

The DAG is built by **value-level hash-consing in one `consteval`** (see [building the
DAG](#cse-building-the-dag)): each distinct sub-expression is emitted once into a
`std::vector<Node>`, deduplicated through an open-addressed hash index, and the C++20
two-step materialises a tight `std::array` — so the build is **O(DAG)** with **no fixed
capacity**. Only the root type materialises an `Ssa`; the O(DAG) intermediate sub-types
instantiate only the cheap `kHash`/`walk`.

This matters because the build cost stays linear in the DAG with no fixed capacity, which is
what lets large traces lower cheaply. The largest traces can still raise the compiler's
`constexpr` step/operation limits (see the flags in `CMakeLists.txt`); the real-SSA unroll
cutoff (`kUnrollMax`, overridable per-TU via `-DNT_UNROLL_MAX=`) trades a little compile time for
the straight-line hot path up to its slot bound, switching to an interpreter loop above it. All
of this is compile-time only — the generated runtime kernel is unaffected.

## Next

You have now seen every lowering layer in isolation. The [worked example](worked-example.md)
walks a complete physical integrand — the quark self-energy — from network to generated kernel.
