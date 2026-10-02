# CSE and Horner lowering

> Headers: `codegen/lower.hpp`, `codegen/real_cse.hpp`, `codegen/gen.hpp` · Namespace: `numtracer::network`

After the [numeric engine](numeric-engine.md) contracts a diagram, the result is one small
**polynomial** — an `Poly`: a sum of monomials, each a complex coefficient times a product of
the frame's scalar symbols (`l·p`, `l²`, …) and any surviving inverse-propagator atoms (`1/k²`).
That polynomial is already *collected* (like terms merged during the contraction). What remains
is to turn it into a fast, flat, straight-line real kernel. This page covers the two passes that do
it, both build-time:

1. **Horner factoring** (`lower.hpp`) — rewrite the monomial set as a nest of
   `pivot * (…) + (…)` so shared factors are computed once.
2. **Real value-numbering / CSE** (`real_cse.hpp`) — accumulate the factored arithmetic into a
   real SSA that dedups every repeated subexpression and folds away trivial ops.

The generator (`gen.hpp`) then prints the SSA as `s0 = …; s1 = …; …` C++.

```{note}
There is no compile-time expression-template evaluation here: the contraction runs numerically in
the generator, and everything below operates on plain `Poly` monomials and an ordinary runtime
builder.
```

## From a Poly to a real monomial list

Lowering starts in `numeric::to_genprog`, which flattens the `Poly` into a list of `LMono`
(`lower.hpp`) — the lowering's own monomial type:

```cpp
// lower.hpp
struct LMono {
  double c = 0;                        // real coefficient
  std::vector<std::pair<int,int>> vp;  // (envId, exponent), sorted by envId
};
```

Two things happen in the flattening:

* **Real-part extraction.** A trace of Dirac/colour/Lorentz tensors over real momenta is real,
  so the kernel only needs `Re`. Each `LMono.c` is set to `c.re`; the imaginary half is dropped.
  `to_genprog` scans the coefficients first and only if some `c.im` genuinely survives does it
  also lower an imaginary root — so the pure-real case (the overwhelming majority) never pays for
  complex arithmetic at all. This is where the ~2× complex overhead a naive evaluator would carry
  simply never appears.
* **Env ids.** Every distinct symbol gets one global environment index (an `f[]` slot, assigned by
  `GlobalEnv` in `gen.hpp`): a scalar product becomes a `var` id, a surviving inverse atom an `inv`
  id, and — for dressed diagrams — a dressing factor a `dress` id. `vp` records `(envId, exponent)`
  per factor. Tiny round-off monomials (analytically zero, ~1e-13 of the largest coefficient) are
  pruned here.

```{note}
CSE and Horner never **collect like terms** or change basis. That is fine because the
[numeric engine](numeric-engine.md) already delivers a *collected* polynomial in the frame's
scalar symbols — the like-term collection happened during the contraction, so lowering only has
to remove redundant work.
```

## Horner: factoring the monomial set

`horner` (`lower.hpp`) rewrites a monomial set into a factored nest, emitting into the real
builder `w` as it goes:

```cpp
// lower.hpp (abridged)
constexpr int horner(rdetail::RBuilder &w, std::vector<LMono> terms) {
  // base case: no variables left -> emit the sum of the constants
  if (all terms are constant) return rconst(w, Σ c);

  auto [pivot, pivotExp] = choose_pivot(terms);              // most-frequent var, lowest power
  auto [with, without]   = partition_pivot(terms, pivot, pivotExp);
  const int withSlot    = horner(w, with);                   // recurse on the factored part
  const int withoutSlot = horner(w, without);                // recurse on the rest
  int pivotPow = pivot^pivotExp;                             // built from rmul(pivotVar, …)
  return radd(w, rmul(w, pivotPow, withSlot), withoutSlot);  // pivot^e * horner(with) + horner(without)
}
```

* **`choose_pivot`** picks the variable that occurs in the most terms (a linear tally, no
  `std::map`, so it stays `constexpr`) and the *lowest* exponent it appears with — the most that
  can be factored out.
* **`partition_pivot`** splits the terms into `with` (those containing the pivot, with
  `pivot^pivotExp` divided out) and `without` (the rest).
* The recursion combines them as `pivot^pivotExp · horner(with) + horner(without)`, so a factor
  shared by many terms is emitted once and reused — the standard per-diagram Horner factoring.

`best_into` (`gen.hpp`) runs `horner` under a few pivot orderings on throwaway builders and keeps
the ordering that emits the fewest instructions.

## Real value-numbering: `RBuilder`

Every `rconst`/`rvar`/`rmul`/`radd`/`rsub`/`rneg` above appends into an `RBuilder`
(`real_cse.hpp`): a growable array of real SSA instructions with **hash-consed value numbering**.
`find_or_add` returns the slot of an identical existing instruction instead of appending a
duplicate, so a subexpression is stored — and later computed — exactly once. The array grows
dynamically (no fixed capacity) and the dedup stays O(1) amortised via an open-addressed hash
index.

```cpp
// real_cse.hpp
enum ROp : int { RCONST, RVAR, RADD, RMUL, RNEG };
struct RInstr { int op = RCONST; int a = -1; int b = -1; double value = 0; };
```

The emit helpers do algebraic folding *as they build*, which is where the arithmetic shrinks to a
compact flop count:

* `rconst(0)` returns the `-1` **structural-zero sentinel**; any op with a `-1` operand
  short-circuits (`rmul` → `-1`, `radd` → the other operand). Zero terms cost nothing.
* `rmul` folds `×1` to the operand and `×(-1)` to a negation — killing the trivial multiplies from
  ±1 tensor entries — and canonicalises operand order so `a*b` and `b*a` dedup.
* `rneg` folds `-(const)` and cancels `-(-x)`.

```cpp
// real_cse.hpp
constexpr int rmul(RBuilder &w, int x, int y) {
  if (x < 0 || y < 0) return -1;                                  // ×0 -> structurally zero
  if (w.ins[x].op == RCONST) { if (w.ins[x].value == 1.0) return y; if (w.ins[x].value == -1.0) return rneg(w, y); }
  if (w.ins[y].op == RCONST) { if (w.ins[y].value == 1.0) return x; if (w.ins[y].value == -1.0) return rneg(w, x); }
  return w.find_or_add({RMUL, x < y ? x : y, x < y ? y : x, 0}); // commutative -> canonical
}
```

The result is a `GenProg` (`gen.hpp`): the flat `RInstr` array plus the root slot (and an
imaginary root only when the diagram was genuinely complex).

## Emitting the kernel

`gen.hpp` prints each `GenProg` as a straight-line trace function — one C++ statement per SSA slot,
with `RVAR` slots reading the shared `f[]` environment:

```cpp
// what the printer emits, per opcode
case RCONST: out << in.value;                    break;
case RVAR:   out << "f[" << in.a << "]";        break;   // shared scalar-symbol env
case RADD:   out << "s" << in.a << "+s" << in.b; break;
case RMUL:   out << "s" << in.a << "*s" << in.b; break;
default:     out << "-s" << in.a;               break;   // RNEG
```

Across a kernel's diagrams the env is **shared** (`GlobalEnv`), so each scalar product `l·p` or
inverse `1/k²` is computed once per call in `fill()` and every trace function reads the same
`f[]`. The emitted kernel is therefore `trN(const double* f)` trace functions + a `fill()` + the
per-diagram assembly — self-contained straight-line real arithmetic at a compact flop count.

## Next

You have now seen the lowering end to end. The [worked example](worked-example.md) walks a
complete physical integrand — the quark self-energy — from DSL network to generated kernel.
