# The expression algebra

> Headers: `cx.hpp`, `expr.hpp` · Namespace: `numtracer::expr`

This is the algebra the kernel **lowering** works in. After the
[numeric engine](numeric-engine.md) reduces a diagram to a polynomial, that polynomial is
turned into fast straight-line code through the expression types defined here. An *expression*
in NumTracer is a C++ type; this page explains what those types are, how building them
automatically simplifies the algebra (the "structural pruning" you will hear about
everywhere), and how a finished expression is turned into a number.

## A complex number you can pass as a template argument

The entries of the tensors NumTracer contracts are complex constants (gamma-matrix
entries are $\pm 1$ and $\pm i$; SU(N) structure constants are real but live
alongside complex generators). To embed a complex *constant* into the type system
we need a complex type that can be a **non-type template parameter (NTTP)** — a
value passed as a template argument.

`std::complex` is not fully usable in `constexpr` contexts before C++23 and is not
a structural type, so NumTracer defines its own minimal complex number, `Cx`:

```cpp
// cx.hpp
struct Cx {
  double re = 0; ///< Real part.
  double im = 0; ///< Imaginary part.
  constexpr Cx operator+(Cx o) const { return {re + o.re, im + o.im}; }
  constexpr Cx operator*(Cx o) const { return {re * o.re - im * o.im, re * o.im + im * o.re}; }
};
```

Two things make this work as a template argument:

* It is a **structural type** — all members are public and of literal type — which
  is exactly the C++20 condition for a class type to be usable as an NTTP.
* Its operations are `constexpr`, so arithmetic on `Cx` values happens at compile
  time.

That is what makes a compile-time complex literal such as `expr::Lit<Cx{0,1}>`
(the imaginary unit $i$) legal. Only `+` and `*` are provided, because that is all
the expression algebra ever needs. The header also has `cabs_real` (a `constexpr`
absolute value) and `approx` (tolerance comparison) used by the tests.

```{note}
`Cx` is deliberately dependency-free — it includes nothing from the Lorentz / Dirac
/ SU(N) sectors. It sits at the very bottom of the layer stack so every other
header can use it.
```

## The six node types

An *expression type* is an empty C++ type carrying no runtime state: the type
itself *is* the expression. There are exactly six node kinds (`expr.hpp`):

```cpp
// expr.hpp — node types (all empty)
struct Zero {};                              // additive identity, multiplicative annihilator (= 0)
struct One {};                               // multiplicative identity (= 1)
template <Cx C> struct Lit {};               // any other compile-time complex constant
template <int Id> struct Var {};             // a runtime value, looked up as env[Id]
template <class L, class R> struct Add {};   // binary sum  L + R
template <class L, class R> struct Mul {};   // binary product  L * R
```

Read these as a small grammar. A *leaf* is `Zero`, `One`, `Lit<C>`, or `Var<Id>`;
an *internal node* is `Add<L,R>` or `Mul<L,R>` whose children are themselves
expression types. For example the type

```cpp
Add<Mul<Lit<Cx{2,0}>, Var<0>>, Mul<Var<0>, Var<1>>>
```

is the expression $2\,v_0 + v_0 v_1$, where $v_0 = \texttt{env[0]}$ and
$v_1 = \texttt{env[1]}$.

Why split out `Zero`, `One`, and `Lit`? Because `Zero` and `One` are the algebra's
identity elements, and giving them dedicated types (rather than `Lit<Cx{0,0}>` and
`Lit<Cx{1,0}>`) lets the combinators recognise them with a trivial type comparison
and apply the simplification laws. That recognition is the whole game.

### Traits over the node types

A handful of compile-time predicates classify these types (`expr.hpp`):

```cpp
template <class T> inline constexpr bool is_zero = std::is_same_v<T, Zero>;
template <class T> inline constexpr bool is_one  = std::is_same_v<T, One>;
```

`const_val<T>` reports whether a type has a known compile-time value and what it
is — specialised to `true` for `Zero` (value $0$), `One` (value $1$) and `Lit<C>`
(value $C$), and `false` otherwise:

```cpp
template <class T> struct const_val { static constexpr bool is = false; };
template <> struct const_val<Zero> { static constexpr bool is = true; static constexpr Cx value{0, 0}; };
template <> struct const_val<One>  { static constexpr bool is = true; static constexpr Cx value{1, 0}; };
template <Cx C> struct const_val<Lit<C>> { static constexpr bool is = true; static constexpr Cx value = C; };
```

`is_expr<T>` marks the six node types and nothing else — its job is to stop the
operator overloads (below) from hijacking ordinary arithmetic like `Cx * Cx`.
Finally `node_count<E>` counts the nodes in a type (leaves count as 1, `Add`/`Mul`
add 1 plus both children) — a direct measure of how much pruning shrank a tree,
and the size bound the CSE pass allocates against.

## The normalizing combinators — where pruning happens

You never write `Mul<A, B>` by hand when building tensors. Instead you call the
combinators `mul(a, b)` and `add(a, b)`, which apply the algebra's rewrite laws as
they construct the result type. They are the single most important functions in the
library.

```cpp
// expr.hpp
template <class A, class B> constexpr auto mul(A, B) {
  if constexpr (is_zero<A> || is_zero<B>)
    return Zero{};                                       // 0*x = 0,  x*0 = 0
  else if constexpr (is_one<A>)
    return B{};                                          // 1*x = x
  else if constexpr (is_one<B>)
    return A{};                                          // x*1 = x
  else if constexpr (const_val<A>::is && const_val<B>::is) // Lit*Lit -> fold the constants
    return lit_norm<const_val<A>::value * const_val<B>::value>();
  else
    return Mul<A, B>{};                                  // otherwise: keep a real product node
}
```

`add` mirrors it:

```cpp
template <class A, class B> constexpr auto add(A, B) {
  if constexpr (is_zero<A>)
    return B{};                                          // 0 + x = x
  else if constexpr (is_zero<B>)
    return A{};                                          // x + 0 = x
  else if constexpr (const_val<A>::is && const_val<B>::is) // const + const -> fold
    return lit_norm<const_val<A>::value + const_val<B>::value>();
  else
    return Add<A, B>{};
}
```

Each branch is an `if constexpr`, so the compiler picks exactly one and the others
do not exist in the generated code. The return *type* differs per branch — that is
the point: `mul(Zero{}, anything)` has return type `Zero`, while `mul(Var<0>{},
Var<1>{})` has return type `Mul<Var<0>, Var<1>>`.

The constant-folding branches call `lit_norm`, which maps a folded `Cx` value back
to the canonical node type so the identity laws keep firing downstream:

```cpp
template <Cx C> constexpr auto lit_norm() {
  if constexpr (C.re == 0.0 && C.im == 0.0) return Zero{}; // 0 -> Zero
  else if constexpr (C.re == 1.0 && C.im == 0.0) return One{}; // 1 -> One
  else return Lit<C>{};
}
```

So `Lit<Cx{2,0}> * Lit<Cx{0.5,0}>` folds to `Cx{1,0}`, which `lit_norm` turns into
`One` — and then a later `One * x` collapses to `x`. The normalisation is what
keeps simplifications cascading.

```{important}
**Structural pruning, precisely.** When a contraction sums products of entries, any
product that touches a `Zero` becomes `Zero` (annihilator), and `add` then drops
that `Zero` from the sum. The dead branch is absent from the final type — not
stored, not evaluated. This is why a contraction over largely-sparse gamma matrices
collapses to a handful of surviving terms instead of $4^k$ of them.
```

For convenience there are `operator*` and `operator+` overloads that forward to
`mul`/`add`, constrained by `is_expr` so they only apply to expression nodes:

```cpp
template <class A, class B> requires(is_expr<A> && is_expr<B>)
constexpr auto operator*(A a, B b) { return mul(a, b); }
```

This is why a contraction can literally accumulate `expr::Zero{} + ... + expr::mul(...)` and
have the fold normalize itself — dead branches drop out as the type is built.

## Evaluation: plugging in the runtime values

Building and pruning produce a type. To get a number you call `evaluate<E>(env)`,
which walks the type and substitutes each `Var<Id>` with `env[Id]`:

```cpp
// expr.hpp — one overload per node kind
constexpr Cx eval(Zero, const Cx *) { return {0, 0}; }
constexpr Cx eval(One,  const Cx *) { return {1, 0}; }
template <Cx C>   constexpr Cx eval(Lit<C>, const Cx *)      { return C; }
template <int Id> constexpr Cx eval(Var<Id>, const Cx *env)  { return env[Id]; }
template <class L, class R> constexpr Cx eval(Add<L, R>, const Cx *env) { return eval(L{}, env) + eval(R{}, env); }
template <class L, class R> constexpr Cx eval(Mul<L, R>, const Cx *env) { return eval(L{}, env) * eval(R{}, env); }

template <class E> constexpr Cx evaluate(const Cx *env) { return eval(E{}, env); }
```

`env` is just a `const Cx *` — a flat array indexed by `Var`'s `Id`. For a
`Var`-free expression (a pure colour factor or gamma trace, which pruned down to a
single `Lit`/`One`/`Zero`) you can pass `nullptr`, because no `Var<Id>` is ever
reached.

```{warning}
`evaluate` is the *eager tree walk*: a sub-expression shared between two branches is
recomputed at each occurrence. That is correct but slow on hot paths, and is the
reason the CSE DAG and real-arithmetic lowering exist — see
[cse-and-lowering](cse-and-lowering.md). Use `evaluate` for compile-time constants
and tests; use `eval_real` in the numerical inner loop.
```

## Worked example: $2v_0 + v_0 v_1$

Build the expression with the combinators, then evaluate it:

```cpp
using namespace numtracer;
using namespace numtracer::expr;

using X = Var<0>;   // env[0]
using Y = Var<1>;   // env[1]

// Build the type: (2 * X) + (X * Y)
using E = decltype(add(mul(Lit<Cx{2, 0}>{}, X{}),
                       mul(X{}, Y{})));

Cx env[2] = {Cx{3.0, 0.0}, Cx{-2.0, 1.0}}; // v0 = 3,  v1 = -2 + i
Cx got = evaluate<E>(env);                  // 2*3 + 3*(-2+i) = 6 + (-6+3i) = 0 + 3i
```

The math: $2v_0 + v_0 v_1 = 2(3) + 3(-2+i) = 6 + (-6 + 3i) = 3i$, so
`got == Cx{0, 3}`.

And here is the pruning in action — these are real `static_assert`s from the test,
checking the *types* that the combinators produce:

```cpp
static_assert(std::is_same_v<decltype(mul(Zero{}, X{})), Zero>, "0*x = 0");
static_assert(std::is_same_v<decltype(mul(One{},  X{})), X>,    "1*x = x");
static_assert(std::is_same_v<decltype(add(Zero{}, X{})), X>,    "0+x = x");

// constant folding normalizes back into Zero/One:
static_assert(std::is_same_v<decltype(mul(Lit<Cx{2,0}>{}, Lit<Cx{0.5,0}>{})), One>,         "2*0.5 -> One");
static_assert(std::is_same_v<decltype(mul(Lit<Cx{0,1}>{}, Lit<Cx{0,1}>{})),   Lit<Cx{-1,0}>>, "i*i -> -1");

// a structural zero kills a whole product branch, then the sum drops it:
//   x*Zero + One*y  ->  Zero + y  ->  y
static_assert(std::is_same_v<decltype(add(mul(X{}, Zero{}), mul(One{}, Y{}))), Y>, "branch pruning");
```

The last assertion is the canonical illustration: a product touching `Zero` becomes
`Zero`, and the surrounding sum drops it, leaving just `Y`. No runtime branch, no
stored zero.

## Where to go next

You now understand the leaf-level algebra. [CSE and real lowering](cse-and-lowering.md) is what
turns a `Var`-carrying expression — the polynomial the [numeric engine](numeric-engine.md)
produces — into a fast, straight-line real kernel.
```
