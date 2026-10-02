# Cx: the compile-time complex type

> Headers: `core/cx.hpp`, `core/lit.hpp` · Namespace: `numtracer`

The entries of the tensors NumTracer contracts are complex constants — gamma-matrix entries are
$\pm 1$ and $\pm i$, and the imaginary unit shows up in projectors and non-abelian colour
factors. NumTracer represents them with its own minimal complex type, `Cx`, rather than
`std::complex`, for two reasons: `std::complex` is not fully usable in `constexpr` contexts before
C++23, and it is not a *structural type*, so it cannot be a non-type template parameter (NTTP).

```cpp
// core/cx.hpp
struct Cx {
  double re = 0; ///< Real part.
  double im = 0; ///< Imaginary part.
  constexpr Cx operator+(Cx o) const { return {re + o.re, im + o.im}; }
  constexpr Cx operator*(Cx o) const { return {re * o.re - im * o.im, re * o.im + im * o.re}; }
};
```

Two properties make `Cx` work as a template argument:

* It is a **structural type** — all members public and of literal type — which is exactly the
  C++20 condition for a class type to be usable as an NTTP.
* Its operations are `constexpr`, so arithmetic on `Cx` values happens at compile time.

`Cx` is deliberately dependency-free: it includes nothing from the Lorentz / Dirac / SU(N)
sectors, and sits at the bottom of the layer stack so every other header can use it. It is the
element type of the dense trace and the numeric `Poly`, and the coefficient type the
[numeric engine](numeric-engine.md) carries through a contraction.

## Lit: a compile-time constant carrier

Because `Cx` is an NTTP, a complex *constant* can be embedded into a type:

```cpp
// core/lit.hpp
template <Cx C> struct Lit {};   // e.g. Lit<Cx{0,1}> is the imaginary unit i
```

`Lit<C>` is an empty type whose only job is to carry a `Cx` value as a template argument. The
code generator uses it to pass scalar coefficients into the emitted helper templates: a diagram's
numeric coefficient is written as `numtracer::Lit<Cx{…}>` and recovered by a small trait
specialisation in the generator program (`litco`/`sc<L>`, emitted by `mathematica/CodegenGenerator.m`).

```{note}
`Lit` is *only* a constant carrier — it is not an expression algebra. The contraction itself runs
numerically in the generator (see [CSE and Horner lowering](cse-and-lowering.md)), so the only
compile-time role left for a type-encoded value is to hand a scalar coefficient to the emitted
helpers.
```

## Where to go next

[CSE and Horner lowering](cse-and-lowering.md) shows how the polynomial the
[numeric engine](numeric-engine.md) produces — with `Cx` coefficients — becomes a fast,
straight-line real kernel.
