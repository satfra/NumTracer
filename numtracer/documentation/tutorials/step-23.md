# step-23: Trace tables — fill once, evaluate at every scale

*Builds on: [step-09](step-09.md), [step-15](step-15.md) · Built on by: — · Tags: `codegen`,
`performance`, `diffrg` · **Tier B** (a Wolfram kernel and FunKit)*

## Introduction

An fRG integrator evaluates the integrand on a fixed set of quadrature nodes, again and again: once
per external point, once per RG scale $k$, once per right-hand-side evaluation. Most of the
kinematics never changes between those calls. The angles sit on the same nodes all run long; only
the momentum *scales* move — the loop momentum $l = k\,\hat q$ and the external momentum $p$. Yet an
inline trace redoes all of its angle-dependent arithmetic on every call.

A trace table splits the two. Every trace is homogeneous in the momenta, so with $c = \cos\theta$
and the remaining denominators $\mathrm{inv}_j$ (atoms like $1/(l-p)^2$, step-04) set aside, it is

$$
T = N(l, p)\,\prod_j \mathrm{inv}_j^{\,m_j}, \qquad
N(l, p) = \sum_{n=n_0}^{n_1} c_n(c)\; l^{\,n}\, p^{\,D-n},
$$

a short polynomial of total degree $D$ in the momenta, whose coefficients $c_n$ depend on the angles
only. The **table** holds those coefficients at every angle node. It is filled once per run;
afterwards each evaluation reads a few numbers and rebuilds the trace for whatever $l$ and $p$ the
call brings.

The rebuild never forms $l^n$ or $p^{D-n}$ on its own. With $s = \max(l, p)$ and
$r = \min(l, p)/s \le 1$,

$$
N = s^D \sum_n c_n\, r^{\,n} \quad (l \le p), \qquad N = s^D \sum_n c_n\, r^{\,D-n} \quad (l > p),
$$

evaluated by Horner's rule in $r$. Each atom is multiplied in as the dimensionless
$s^{w_j}\,\mathrm{inv}_j$, where $w_j$ is the degree of its denominator (2 for $1/(l-p)^2$). The
leftover power of $s$ comes last, one factor at a time, so the running product overflows only if
the result does; a `float` kernel survives scales where $s^D$ alone would overflow. The Results
section shows a second consequence: the rebuilt trace is also *more accurate* than the inline one
when $l \ll p$.

### What you have to provide

Nothing new in the network. Tables are an emission option, `"TraceTables"`, and they work for any
flow that satisfies three conditions:

* **Vacuum, real frame.** A finite-temperature flow (`"MatsubaraVar"` set) or a complex frame is not
  tabulated.
* **The argument order.** The first kernel argument is the radial loop variable. The external
  coordinates (by default the argument named `p`; set `"CoordinateArgs"` otherwise) follow the
  integration variables directly. That is the order DiFfRG uses anyway.
* **At most one external scale.** NumTracer finds it by checking which external coordinate makes
  every frame vector homogeneous together with the radial variable. `<|"Scales" -> {"S0"}|>` names
  it explicitly, `<|"Scales" -> {}|>` declares that there is none. Two dimensionful external
  coordinates would need two-variable tables, which do not exist yet.

The table axes are everything else: the angles, and any external coordinate that is not a scale
(a shape parameter such as an opening angle between external legs).

## The commented program

`Tutorials/step-23-trace-tables/trace_tables.wls` and `check_trace_tables.cpp`.

### The frame and the network

```{literalinclude} ../../../Tutorials/step-23-trace-tables/trace_tables.wls
:language: mathematica
:start-after: "@snip begin: frame"
:end-before: "@snip end: frame"
```

```{literalinclude} ../../../Tutorials/step-23-trace-tables/trace_tables.wls
:language: mathematica
:start-after: "@snip begin: net"
:end-before: "@snip end: net"
```

### The option

```{literalinclude} ../../../Tutorials/step-23-trace-tables/trace_tables.wls
:language: mathematica
:start-after: "@snip begin: emit"
:end-before: "@snip end: emit"
```

`"TraceTables"` takes four forms:

| Value | Meaning |
|---|---|
| `False` | No table (the default of `MakeNTKernel`). |
| `True` | Tabulate every trace; stop with `MakeNTKernel::ttoff` if the flow cannot be tabulated. |
| `Automatic` | Tabulate a trace only where that pays, and keep the flow inline if it cannot be tabulated (with `$NumTracerVerbose = True` the log says why). The default of `MakeNTKernelDiFfRG`. |
| `<\|"Balance" -> B, "Scales" -> {…}\|>` | Explicit settings. A trace is tabulated when evaluating it inline costs more than `B` flops per byte of its table entries (`Automatic` uses `B = 1`; `True` is `B = 0`). |

The balance is decided at generation time, per trace. A trace that is not tabulated keeps its inline
code, and the table version of the kernel simply calls it.

### What the generator emits

The traces header gains, for each tabulated trace, a fill and a rebuild. This is the projector
diagram (abbreviated from `build/gen/trace_tables/Tab_kernels.hh`):

```cpp
static inline void tab_fill1(double *out, const std::size_t stride, double l1, double cos1, double p) {
  double f[1];
  f[0] = cos1;
  const double s0 = f[0];
  const double s3 = NT_FMAD(s0, s0, (-1));
  const double s4 = s0*s3;
  const double s5 = -s3;
  out[0 * stride] = s4;          // c_1 = -cos1 (1 - cos1^2)
  out[1 * stride] = s5;          // c_2 =        1 - cos1^2
}

static inline double ttr1(const double *row, const std::size_t st, const double *f) {
  const double l_ = f[0];        // l1
  const double e_ = f[2];        // p
  const bool le_ = l_ <= e_;
  const double s_ = le_ ? e_ : l_;
  const double r_ = (le_ ? l_ : e_) / s_;
  double h_;
  if (le_) {
    h_ = row[3 * st];
    h_ = NT_FMAD(h_, r_, row[2 * st]);
    h_ *= powr<1>(r_);
  } else {
    h_ = row[2 * st];
    h_ = NT_FMAD(h_, r_, row[3 * st]);
    h_ *= powr<2>(r_);
  }
  h_ *= (s_ * s_ * f[3]);        // f[3] = 1/(l - p)^2, the atom, still computed per call
  h_ = tab_scale<2>(h_, s_);
  return h_;
}
```

The traces share one table row: `tab_fill` writes the Dirac trace's two entries first and calls
`tab_fill1(out + 2 * stride, …)`, so `ttr1` finds its own as entries 2 and 3.

Note what the table does *not* contain. The fill reads only `cos1`. The transverse projector on $l$
came with a $1/l^2$, which the contraction cancelled exactly, so only $1/(l-p)^2$ is left; it is
an atom of the per-call `fill()`, outside the table.

The kernel class gains four members next to `kernel()`:

| Member | What it is |
|---|---|
| `has_tab`, `tab_size` | Whether anything is tabulated, and the number of entries per table node. |
| `tab_box_mask`, `tab_ext_mask` | One flag per integration variable and per external coordinate: `true` marks a table axis. Here `{false, true}` for `(l1, cos1)` and `{false}` for `p`. |
| `tab(out, stride, <kernel() args>)` | Fills one node: entry `e` goes to `out[e * stride]`. It reads only the axes, so the other arguments may hold anything. |
| `kernel_tab(row, stride, <kernel() args>)` | Returns `kernel(args)`, reading the traces from the node's entries `row[e * stride]`. |

A kernel with hoisted lookups (DiFfRG's `post()`) also gets `post_tab()`, in the same way.

### Filling and reading the table by hand

This is the part an integrator normally does for you; done by hand it shows the whole contract.

```{literalinclude} ../../../Tutorials/step-23-trace-tables/check_trace_tables.cpp
:language: cpp
:start-after: "@snip begin: traits"
:end-before: "@snip end: traits"
```

```{literalinclude} ../../../Tutorials/step-23-trace-tables/check_trace_tables.cpp
:language: cpp
:start-after: "@snip begin: fill"
:end-before: "@snip end: fill"
```

The stride is what lets the consumer choose the layout. Here the node index runs fastest, which is
also DiFfRG's choice on the GPU: neighbouring threads handle neighbouring nodes and read neighbouring
addresses.

```{literalinclude} ../../../Tutorials/step-23-trace-tables/check_trace_tables.cpp
:language: cpp
:start-after: "@snip begin: use"
:end-before: "@snip end: use"
```

```{literalinclude} ../../../Tutorials/step-23-trace-tables/check_trace_tables.cpp
:language: cpp
:start-after: "@snip begin: auto"
:end-before: "@snip end: auto"
```

## Results

```bash
cmake --build build --target trace_tables && ./build/trace_tables
```

```text
Tab_kernel: tab_size = 4 entries per node, axes: l1 0, cos1 1, p 0
the entries of the node cos1 = 0.5: -2 -2 -0.375 0.75
relative deviation from the closed form, worst over the 10 nodes
       p     l1/p   kernel_tab       kernel
     0.1    0.001      2.1e-16      4.5e-13
     0.1      0.1      2.6e-16      3.6e-15
     0.1        1      1.9e-16      6.8e-16
     0.1       10      2.2e-16      2.2e-16
     0.1     1000      2.4e-16      3.6e-16
       1    0.001      1.1e-16      6.0e-14
       1      0.1      1.4e-16      4.9e-16
       1        1      1.8e-17      1.4e-16
       1       10      5.9e-16      3.0e-16
       1     1000      3.0e-16      1.5e-16
      30    0.001      4.7e-16      4.7e-16
      30      0.1      1.2e-16      1.2e-16
      30        1      1.6e-16      8.2e-17
      30       10      3.4e-16      3.4e-16
      30     1000      3.4e-16      3.4e-16
checks
  the fill does not depend on l1 or p                  ok
  kernel_tab == closed form (1e-14)                    ok
  kernel == closed form (1e-11)                        ok
  Automatic: nothing tabulated, kernel_tab == kernel   ok
ALL TESTS PASSED
```

The roundoff digits depend on the compiler and the machine; the pattern does not.

**The entries are the closed form's coefficients.** At $c = 1/2$ the four numbers are, in order,

* the Dirac trace, $4 l^2 p^2 (2c^2 - 1) - 4 l p^3 c$: $c_1 = -4c = -2$ and $c_2 = 4(2c^2-1) = -2$;
* the projector diagram, $l p^2 (1-c^2)(l - c p)$ over $(l-p)^2$: $c_1 = -c(1-c^2) = -0.375$ and
  $c_2 = 1 - c^2 = 0.75$.

Four numbers per node, against two inline traces that each redo their arithmetic on every call.

**The table is the more accurate form.** `kernel_tab` stays at roundoff for every $l/p$. The inline
kernel loses digits for $l \ll p$: the projector trace there subtracts terms of order $p^2$ to leave
a result of order $l p$, which costs about $\log_{10}(p/l)$ digits. The rebuild never forms that
difference, because each power of $l/p$ is an explicit factor of $r$.

**`Automatic` declined both traces.** They are a handful of flops each, cheaper to evaluate than to
read from memory, so `Tab_auto_kernel` has `has_tab == false` and its `kernel_tab()` ignores the
row. That is the intended outcome for small traces. The traces of a real vertex flow are thousands
of flops each, and there `Automatic` tabulates most of them.

## Trace tables in DiFfRG

With [step-15](step-15.md)'s `MakeNTKernelDiFfRG` (and `MakeNTKernelGroupDiFfRG`) you do none of
the work above:

* `"TraceTables"` defaults to `Automatic` there. Pass `"TraceTables" -> False` to compare against the
  inline kernel.
* DiFfRG's quadrature integrators recognise a kernel with `has_tab` (the `kernel_has_tab` concept in
  `DiFfRG/physics/integration/tab.hh`). They build the table on first use, over the grid of the
  masked variables, and keep it while those grids and the masked external values stay the same. It
  never depends on $k$ or on the interpolated dressings.
* `map()` and `map_points()` then call `kernel_tab()` (or `post_tab()`). `get()` uses the table only
  if `tab_ext_mask` is all `false`, because otherwise each new external point would need a new table.
* Supported integrators: the vacuum momentum integrators `integrator_p2`, `integrator_p2_1ang`,
  `integrator_p2_4d_2ang` and `integrator_p2_4d_3ang`.

Measured per call, tables against the inline kernel:

| Flow | GPU | Speedup |
|---|---|--:|
| Nf2 four-gluon tadpole, `map()` | RTX 4070 / A100 / H100 | 8.0× / 4.9–5.2× / 6.4–6.9× |
| Yang–Mills full flow, per right-hand side | RTX 4070 / A100 / H100 | 1.50× / 1.46× / 1.57× |
| Nf2 quark–gluon vertex group, `map_points()` | RTX 4070 | 1.38× |
| Yang–Mills three-gluon vertex, `map_points()` | RTX 4070 | 1.13× |

The gain follows the share of the kernel's time spent in traces: in the tadpole the traces are
nearly all of it, while the three-gluon vertex spends most of its time in dressing lookups.

## Possibilities for extensions

1. **Move the balance.** Emit the `Automatic` variant with `"TraceTables" -> <|"Balance" -> B|>` for
   decreasing `B` and watch the `[tables]` log line. Find the `B` at which each trace switches to the
   table. `B = 0` gives the same kernel as `True`.

2. **Break the argument order on purpose.** Change `args` to `{cos1, l1, p}`. The first argument is
   now taken as the radial variable, and the explicit request fails with

   ```text
   MakeNTKernel::ttoff: "TraceTables" was requested, but this flow cannot be tabulated: no choice of
   external scale among {p} makes the frame vectors homogeneous together with cos1.
   ```

   With `Automatic` the same flow is generated inline; with `$NumTracerVerbose = True` the log says
   why.

3. **Two table axes.** Emit a three-point network in `sp3Frame` (step-07) with `"TraceTables" -> True`.
   `tab_box_mask` becomes `{false, true, true}` for `(l1, cos1, cos2)`. Fill a product grid of
   `cos1` and `cos2` nodes, one row per pair, and check `kernel_tab()` against `kernel()` again.

4. **Time it.** Call `kernel()` and `kernel_tab()` a few million times each in a loop, compiled
   with `-march=native`. For this network the table loses (on one laptop, 31 ms against 20 ms for
   2·10⁷ calls), the same verdict `Automatic` reached. Then drop `-march=native`: on x86 every
   `fma()` becomes a library call, and the inline kernel, which has more of them, falls behind
   (242 ms against 111 ms). Time with the flags you build with. Repeat with the four-point network of
   [step-14](step-14.md), where the traces are much larger.

## The plain program

```{literalinclude} ../../../Tutorials/step-23-trace-tables/trace_tables.wls
:language: mathematica
```

```{literalinclude} ../../../Tutorials/step-23-trace-tables/check_trace_tables.cpp
:language: cpp
```
