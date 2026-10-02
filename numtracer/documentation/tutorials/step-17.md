# step-17: Per-flavour and per-component dressings

*Builds on: [step-02](step-02.md) · Built on by: — · Tags: `sun`, `dressing` ·
**Tier A** (a C++20 compiler, nothing else)*

## Introduction

[step-02](step-02.md) folded a closed SU($N$) loop to a number. That is the right answer when every
component of the multiplet is equivalent — and it is exactly what makes it useless when they are
not.

### The problem

A closed flavour loop collapses to $N_f$, the dimension of the multiplet. That single number can
carry a single dressing: whatever propagator dressing the quark has, applied to all flavours alike.
But physics regularly wants otherwise:

* **Broken isospin.** The up and down quarks have different masses, so $Z_u \neq Z_d$.
* **A gluon condensate.** A background field along the Cartan directions of SU(3) dresses the
  $\lambda_3$ and $\lambda_8$ gluons differently from the other six.
* **Any partial-symmetry-breaking pattern** where the group index is no longer a spectator.

The obvious fix is to split the diagram: one copy per component, each with its own dressing. That
works and it is expensive — the Dirac and Lorentz trace, which is the costly part and is *identical*
in every copy, gets recomputed $N$ times.

### The mechanism

NumTracer instead makes the loop's Kronecker delta **group-diagonal**: instead of $\delta_{ij}$,
$\operatorname{diag}(D_1,\dots,D_N)$ with a distinct dressing per component. The same fold then
gives, not a number, but a **polynomial** over independently named runtime dressings:

$$
\delta_{ij}\ \longrightarrow\ \operatorname{diag}(D_a)
\qquad\Longrightarrow\qquad
\sum_i \delta_{ii} = N \;\longrightarrow\; \sum_a D_a .
$$

**No diagram split.** The Dirac/Lorentz trace multiplying it is computed once, exactly as before.
The only thing that changed is what the colour half folds *to*.

This is the SU($N$)-side analogue of [step-18](step-18.md)'s dressed Dirac numerators, and the two
use the same idea: keep the sum eager, carry the dressings as opaque atoms, and fold to a
polynomial instead of a scalar.

### The two entry points

| Function | Returns | Use |
|---|---|---|
| `su.value(net)`, `nt::sun_value(net)` | one `Cx` | a fully contracted net with no diagonal factor ([step-02](step-02.md)) |
| `nt::sun_value_dressed(net)` | a `SUNPoly` | a net containing a `su.diag(…)` factor |

A `SUNPoly` is a list of `SUNTerm{coeff, dress}`: a coefficient times the product of the
dressing-ids in that monomial. A net with **no** diagonal factor comes back as a single constant
term equal to what `sun_value` would have given — the two agree by construction.

### The `comp2dr` map

A diagonal factor, `su.diag(i, j, comp2dr)` (fundamental) or `su.diag(a, b, comp2dr)` (adjoint), is
a delta tagged with a component-to-dressing map. `comp2dr[v]` is the
dressing-id for component `v`, **or `-1` to drop that component entirely**. Dropping is not the same
as dressing with zero at runtime: a dropped component contributes no term at all, so there is no
dead arithmetic in the emitted kernel.

```{admonition} Components are 0-based here
:class: important
The physics numbers fundamental components $1..N$ and adjoint components $1..N^2-1$. The C++ index
is one less. So the SU(3) Cartan directions $\lambda_3$ and $\lambda_8$ are `comp2dr[2]` and
`comp2dr[7]`. The Mathematica DSL heads (`ntSUNDiagFund`, `ntSUNDiagAdj`) take **1-based** component
indices in their `spec` rules-list, which is the more natural spelling for physics — do not carry
one convention into the other.
```

## The commented program

`Tutorials/step-17-dressed-flavour/dressed_flavour.cpp` validates the two invariants the physics
rests on.

### A. Fundamental: an isospin doublet

```{literalinclude} ../../../Tutorials/step-17-dressed-flavour/dressed_flavour.cpp
:language: cpp
:start-after: "@snip begin: fund"
:end-before: "@snip end: fund"
```

`su2.diag(i, j, {0, 1})` is the flavour delta with component 0 (the up quark) carrying dressing-id 0
and component 1 (the down quark) carrying id 1; `su2.delta(j, i)` closes the loop. The fold gives
the two-term `SUNPoly` $D_u + D_d$. The labels `i, j` are `FundIndex` from `su2.fundamental<2>()`,
so `diag` picks the fundamental version; the adjoint labels of section B pick the adjoint one.

Then the two invariants:

* **COLLAPSE.** Setting all dressings equal to 1 must recover the flavour-blind number the plain
  fold gives — here $N_f = 2$. If it does not, the diagonal fold is not a refinement of the ordinary
  one and everything built on it is suspect.
* **DROP.** `{0, -1}` dresses only the up quark. The result must be $D_u$ alone.

### B. Adjoint: a gluon condensate on the Cartan

```{literalinclude} ../../../Tutorials/step-17-dressed-flavour/dressed_flavour.cpp
:language: cpp
:start-after: "@snip begin: adj"
:end-before: "@snip end: adj"
```

Eight adjoint components, six of them dropped. Note the third check: the Cartan-only polynomial
equals the *fully* dressed polynomial evaluated with the other six dressings set to zero. That is
the statement that **dropping a component is exactly zero-dressing it** — the optimisation is
semantically transparent, which is what licenses using it.

## Results

```bash
cmake --build build --target dressed_flavour && ./build/dressed_flavour
```

```{literalinclude} ../../../Tutorials/step-17-dressed-flavour/dressed_flavour.expected.txt
:language: text
```

Read line 3 of block A and line 2 of block B together: `2` and `8` are $N_f$ and $N_c^2-1$, the
dimensions [step-02](step-02.md) would have produced. The diagonal fold **contains** the blind fold
as its all-equal special case, and that is the property that makes it safe to switch on
unconditionally.

The last line of block B is the one that costs nothing to check and would be expensive to get wrong.

## The front-end spelling

```mathematica
ntSUNDiagFund[N, i, j, spec]
ntSUNDiagAdj [N, a, b, spec]
```

`spec` is a **rules list** with 1-based component indices:

```mathematica
(* u and d dressed differently inside one SU(2) flavour trace *)
ntSUNDiagFund[2, i, j, {1 -> Zu[scale], 2 -> Zd[scale]}]

(* a gluon condensed along the SU(3) Cartan; the other six colours drop out with no dead terms *)
ntSUNDiagAdj[3, a, b, {3 -> A03[scale], 8 -> A08[scale]}]
```

Three details:

* Each `expri` is a **complete scalar dressing expression**, applied to its own kinematics by the
  caller — the same opaque-callable mechanism as [step-08](step-08.md). Components may use different
  scales. Kinematic heads (`ntVec`/`ntSP`/`ntSPS`) inside an `expri` are frame-resolved at kernel
  generation, so they expand per diagram.
* A `Default -> defExpr` rule dresses every unnamed component. A pure PROJECTOR component — one
  that pins an index rather than carrying physics — is simply `c -> 1`.
* Components with **neither a rule nor a `Default` are dropped**. That is a convenience and a trap:
  a typo in a component index silently removes that component instead of erroring.

## Possibilities for extensions

1. **Two groups, one dressed.** Put a diagonal SU(2) flavour delta and a plain SU(3) colour loop in
   the same network. The result should be $3 \times (D_u + D_d)$ — the colour half folds to a
   number and multiplies the flavour polynomial.

2. **Two diagonal factors.** Give both the colour and the flavour delta a `comp2dr` map. The fold
   should produce *products* of dressing ids — a `SUNTerm` whose `dress` list has two entries.
   Confirm the term count is $N_c \times N_f$ and that the collapse still gives $N_c N_f$.

3. **Drop everything.** Pass `comp2dr` all `-1`. Predict the result before running it: an empty
   `SUNPoly`, not a polynomial that evaluates to zero. Check `poly.size()`.

4. **The off-diagonal question.** A `diagFund` is diagonal by construction. Work out what you would
   need if a symmetry-breaking pattern mixed components — and convince yourself that is a *different*
   object requiring a full matrix, not a knob on this one.

5. **Measure the saving.** Build the same trace as $N$ separate diagrams (the naive split) and
   compare the total Dirac/Lorentz work against the single diagonal fold. The ratio is $N$, and it
   is why this mechanism exists.

## The plain program

```{literalinclude} ../../../Tutorials/step-17-dressed-flavour/dressed_flavour.cpp
:language: cpp
```

Next, [step-18](step-18.md): the same "keep it eager, fold to a polynomial" idea, on the Dirac side.
