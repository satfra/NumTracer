# Glossary

Two vocabularies show up in this guide: the **engine's own terms** (which you need to use
NumTracer at all) and the **physics terms of the running example** (which are just the vocabulary
of *one* application — QCD/fRG — not concepts the engine knows about). They are separated below so
a reader from another domain can skip the second list.

## Engine terms

Network
: A product of tensor *heads* with index labels. Repeated labels are summed (Einstein); a label
  used once stays free. A network with no free labels is a scalar.

Index label
: The name of an index slot. In C++ it is typed — `LorentzIndex`, `AdjIndex`, `FundIndex` — and
  handed out by the frame (`F.indices<n>()`) or the SU($N$) group (`su.adjoint<n>()`,
  `su.fundamental<n>()`).

Head / builder
: One factor of a network — a metric, vector, projector, gamma, SU($N$) object. `nt…` in the
  Mathematica DSL; `nt::metric`, `nt::vec`, `nt::projT`, `nt::gamma`, `nt::slash`, `su.T`, … in C++
  (see the [DSL↔C++ dictionary](bring-your-own-network.md)).

Frame
: A choice of reference components for every vector in the network — e.g. one vector along an axis,
  another at a fixed angle. It fixes how each scalar product is written in terms of the kernel's
  runtime arguments. In C++, `nt::Frame` holds the symbols, the momenta and the projector
  denominators, and does the contraction (`F.trace`) and evaluation (`F.eval`). Full definition
  under [Key concepts](concepts.md#runtime-numbers-come-from-the-frame) and
  [internals Terminology](../internals/index.md#terminology).

Scalar symbol
: A scalar product (`l·p`, `l²`, …) kept as a named runtime quantity. The generated kernel computes
  each once per call; the lowered arithmetic is a polynomial in them.

Atom
: A surviving inverse propagator $1/k^2$ carried alongside a polynomial's monomials and evaluated at
  evaluation time. Named because it is an indivisible reciprocal the contraction tracks rather than
  expands. The frame registers one per projector momentum and knows its denominator.

Diagram
: One term of the network at the top level — a scalar coefficient times one contraction. A kernel
  is a sum over diagrams. (In the running example each diagram is a Feynman diagram, but the engine
  just sees "a coefficient times a contraction.")

Poly
: The engine's multivariate polynomial in the frame's scalar symbols, carrying inverse atoms. Every
  contraction returns one. (`DPoly` is the same with runtime dressings as extra variables.)

SUN
: An SU($N$) group object, `nt::SUN su(N)`: it hands out the group's typed labels, builds its
  factors (`f`, `T`, `delta`, `diag`) and folds a closed network of them to an exact number
  (`su.value`). One object per group — colour, flavour, or any other SU($N$).

Lowering (CSE + Horner)
: The build-time passes that turn a `Poly` into flat, straight-line real arithmetic — Horner
  factoring plus common-subexpression elimination. See [CSE and Horner lowering](../internals/cse-and-lowering.md).

## Physics terms, as used in the running example

These describe *what the example computes*, not the engine. If you are not doing QFT you can treat
them as opaque names.

Dressing
: A scalar function (of momentum/scale) multiplying a propagator or vertex — a runtime coefficient
  the kernel takes as an argument. To the engine it is just an opaque scalar.

Regulator ($\partial_t R$)
: The fRG infrared cutoff insertion; here, another scalar coefficient on a diagram.

fRG / Wetterich equation
: The functional Renormalization Group and its flow equation, the framework the example diagrams
  come from. NumTracer traces the diagrams; it does not know the flow equation.

Propagator
: A field's two-point line in a diagram; its numerator is a tensor the engine traces, its scalar
  denominator a coefficient the user supplies.

Colour / flavour
: The two SU($N$) sectors of the example (SU($N_c$) gauge colour, SU($N_f$) flavour). To the engine
  they are just two SU($N$) groups, one `nt::SUN` object each.

Fundamental / adjoint
: The two SU($N$) representations the engine handles — the $N$-dimensional (quark-like) and
  $(N^2-1)$-dimensional (gluon-like) index spaces.

Transverse / longitudinal / electric / magnetic projector
: Lorentz projectors onto subspaces relative to a momentum: transverse
  $P^T = \delta - kk/k^2$, longitudinal $P^L = kk/k^2$, and (at finite $T$) the electric/magnetic
  split of $P^T$. See [step-16](../tutorials/step-16.md).

Spinor / Dirac trace / Clifford algebra
: The spin-½ index space (dimension 4), a closed loop of gamma matrices over it, and the relation
  $\{\gamma^\mu,\gamma^\nu\}=2\delta^{\mu\nu}$ they satisfy.

Matsubara / finite temperature / heat bath
: The imaginary-time thermal formalism, in which a rest frame singles out the temporal component of
  each momentum. Handled by the finite-$T$ frames and projectors, not by any change to the engine.

Condensate
: A field acquiring a direction-selective background — modelled here by dressing selected SU($N$)
  components differently (see [step-17](../tutorials/step-17.md)).
