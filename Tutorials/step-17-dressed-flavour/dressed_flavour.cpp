// step-17 — Per-flavour and per-component dressings: the group-diagonal fold.
//
// The flows so far dress every quark the same way (one propagator dressing for the whole
// colour/flavour multiplet). Physics often needs more: the u and d quarks dressed differently
// (broken isospin), or a single colour direction dressed on its own (a gluon condensate).
//
// The mechanism is the SU(N) fold. A closed SU(N) loop normally collapses to one
// flavour-BLIND number (its dimension N) — which can carry only ONE dressing. Making the loop's
// delta a *group-diagonal dressing* — diag(D_1, ..., D_N) instead of a plain delta — folds the
// same loop to a POLYNOMIAL  Sum_a D_a  over independently-named runtime dressings, WITHOUT
// splitting the trace into one diagram per component. The Dirac/Lorentz trace multiplying it is
// still computed once.
//
// Two entry points (network/sun_net.hpp):
//   sun_value(net)         -- plain fold: a fully-contracted net -> one number (Cx).
//   sun_value_dressed(net) -- diagonal fold: a net with a diag(...) factor -> a SUNPoly,
//                             i.e. Sum_t coeff_t * Prod D^dr. Each term is SUNTerm{coeff, dress},
//                             where `dress` lists the dressing-ids in that monomial. A net with no
//                             diagonal factor comes back as a single constant term == sun_value.
//
// A diag factor is a delta tagged with a comp2dr map: comp2dr[v] is the dressing-id for component
// v, or -1 to DROP that component (it contributes nothing — no dead terms). Components are 0-based
// here (1..N fundamental, 1..N^2-1 adjoint in the physics, minus one for the index).
//
// We validate the two invariants the physics rests on: the COLLAPSE (all dressings equal recovers
// the flavour-blind number the plain fold gives) and the DROP (a -1 component == dressing it with
// a zero function).
#include <numtracer.hpp> // the whole NumTracer API

#include <cmath>
#include <cstdio>
#include <vector>

namespace nt = numtracer;
using nt::Cx;

// Evaluate a SUNPoly at a dressing-id -> value assignment D:  Sum_t coeff_t * Prod_{id in dress} D(id).
static Cx eval_poly(const nt::SUNPoly &p, double (*D)(int)) {
  Cx s{0, 0};
  for (const nt::SUNTerm &t : p) {
    Cx c = t.coeff;
    for (int id : t.dress) c = c * Cx{D(id), 0.0};
    s = s + c;
  }
  return s;
}

int main() {
  bool ok = true;

  // One group object per SU(N): it hands out the labels of the closed delta loop.
  nt::SUN su2(2); // SU(2) flavour doublet (section A)
  nt::SUN su3(3); // SU(3) colour adjoint (section B)
  auto [i, j] = su2.fundamental<2>();
  auto [a, b] = su3.adjoint<2>();

  // ---- A. Fundamental: a u/d isospin doublet (SU(2) flavour) --------------------------------
  //
  // A closed flavour delta-loop over the doublet, its delta made flavour-DIAGONAL: component 0 (u)
  // carries dressing-id 0, component 1 (d) carries id 1. su2.diag(i, j, comp2dr) is that delta;
  // su2.delta(j, i) closes the loop. The fold gives the SUNPoly  D_u + D_d.
  // @snip begin: fund
  const nt::SUNPoly ud = nt::sun_value_dressed(su2.diag(i, j, {0, 1}) * su2.delta(j, i));

  auto broken = [](int id) { return id == 0 ? 2.0 : 5.0; }; // D_u = 2, D_d = 5 (broken isospin)
  auto ones = [](int) { return 1.0; };                      // all dressings equal to 1

  const Cx ud_broken = eval_poly(ud, broken); // 2 + 5 = 7
  const Cx ud_blind = eval_poly(ud, ones);    // 1 + 1 = 2 = N_f  (the COLLAPSE)

  // DROP: the same loop dressing ONLY component 0 (id -1 drops component 1) folds to just D_u.
  const nt::SUNPoly u_only = nt::sun_value_dressed(su2.diag(i, j, {0, -1}) * su2.delta(j, i));
  const Cx u_drop = eval_poly(u_only, broken); // D_u = 2
  // @snip end: fund

  std::printf("A. fundamental u/d doublet (SU(2) flavour)\n");
  std::printf("   D_u + D_d           = %g   (D_u=2, D_d=5 -> 7)\n", ud_broken.re);
  std::printf("   collapse D_u=D_d=1  = %g   (-> flavour-blind N_f = 2)\n", ud_blind.re);
  std::printf("   drop d (only u)     = %g   (-> D_u = 2)\n", u_drop.re);
  ok = ok && nt::approx(ud_broken, Cx{7, 0}) && nt::approx(ud_blind, Cx{2, 0}) && nt::approx(u_drop, Cx{2, 0});

  // ---- B. Adjoint: a gluon condensate on the Cartan directions (SU(3) colour) ---------------
  //
  // The adjoint has N^2-1 = 8 components. su3.diag(a, b, comp2dr) dresses them individually. A
  // condensate lives in the CARTAN directions -- lambda_3, lambda_8 for SU(3), i.e. components 3
  // and 8 (0-based indices 2 and 7). We dress those two and DROP the other six.
  // @snip begin: adj
  std::vector<int> cartan(8, -1); // start with everything dropped
  cartan[2] = 0;                  // lambda_3 -> dressing-id 0  (Z_3)
  cartan[7] = 1;                  // lambda_8 -> dressing-id 1  (Z_8)
  const nt::SUNPoly cond = nt::sun_value_dressed(su3.diag(a, b, cartan) * su3.delta(b, a));
  const Cx cond_blind = eval_poly(cond, ones); // Z_3 + Z_8 at 1 -> 2

  // The colour-BLIND adjoint loop, every component dressed by one id, is the collapse baseline:
  // all-equal -> N^2-1 = 8 (what the plain sun_value delta-loop would give).
  std::vector<int> all8(8);
  for (int c = 0; c < 8; ++c) all8[c] = c;
  const nt::SUNPoly full = nt::sun_value_dressed(su3.diag(a, b, all8) * su3.delta(b, a));
  const Cx full_blind = eval_poly(full, ones); // 8

  // DROP == zero-dressing: the Cartan-only poly equals the full poly with the other six dressings
  // set to zero. Here: full evaluated with D(id)=1 only for the Cartan ids {2,7}, else 0 -> 2.
  auto cartan_mask = [](int id) { return (id == 2 || id == 7) ? 1.0 : 0.0; };
  const Cx full_masked = eval_poly(full, cartan_mask); // 2, matching cond_blind
  // @snip end: adj

  std::printf("B. adjoint gluon condensate on the Cartan (SU(3) colour)\n");
  std::printf("   full loop, all Z=1  = %g   (-> N^2-1 = 8)\n", full_blind.re);
  std::printf("   Cartan {3,8}, Z=1   = %g   (only lambda_3, lambda_8 survive -> 2)\n", cond_blind.re);
  std::printf("   full with rest=0    = %g   (drop == zero-dressing -> matches Cartan)\n", full_masked.re);
  ok = ok && nt::approx(full_blind, Cx{8, 0}) && nt::approx(cond_blind, Cx{2, 0}) &&
       nt::approx(full_masked, cond_blind);

  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
