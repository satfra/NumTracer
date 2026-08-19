/// @file numeric_contract.hpp
/// @brief The numeric (matrix-product) contraction engine: fold a diagram's Dirac trace
///        by 4×4 spinor matrix **products** over @ref MPoly entries (γ numeric, momenta symbolic),
///        then contract the surviving free gluon legs against the pure-Lorentz network (projectors /
///        metrics / vectors / Levi-Civita) by a bounded index sum. Neither step suffers the
///        `(2n−1)!!` Wick blowup or the `2^np` projector-mask blowup of the symbolic `reduce` path.
///
/// Inputs per diagram:
///   - a closed Dirac chain @ref network::DiracNet (free legs `dgamma(mu)`, slashed propagators
///     `dslash(vlc)`); reused verbatim from the inv backend so the Mathematica front-end is unchanged.
///   - a pure-Lorentz network @ref NNet (metrics / vectors / transverse projectors / Levi-Civita)
///     over the same Lorentz ids.
///   - a **component table** `comp[vid]` = the 4 @ref MPoly components of fundamental momentum `vid`
///     (user-supplied, partially numeric / partially symbolic), and the projector denominators
///     `atomDen[aid] = k²` for monomial-cancellation bookkeeping.
///
/// Output: the diagram's scalar trace as one @ref MPoly (surviving `1/k²` atoms in its monomials),
/// ready for @ref to_genprog → the shared `gdetail::best_into` CSE/Horner emission.
#pragma once

#include "numtracer/core/export.hpp"   // NUMTRACER_FUNC / NUMTRACER_DEFINE_BODIES (compiled vs header-only)
#include "numtracer/core/envvar.hpp" // env_flag / env_int — the single truth test for NT_* switches
#include "numtracer/core/config.hpp"    // NT_THROW (exception-optional guard for -fno-exceptions builds)
#include "numtracer/codegen/gen.hpp"   // network::GlobalEnv / GenProg / LMono / gdetail::best_into
#include "numtracer/network/dirac.hpp" // network::DiracNet / DFac (reused chain representation)
#include "numtracer/numeric/dpoly.hpp"      // DPoly / DMono (dressing-atom layer)
#include "numtracer/numeric/mpoly.hpp"      // MPoly (the polynomial type)
#include "numtracer/numeric/spinor_mat.hpp" // Mat4 / matmul / mtrace / gammaC / slashC

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string> // the open-index guard's diagnostic (assert_no_open_ids)
#include <vector>

namespace numtracer::numeric
{


  /// @brief A pure-Lorentz network factor (numeric counterpart of @ref network::Elem), tagged by
  ///        @ref NElem::Kind and built by the nXxx() helpers below:
  ///   - `Metric`  `nmet(a,b)`                     — δ_{ab}
  ///   - `Vector`  `nvec(a, vlc)`                  — vector leg `Σ coeff·comp(vid)` on Lorentz id `a`
  ///   - `Epsilon` `neps(a,b,c,d)`                 — Levi-Civita ε_{abcd}
  ///   - `ProjT`   `nprojT(a,b, vlc, atom)`        — transverse `P_T(k)_{ab}=δ_{ab}−k_a k_b·INV(k)`
  ///   - `ProjL`   `nprojL(a,b, vlc, atom)`        — longitudinal `P_L(k)_{ab}=k_a k_b·INV(k)`
  ///   - `ProjE`   `nprojE(a,b, vlc, atom, atomS)` — finite-T electric `P_E = P_T − P_M`
  ///   - `ProjM`   `nprojM(a,b, vlc, atomS)`       — finite-T magnetic `P_M_{ij}=δ_{ij}−k_i k_j·INVS(k)`
  ///                                                 (i,j spatial; the temporal row/col 0 vanishes)
  /// where `k = Σ coeff·comp(vid)`, `INV(k)=1/k²` is inverse atom `atom`, `INVS(k)=1/|k⃗|²` is `atomS`.
  ///
  /// Members are non-const so NElem stays movable in its std::vector (const members would delete
  /// move-assignment and copy the vlc vector on every reallocation); the nXxx() builders are the only
  /// constructors, so the values are still effectively immutable in practice.
  struct NElem {
    enum Kind { Metric, Vector, Epsilon, ProjT, ProjL, ProjE, ProjM };
    Kind kind = Metric;
    /// Lorentz index ids, named by the variant that uses them (never by enum ordinal):
    ///   Metric      δ_{a b}                 — a, b
    ///   Vector      q_a                     — a
    ///   ProjT/L/E/M P_{a b}(k)              — a, b
    ///   Epsilon     ε_{a b c d}             — a, b, c, d (the only kind using c and d)
    /// Two factors contract exactly when they share an id, so these are the network's edges.
    int a = 0, b = 0, c = 0, d = 0;
    std::vector<std::pair<double, int>> vlc; ///< momentum (Vector leg, or any projector's `k`)
    int atom = -1;                           ///< full inverse-atom id `1/k²` (ProjT / ProjL / ProjE)
    int atomS = -1;                          ///< spatial inverse-atom id `1/|k⃗|²` (ProjE / ProjM)
  };
  struct NTerm {
    Cx coeff{1, 0};
    std::vector<NElem> e;
  };
  using NNet = std::vector<NTerm>;

  inline NElem nmet(int a, int b) { return {NElem::Metric, a, b, 0, 0, {}, -1}; }
  inline NElem nvec(int a, std::vector<std::pair<double, int>> vlc)
  {
    return {NElem::Vector, a, 0, 0, 0, std::move(vlc), -1};
  }
  inline NElem nprojT(int a, int b, std::vector<std::pair<double, int>> vlc, int atom)
  {
    return {NElem::ProjT, a, b, 0, 0, std::move(vlc), atom};
  }
  inline NElem nprojL(int a, int b, std::vector<std::pair<double, int>> vlc, int atom)
  {
    return {NElem::ProjL, a, b, 0, 0, std::move(vlc), atom};
  }
  /// Finite-T electric projector `P_E = P_T − P_M`: `atom`=1/k², `atomS`=1/|k⃗|².
  inline NElem nprojE(int a, int b, std::vector<std::pair<double, int>> vlc, int atom, int atomS)
  {
    return {NElem::ProjE, a, b, 0, 0, std::move(vlc), atom, atomS};
  }
  /// Finite-T magnetic projector `P_M_{ij}=δ_{ij}−k_i k_j/|k⃗|²` (spatial only): `atomS`=1/|k⃗|².
  inline NElem nprojM(int a, int b, std::vector<std::pair<double, int>> vlc, int atomS)
  {
    return {NElem::ProjM, a, b, 0, 0, std::move(vlc), -1, atomS};
  }
  inline NElem neps(int a, int b, int c, int d) { return {NElem::Epsilon, a, b, c, d, {}, -1}; }

  namespace ndetail
  {

    /// The 4 components of momentum `Σ coeff·comp(vid)`.
    inline std::array<MPoly, 4> mom_components(int nsym, const std::vector<std::pair<double, int>> &vlc,
                                               const std::vector<std::array<MPoly, 4>> &comp)
    {
      std::array<MPoly, 4> r = {MPolyFactory::zero(nsym), MPolyFactory::zero(nsym), MPolyFactory::zero(nsym), MPolyFactory::zero(nsym)};
      // `scaled`, not `constant(coeff) * cv[mu]`: same coefficient product in the same operand order,
      // same monomials, no scratch and no sort. Four of these per `vlc` entry, and this runs for every
      // Slash token and every projector/vector `elem_factor` — one of the hottest sites in the engine.
      // The accumulate moves too: on the first entry `r[mu]` is empty, which is exactly the empty-side
      // deep copy the rvalue `operator+` lever exists to kill.
      for (const auto &[coeff, vid] : vlc) {
        const auto &cv = comp[vid];
        for (int mu = 0; mu < 4; ++mu)
          r[mu] = std::move(r[mu]) + MPolyFactory::scaled(nsym, cv[mu], Cx{coeff, 0});
      }
      return r;
    }

    /// Totally antisymmetric ε_{abcd}, ε_{0123}=+1.
    inline double levi(int a, int b, int c, int d)
    {
      int idx[4] = {a, b, c, d};
      double sgn = 1.0;
      for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j) {
          if (idx[i] == idx[j]) return 0.0;
          if (idx[i] > idx[j]) sgn = -sgn;
        }
      return sgn;
    }

    /// A dense Lorentz factor: a tensor over `ids` (each extent 4), row-major flat `v` of MPoly.
    struct Factor {
      std::vector<int> ids;
      std::vector<MPoly> entries; ///< the dense tensor, flattened: size 4^ids.size()
    };

  } // namespace ndetail

  // Public engine entry points: declared always, defined once (in the library TU, or inline in a
  // header-only build). See core/export.hpp. Their heavy bodies live in the NUMTRACER_DEFINE_BODIES
  // regions below so a normal consumer/generator TU only parses the declarations and links the lib.
  NUMTRACER_FUNC ndetail::Factor numeric_dirac(int nsym, const network::DiracNet &chain,
                                               const std::vector<std::array<MPoly, 4>> &comp);

  /// @brief Fold the closed Dirac chain into a tensor over its free gluon legs (the `dgamma` ids).
  ///        For each assignment of the free legs to concrete indices 0..3, build the slashed/free
  ///        γ chain as 4×4 @ref MPoly matrices and take the trace. Returns the free-leg ids and the
  ///        row-major tensor of trace polynomials (one entry per `4^f` assignment).
#if NUMTRACER_DEFINE_BODIES
  NUMTRACER_FUNC ndetail::Factor numeric_dirac(int nsym, const network::DiracNet &chain,
                                               const std::vector<std::array<MPoly, 4>> &comp)
  {
    // Algorithm:
    //   1. Walk the chain once: record which tokens are FREE γ legs (open ids, summed below), build the
    //      fixed slash matrices, and split the commutator legs into free/slash. Precompute every factor's
    //      two Weyl 2×2 blocks (P = upper-right, Q = lower-left); each γ^μ and slash is block-antidiagonal.
    //   2. If the count of block-antidiagonal factors (γ + slash) is odd, the product is antidiagonal and
    //      the trace vanishes — return an all-zero tensor.
    //   3. For each of the 4^f assignments of the free legs to concrete indices (`legComp`), multiply the
    //      running product as two 2×2 blocks (`m0`/`m1`) whose upper/lower roles alternate with each
    //      antidiagonal factor (`antidiag` flag). γ5 and bare commutators are block-DIAGONAL: they never
    //      flip the parity and cost no full multiply (γ5 just negates one block). The trace of the (even)
    //      block-diagonal product is tr(m0)+tr(m1).
    // collect free-leg ids (in chain order) and precompute the per-token matrices that don't vary.
    std::vector<int> freeLegs;
    std::vector<int> tokenFree(chain.size(), -1); // token index → free-leg slot, or -1
    std::vector<int> tokenFree2(chain.size(),
                                -1); // kind 3 (commutator): leg-B free-leg slot (or -1 if leg-B is a slash)
    std::vector<Mat4> slashMat;      // precomputed Slash matrices, indexed by token
    slashMat.reserve(chain.size());
    // kind 3 commutator [A,B]: each leg is FREE (a looped component) or SLASH (a fixed momentum matrix).
    // For a slashed leg we stash its slash matrix here and turn it into 2×2 blocks below (next to sP/sQ).
    std::vector<Mat4> commA(chain.size(), Mat4(nsym)), commB(chain.size(), Mat4(nsym));
    std::vector<char> commAslash(chain.size(), 0), commBslash(chain.size(), 0);
    for (std::size_t i = 0; i < chain.size(); ++i) {
      const network::DFac &d = chain[i];
      if (d.kind == network::DFac::Gamma) {
        tokenFree[i] = freeLegs.size();
        freeLegs.push_back(d.mu);
        slashMat.emplace_back(nsym); // placeholder
      } else if (d.kind == network::DFac::Slash) {
        slashMat.push_back(slashC(nsym, ndetail::mom_components(nsym, d.vlc, comp)));
      } else if (d.kind ==
                 network::DFac::Comm) { // commutator [A,B]: leg-A then leg-B, each free (open id) or slash (momentum)
        if (d.mu >= 0) {
          tokenFree[i] = freeLegs.size();
          freeLegs.push_back(d.mu);
        } else {
          commAslash[i] = 1;
          commA[i] = slashC(nsym, ndetail::mom_components(nsym, d.vlc, comp));
        }
        if (d.nu >= 0) {
          tokenFree2[i] = freeLegs.size();
          freeLegs.push_back(d.nu);
        } else {
          commBslash[i] = 1;
          commB[i] = slashC(nsym, ndetail::mom_components(nsym, d.vlc2, comp));
        }
        slashMat.emplace_back(nsym); // placeholder
      } else {                       // kind 2: γ5 — block-diagonal in the Weyl basis, folded in below without a matrix.
        slashMat.emplace_back(nsym); // placeholder (keeps slashMat aligned with the token index)
      }
    }
    const int f = freeLegs.size();
    NT_STAT_ADD(nd_calls, 1);
    NT_STAT_ADD(nd_tokens, chain.size());
    // Trace parity is set by the BLOCK-ANTIDIAGONAL factors only — Gamma and Slash. Gamma5 and Comm
    // (σ) are block-DIAGONAL, so they do not flip it: a Comm is two gammas and a LoopSep is none,
    // both 0 mod 2. An odd antidiagonal count ⇒ the final product is antidiagonal ⇒ tr = 0.
    // Keep this count in step with its sibling in network/dirac.hpp; the two engines must agree.
    std::size_t nAntidiag = 0;
    for (const network::DFac &d : chain)
      if (d.kind == network::DFac::Gamma || d.kind == network::DFac::Slash) ++nAntidiag;
    ndetail::Factor F;
    F.ids = freeLegs;
    // The free-leg tensor has 4^f entries; `total` is a 32-bit int, so f >= 16 (4^16 = 2^32)
    // would silently overflow and under-size F.entries, corrupting the later index writes. Such a chain
    // is astronomically large anyway — refuse it loudly rather than miscompute.
    if (f > 15)
      NT_THROW(std::runtime_error, "numeric_dirac: Dirac chain has too many free Lorentz legs "
                                   "(4^f index space overflows a 32-bit int at f>=16)");
    int total = 1;
    for (int k = 0; k < f; ++k)
      total *= 4;
    F.entries.assign(total, MPolyFactory::zero(nsym));
    if (nAntidiag % 2 == 1) {
      NT_STAT_ADD(nd_odd_skip, 1);
      return F; // odd → all zero
    }
    NT_STAT_ADD(nd_assign, total);

    // Precompute each factor's Weyl 2×2 blocks (P = upper-right, Q = lower-left) once. The fold carries
    // only these two blocks rather than the full 4×4 (4× fewer MPoly multiplies); see the algorithm note.
    using B2 = std::array<MPoly, 4>; // row-major 2×2
    auto blocksOf = [&](const Mat4 &M, B2 &P, B2 &Q) {
      for (int rr = 0; rr < 2; ++rr)
        for (int cc = 0; cc < 2; ++cc) {
          P[rr * 2 + cc] = M.entries[rr][cc + 2];
          Q[rr * 2 + cc] = M.entries[rr + 2][cc];
        }
    };
    auto mul2 = [&](const B2 &x, const B2 &y) {
      NT_STAT_ADD(mul2_calls, 1);
      return B2{x[0] * y[0] + x[1] * y[2], x[0] * y[1] + x[1] * y[3], x[2] * y[0] + x[3] * y[2],
                x[2] * y[1] + x[3] * y[3]};
    };
    std::array<B2, 4> gP, gQ;
    for (int mu = 0; mu < 4; ++mu) {
      Mat4 g = gammaC(nsym, mu);
      blocksOf(g, gP[mu], gQ[mu]);
    }
    std::vector<B2> sP(chain.size()), sQ(chain.size());
    for (std::size_t i = 0; i < chain.size(); ++i)
      if (chain[i].kind == network::DFac::Slash) blocksOf(slashMat[i], sP[i], sQ[i]);
    // kind 3 commutator: 2×2 blocks of any SLASHED leg (the FREE legs index gP/gQ per assignment below).
    std::vector<B2> cAP(chain.size()), cAQ(chain.size()), cBP(chain.size()), cBQ(chain.size());
    for (std::size_t i = 0; i < chain.size(); ++i)
      if (chain[i].kind == network::DFac::Comm) {
        if (commAslash[i]) blocksOf(commA[i], cAP[i], cAQ[i]);
        if (commBslash[i]) blocksOf(commB[i], cBP[i], cBQ[i]);
      }
    // γ5 = diag(+I,−I) (Weyl): block-diagonal, so it costs no multiply and never flips the parity — it just
    // negates one running block (upper when the product is antidiagonal, lower when diagonal). A leading γ5
    // seeds the product as diag(+I,−I).
    const B2 id2 = {MPolyFactory::constant(nsym, Cx{1, 0}), MPolyFactory::zero(nsym), MPolyFactory::zero(nsym), MPolyFactory::constant(nsym, Cx{1, 0})};
    auto neg2 = [&](const B2 &x) {
      return B2{MPolyFactory::zero(nsym) - x[0], MPolyFactory::zero(nsym) - x[1], MPolyFactory::zero(nsym) - x[2], MPolyFactory::zero(nsym) - x[3]};
    };
    // BARE commutator [A,B] = A·B − B·A is block-diagonal (A,B antidiagonal ⇒ A·B diagonal), like γ5 so it
    // never flips the parity: upper = P_a Q_b − P_b Q_a, lower = Q_a P_b − Q_b P_a, from the (upper-right,
    // lower-left) blocks of legs A,B. (The σ^{μν}=(i/2)[γ^μ,γ^ν] normalization is applied in the scalar.)
    auto subB2 = [&](const B2 &x, const B2 &y) { return B2{x[0] - y[0], x[1] - y[1], x[2] - y[2], x[3] - y[3]}; };
    auto commBlocks = [&](const B2 &Pa, const B2 &Qa, const B2 &Pb, const B2 &Qb, B2 &Su, B2 &Sl) {
      Su = subB2(mul2(Pa, Qb), mul2(Pb, Qa));
      Sl = subB2(mul2(Qa, Pb), mul2(Qb, Pa));
    };
    // DFS over the free legs IN CHAIN ORDER, carrying the partial two-block product down the tree.
    // The running product is two Weyl blocks (m0,m1) with alternating roles (`antidiag`); the
    // antidiagonal-factor count is even, so the product is block-diagonal and each leaf's trace is
    // tr(m0)+tr(m1).
    //
    // Assignments sharing a free-leg prefix share the identical prefix product. The obvious
    // alternative — decode each of the 4^f assignments and re-fold the whole chain from scratch —
    // recomputes every such prefix 4^(remaining legs) times; it lived here behind an NT_DIRAC_FLAT=1
    // hatch as a bit-identical A/B control, measured at 3.64 B mul2 calls against the DFS's 522 M on
    // a production flow. Nothing ever referenced the hatch, so it and the flat fold are gone.
    //
    // The token→mul2 sequence per leaf is unchanged by the sharing, so every tensor entry is
    // bit-identical to the flat fold — only the redundant recomputation is gone (mul2 count drops
    // from ~4^f·L toward the 4/3·4^f tree sum). Leaf visit order is ascending flat index: the first
    // free leg in chain order is the most significant digit, matching the `idx = idx*4 + val` read
    // convention. Live memory is one (m0,m1) pair per recursion level — chain-depth bounded, no
    // cache, no eviction.
    auto walk = [&](auto &&self, std::size_t i, int flat, const B2 &m0, const B2 &m1, bool started,
                    bool antidiag) -> void {
      if (i == chain.size()) {
        F.entries[static_cast<std::size_t>(flat)] = m0[0] + m0[3] + m1[0] + m1[3];
        return;
      }
      const network::DFac &d = chain[i];
      if (d.kind == network::DFac::Gamma5) { // block-diagonal: negate one block, parity unchanged
        if (!started) {
          self(self, i + 1, flat, id2, neg2(id2), true, false);
          return;
        }
        if (antidiag)
          self(self, i + 1, flat, neg2(m0), m1, started, antidiag);
        else
          self(self, i + 1, flat, m0, neg2(m1), started, antidiag);
        return;
      }
      if (d.kind == network::DFac::Comm) { // block-diagonal diag(Su,Sl); free legs branch 4-way each
        const bool aFree = !commAslash[i], bFree = !commBslash[i];
        for (int muA = 0; muA < (aFree ? 4 : 1); ++muA)
          for (int muB = 0; muB < (bFree ? 4 : 1); ++muB) {
            const B2 &Pa = aFree ? gP[muA] : cAP[i];
            const B2 &Qa = aFree ? gQ[muA] : cAQ[i];
            const B2 &Pb = bFree ? gP[muB] : cBP[i];
            const B2 &Qb = bFree ? gQ[muB] : cBQ[i];
            B2 Su, Sl;
            commBlocks(Pa, Qa, Pb, Qb, Su, Sl);
            int nf = flat;
            if (aFree) nf = nf * 4 + muA;
            if (bFree) nf = nf * 4 + muB;
            if (!started)
              self(self, i + 1, nf, Su, Sl, true, false);
            else if (antidiag)
              self(self, i + 1, nf, mul2(m0, Sl), mul2(m1, Su), started, antidiag);
            else
              self(self, i + 1, nf, mul2(m0, Su), mul2(m1, Sl), started, antidiag);
          }
        return;
      }
      if (d.kind == network::DFac::Gamma) { // free leg: 4-way branch on the concrete component
        for (int mu = 0; mu < 4; ++mu) {
          const int nf = flat * 4 + mu;
          if (!started)
            self(self, i + 1, nf, gP[mu], gQ[mu], true, true);
          else if (antidiag)
            self(self, i + 1, nf, mul2(m0, gQ[mu]), mul2(m1, gP[mu]), started, !antidiag);
          else
            self(self, i + 1, nf, mul2(m0, gP[mu]), mul2(m1, gQ[mu]), started, !antidiag);
        }
        return;
      }
      // Slash: fixed antidiagonal factor
      if (!started) {
        self(self, i + 1, flat, sP[i], sQ[i], true, true);
        return;
      }
      if (antidiag)
        self(self, i + 1, flat, mul2(m0, sQ[i]), mul2(m1, sP[i]), started, !antidiag);
      else
        self(self, i + 1, flat, mul2(m0, sP[i]), mul2(m1, sQ[i]), started, !antidiag);
    };
    const B2 seed0{}, seed1{}; // pre-start state: empty blocks (the fold's `started` flag seeds them)
    walk(walk, 0, 0, seed0, seed1, false, true);
    return F;
  }
#endif // NUMTRACER_DEFINE_BODIES

  namespace ndetail
  {
    /// Externally referenced helper (test_projector_fusion): declared always, defined below.
    NUMTRACER_FUNC void fuse_projectors(std::vector<NElem> &e, Cx &coeff);
  } // namespace ndetail

#if NUMTRACER_DEFINE_BODIES
  namespace ndetail
  {

    /// @brief Build the dense `4^rank` tensor for one @ref NElem, entry by entry.
    ///
    /// One branch per variant; the formula each writes, with `k` = the element's momentum resolved
    /// through the component table and `at`/`atS` the inverse atoms `1/k²` and `1/|k⃗|²`:
    ///
    ///   Metric    δ_{ab}                                     rank 2
    ///   Vector    k_a                                        rank 1
    ///   ProjT     δ_{ab} − k_a k_b · at                       rank 2  (transverse)
    ///   ProjL     k_a k_b · at                                rank 2  (longitudinal)
    ///   ProjM     δ_{ab}(spatial) − k_a k_b · atS             rank 2  (magnetic; row/col 0 vanish)
    ///   ProjE     P_T − P_M, written out                      rank 2  (electric)
    ///   Epsilon   ε_{abcd}                                    rank 4
    ///
    /// This is the DENSE form. Callers go through @ref push_elem_factors, which splits an electric
    /// projector into two rank-1 factors where it can — that is what lets the elimination cut the
    /// network there — and falls back to this builder for everything else.
    inline Factor elem_factor(int nsym, const NElem &el, const std::vector<std::array<MPoly, 4>> &comp)
    {
      Factor F;
      if (el.kind == NElem::Metric) { // metric δ_{a,b}
        F.ids = {el.a, el.b};
        F.entries.assign(16, MPolyFactory::zero(nsym));
        for (int i = 0; i < 4; ++i)
          F.entries[i * 4 + i] = MPolyFactory::constant(nsym, Cx{1, 0});
      } else if (el.kind == NElem::Vector) { // vector on id a
        F.ids = {el.a};
        auto cc = mom_components(nsym, el.vlc, comp);
        F.entries = {cc[0], cc[1], cc[2], cc[3]};
      } else if (el.kind == NElem::ProjT) { // transverse projector P_T(k)_{a,b} = δ − k_a k_b INV(k)
        F.ids = {el.a, el.b};
        F.entries.assign(16, MPolyFactory::zero(nsym));
        auto k = mom_components(nsym, el.vlc, comp);
        const MPoly at = MPolyFactory::atom(nsym, el.atom);
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j) {
            MPoly e = (i == j) ? MPolyFactory::constant(nsym, Cx{1, 0}) : MPolyFactory::zero(nsym);
            e = e - (k[i] * k[j]) * at;
            F.entries[i * 4 + j] = std::move(e);
          }
      } else if (el.kind == NElem::ProjL) { // longitudinal projector P_L(k)_{a,b} = k_a k_b INV(k)
        F.ids = {el.a, el.b};
        F.entries.assign(16, MPolyFactory::zero(nsym));
        auto k = mom_components(nsym, el.vlc, comp);
        const MPoly at = MPolyFactory::atom(nsym, el.atom);
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
            F.entries[i * 4 + j] = (k[i] * k[j]) * at;
      } else if (el.kind == NElem::ProjM) { // magnetic projector P_M_{i,j}=δ_{ij}−k_i k_j INVS(k) (spatial)
        F.ids = {el.a, el.b};
        F.entries.assign(16, MPolyFactory::zero(nsym));
        auto k = mom_components(nsym, el.vlc, comp); // spatial part: zero the temporal (slot 0) component
        const MPoly atS = MPolyFactory::atom(nsym, el.atomS);
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j) {
            MPoly e = (i == j && i > 0) ? MPolyFactory::constant(nsym, Cx{1, 0}) : MPolyFactory::zero(nsym);
            if (i > 0 && j > 0) e = e - (k[i] * k[j]) * atS;
            F.entries[i * 4 + j] = std::move(e);
          }
      } else if (el.kind == NElem::ProjE) { // electric projector P_E = P_T − P_M
        F.ids = {el.a, el.b};
        F.entries.assign(16, MPolyFactory::zero(nsym));
        auto k = mom_components(nsym, el.vlc, comp);
        const MPoly at = MPolyFactory::atom(nsym, el.atom);
        const MPoly atS = MPolyFactory::atom(nsym, el.atomS);
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j) {
            MPoly full = (i == j) ? MPolyFactory::constant(nsym, Cx{1, 0}) : MPolyFactory::zero(nsym);
            full = full - (k[i] * k[j]) * at; // P_T entry
            MPoly mag = (i == j && i > 0) ? MPolyFactory::constant(nsym, Cx{1, 0}) : MPolyFactory::zero(nsym);
            if (i > 0 && j > 0) mag = mag - (k[i] * k[j]) * atS; // P_M entry
            F.entries[i * 4 + j] = full - mag;
          }
      } else { // Levi-Civita ε_{a,b,c,d}
        F.ids = {el.a, el.b, el.c, el.d};
        F.entries.assign(256, MPolyFactory::zero(nsym));
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
            for (int p = 0; p < 4; ++p)
              for (int q = 0; q < 4; ++q) {
                const double s = levi(i, j, p, q);
                if (s != 0.0) F.entries[((i * 4 + j) * 4 + p) * 4 + q] = MPolyFactory::constant(nsym, Cx{s, 0});
              }
      }
      return F;
    }

    /// @brief Push the dense factor(s) for one @ref NElem onto @p out.
    ///
    /// Everything except the finite-T ELECTRIC projector contributes exactly one factor, i.e. this is
    /// `out.push_back(elem_factor(...))`. `ProjE` is the exception, and the reason this wrapper
    /// exists: it has RANK 1 (the engine already knows — `projTrace4D` returns 1 for it), so instead
    /// of one dense 4x4 factor carrying BOTH Lorentz indices it can be emitted as an outer product of
    /// two 4-vectors, one per index:
    ///
    ///     P_E = v ⊗ v · INV(k) · INVS(k),   v_0 = |k⃗|²,   v_i = −k_0 k_i  (i = 1,2,3)
    ///
    /// Check the two corners against `P_E = P_T − P_M`: at (0,0) it gives `|k⃗|⁴/(k²|k⃗|²) = |k⃗|²/k²`,
    /// which is `1 − k_0²/k²`; at (i,j) it gives `k_0² k_i k_j/(k²|k⃗|²)`, which is
    /// `k_i k_j (1/|k⃗|² − 1/k²)`. Exact, not an approximation.
    ///
    /// WHY IT PAYS. `contract_factors` is greedy min-width variable elimination: its cost is
    /// exponential in the size of the id-union of each factor group. A dense `{a,b}` factor forces
    /// those two indices into a common group forever; two rank-1 factors `{a}` and `{b}` do not, so
    /// the elimination can cut the network there. That is the classic tensor-network separator win,
    /// and it is worth far more than halving one factor. Nothing in the elimination itself changes —
    /// it simply finds cheaper orderings once the artificial coupling is gone.
    ///
    /// `v_0` is `atomDen[el.atomS]`: the SPATIAL denominator |k⃗|², which the caller already holds, so
    /// no `k²` polynomial is needed anywhere (`k² − k_0² = |k⃗|²` is what makes that work).
    ///
    /// The split is unconditional. It used to sit behind an NT_NO_RANK1_PROJE=1 A/B control that
    /// nothing ever set — `test_rank1_proje` grades the two forms by calling both builders directly,
    /// precisely so it does not depend on the hatch. The dense builder is still reached, by every
    /// element the guard below rejects.
    inline void push_elem_factors(std::vector<Factor> &out, int nsym, const NElem &el,
                                  const std::vector<std::array<MPoly, 4>> &comp,
                                  const std::vector<MPoly> &atomDen)
    {
      // atomS must be a real id with a filled denominator; a malformed net falls back to dense
      // rather than silently building a wrong factor. A self-contracted P_E (a == b) is a trace,
      // not a separable pair, so it goes dense too.
      const bool canSplit = el.kind == NElem::ProjE && el.a != el.b && el.atom >= 0 &&
                            el.atomS >= 0 && static_cast<std::size_t>(el.atomS) < atomDen.size();
      if (!canSplit) {
        out.push_back(elem_factor(nsym, el, comp));
        return;
      }
      const auto k = mom_components(nsym, el.vlc, comp);
      const MPoly at = MPolyFactory::atom(nsym, el.atom);   // 1/k²
      const MPoly atS = MPolyFactory::atom(nsym, el.atomS); // 1/|k⃗|²
      std::array<MPoly, 4> v{atomDen[static_cast<std::size_t>(el.atomS)], MPolyFactory::zero(nsym),
                             MPolyFactory::zero(nsym), MPolyFactory::zero(nsym)};
      for (int i = 1; i < 4; ++i)
        v[static_cast<std::size_t>(i)] = MPolyFactory::zero(nsym) - k[0] * k[static_cast<std::size_t>(i)];
      // The two scalar atoms ride on ONE leg, so the product over the pair is v_a v_b · at · atS.
      Factor A, B;
      A.ids = {el.a};
      B.ids = {el.b};
      A.entries.reserve(4);
      B.entries.reserve(4);
      for (int i = 0; i < 4; ++i) {
        A.entries.push_back(v[static_cast<std::size_t>(i)] * at * atS);
        B.entries.push_back(v[static_cast<std::size_t>(i)]);
      }
      out.push_back(std::move(A));
      out.push_back(std::move(B));
    }

    /// First-seen-ordered union of all Lorentz ids carried by a factor group (duplicates dropped). The
    /// order matters: it fixes the union-index layout the elimination decodes against.
    inline std::vector<int> first_seen_union(const std::vector<Factor> &group)
    {
      std::vector<int> ids;
      for (const Factor &F : group)
        for (int id : F.ids)
          if (std::find(ids.begin(), ids.end(), id) == ids.end()) ids.push_back(id);
      return ids;
    }

    /// Merge the factors @p group (all sharing the Lorentz id @p elim) into one, then SUM over @p elim
    /// (0..3) — one step of variable elimination. The result carries the union of the group's ids minus
    /// @p elim. Each factor maps its own slots onto the union positions, so an id repeated within a
    /// factor (a self-trace) is handled by sharing the same union slot.
    inline Factor eliminate(int nsym, const std::vector<Factor> &group, int elim, const std::vector<MPoly> &atomDen,
                            const std::vector<std::vector<int>> &units)
    {
      // Algorithm (one variable-elimination step):
      //   1. `unionIds` = all Lorentz ids carried by the group, first-seen order; `elim` is one of them.
      //   2. For each output assignment of the surviving ids (= union minus `elim`), sum over the 4
      //      values of `elim`: at each value, multiply every group factor's entry at the matching slot.
      //   3. Both the output index `outFlat` and each factor's read index are decoded/encoded MSB-first
      //      (id[0] most significant) so a rank>=2 asymmetric result is not transposed.
      //   4. Reduce the result tensor in place (unit-constraint + monomial cancellation) so projector
      //      `k⊗k·INV` factors don't bloat the running tensor across a long chain.
      // union of ids (in first-seen order), and `elim`'s position within it.
      const std::vector<int> unionIds = first_seen_union(group);
      int elimPos = -1;
      for (int p = 0; p < (int)unionIds.size(); ++p)
        if (unionIds[p] == elim) {
          elimPos = p;
          break;
        }
      const int unionRank = unionIds.size();
      // per factor, the union-position of each of its slots.
      std::vector<std::vector<int>> slotPos(group.size());
      for (std::size_t fIdx = 0; fIdx < group.size(); ++fIdx)
        for (int id : group[fIdx].ids)
          for (int p = 0; p < unionRank; ++p)
            if (unionIds[p] == id) {
              slotPos[fIdx].push_back(p);
              break;
            }
      // output ids = union minus elim.
      Factor out;
      for (int id : unionIds)
        if (id != elim) out.ids.push_back(id);
      const int outRank = out.ids.size();
      int outTotal = 1;
      for (int k = 0; k < outRank; ++k)
        outTotal *= 4;
      out.entries.assign(outTotal, MPolyFactory::zero(nsym));
      // map output position -> union position (output is union with elimPos removed).
      std::vector<int> outToUnion;
      for (int p = 0; p < unionRank; ++p)
        if (p != elimPos) outToUnion.push_back(p);
      std::vector<int> idxVal(unionRank, 0); // current value at each union position
      for (int outFlat = 0; outFlat < outTotal; ++outFlat) {
        // decode MSB-first (out.ids[0] most significant) to match the read convention
        // `idx = idx*4 + idxVal[pos]` used everywhere — else a rank>=2 asymmetric output is transposed.
        int r = outFlat;
        for (int k = outRank - 1; k >= 0; --k) {
          idxVal[outToUnion[k]] = r % 4;
          r /= 4;
        }
        MPoly acc = MPolyFactory::zero(nsym);
        for (int elimVal = 0; elimVal < 4; ++elimVal) {
          idxVal[elimPos] = elimVal;
          // Seed from the FIRST factor rather than from `constant(1)`: the old form put a one-term
          // polynomial through the full `operator*` on the first iteration, paying an |e|-entry
          // scratch and a `from_scratch` sort to compute what is just `1·e`. With 2-4 factors per
          // group that identity multiply was a quarter to a half of all elimination multiplies, and
          // this is the innermost loop of the innermost operation. `scaled(e, 1)` applies the very
          // same coefficient product (`Cx{1,0} * c`) to the same monomials in the same order, so it is
          // bit-identical — it just skips the scratch and the sort.
          MPoly prod = MPolyFactory::zero(nsym);
          bool seeded = false;
          bool zero = false;
          for (std::size_t fIdx = 0; fIdx < group.size() && !zero; ++fIdx) {
            int idx = 0;
            for (int p : slotPos[fIdx])
              idx = idx * 4 + idxVal[p];
            const MPoly &e = group[fIdx].entries[idx];
            if (e.empty()) {
              zero = true;
              break;
            }
            if (!seeded) {
              prod = MPolyFactory::scaled(nsym, e, Cx{1, 0});
              seeded = true;
            } else
              prod = prod * e;
          }
          // an EMPTY group is the empty product = 1, exactly what the old constant(1) seed gave
          if (!zero && !seeded) prod = MPolyFactory::constant(nsym, Cx{1, 0});
          if (!zero) acc = std::move(acc) + std::move(prod); // rvalue +: first iteration MOVES prod
        }
        // Reduce the intermediate IN PLACE: cancel bare-loop atoms (sin²→1-cos² makes k²=l1² a monomial)
        // and collapse sin² powers, so projector `k⊗k·INV` factors don't bloat the running tensor across
        // a long pure-gauge chain (the ZA4 monster nets) — only cleaning at the end is too late.
        if (!units.empty()) {
          NT_STAT_TIMER(t_cf_reduce);
          // Both reductions pass through on the large majority of calls, and on a `const MPoly&` a
          // pass-through `return p` is a full DEEP COPY — two of them, per outFlat, per elimination
          // step. Moving picks the rvalue overloads, which move instead. `acc` is overwritten by the
          // very assignment whose RHS moves from it, so nothing reads the moved-from value.
          acc = divThroughMonomialAtoms(reduce_units(std::move(acc), units), atomDen);
        }
        out.entries[outFlat] = std::move(acc);
      }
      return out;
    }

    /// @brief Reject an OPEN (uncontracted) Lorentz index before it is silently summed away.
    ///
    /// `eliminate` sums EVERY live id over 0..3, whether or not a second slot carries it. An id
    /// appearing exactly once across the factor list is therefore not "left free": it is contracted
    /// against `(1,1,1,1)`, which is not a tensor operation at all. The result is a plausible-looking
    /// but meaningless scalar, with no diagnostic anywhere — the worst failure mode the engine can
    /// have. A finished network is a scalar (see documentation/getting_started/scope-and-conventions.md),
    /// so a once-occurring id is always a caller error: an unclosed hand-built net, a mis-split
    /// diagram, or a projector/basis element left off.
    ///
    /// Only count 1 is rejected. Count 2 is the ordinary contraction (including a self-trace whose two
    /// slots sit on the SAME factor, e.g. `g^{μμ}`); counts ≥ 3 are caught upstream by the front-end
    /// (`NumTrace::badlabel`) and are left alone here so this guard cannot change what a valid net does.
    inline void assert_no_open_ids(const std::vector<Factor> &facs)
    {
      // Flat parallel arrays, not a map: id lists are tiny (≤ ~14 distinct) and this runs per Lorentz
      // term of every diagram, so a linear scan beats any allocation.
      gch::small_vector<int, 16> ids, cnt;
      for (const Factor &F : facs)
        for (int id : F.ids) {
          std::size_t p = 0;
          for (; p < ids.size(); ++p)
            if (ids[p] == id) break;
          if (p == ids.size()) {
            ids.push_back(id);
            cnt.push_back(1);
          } else
            ++cnt[p];
        }
      for (std::size_t p = 0; p < ids.size(); ++p)
        if (cnt[p] == 1) {
          const std::string msg = "numtracer: Lorentz index id " + std::to_string(ids[p]) +
                                  " is OPEN (occurs once) — the network does not close to a scalar. "
                                  "Every Lorentz index must be contracted: tie the leg off with a "
                                  "metric/vector/projector, or pin it to a fixed component. (Contracting "
                                  "an open index would silently sum it over 0..3 and return a "
                                  "meaningless number.)";
          NT_THROW(std::runtime_error, msg.c_str());
        }
    }

    /// @brief Size of the union of the Lorentz ids carried by every factor incident to @p cand.
    ///
    /// This is the min-width score: eliminating `cand` merges exactly those factors into one
    /// intermediate, whose rank is this union's size. Linear scans (not a set) on purpose — a factor
    /// carries a handful of ids and the union is bounded by the network's treewidth, so the constant
    /// factor of a hash set costs more than it saves.
    inline std::size_t incident_union_size(const std::vector<Factor> &facs, int cand)
    {
      std::vector<int> unionIds;
      for (const Factor &F : facs) {
        const bool incident = std::find(F.ids.begin(), F.ids.end(), cand) != F.ids.end();
        if (!incident) continue;
        for (int id : F.ids)
          if (std::find(unionIds.begin(), unionIds.end(), id) == unionIds.end()) unionIds.push_back(id);
      }
      return unionIds.size();
    }

    /// Contract a set of dense factors over their shared Lorentz ids into one scalar MPoly, by GREEDY
    /// VARIABLE ELIMINATION: repeatedly pick the id whose incident factors have the smallest combined
    /// id-set (the MIN-WIDTH heuristic — the score is the size of the incident-id UNION, i.e. the rank
    /// of the intermediate the elimination would create, NOT the number of incident factors, which is
    /// what min-degree would count), merge those factors, and sum that id out. This keeps every
    /// intermediate tensor treewidth-bounded — a single global 4^G sum (G = all ids) is exponential and
    /// blows up for the quark triangle (G≈14); elimination stays at 4^(few) per step.
    ///
    /// Throws (via @ref assert_no_open_ids) if the factor list does not close — see there for why an
    /// open index cannot simply be left free.
    inline MPoly contract_factors(int nsym, std::vector<Factor> facs, const std::vector<MPoly> &atomDen = {},
                                  const std::vector<std::vector<int>> &units = {})
    {
      assert_no_open_ids(facs);
      // Algorithm (greedy min-width elimination), repeated until no shared ids remain:
      //   1. collect the distinct live Lorentz ids;
      //   2. score each id by the size of the combined id-set of its incident factors (the rank of
      //      the intermediate that eliminating it would create) and pick the smallest;
      //   3. partition factors into those touching that id (`group`) and the rest, merge+sum the group
      //      over that id via `eliminate`, push the result back, repeat;
      //   4. once nothing is shared, the leftovers are scalars — multiply their single entries.
      for (;;) {
        int bestId = -1;
        {
          NT_STAT_TIMER(t_cf_score);
          // distinct ids still live, and which factors carry each.
          std::vector<int> ids;
          for (const Factor &F : facs)
            for (int id : F.ids) {
              bool seen = false;
              for (int u : ids)
                if (u == id) {
                  seen = true;
                  break;
                }
              if (!seen) ids.push_back(id);
            }
          if (ids.empty()) break;
          // choose the id with the smallest incident-factor union (min width).
          std::size_t bestUnionSize = SIZE_MAX;
          for (int cand : ids) {
            const std::size_t width = incident_union_size(facs, cand);
            if (width < bestUnionSize) {
              bestUnionSize = width;
              bestId = cand;
            }
          }
        }
        NT_STAT_ADD(cf_steps, 1);
        // partition: factors incident to bestId vs the rest.
        std::vector<Factor> group, rest;
        for (Factor &F : facs) {
          bool has = false;
          for (int id : F.ids)
            if (id == bestId) {
              has = true;
              break;
            }
          (has ? group : rest).push_back(std::move(F));
        }
        rest.push_back(eliminate(nsym, group, bestId, atomDen, units));
        facs = std::move(rest);
      }
      // remaining factors are scalars (no ids): multiply their single entries. Seeded from the first
      // factor for the same reason as in `eliminate` above — `constant(1) * e` pays a scratch + sort
      // for an identity. A zero factor absorbs, so stop there rather than multiplying zero through.
      MPoly prod = MPolyFactory::zero(nsym);
      bool seeded = false;
      for (const Factor &F : facs) {
        if (F.entries.empty()) {
          prod = MPolyFactory::zero(nsym);
          seeded = true;
          break;
        }
        if (!seeded) {
          prod = MPolyFactory::scaled(nsym, F.entries[0], Cx{1, 0});
          seeded = true;
        } else
          prod = prod * F.entries[0];
      }
      if (!seeded) prod = MPolyFactory::constant(nsym, Cx{1, 0}); // empty product = 1, as before
      return prod;
    }

    /// @brief Reduce a Dirac-trace @ref Factor to a scalar @ref MPoly when there is NO surrounding
    ///        Lorentz net. Three cases: (a) no free legs → the scalar entry; (b) the trace is
    ///        structurally zero (odd antidiagonal count — `numeric_dirac` zeroed every entry) → 0; (c)
    ///        every free-leg id is PAIRED (a γ^μ…γ^μ self-contraction, e.g. a dressed numerator whose
    ///        IDENT/SLASH choice left two free legs on the same id) → sum those ids via
    ///        @ref contract_factors. An UNPAIRED free leg with no Lorentz partner is a genuinely
    ///        uncontracted Lorentz index (an ill-formed/mis-split diagram) and still throws.
    ///        This is what lets the dressed-collection combination loop produce odd / self-paired
    ///        chains (which the front-end odd-trace filter never sees) without a spurious abort.
    inline MPoly close_free_legs(int nsym, const Factor &T, const std::vector<MPoly> &atomDen = {},
                                 const std::vector<std::vector<int>> &units = {})
    {
      if (T.ids.empty()) return T.entries.empty() ? MPolyFactory::constant(nsym, Cx{1, 0}) : T.entries[0];
      bool allZero = true;
      for (const MPoly &v : T.entries)
        if (!v.empty()) {
          allZero = false;
          break;
        }
      if (allZero) return MPolyFactory::zero(nsym); // odd antidiagonal count: trace ≡ 0
      std::vector<int> s = T.ids;
      std::sort(s.begin(), s.end());
      bool paired = true;
      for (std::size_t i = 0; i < s.size();) {
        std::size_t j = i;
        while (j < s.size() && s[j] == s[i])
          ++j;
        if ((j - i) % 2 != 0) {
          paired = false;
          break;
        }
        i = j;
      }
      if (paired) return contract_factors(nsym, {T}, atomDen, units); // self-contract paired legs
      NT_THROW(std::runtime_error, "numeric_value: Dirac chain has an UNPAIRED free leg but no Lorentz net "
                                   "(uncontracted Lorentz index — mis-split diagram)");
    }

    /// @brief Split a Dirac chain at @ref network::DFac::LoopSep markers into its independent closed
    ///        spinor loops. A chain with no separator yields one segment (the whole chain) — so
    ///        single-loop diagrams are unaffected. Each segment is traced separately and the resulting
    ///        Lorentz tensors are multiplied / contracted by @ref contract_factors.
    inline std::vector<network::DiracNet> split_loops(const network::DiracNet &chain)
    {
      std::vector<network::DiracNet> segs;
      network::DiracNet cur;
      for (const network::DFac &d : chain) {
        if (d.kind == network::DFac::LoopSep) {
          if (!cur.empty()) segs.push_back(std::move(cur));
          cur.clear();
        } else
          cur.push_back(d);
      }
      if (!cur.empty()) segs.push_back(std::move(cur));
      return segs;
    }

    /// @brief Trace each spinor loop of @p dirac into a Lorentz @ref Factor (one per loop). Each loop is
    ///        an independent γ-trace; the gluon legs they share are contracted later via the Lorentz net
    ///        (or, with no Lorentz net, among themselves). Empty segments are skipped.
    inline std::vector<Factor> dirac_loop_factors(int nsym, const network::DiracNet &dirac,
                                                  const std::vector<std::array<MPoly, 4>> &comp)
    {
      std::vector<Factor> fs;
      for (const network::DiracNet &seg : split_loops(dirac))
        fs.push_back(numeric_dirac(nsym, seg, comp));
      return fs;
    }

    /// @brief Close a set of Dirac-loop Factors with NO surrounding Lorentz net: contract them together
    ///        (gluon legs shared between loops, or a loop's self-paired legs, are summed). One loop with
    ///        no free legs is its scalar; one loop with free legs defers to @ref close_free_legs (odd→0 /
    ///        paired self-contract / unpaired→throw).
    inline MPoly close_loops(int nsym, const std::vector<Factor> &fs, const std::vector<MPoly> &atomDen = {},
                             const std::vector<std::vector<int>> &units = {})
    {
      if (fs.empty()) return MPolyFactory::constant(nsym, Cx{1, 0});
      if (fs.size() == 1) return close_free_legs(nsym, fs[0], atomDen, units);
      return contract_factors(nsym, fs, atomDen, units); // multiple loops: contract their shared legs
    }

    // ──────────────────────────── projector-algebra fusion (T/L/E/M) ────────────────────────────

    inline bool isProjKind(NElem::Kind k)
    {
      return k == NElem::ProjT || k == NElem::ProjL || k == NElem::ProjE || k == NElem::ProjM;
    }
    /// `P_L ⟂ {P_T,P_E,P_M}` and `P_E ⟂ P_M` (on the same momentum) ⇒ their product is the zero tensor.
    inline bool projOrthogonal(NElem::Kind A, NElem::Kind B)
    {
      auto lOrtho = [](NElem::Kind x) { return x == NElem::ProjT || x == NElem::ProjE || x == NElem::ProjM; };
      if (A == NElem::ProjL && lOrtho(B)) return true;
      if (B == NElem::ProjL && lOrtho(A)) return true;
      if ((A == NElem::ProjE && B == NElem::ProjM) || (A == NElem::ProjM && B == NElem::ProjE)) return true;
      return false;
    }
    /// `tr P` in 4D: `tr P_T=3, tr P_L=1, tr P_E=1, tr P_M=2` (cf. test_projector_identity).
    inline int projTrace4D(NElem::Kind k)
    {
      switch (k) {
      case NElem::ProjT: return 3;
      case NElem::ProjL: return 1;
      case NElem::ProjE: return 1;
      case NElem::ProjM: return 2;
      default: return 0;
      }
    }
    inline NElem makeProj(NElem::Kind k, int a, int b, const std::vector<std::pair<double, int>> &vlc, int atom,
                          int atomS)
    {
      switch (k) {
      case NElem::ProjT: return nprojT(a, b, vlc, atom);
      case NElem::ProjL: return nprojL(a, b, vlc, atom);
      case NElem::ProjE: return nprojE(a, b, vlc, atom, atomS);
      case NElem::ProjM: return nprojM(a, b, vlc, atomS);
      default: return nmet(a, b); // unreachable
      }
    }

    /// @brief Projector-algebra fusion on one term's element list, applied BEFORE the projectors are
    ///        expanded into `δ − kk/k²`. Two projectors on the **same momentum** contracting over a
    ///        shared dummy Lorentz index collapse by the projector algebra:
    ///          - idempotency `P·P → P` (one shared index) / `tr P` (both shared: T→3,L→1,E→1,M→2);
    ///          - orthogonality `L⟂{T,E,M}`, `E⟂M` ⇒ the contraction is the zero tensor ⇒ the term dies.
    ///        `T·E`/`T·M` are left untouched (they need `Pᵀ=Pᴱ+Pᴹ`; omitting an identity is value-safe).
    ///        Value-preserving and a strict **no-op unless such a pair exists**, so non-redundant nets
    ///        (≤ one projector per line) are unchanged. This is what keeps a transverse-vertex
    ///        (`…ClassTrans`) derivation as compact as the bare one: the backend cannot otherwise
    ///        re-cancel the redundant `kk/k²` of shifted lines (their multi-term `k²` denominator is not a
    ///        single monomial — see @ref divThroughMonomialAtoms). Mutates `e` (rebuilt) and `coeff`.
    NUMTRACER_FUNC void fuse_projectors(std::vector<NElem> &e, Cx &coeff)
    {
      // cheap early-exit: nothing to fuse without ≥2 projector elements.
      int nproj = 0, maxIdx = -1;
      for (const NElem &el : e) {
        if (isProjKind(el.kind)) ++nproj;
        maxIdx = std::max(maxIdx, std::max(std::max(el.a, el.b), std::max(el.c, el.d)));
      }
      if (nproj < 2) return;

      bool changed = true;
      while (changed) {
        changed = false;
        if (coeff.re == 0 && coeff.im == 0) return;

        // Lorentz-index multiplicity over the whole element list (a shared index is a dummy iff == 2).
        std::vector<int> mult(maxIdx + 1, 0);
        auto bump = [&](int idx) { if (idx >= 0 && idx <= maxIdx) ++mult[idx]; };
        for (const NElem &el : e) {
          if (el.kind == NElem::Epsilon) {
            bump(el.a); bump(el.b); bump(el.c); bump(el.d);
          } else if (el.kind == NElem::Vector) {
            bump(el.a);
          } else { // Metric / ProjT / ProjL / ProjE / ProjM
            bump(el.a); bump(el.b);
          }
        }

        for (std::size_t i = 0; i < e.size() && !changed; ++i) {
          if (!isProjKind(e[i].kind) || e[i].a == e[i].b) continue;
          for (std::size_t j = i + 1; j < e.size(); ++j) {
            if (!isProjKind(e[j].kind) || e[j].a == e[j].b) continue;
            if (e[i].vlc != e[j].vlc) continue; // same momentum only
            const NElem::Kind kA = e[i].kind, kB = e[j].kind;
            // same projector (full field identity) ⇒ idempotent; orthogonal pair ⇒ vanishes.
            const bool same = (kA == kB) && e[i].atom == e[j].atom && e[i].atomS == e[j].atomS;
            const bool ortho = projOrthogonal(kA, kB);
            if (!same && !ortho) continue; // T·E / T·M: leave untouched

            // shared indices among {a,b} that are CLEAN dummies (multiplicity exactly 2).
            int sharedCount = 0, sharedIdx = -1;
            for (int x : {e[i].a, e[i].b})
              if ((x == e[j].a || x == e[j].b) && mult[x] == 2) {
                ++sharedCount;
                sharedIdx = x;
              }
            if (sharedCount == 0) continue;

            std::vector<NElem> out;
            out.reserve(e.size());
            for (std::size_t t = 0; t < e.size(); ++t)
              if (t != i && t != j) out.push_back(e[t]);
            if (sharedCount == 2) { // full contraction → scalar trace (orthogonal ⇒ 0)
              const int tr = same ? projTrace4D(kA) : 0;
              coeff = Cx{coeff.re * tr, coeff.im * tr};
            } else if (same) { // single shared index → fused projector on the outer indices
              const int outA = (e[i].a == sharedIdx) ? e[i].b : e[i].a;
              const int outB = (e[j].a == sharedIdx) ? e[j].b : e[j].a;
              out.push_back(makeProj(kA, outA, outB, e[i].vlc, e[i].atom, e[i].atomS));
            } else { // orthogonal single contraction → zero tensor → whole term vanishes
              coeff = Cx{0, 0};
            }
            e.swap(out);
            changed = true;
            break;
          }
        }
      }
    }

  } // namespace ndetail
#endif // NUMTRACER_DEFINE_BODIES

  // Public entry points (numeric_value / numeric_value_netval); the dressed variants, to_genprog and
  // collect_atom_denoms are declared further below, after the dressing types they reference.
  NUMTRACER_FUNC MPoly numeric_value(int nsym, const network::DiracNet &dirac, const NNet &lorentz,
                                     const std::vector<std::array<MPoly, 4>> &comp,
                                     const std::vector<MPoly> &atomDen);
  NUMTRACER_FUNC MPoly numeric_value_netval(int nsym, const network::DiracNet &dirac,
                                            const network::NetVal &lor,
                                            const std::vector<std::array<MPoly, 4>> &comp,
                                            const std::vector<MPoly> &atomDen,
                                            const std::vector<std::vector<int>> &units = {});

#if NUMTRACER_DEFINE_BODIES
  /// @brief Contract a diagram (Dirac chain ⊗ Lorentz network) to its scalar trace polynomial.
  /// @param nsym     number of user symbols (MPoly variable count)
  /// @param dirac    the closed Dirac chain (may be empty for a pure-gauge diagram)
  /// @param lorentz  the pure-Lorentz network (metrics / vectors / projectors / Levi-Civita)
  /// @param comp     component table `comp[vid]` = 4 MPoly components of fundamental momentum `vid`
  /// @param atomDen  projector denominators `atomDen[aid] = k²` (for monomial cancellation)
  NUMTRACER_FUNC MPoly numeric_value(int nsym, const network::DiracNet &dirac, const NNet &lorentz,
                                     const std::vector<std::array<MPoly, 4>> &comp, const std::vector<MPoly> &atomDen)
  {
    // a component may hold SEVERAL independent spinor loops (e.g. a quark loop + the projection-closed
    // external line, tied only by gluon propagators) — trace each into its own Lorentz tensor; the
    // gluon legs they share contract via the Lorentz net below. tr(1)=4 rides each loop's trace.
    std::vector<ndetail::Factor> loops = ndetail::dirac_loop_factors(nsym, dirac, comp);
    MPoly result = MPolyFactory::zero(nsym);
    for (const NTerm &nt : lorentz) {
      // fold same-momentum projector chains (P·P→P, orthogonal→0) before expansion; no-op unless present.
      std::vector<NElem> elems = nt.e;
      Cx co = nt.coeff;
      ndetail::fuse_projectors(elems, co);
      if (co.re == 0 && co.im == 0) continue;
      std::vector<ndetail::Factor> facs = loops;
      facs.reserve(loops.size() + elems.size());
      for (const NElem &el : elems)
        ndetail::push_elem_factors(facs, nsym, el, comp, atomDen);
      // contract_factors consumes its `facs` argument by value — move so the per-term
      // factor list is built once, not copied again into the call.
      MPoly term = ndetail::contract_factors(nsym, std::move(facs));
      // `scaled`, not `* constant(co)`: same coefficient product (Cx multiply is commutative
      // bit-for-bit), same monomials in the same order, but without the |term|-entry scratch and
      // the `from_scratch` sort that multiplying by a one-term polynomial otherwise pays.
      term = MPolyFactory::scaled(nsym, term, co);
      result = std::move(result) + std::move(term);
    }
    if (lorentz.empty())
      result = ndetail::close_loops(nsym, loops, atomDen);
    return divThroughMonomialAtoms(std::move(result), atomDen);
  }

  /// @brief Map one inv-backend @ref network::Elem to a numeric @ref NElem. The projector's loop momentum
  ///        is stored as a single vector id (`vid`) and its `1/k²` env id (`inv`) becomes the atom id.
  inline NElem elem_to_nelem(const network::Elem &e)
  {
    switch (e.kind) {
    case network::Elem::Metric:
      return nmet(e.a, e.b);
    case network::Elem::Vector:
      return nvec(e.a, e.vlc);
    case network::Elem::ProjT:
      return nprojT(e.a, e.b, {{1.0, e.vid}}, e.inv);
    case network::Elem::ProjL:
      return nprojL(e.a, e.b, {{1.0, e.vid}}, e.inv);
    case network::Elem::ProjE:
      return nprojE(e.a, e.b, {{1.0, e.vid}}, e.inv, e.invS);
    case network::Elem::ProjM:
      return nprojM(e.a, e.b, {{1.0, e.vid}}, e.invS);
    case network::Elem::Epsilon:
      return neps(e.a, e.b, e.c, e.d);
    }
    return nmet(e.a, e.b); // unreachable; silences -Wreturn-type
  }

  /// @brief Append a collected-slot combination's open-leg net factors (@ref DSlotOpt::netFacs — a
  ///        vector `p^μ`, a metric `g^{μν}`, …) to EVERY term of a Lorentz net, so the surrounding net
  ///        closes their legs. This is the exact Lorentz structure the distributed diagram emits, which
  ///        is why it can reuse all the validated contraction machinery unchanged.
  ///
  /// An empty @p facs returns @p lor untouched — that is the case for every option whose structure is
  /// purely Dirac-side, and it keeps the byte-identical fast path of the pre-Stage-4 δ/slash collection.
  /// Both overloads exist because the two dressed entry points read different net representations
  /// (@ref network::NetVal for the codegen path, @ref NNet for the numeric one).
  inline network::NetVal with_slot_facs(const network::NetVal &lor, const std::vector<network::Elem> &facs)
  {
    if (facs.empty()) return lor;
    network::NetVal out = lor;
    for (network::PTerm &pt : out)
      pt.e.insert(pt.e.end(), facs.begin(), facs.end());
    return out;
  }
  inline NNet with_slot_facs(const NNet &lor, const std::vector<network::Elem> &facs)
  {
    if (facs.empty()) return lor;
    NNet out = lor;
    for (NTerm &t : out)
      for (const network::Elem &e : facs)
        t.e.push_back(elem_to_nelem(e));
    return out;
  }

  /// Multi-term denominator cancellation (@ref divThroughPolyAtoms) — on by default; `NT_GEN_NO_POLYDIV=1`
  /// disables it, which is how it was A/B'd. See that function for the measurement.
  inline bool polydiv_enabled()
  {
    static const bool on = !env_flag("NT_GEN_NO_POLYDIV");
    return on;
  }

  /// @brief Contract a diagram given the Lorentz part as an inv-backend @ref network::NetVal (so the
  ///        generator reuses the existing net-string emission: `proj`/`met`/`vec`/`epsilon` builders).
  ///        Equivalent to @ref numeric_value but reading `network::Elem` instead of @ref NElem.
  NUMTRACER_FUNC MPoly numeric_value_netval(int nsym, const network::DiracNet &dirac, const network::NetVal &lor,
                                            const std::vector<std::array<MPoly, 4>> &comp,
                                            const std::vector<MPoly> &atomDen,
                                            const std::vector<std::vector<int>> &units)
  {
    NT_STAT_ADD(traces, 1);
    // trace each independent spinor loop of the component (see numeric_value); the shared gluon legs
    // contract via the Lorentz net below.
    std::vector<ndetail::Factor> loops;
    {
      NT_STAT_TIMER(t_dirac);
      loops = ndetail::dirac_loop_factors(nsym, dirac, comp);
    }
    // pre-reduce the atom denominators (idempotent if the caller already did) so monomial-cancellation
    // detection works during the per-step intermediate reduction inside contract_factors.
    //
    // The caller — the generated driver — already does exactly this once at setup, so on every trace
    // the reduction is a NO-OP. It nevertheless used to cost a deep copy of the whole table plus one
    // more deep copy per entry (reduce_units' pass-through returned `p` by value), on 10^5 traces ×
    // combinations. So: test the pass-through predicate first and only materialise a reduced copy when
    // an entry genuinely needs rewriting; otherwise alias the caller's table. When every entry passes
    // through, `reduce_units` would have returned it unchanged, so aliasing is value-identical.
    std::vector<MPoly> adenOwned;
    bool needReduce = false;
    for (const MPoly &a : atomDen)
      if (!reduceUnitsIsNoop(a, units)) {
        needReduce = true;
        break;
      }
    if (needReduce) {
      adenOwned = atomDen;
      for (MPoly &a : adenOwned)
        a = reduce_units(std::move(a), units);
    }
    const std::vector<MPoly> &aden = needReduce ? adenOwned : atomDen;
    MPoly result = MPolyFactory::zero(nsym);
    for (const network::PTerm &pt : lor) {
      // build the numeric element list, then fold same-momentum projector chains before expansion.
      Cx co = pt.coeff;
      std::vector<ndetail::Factor> facs;
      {
        NT_STAT_TIMER(t_elem);
        std::vector<NElem> elems;
        elems.reserve(pt.e.size());
        for (const network::Elem &el : pt.e)
          elems.push_back(elem_to_nelem(el));
        ndetail::fuse_projectors(elems, co);
        if (co.re == 0 && co.im == 0) continue;
        facs = loops;
        facs.reserve(loops.size() + elems.size());
        for (const NElem &el : elems)
          ndetail::push_elem_factors(facs, nsym, el, comp, aden);
      }
      // move the per-term factor list into contract_factors (consumed by value) — avoids
      // a redundant deep copy of every Factor's MPoly entries.
      MPoly term;
      {
        NT_STAT_TIMER(t_contract);
        term = ndetail::contract_factors(nsym, std::move(facs), aden, units);
      }
      // `scaled` instead of `* constant(co)`: bit-identical (see MPolyFactory::scaled) without the
      // scratch + sort. Once per Lorentz PTerm per trace, on the FULLY CONTRACTED term.
      term = MPolyFactory::scaled(nsym, term, co);
      result = std::move(result) + std::move(term);
    }
    if (lor.empty())
      result = ndetail::close_loops(nsym, loops, aden, units);
    // sin^2 -> 1 - cos^2 BEFORE cancellation: collapses bare-loop k²=l1²(cos²+sin²) to the monomial l1²
    // (so its atom cancels) and shrinks the polynomial to the FORM angular basis.
    {
      NT_STAT_TIMER(t_reduce);
      result = reduce_units(std::move(result), units); // pass-through MOVES the whole trace polynomial
    }
    {
      NT_STAT_TIMER(t_divmono);
      result = divThroughMonomialAtoms(std::move(result), aden);
    }
    // Then cancel the MULTI-TERM (shifted-line) denominators by exact polynomial division — the case
    // divThroughMonomialAtoms structurally cannot reach. Off via NT_GEN_NO_POLYDIV=1 for A/B.
    // NOTE `aden`, not `atomDen`: the trial division must see the SAME unit-reduced denominators the
    // intermediate reductions used, or a numerator reduced mod ΣU²=1 will not divide by an unreduced D.
    if (!polydiv_enabled()) return result;
    NT_STAT_TIMER(t_divpoly);
    return divThroughPolyAtoms(result, aden);
  }
#endif // NUMTRACER_DEFINE_BODIES

  // ───────────────────────────── dressed structure sums (symbolic dressing collection) ────────────
  //
  // A dressed propagator numerator like `Mq·δ + Z(p)·(γ·p)` is a SUM of Dirac structures whose
  // per-structure coefficients are *runtime* dressings. The front-end keeps such a sum eager (it no
  // longer distributes the diagram into one copy per structure), and the generator hands it here as a
  // Dirac chain with SLOT tokens: each slot is the structure sum of one dressed numerator. We collect
  // it WITHOUT the `2^D`-diagram blowup: enumerate the structure choices, contract each concrete chain
  // with the EXISTING @ref numeric_value_netval (so all the validated Dirac/Lorentz machinery is reused
  // verbatim), and accumulate the results into ONE @ref DPoly keyed by the dressing monomial. The whole
  // diagram then lowers to ONE trace function whose dressing factors the shared CSE/Horner collects.

  /// @brief One structure option of a collected Dirac slot (Stage 4, general form): a numeric
  ///        coefficient × a product of dressing atoms (`dress`) × a Dirac structure, where the
  ///        structure is
  ///          - a **Dirac-token chain** @ref toks spliced in place of the slot (spinor `din→dout`);
  ///            its FREE-Lorentz tokens (a `dgamma(μ)` / an open `dcomm` leg) are open legs on the
  ///            Dirac side, and
  ///          - a set of **Lorentz-net factors** @ref netFacs (a vector `p^μ`, a metric `g^{μν}`, …)
  ///            carrying any remaining open legs, appended to the surrounding net.
  ///
  /// The union of the two sets' free Lorentz ids is the slot's shared open-leg set `{μ₁,…,μ_k}`, closed
  /// by the surrounding net — the SAME set for every option (so the net contracts a fixed leg set
  /// regardless of structure choice). This one form covers every case, any leg count `k ≥ 0`:
  ///  - propagator numerator `δ` ⇒ `toks={}, netFacs={}`; slash `γ·p̸` ⇒ `toks={dslash(vlc)}` (k=0);
  ///  - single-gluon vertex T1 `γ^μ` ⇒ `toks={dgamma(μ)}`; T4 `p̸₁γ^μ` ⇒ `toks={dslash(p1),dgamma(μ)}`;
  ///    T7 `σ^{μν}p̸_ν` ⇒ `toks={dcomm_fs(μ,vlc)}` (k=1);
  ///  - open leg on a vector `p^μ·δ` ⇒ `netFacs={Vector(μ,p)}`; two-gluon `g^{μν}·δ` ⇒
  ///    `netFacs={Metric(μ,ν)}`; two open γ's ⇒ `toks={dgamma(μ),dgamma(ν)}` (k=2), etc.
  struct DSlotOpt {
    Cx coeff{1, 0};
    std::vector<int> dress; ///< dressing-atom ids (need not be pre-sorted; merged sorted on use)
    std::vector<network::DFac> toks; ///< Dirac-token chain spliced for this option (free-Lorentz tokens = open legs)
    std::vector<network::Elem> netFacs; ///< extra Lorentz-net factors for this option (open legs on vectors/metrics/…)
  };
  /// @brief A dressed numerator = the sum of its structure options.
  using DSlot = std::vector<DSlotOpt>;
  /// @brief One token of a dressed Dirac chain: either a FIXED @ref network::DFac (γ / slash / γ5 /
  ///        commutator) or a SLOT reference (an index into the diagram's slot list).
  struct DChainTok {
    bool isSlot = false;
    network::DFac fac{network::DFac::Gamma, -1, {}, -1, {}};
    int slot = -1;
  };
  inline DChainTok dtfix(network::DFac f) { return {false, std::move(f), -1}; }
  inline DChainTok dtslot(int s) { return {true, {network::DFac::Gamma, -1, {}, -1, {}}, s}; }

  // Public entry points that reference the dressing types above.
  NUMTRACER_FUNC DPoly numeric_value_dressed_netval(int nsym, const std::vector<DChainTok> &chain,
                                                    const std::vector<DSlot> &slots, const network::NetVal &lor,
                                                    const std::vector<std::array<MPoly, 4>> &comp,
                                                    const std::vector<MPoly> &atomDen,
                                                    const std::vector<std::vector<int>> &units = {});
  /// @brief STRUCTURAL variant of @ref numeric_value_dressed_netval: contract the collected slots into a
  ///        PLAIN @ref MPoly, DISCARDING the dressing dimension (the dressing monomial keys are summed
  ///        away). This is lever (b): the generator strips each option's dressing (`coeff`→1, `dress`→{})
  ///        into a per-sub-term scalar/monomial and feeds the dressing-free structure here, so the trace
  ///        table dedups on structure alone (dressing variants of one concrete trace collapse) and each
  ///        entry is a plain `MPoly` with no dressing dimension. The 4^(#collapsed loop) tr(1)=4 factor
  ///        is STRUCTURAL (it depends only on the concrete chain) and is kept in the returned trace.
  ///        Only meaningful when every option's `dress` is empty (the generator guarantees this); if fed
  ///        genuine dressings it would silently sum them, which is why this is not the collection path.
  NUMTRACER_FUNC MPoly numeric_value_dressed_netval_mp(int nsym, const std::vector<DChainTok> &chain,
                                                       const std::vector<DSlot> &slots, const network::NetVal &lor,
                                                       const std::vector<std::array<MPoly, 4>> &comp,
                                                       const std::vector<MPoly> &atomDen,
                                                       const std::vector<std::vector<int>> &units = {});
  NUMTRACER_FUNC DPoly numeric_value_dressed(int nsym, const std::vector<DChainTok> &chain,
                                             const std::vector<DSlot> &slots, const NNet &lorentz,
                                             const std::vector<std::array<MPoly, 4>> &comp,
                                             const std::vector<MPoly> &atomDen);
  NUMTRACER_FUNC std::vector<MPoly> collect_atom_denoms(int nsym, const std::vector<network::NetVal> &lors,
                                                        const std::vector<std::array<MPoly, 4>> &comp);

#if NUMTRACER_DEFINE_BODIES
  namespace ndetail
  {
    /// Enumerate the Cartesian product of slot-structure choices, contract each concrete Dirac chain
    /// via @p contract (an `MPoly`-returning closure that runs the full Dirac+Lorentz contraction), and
    /// collect the results into one @ref DPoly keyed by the dressing monomial. Each slot is referenced
    /// at most once in @p chain (a numerator occupies one chain position). @p contract receives the
    /// concrete Dirac chain AND this combination's extra Lorentz-net factors (the chosen options'
    /// @ref DSlotOpt::netFacs — open legs on vectors/metrics/…); empty for every option whose structure
    /// is purely Dirac-side ⇒ the caller takes a byte-identical fast path.
    ///
    /// The enumeration is factored out of @ref dress_collect so the two accumulation policies share it:
    /// the DPoly collection (`emit(dmono, mp)` → `out.add(dmono, mp)`) and the STRUCTURAL MPoly reduction
    /// (@ref dress_collect_mp: `emit(_, mp)` → `out += mp`, dropping the dressing key). @p emit is called
    /// once per surviving combination with the sorted dressing monomial and the (coeff·trace) `MPoly`.
    template <class ContractFn, class Emit>
    inline void dress_enumerate(int nsym, const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                                ContractFn &&contract, Emit &&emit)
    {
      const int nSlots = static_cast<int>(slots.size());
      // The dressed chain is a quark line with ≥1 closed spinor loop ((LoopSep markers)+1). A combination
      // that collapses a WHOLE loop to the identity (every slot on it takes the δ option, no fixed γ) leaves
      // an empty segment: split_loops drops it and close_loops returns 1, silently losing that loop's
      // tr(1)=4 — a 4× undercount (e.g. the all-mass channel of a Yukawa self-energy quark loop). Count the
      // loops up front so each combination can restore 4^(#collapsed loops). Flows whose every loop keeps a
      // fixed γ (Zq/ZAqbq{1,4,7}, …) never collapse ⇒ nEmpty==0 ⇒ output byte-identical.
      int nloops = 1;
      for (const DChainTok &tok : chain)
        if (!tok.isSlot && tok.fac.kind == network::DFac::LoopSep) ++nloops;
      std::vector<int> choice(nSlots, 0);
      bool done = false;
      while (!done) {
        // build the concrete chain + accumulate this combination's numeric coeff and dressing monomial
        Cx combCoeff{1, 0};
        DMono dressMono;
        network::DiracNet concrete;
        concrete.reserve(chain.size());
        std::vector<network::Elem> extraNet; // this combination's options' net factors (usually empty)
        for (const DChainTok &tok : chain) {
          if (!tok.isSlot) {
            concrete.push_back(tok.fac);
            continue;
          }
          const DSlotOpt &opt = slots[tok.slot][choice[tok.slot]];
          combCoeff = combCoeff * opt.coeff;
          dressMono.insert(dressMono.end(), opt.dress.begin(), opt.dress.end());
          // General collected Dirac slot (Stage 4): splice this option's Dirac-token chain in place of
          // the slot (its free-Lorentz tokens are open legs), and collect its Lorentz-net factors — the
          // surrounding net closes every open leg for all structure choices alike. An option with an
          // empty `toks` is the spinor identity δ (tr(1)=4 restored by nCollapsed below when a whole
          // loop collapses); empty `netFacs` ⇒ no net change (the byte-identical fast path).
          concrete.insert(concrete.end(), opt.toks.begin(), opt.toks.end());
          extraNet.insert(extraNet.end(), opt.netFacs.begin(), opt.netFacs.end());
        }
        if (!(combCoeff.re == 0 && combCoeff.im == 0)) {
          // restore tr(1)=4 for every spinor loop that collapsed to the identity in this combination
          const int nCollapsed = nloops - static_cast<int>(split_loops(concrete).size());
          for (int c = 0; c < nCollapsed; ++c) combCoeff = combCoeff * Cx{4, 0};
          MPoly mp = contract(concrete, extraNet);
          if (!mp.empty()) {
            mp = mp * MPolyFactory::constant(nsym, combCoeff);
            emit(dmono_sorted(std::move(dressMono)), std::move(mp));
          }
        }
        // odometer over the slot choices: advance the lowest slot, carrying into the next on wrap.
        if (nSlots == 0) break;
        int k = 0;
        for (; k < nSlots; ++k) {
          if (++choice[k] < static_cast<int>(slots[k].size())) break;
          choice[k] = 0;
        }
        if (k == nSlots) done = true;
      }
    }

    /// @brief The DPoly collection: accumulate each combination into its dressing-monomial channel.
    template <class ContractFn>
    inline DPoly dress_collect(int nsym, const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                               ContractFn &&contract)
    {
      DPoly out = DPolyFactory::zero(nsym);
      dress_enumerate(nsym, chain, slots, std::forward<ContractFn>(contract),
                      [&](const DMono &d, MPoly &&mp) { out.add(d, mp); });
      return out;
    }

    /// @brief The STRUCTURAL MPoly reduction (lever (b)): sum every combination into ONE plain MPoly,
    ///        discarding the dressing monomial. Correct only when the slots carry no dressing (the
    ///        generator strips it out first); see @ref numeric_value_dressed_netval_mp. The debug
    ///        assert below makes that precondition checkable instead of comment-only: fed genuine
    ///        dressings, this would silently sum structures that belong in different channels.
    ///
    /// The left-to-right accumulation is deliberate and is NOT a quadratic hazard, despite looking
    /// like one. Under lever (b) the generator expands each structure×dressing combination into its
    /// own SINGLE-option sub-term at codegen time (`Codegen.m`, `sdsl[k].push_back(DSlot{optp[oi]})`),
    /// so every production call has exactly one combination and this loop body runs once. Folding the
    /// sum as a balanced tree instead would reassociate the like-term coefficient sums (≤ 1 ulp) and
    /// so could shift the emitted kernel's literals — a real cost for a case no caller reaches. Every
    /// caller today (the generator, and test_dpoly.cpp case K, which mirrors it) passes single-option
    /// dressing-free slots; the general loop is kept only so the entry point stays total.
    template <class ContractFn>
    inline MPoly dress_collect_mp(int nsym, const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                                  ContractFn &&contract)
    {
#ifndef NDEBUG
      for (const DSlot &s : slots)
        for (const DSlotOpt &o : s)
          assert(o.dress.empty() && "dress_collect_mp requires dressing-free slots (lever (b))");
#endif
      MPoly out = MPolyFactory::zero(nsym);
      dress_enumerate(nsym, chain, slots, std::forward<ContractFn>(contract),
                      [&](const DMono &, MPoly &&mp) { out = std::move(out) + std::move(mp); });
      return out;
    }
  } // namespace ndetail

  /// @brief Dressed analogue of @ref numeric_value_netval: contract a diagram whose Dirac chain carries
  ///        dressed-numerator SLOTS, returning the collected @ref DPoly (one term per dressing monomial).
  ///        Reuses @ref numeric_value_netval per structure combination, so the Dirac/Lorentz contraction
  ///        is the exact same validated code path; only the collection over dressing atoms is new.
  NUMTRACER_FUNC DPoly numeric_value_dressed_netval(int nsym, const std::vector<DChainTok> &chain,
                                                    const std::vector<DSlot> &slots, const network::NetVal &lor,
                                                    const std::vector<std::array<MPoly, 4>> &comp,
                                                    const std::vector<MPoly> &atomDen,
                                                    const std::vector<std::vector<int>> &units)
  {
    return ndetail::dress_collect(
        nsym, chain, slots, [&](const network::DiracNet &d, const std::vector<network::Elem> &slotFacs) {
          return numeric_value_netval(nsym, d, with_slot_facs(lor, slotFacs), comp, atomDen, units);
        });
  }

  /// @brief STRUCTURAL MPoly reduction of a collected diagram (lever (b)). Same contraction machinery as
  ///        @ref numeric_value_dressed_netval, but the dressing dimension is summed away, yielding a plain
  ///        @ref MPoly. The generator feeds this dressing-free slots (each option's `coeff`=1, `dress`={})
  ///        so the returned MPoly is the concrete structural trace (including its tr(1)=4 collapse
  ///        factors); the dressing rides the sub-term scalar/monomial the generator carries alongside.
  NUMTRACER_FUNC MPoly numeric_value_dressed_netval_mp(int nsym, const std::vector<DChainTok> &chain,
                                                       const std::vector<DSlot> &slots, const network::NetVal &lor,
                                                       const std::vector<std::array<MPoly, 4>> &comp,
                                                       const std::vector<MPoly> &atomDen,
                                                       const std::vector<std::vector<int>> &units)
  {
    return ndetail::dress_collect_mp(
        nsym, chain, slots, [&](const network::DiracNet &d, const std::vector<network::Elem> &slotFacs) {
          return numeric_value_netval(nsym, d, with_slot_facs(lor, slotFacs), comp, atomDen, units);
        });
  }

  /// @brief Dressed analogue of @ref numeric_value (reading the numeric @ref NNet Lorentz network).
  NUMTRACER_FUNC DPoly numeric_value_dressed(int nsym, const std::vector<DChainTok> &chain,
                                             const std::vector<DSlot> &slots, const NNet &lorentz,
                                             const std::vector<std::array<MPoly, 4>> &comp,
                                             const std::vector<MPoly> &atomDen)
  {
    return ndetail::dress_collect(
        nsym, chain, slots, [&](const network::DiracNet &d, const std::vector<network::Elem> &slotFacs) {
          return numeric_value(nsym, d, with_slot_facs(lorentz, slotFacs), comp, atomDen);
        });
  }

  /// @brief Build the projector inverse-atom denominators by scanning the Lorentz nets for every
  ///        projector kind. Each carries its loop momentum's fundamental id (`vid`); a transverse /
  ///        longitudinal / electric projector fills its full atom `atomDen[inv] = k² = Σ_μ comp[μ]²`,
  ///        and an electric / magnetic projector fills its spatial atom
  ///        `atomDen[invS] = |k⃗|² = Σ_{μ=1..3} comp[μ]²` (component 0 = temporal). The result is sized
  ///        to hold every `inv`/`invS` id seen (others are unused all-zero MPolys).
  NUMTRACER_FUNC std::vector<MPoly> collect_atom_denoms(int nsym, const std::vector<network::NetVal> &lors,
                                                        const std::vector<std::array<MPoly, 4>> &comp)
  {
    const auto isProj = [](const network::Elem &e) {
      return e.kind == network::Elem::ProjT || e.kind == network::Elem::ProjL ||
             e.kind == network::Elem::ProjE || e.kind == network::Elem::ProjM;
    };
    int maxId = -1;
    for (const network::NetVal &nv : lors)
      for (const network::PTerm &pt : nv)
        for (const network::Elem &e : pt.e)
          if (isProj(e)) {
            if (e.inv > maxId) maxId = e.inv;   // full 1/k² atom (ProjT/ProjL/ProjE)
            if (e.invS > maxId) maxId = e.invS; // spatial 1/|k⃗|² atom (ProjE/ProjM)
          }
    std::vector<MPoly> atomDen(maxId + 1, MPolyFactory::zero(nsym));
    // An atom id names ONE denominator. The front end is supposed to guarantee that, but the writes
    // below were unconditional `atomDen[id] = ...`, so two projectors sharing an id while carrying
    // different momenta silently gave the second one's k² to both — every cancellation against that
    // id then divides by the wrong polynomial. Record which ids have been written so a genuine
    // rewrite (same id, DIFFERENT denominator) is caught; a repeat of the identical denominator is
    // the normal case (the same projector appearing in many terms) and must stay silent.
    std::vector<char> written(static_cast<std::size_t>(maxId + 1), 0);
    // Exact term-wise equality is the right test here (there is no MPoly::operator==): both sides
    // are built by the same deterministic Σ_μ comp[μ]² over the same component table, so equal
    // momenta give bit-equal polynomials and any difference means genuinely different momenta.
    const auto sameDen = [](const MPoly &x, const MPoly &y) {
      if (x.terms.size() != y.terms.size()) return false;
      for (std::size_t i = 0; i < x.terms.size(); ++i)
        if (!(x.terms[i].first == y.terms[i].first) || x.terms[i].second.re != y.terms[i].second.re ||
            x.terms[i].second.im != y.terms[i].second.im)
          return false;
      return true;
    };
    const auto claim = [&](int id, MPoly &&den, const char *what) {
      auto &slot = atomDen[static_cast<std::size_t>(id)];
      if (written[static_cast<std::size_t>(id)]) {
        if (!sameDen(slot, den)) NT_THROW(std::runtime_error, what);
        return;
      }
      written[static_cast<std::size_t>(id)] = 1;
      slot = std::move(den);
    };
    for (const network::NetVal &nv : lors)
      for (const network::PTerm &pt : nv)
        for (const network::Elem &e : pt.e)
          if (isProj(e)) {
            const auto &cv = comp[e.vid];
            if (e.inv >= 0) { // full k² = Σ_μ comp[μ]²
              MPoly k2 = MPolyFactory::zero(nsym);
              for (int mu = 0; mu < 4; ++mu)
                k2 = k2 + cv[mu] * cv[mu];
              claim(e.inv, std::move(k2),
                    "collect_atom_denoms: two projectors share an `inv` atom id but carry different "
                    "momenta — one denominator would silently overwrite the other");
            }
            if (e.invS >= 0) { // spatial |k⃗|² = Σ_{μ=1..3} comp[μ]² (component 0 = temporal)
              MPoly ks2 = MPolyFactory::zero(nsym);
              for (int mu = 1; mu < 4; ++mu)
                ks2 = ks2 + cv[mu] * cv[mu];
              claim(e.invS, std::move(ks2),
                    "collect_atom_denoms: two projectors share an `invS` spatial atom id but carry "
                    "different momenta — one denominator would silently overwrite the other");
            }
          }
    return atomDen;
  }
#endif // NUMTRACER_DEFINE_BODIES

  /// @brief Relative noise-prune tolerance: a monomial whose |coefficient| is below this fraction of
  ///        the largest coefficient is round-off from the numeric frame (a ~10-order gap separates it
  ///        from physics), so @ref to_genprog drops it. Tuning this changes which monomials survive —
  ///        validate against the INTEGRATED numeric-vs-FORM error, not a pointwise round-off floor.
  inline constexpr double kNoisePruneRelTol = 1e-9;

  /// @brief Significant decimal digits kept when snapping a lowered coefficient (@ref snap_coeff).
  ///
  /// The numeric backend contracts components in `double`, so one correctly-rounded frame literal
  /// (`Sqrt[3]/2 -> 0.86602540378443865`, from the 120-degree symmetric point) is sprayed by ~1e4
  /// double operations into thousands of near-duplicates of the SAME mathematical constant:
  /// `16.000000000000007` and `16.000000000000004` both mean 16, reached by different association
  /// orders. That matters because the SSA CSE keys constants on RAW BITS
  /// (`codegen/real_cse.hpp` ihash/ieq), so each variant becomes its own RCONST, every downstream
  /// RMUL that consumes it also differs, and one dirty ulp at a leaf duplicates an entire subtree.
  /// The bit-exact folds (`k == 1.0`, `k == 0.0`) miss for the same reason, emitting free multiplies.
  ///
  /// MEASURED tradeoff on ZAqbq1_147 (Mq in), 200k points, vs the FormTracer oracle. Note the trace
  /// sums cancel heavily, so a coefficient perturbation is amplified ~1e6 in the result — do NOT
  /// pick this from the literal-collapse count alone:
  ///
  ///   digits | multiplies | distinct lits | ns/eval | NT/FORM | rel-err vs FORM
  ///   -------|------------|---------------|---------|---------|----------------
  ///   off    |     29,395 |          2931 |    3508 |  2.90x  | 4.26e-09
  ///   15     |     27,049 |           740 |    3371 |  2.75x  | 2.98e-09
  ///   14     |     26,096 |           378 |    3283 |  2.64x  | 6.43e-09   <- default
  ///   12     |     25,988 |           338 |    3220 |  2.63x  | 4.55e-06   <- 1000x accuracy loss
  ///
  /// 14 captures essentially the whole speed win at no accuracy cost; 12 buys a further 2% for three
  /// orders of magnitude of accuracy, which is a bad trade. Override with NT_GEN_SNAP_DIGITS
  /// (0 = disable). Re-measure this table before changing the default.
  /// @brief `NT_GEN_POLYSTATS` verbosity: 0 off, 1 the summary counts, 2 the per-monomial key dump.
  ///
  /// Level 2 is a DIFFERENT report, not a superset of level 1: the dump is meant to be piped through
  /// `sort -u` to count distinct monomials across traces, and the summary lines would corrupt that.
  inline int polystats_level()
  {
    static const int lvl = static_cast<int>(env_int("NT_GEN_POLYSTATS", 0));
    return lvl;
  }

  inline constexpr int kCoeffSnapDigits = 14;

  /// @brief Round @p v to @ref kCoeffSnapDigits significant decimal digits.
  ///
  /// Applied at the single point where a polynomial coefficient becomes a codegen `LMono::c` — i.e.
  /// AFTER all polynomial arithmetic and after the noise-prune above — so it cannot perturb a
  /// cancellation, only canonicalise what survived. Perturbation is ~1e-12 relative, three orders
  /// below the numeric-vs-FORM correctness gate.
  ///
  /// It is also the chokepoint where a NON-FINITE coefficient is caught. This function used to pass
  /// NaN/Inf through untouched, and nothing downstream stopped them: `RBuilder::ieq`
  /// (`codegen/real_cse.hpp`) compares constants with `k == k`, so a NaN never compares equal to
  /// itself and `find_or_add` appends a fresh SSA slot on every single call — unbounded duplicate
  /// instructions — and the printer then writes a literal `nan` into a COMMITTED kernel header. A
  /// non-finite coefficient is always an upstream bug (a division by an identically-zero denominator,
  /// an uninitialised component), never something to round, so it aborts here where the trace is
  /// still identifiable rather than becoming a mysterious `nan` in generated source.
  inline double snap_coeff(double v)
  {
    static const int digits = [] {
      // 0 is a meaningful value here (snapping off), so "unset" is the empty/absent case, not 0;
      // a value outside [0,17] is a caller error and falls back to the default rather than
      // silently truncating to nothing.
      const long d = env_int("NT_GEN_SNAP_DIGITS", kCoeffSnapDigits);
      return (d >= 0 && d <= 17) ? static_cast<int>(d) : kCoeffSnapDigits;
    }();
    if (!std::isfinite(v))
      NT_THROW(std::runtime_error,
               "snap_coeff: non-finite (NaN/Inf) polynomial coefficient reached the lowering — it "
               "would be emitted as a literal `nan` and would defeat the SSA constant CSE");
    if (digits == 0 || v == 0.0) return v;
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.*e", digits - 1, v);
    return std::strtod(buf, nullptr);
  }

  NUMTRACER_FUNC network::GenProg to_genprog(const MPoly &p, network::GlobalEnv &g, bool realOnly = false);
  NUMTRACER_FUNC network::GenProg to_genprog(const DPoly &p, network::GlobalEnv &g, bool realOnly = false);
  NUMTRACER_FUNC std::vector<network::FusedProg> to_genprog_fused(const std::vector<MPoly> &ps,
                                                                  network::GlobalEnv &g,
                                                                  const std::vector<int> &realOnly);
  NUMTRACER_FUNC std::vector<network::FusedProg> to_genprog_fused(const std::vector<DPoly> &ps,
                                                                  network::GlobalEnv &g,
                                                                  const std::vector<int> &realOnly);

  /// @brief Incremental form of @ref to_genprog_fused: lower one trace polynomial at a time.
  ///
  /// Same output as handing the whole `std::vector<Poly>` over at once — @ref to_genprog_fused is
  /// literally a loop over `add` — but the caller never has to build that vector. That is what lets
  /// the generator's streaming phase B (`trace_fold.hpp`'s `fold_groups_streaming`) lower each trace
  /// group the moment it is folded and then drop it, instead of accumulating every group's polynomial
  /// first. Traces are added in ascending order and flushed into a new chunk every
  /// `NT_GEN_CC_CHUNK` traces, so the chunk boundaries (and hence `FusedProg::offset`) are identical
  /// to the batch version's `[base, hi)` slices.
  ///
  /// `add` is deliberately two overloads rather than a template: the definitions live in the library
  /// TU alongside `lower_into`, and MPoly/DPoly are the only two backends — the same split every other
  /// entry point here uses (see core/export.hpp).
  class FusedStream
  {
  public:
    NUMTRACER_FUNC FusedStream(network::GlobalEnv &g, const std::vector<int> &realOnly);
    NUMTRACER_FUNC void add(const MPoly &p);
    NUMTRACER_FUNC void add(const DPoly &p);
    /// Flushes the tail (if any). An empty stream yields an empty vector, not one empty chunk —
    /// `emit_cpp_fused` would otherwise emit a bodiless `trace_all_c0` that stores nothing.
    NUMTRACER_FUNC std::vector<network::FusedProg> finish();

  private:
    template <class Poly> void add_(const Poly &p);
    void begin_();
    void flush_();

    network::GlobalEnv *g_;
    const std::vector<int> *ro_;
    std::size_t step_;
    std::size_t n_ = 0;     ///< global trace index of the next add()
    std::size_t total_ = 0; ///< summed SSA instruction count, for the polystats line
    network::rdetail::RBuilder w_;
    network::FusedProg fp_;
    std::vector<network::FusedProg> out_;
  };

#if NUMTRACER_DEFINE_BODIES
  /// @brief Lower a contracted diagram polynomial into the shared env via the existing CSE + Horner
  ///        back-end (`gdetail::best_into`). User symbols intern as env kind 3 (`var`), surviving
  ///        inverse atoms as kind 1 (`inv`). Returns an @ref network::GenProg (real, or complex when the
  ///        polynomial carries an imaginary coefficient).
  /// @param realOnly when true, emit a REAL program even if the polynomial has imaginary coefficients
  ///        — the caller has proven only `Re(this trace)` is consumed (its assembly coefficient is
  ///        real). `Re(Σ c·mono) = Σ Re(c)·mono` for real monomials, so the imaginary half is dead;
  ///        skipping it avoids computing+returning a `std::complex` whose `.imag()` nobody reads.
  /// @brief Lower ONE polynomial into an EXISTING CSE builder; returns `{reRoot, imRoot}` with
  ///        `imRoot == network::kRealProgram` when the polynomial is real.
  ///
  /// Split out of @ref to_genprog so the same lowering can be replayed into a builder that already
  /// holds other traces (@ref to_genprog_fused / `CrossTraceCSE`). `to_genprog` is now this function
  /// plus a fresh builder, so single-trace output is byte-identical to before the split.
  NUMTRACER_FUNC std::pair<int, int> lower_into(const MPoly &p, network::GlobalEnv &g,
                                                network::rdetail::RBuilder &builder, bool realOnly)
  {
    // Prune numerically-zero monomials. A NUMERIC frame fixes the external momenta to concrete
    // components, so exact (analytic) cancellations in the trace surface as tiny residual coefficients
    // (~1e-12 … 1e-30 against real coefficients of O(10²)) — pure round-off, not physics. They don't
    // affect the value (the kernel matches FORM to ~1e-13) but each spurious monomial costs runtime
    // arithmetic; on the dense 1/4/7 trace ~half the monomials are such noise. Drop |c| below a RELATIVE
    // tolerance vs the largest coefficient (a clean ~10-order gap separates noise from real terms).
    double maxabs = 0.0;
    for (const auto &[m, c] : p.terms)
      maxabs = std::max(maxabs, std::max(std::fabs(c.re), std::fabs(c.im)));
    const double tol = kNoisePruneRelTol * maxabs;
    auto keep = [&](const Cx &c) { return std::max(std::fabs(c.re), std::fabs(c.im)) >= tol; };
    bool cplx = false;
    if (!realOnly)
      for (const auto &[m, c] : p.terms)
        if (keep(c) && (c.im > tol || c.im < -tol)) {
          cplx = true;
          break;
        }
    // ONE walk builds BOTH halves of a complex trace. The old shape — a `build(imag)` lambda called
    // once for re and again for im — re-derived and re-sorted every monomial's vp (env interning,
    // atom run-lengths, snap) twice for complex traces, doubling the flatten. The im list must
    // contain an entry for EVERY kept monomial (with snap(c.im), possibly 0.0), exactly as the
    // second pass produced; and the env ids intern in the same first-seen sequence as the old re
    // pass, so the emitted kernel is byte-identical.
    std::vector<network::LMono> monosRe, monosIm;
    {
      monosRe.reserve(p.terms.size());
      if (cplx) monosIm.reserve(p.terms.size());
      for (const auto &[m, c] : p.terms) {
        if (!keep(c)) continue;
        std::vector<std::pair<int, int>> vp;
        for (int k = 0; k < p.nsym; ++k)
          if (m.e[k] > 0) vp.push_back({g.var_id(k), m.e[k]});
        for (int i = 0; i < (int)m.atoms.size();) {
          int j = i;
          while (j < (int)m.atoms.size() && m.atoms[j] == m.atoms[i])
            ++j;
          vp.push_back({g.inv_id(m.atoms[i]), j - i});
          i = j;
        }
        std::sort(vp.begin(), vp.end());
        if (cplx) {
          network::LMono lmIm;
          lmIm.c = snap_coeff(c.im);
          lmIm.vp = vp; // copy — the re half moves it below
          monosIm.push_back(std::move(lmIm));
        }
        network::LMono lm;
        lm.c = snap_coeff(c.re);
        lm.vp = std::move(vp);
        monosRe.push_back(std::move(lm));
      }
    }
    {
      const auto &monos = monosRe;
      // NT_GEN_POLYSTATS=1: report the MONOMIAL count handed to the Horner lowering. Monomial count
      // is canonical (independent of any factorisation strategy), so comparing it against the
      // emitted multiply count says whether a large op count comes from the ALGEBRA (many monomials)
      // or from weak LOWERING (few monomials, many ops).
      if (polystats_level() == 1)
        std::fprintf(stderr, "[polystats] mpoly terms=%zu kept=%zu nsym=%d\n", p.terms.size(),
                     monos.size(), p.nsym);
      {
        // NT_GEN_POLYSTATS=2: additionally dump each monomial's KEY (variable powers + inv/dressing
        // atoms). Piping through `sort -u` then answers: how many of the monomials summed over all
        // traces are actually DISTINCT? FormTracer sums every diagram into ONE polynomial before
        // expanding, so identical monomials from different diagrams collect; NumTracer keeps one
        // polynomial per trace, so they cannot. That is term COLLECTION, not CSE -- the compiler
        // can never do it (it would be a floating-point reassociation across function boundaries).
        if (polystats_level() == 2)
          for (const auto &lm : monos) {
            std::string key;
            for (const auto &[id, pw] : lm.vp)
              key += std::to_string(id) + "^" + std::to_string(pw) + " ";
            std::fprintf(stderr, "[mono] %s\n", key.c_str());
          }
      }
    }
    const int reRoot = network::gdetail::best_into(std::move(monosRe), builder);
    if (!cplx) return {reRoot, network::kRealProgram};
    return {reRoot, network::gdetail::best_into(std::move(monosIm), builder)};
  }

  NUMTRACER_FUNC network::GenProg to_genprog(const MPoly &p, network::GlobalEnv &g, bool realOnly)
  {
    network::rdetail::RBuilder builder;
    const auto [reRoot, imRoot] = lower_into(p, g, builder, realOnly);
    network::GenProg gp{std::move(builder.ins), reRoot};
    gp.rootIm = imRoot;
    return gp;
  }

  /// @brief Lower a dressed-diagram @ref DPoly into the shared env. Each kinematic monomial is emitted
  ///        exactly as in the @ref MPoly overload (user symbols → kind-3 `var`, surviving inverse atoms
  ///        → kind-1 `inv`) and additionally carries the dressing monomial's atoms as kind-2 `dress`
  ///        leaves (@ref network::GlobalEnv::dr_id). The whole `DPoly` is flattened into ONE monomial
  ///        list so the shared CSE/Horner (`gdetail::best_into`) collects the dressing factors across
  ///        monomials — FormTracer-parity collection in one trace function. A `DPoly` with a single
  ///        empty dressing monomial reduces to exactly the @ref MPoly path (no `dress` leaves).
  /// @brief @ref DPoly counterpart of the @ref MPoly `lower_into` — lower into an existing builder.
  NUMTRACER_FUNC std::pair<int, int> lower_into(const DPoly &p, network::GlobalEnv &g,
                                                network::rdetail::RBuilder &builder, bool realOnly)
  {
    // PER-DRESSING-CHANNEL noise prune (not one global tolerance). A DPoly is Σ_d (dressing_d)·(kinematic_d);
    // each channel `d` is reweighted at runtime by its dressing product, which can swing by many orders across
    // the loop domain. A single global tolerance, taken from the largest-coefficient channel, would delete a
    // small channel's GENUINE terms — fine while that channel is runtime-small, but wrong the moment a dressing
    // suppresses the large one (dr→0) and the small one should dominate. So prune each channel against ITS OWN
    // max, exactly as the @ref MPoly overload does and as the distributed (collection-off) path already does
    // per trace. `dmonoTol(mp) ≤` any global tol, so this only ever KEEPS more — strictly safer, and identical
    // when all channels share a scale.
    auto dmonoTol = [](const MPoly &mp) {
      double maxabs = 0.0;
      for (const auto &[m, c] : mp.terms)
        maxabs = std::max(maxabs, std::max(std::fabs(c.re), std::fabs(c.im)));
      return kNoisePruneRelTol * maxabs;
    };
    auto keepAt = [](const Cx &c, double tol) { return std::max(std::fabs(c.re), std::fabs(c.im)) >= tol; };
    bool cplx = false;
    if (!realOnly)
      for (const auto &[d, mp] : p.terms) {
        const double tol = dmonoTol(mp);
        for (const auto &[m, c] : mp.terms)
          if (keepAt(c, tol) && (c.im > tol || c.im < -tol)) {
            cplx = true;
            break;
          }
        if (cplx) break;
      }
    // ONE walk builds BOTH halves — see the MPoly overload: the old build(imag) lambda re-derived
    // every monomial's vp (dress/var/inv interning + sort) twice for complex traces. Interning
    // order and im-list contents (an entry per kept monomial, snap(c.im), possibly 0.0) reproduce
    // the old two-pass shape exactly, so the emitted kernel is byte-identical.
    std::vector<network::LMono> monosRe, monosIm;
    {
      for (const auto &[d, mp] : p.terms) {
        const double tol = dmonoTol(mp);
        // the dressing monomial's atoms become kind-2 `dress` env leaves shared by every monomial of mp
        std::vector<std::pair<int, int>> drvp;
        for (int i = 0; i < (int)d.size();) {
          int j = i;
          while (j < (int)d.size() && d[j] == d[i])
            ++j;
          drvp.push_back({g.dr_id(d[i]), j - i});
          i = j;
        }
        for (const auto &[m, c] : mp.terms) {
          if (!keepAt(c, tol)) continue;
          std::vector<std::pair<int, int>> vp = drvp;
          for (int k = 0; k < mp.nsym; ++k)
            if (m.e[k] > 0) vp.push_back({g.var_id(k), m.e[k]});
          for (int i = 0; i < (int)m.atoms.size();) {
            int j = i;
            while (j < (int)m.atoms.size() && m.atoms[j] == m.atoms[i])
              ++j;
            vp.push_back({g.inv_id(m.atoms[i]), j - i});
            i = j;
          }
          std::sort(vp.begin(), vp.end());
          if (cplx) {
            network::LMono lmIm;
            lmIm.c = snap_coeff(c.im);
            lmIm.vp = vp; // copy — the re half moves it below
            monosIm.push_back(std::move(lmIm));
          }
          network::LMono lm;
          lm.c = snap_coeff(c.re);
          lm.vp = std::move(vp);
          monosRe.push_back(std::move(lm));
        }
      }
      // see the MPoly overload: monomial count is the optimizer-independent baseline against which
      // the emitted instruction count is judged.
      {
        if (polystats_level() == 1) std::fprintf(stderr, "[polystats] dpoly monos=%zu\n", monosRe.size());
        if (polystats_level() == 2)
          for (const auto &lm : monosRe) {
            std::string key;
            for (const auto &[id, pw] : lm.vp)
              key += std::to_string(id) + "^" + std::to_string(pw) + " ";
            std::fprintf(stderr, "[mono] %s\n", key.c_str());
          }
      }
    }
    const int reRoot = network::gdetail::best_into(std::move(monosRe), builder);
    if (!cplx) {
      if (polystats_level() == 1) std::fprintf(stderr, "[polystats] ssa instrs=%zu\n", builder.ins.size());
      return {reRoot, network::kRealProgram};
    }
    const int imRoot = network::gdetail::best_into(std::move(monosIm), builder);
    if (polystats_level() == 1) std::fprintf(stderr, "[polystats] ssa instrs=%zu (re+im)\n", builder.ins.size());
    return {reRoot, imRoot};
  }

  NUMTRACER_FUNC network::GenProg to_genprog(const DPoly &p, network::GlobalEnv &g, bool realOnly)
  {
    network::rdetail::RBuilder builder;
    const auto [reRoot, imRoot] = lower_into(p, g, builder, realOnly);
    network::GenProg gp{std::move(builder.ins), reRoot};
    gp.rootIm = imRoot;
    return gp;
  }

  /// @brief Lower EVERY trace group into one shared CSE builder (`CrossTraceCSE`), so a subexpression
  ///        reached by several traces is emitted once. See @ref network::FusedProg.
  ///
  /// Measured on ZAqbq1_147 Mq-in (108 traces, 54 complex): the shared stream is **30,547** SSA
  /// instructions against **47,558** for the same traces lowered independently — 0.64x, at unchanged
  /// lowering cost (0.30 s either way). The saving is bounded by the fact that each of the 28,856
  /// monomial OCCURRENCES still needs its own accumulate into its own trace; only the products are
  /// shared. So the duplication factor (3.54x distinct-vs-total monomials) is an upper bound that is
  /// nowhere near attainable — do not quote it as the expected speedup.
  ///
  /// NOTE `network::gdetail::best_into` costs its candidate Horner orderings on scratch builders that
  /// start EMPTY, so it cannot see CSE hits against the already-populated `builder`: trace k's ordering is
  /// chosen as if traces 0..k-1 did not exist. That caps the achievable sharing (mostly moot — the
  /// sweep collapses to a single ordering above 2000 monomials).
  /// Traces per fused program. 0 (the default) = fuse everything into one. `NT_GEN_CC_CHUNK` overrides.
  ///
  /// Fusing all traces maximises sharing but builds one enormous basic block that SPILLS: on
  /// ZAqbq1_147 Mq-in the fully-fused kernel executes ~8,200 instructions/eval more than its own
  /// arithmetic op count, against ~2,900 for the unfused baseline — i.e. fusion bought ~17k ops and
  /// handed back ~5,400 instructions of spill traffic, which is the measured IPC drop (2.405 -> 2.278).
  /// Chunking trades a little cross-trace sharing back for register pressure.
  inline int cc_chunk_size()
  {
    static const int n = [] {
      const long v = env_int("NT_GEN_CC_CHUNK", 0);
      return v > 0 ? static_cast<int>(v) : 0; // 0 = fuse every trace into one program
    }();
    return n;
  }

  // The batch entry point is a thin loop over FusedStream::add, so the streaming and non-streaming
  // lowerings cannot drift apart: they are the same code.
  template <class Poly>
  NUMTRACER_FUNC std::vector<network::FusedProg> to_genprog_fused_impl(const std::vector<Poly> &ps,
                                                                       network::GlobalEnv &g,
                                                                       const std::vector<int> &realOnly)
  {
    FusedStream fs(g, realOnly);
    for (const auto &p : ps)
      fs.add(p);
    return fs.finish();
  }

  FusedStream::FusedStream(network::GlobalEnv &g, const std::vector<int> &realOnly)
      : g_(&g), ro_(&realOnly),
        step_(cc_chunk_size() > 0 ? static_cast<std::size_t>(cc_chunk_size())
                                  : std::numeric_limits<std::size_t>::max())
  {
    begin_();
  }

  template <class Poly> void FusedStream::add_(const Poly &p)
  {
    const bool ro = n_ < ro_->size() && (*ro_)[n_] != 0;
    const auto [re, im] = lower_into(p, *g_, w_, ro);
    fp_.root.push_back(re);
    fp_.rootIm.push_back(im);
    ++n_;
    if (fp_.root.size() == step_) {
      flush_();
      begin_();
    }
  }
  void FusedStream::add(const MPoly &p) { add_(p); }
  void FusedStream::add(const DPoly &p) { add_(p); }

  std::vector<network::FusedProg> FusedStream::finish()
  {
    if (!fp_.root.empty()) flush_();
    if (polystats_level() == 1)
      std::fprintf(stderr, "[polystats] FUSED ssa instrs=%zu over %zu traces in %zu chunk(s)\n", total_,
                     n_, out_.size());
    return std::move(out_);
  }

  void FusedStream::begin_()
  {
    w_ = network::rdetail::RBuilder{};
    fp_ = network::FusedProg{};
    fp_.offset = static_cast<int>(n_);
  }
  void FusedStream::flush_()
  {
    total_ += w_.ins.size();
    fp_.ins = std::move(w_.ins);
    out_.push_back(std::move(fp_));
  }

  NUMTRACER_FUNC std::vector<network::FusedProg> to_genprog_fused(const std::vector<MPoly> &ps,
                                                                  network::GlobalEnv &g,
                                                                  const std::vector<int> &realOnly)
  {
    return to_genprog_fused_impl(ps, g, realOnly);
  }
  NUMTRACER_FUNC std::vector<network::FusedProg> to_genprog_fused(const std::vector<DPoly> &ps,
                                                                  network::GlobalEnv &g,
                                                                  const std::vector<int> &realOnly)
  {
    return to_genprog_fused_impl(ps, g, realOnly);
  }
#endif // NUMTRACER_DEFINE_BODIES

} // namespace numtracer::numeric
