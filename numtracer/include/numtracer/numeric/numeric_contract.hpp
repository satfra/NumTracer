/// @file numeric_contract.hpp
/// @brief The numeric (matrix-product) contraction engine (task #22): fold a diagram's Dirac trace
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
#include <stdexcept>
#include <vector>

namespace numtracer::numeric
{


  /// @brief A pure-Lorentz network factor (numeric counterpart of @ref network::Elem).
  ///   - kind 0 `nmet(a,b)`            : metric δ_{ab}
  ///   - kind 1 `nvec(a, vlc)`         : vector leg `Σ coeff·comp(vid)` on Lorentz id `a`
  ///   - kind 2 `nprojT(a,b, vlc, atom)`: transverse projector `P_T(k)_{ab}=δ_{ab}−k_a k_b·INV(k)`,
  ///                                     `k=Σ coeff·comp(vid)`, `INV(k)` = inverse atom `atom`
  ///   - kind 3 `neps(a,b,c,d)`        : Levi-Civita ε_{abcd}
  ///   - kind 4 `nprojL(a,b, vlc, atom)`: longitudinal projector `P_L(k)_{ab}=k_a k_b·INV(k)`
  ///   - kind 5 `nprojE(a,b, vlc, atom, atomS)`: finite-T electric projector `P_E = P_T − P_M`
  ///   - kind 6 `nprojM(a,b, vlc, atomS)`: finite-T magnetic projector `P_M_{ij}=δ_{ij}−k_i k_j·INVS(k)`
  ///     (i,j spatial; rows/cols of the temporal component 0 vanish). `INVS(k)=1/|k⃗|²` is atom `atomS`.
  struct NElem {
    enum Kind { Metric, Vector, Epsilon, ProjT, ProjL, ProjE, ProjM };
    const Kind kind = Metric;
    const int a = 0, b = 0, c = 0, d = 0;
    const std::vector<std::pair<double, int>> vlc; ///< momentum (kind 1 vector / kind 2/4/5/6 projector)
    const int atom = -1;                           ///< full inverse-atom id `1/k²` (kind 2/4/5)
    const int atomS = -1;                          ///< spatial inverse-atom id `1/|k⃗|²` (kind 5/6)
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
      std::array<MPoly, 4> r = {MPoly(nsym), MPoly(nsym), MPoly(nsym), MPoly(nsym)};
      for (const auto &[coeff, vid] : vlc) {
        const auto &cv = comp[vid];
        const MPoly s = MPoly::constant(nsym, Cx{coeff, 0});
        for (int mu = 0; mu < 4; ++mu)
          r[mu] = r[mu] + s * cv[mu];
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
      std::vector<MPoly> v; ///< size 4^ids.size()
    };

  } // namespace ndetail

  /// @brief Fold the closed Dirac chain into a tensor over its free gluon legs (the `dgamma` ids).
  ///        For each assignment of the free legs to concrete indices 0..3, build the slashed/free
  ///        γ chain as 4×4 @ref MPoly matrices and take the trace. Returns the free-leg ids and the
  ///        row-major tensor of trace polynomials (one entry per `4^f` assignment).
  inline ndetail::Factor numeric_dirac(int nsym, const network::DiracNet &chain,
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
    std::vector<Mat4> slashMat;      // precomputed slash matrices, indexed by token (kind 1)
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
    // trace parity is set by the BLOCK-ANTIDIAGONAL factors only (γ^μ kind 0, slash kind 1); γ5 (kind 2)
    // and σ (kind 3) are block-DIAGONAL. An odd antidiagonal count → final product antidiagonal → tr=0.
    std::size_t ng = 0;
    for (const network::DFac &d : chain)
      if (d.kind == network::DFac::Gamma || d.kind == network::DFac::Slash) ++ng;
    ndetail::Factor F;
    F.ids = freeLegs;
    int total = 1;
    for (int k = 0; k < f; ++k)
      total *= 4;
    F.v.assign(total, MPoly(nsym));
    if (ng % 2 == 1) return F; // odd → all zero

    // CHIRAL 2×2 BLOCK fold (4× fewer MPoly multiplies than the full 4×4, tr(odd)=0 free): every γ^μ and
    // slash is block-antidiagonal in the Weyl basis, so the running product is two 2×2 blocks whose roles
    // alternate. Precompute each factor's (P=upper-right, Q=lower-left) 2×2 blocks once.
    using B2 = std::array<MPoly, 4>; // row-major 2×2
    auto blocksOf = [&](const Mat4 &M, B2 &P, B2 &Q) {
      for (int rr = 0; rr < 2; ++rr)
        for (int cc = 0; cc < 2; ++cc) {
          P[rr * 2 + cc] = M.a[rr][cc + 2];
          Q[rr * 2 + cc] = M.a[rr + 2][cc];
        }
    };
    auto mul2 = [&](const B2 &x, const B2 &y) {
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
    // γ5 = diag(+I,−I) in the Weyl basis: block-DIAGONAL, so it never toggles the antidiag parity and
    // costs no matrix multiply — it just negates one running block (the upper when the product is
    // antidiagonal, the lower when diagonal; derived from R·diag(+I,−I)). A leading γ5 seeds the
    // product as the diagonal diag(+I,−I) itself.
    const B2 id2 = {MPoly::constant(nsym, Cx{1, 0}), MPoly(nsym), MPoly(nsym), MPoly::constant(nsym, Cx{1, 0})};
    auto neg2 = [&](const B2 &x) {
      return B2{MPoly(nsym) - x[0], MPoly(nsym) - x[1], MPoly(nsym) - x[2], MPoly(nsym) - x[3]};
    };
    // BARE commutator [A,B] = A·B − B·A is block-DIAGONAL in the Weyl basis (A,B antidiagonal ⇒ A·B
    // diagonal): upper = P_a Q_b − P_b Q_a, lower = Q_a P_b − Q_b P_a, where (P_a,Q_a)/(P_b,Q_b) are the
    // (upper-right, lower-left) blocks of leg A / leg B. Like γ5, a diagonal factor — it never toggles the
    // antidiag parity. The σ^{μν}=(i/2)[γ^μ,γ^ν] normalization (i/2 and sign) lives in the emitted SCALAR.
    auto subB2 = [&](const B2 &x, const B2 &y) { return B2{x[0] - y[0], x[1] - y[1], x[2] - y[2], x[3] - y[3]}; };
    auto commBlocks = [&](const B2 &Pa, const B2 &Qa, const B2 &Pb, const B2 &Qb, B2 &Su, B2 &Sl) {
      Su = subB2(mul2(Pa, Qb), mul2(Pb, Qa));
      Sl = subB2(mul2(Qa, Pb), mul2(Qb, Pa));
    };
    std::vector<int> legComp(f, 0); // concrete index 0..3 assigned to each free leg
    for (int flat = 0; flat < total; ++flat) {
      int r = flat;
      for (int k = f - 1; k >= 0; --k) {
        legComp[k] = r % 4;
        r /= 4;
      }
      B2 m0, m1;
      bool started = false, antidiag = true;
      for (std::size_t i = 0; i < chain.size(); ++i) {
        if (chain[i].kind == network::DFac::Gamma5) { // γ5: block-diagonal diag(+I,−I), parity unchanged
          if (!started) {
            m0 = id2;
            m1 = neg2(id2);
            started = true;
            antidiag = false;
            continue;
          }
          if (antidiag)
            m0 = neg2(m0);
          else
            m1 = neg2(m1);
          continue;
        }
        if (chain[i].kind == network::DFac::Comm) { // commutator [A,B]: block-diagonal diag(Su,Sl), parity unchanged
          B2 Su, Sl;
          // leg-A blocks: free ⇒ gP/gQ at this assignment's component; slash ⇒ precomputed cAP/cAQ.
          const B2 &Pa = commAslash[i] ? cAP[i] : gP[legComp[tokenFree[i]]];
          const B2 &Qa = commAslash[i] ? cAQ[i] : gQ[legComp[tokenFree[i]]];
          const B2 &Pb = commBslash[i] ? cBP[i] : gP[legComp[tokenFree2[i]]];
          const B2 &Qb = commBslash[i] ? cBQ[i] : gQ[legComp[tokenFree2[i]]];
          commBlocks(Pa, Qa, Pb, Qb, Su, Sl);
          if (!started) {
            m0 = Su;
            m1 = Sl;
            started = true;
            antidiag = false;
            continue;
          }
          if (antidiag) {
            m0 = mul2(m0, Sl);
            m1 = mul2(m1, Su);
          } // R·diag(Su,Sl), R antidiagonal
          else {
            m0 = mul2(m0, Su);
            m1 = mul2(m1, Sl);
          }
          continue;
        }
        const std::size_t mu = (chain[i].kind == network::DFac::Gamma) ? legComp[tokenFree[i]] : 0;
        const B2 &Pi = (chain[i].kind == network::DFac::Gamma) ? gP[mu] : sP[i];
        const B2 &Qi = (chain[i].kind == network::DFac::Gamma) ? gQ[mu] : sQ[i];
        if (!started) {
          m0 = Pi;
          m1 = Qi;
          started = true;
          antidiag = true;
          continue;
        }
        if (antidiag) {
          m0 = mul2(m0, Qi);
          m1 = mul2(m1, Pi);
        } // [[0,m0],[m1,0]]·[[0,Pi],[Qi,0]]
        else {
          m0 = mul2(m0, Pi);
          m1 = mul2(m1, Qi);
        }
        antidiag = !antidiag;
      }
      // ng (γ-count, γ5 excluded) even → block-diagonal product → trace = tr(m0)+tr(m1).
      F.v[flat] = m0[0] + m0[3] + m1[0] + m1[3];
    }
    return F;
  }

  namespace ndetail
  {

    /// Build the dense factor for one @ref NElem (metric / vector / projector / Levi-Civita).
    inline Factor elem_factor(int nsym, const NElem &el, const std::vector<std::array<MPoly, 4>> &comp)
    {
      Factor F;
      if (el.kind == NElem::Metric) { // metric δ_{a,b}
        F.ids = {el.a, el.b};
        F.v.assign(16, MPoly(nsym));
        for (int i = 0; i < 4; ++i)
          F.v[i * 4 + i] = MPoly::constant(nsym, Cx{1, 0});
      } else if (el.kind == NElem::Vector) { // vector on id a
        F.ids = {el.a};
        auto cc = mom_components(nsym, el.vlc, comp);
        F.v = {cc[0], cc[1], cc[2], cc[3]};
      } else if (el.kind == NElem::ProjT) { // transverse projector P_T(k)_{a,b} = δ − k_a k_b INV(k)
        F.ids = {el.a, el.b};
        F.v.assign(16, MPoly(nsym));
        auto k = mom_components(nsym, el.vlc, comp);
        const MPoly at = MPoly::atom(nsym, el.atom);
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j) {
            MPoly e = (i == j) ? MPoly::constant(nsym, Cx{1, 0}) : MPoly(nsym);
            e = e - (k[i] * k[j]) * at;
            F.v[i * 4 + j] = std::move(e);
          }
      } else if (el.kind == NElem::ProjL) { // longitudinal projector P_L(k)_{a,b} = k_a k_b INV(k)
        F.ids = {el.a, el.b};
        F.v.assign(16, MPoly(nsym));
        auto k = mom_components(nsym, el.vlc, comp);
        const MPoly at = MPoly::atom(nsym, el.atom);
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
            F.v[i * 4 + j] = (k[i] * k[j]) * at;
      } else if (el.kind == NElem::ProjM) { // magnetic projector P_M_{i,j}=δ_{ij}−k_i k_j INVS(k) (spatial)
        F.ids = {el.a, el.b};
        F.v.assign(16, MPoly(nsym));
        auto k = mom_components(nsym, el.vlc, comp); // spatial part: zero the temporal (slot 0) component
        const MPoly atS = MPoly::atom(nsym, el.atomS);
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j) {
            MPoly e = (i == j && i > 0) ? MPoly::constant(nsym, Cx{1, 0}) : MPoly(nsym);
            if (i > 0 && j > 0) e = e - (k[i] * k[j]) * atS;
            F.v[i * 4 + j] = std::move(e);
          }
      } else if (el.kind == NElem::ProjE) { // electric projector P_E = P_T − P_M
        F.ids = {el.a, el.b};
        F.v.assign(16, MPoly(nsym));
        auto k = mom_components(nsym, el.vlc, comp);
        const MPoly at = MPoly::atom(nsym, el.atom);
        const MPoly atS = MPoly::atom(nsym, el.atomS);
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j) {
            MPoly full = (i == j) ? MPoly::constant(nsym, Cx{1, 0}) : MPoly(nsym);
            full = full - (k[i] * k[j]) * at; // P_T entry
            MPoly mag = (i == j && i > 0) ? MPoly::constant(nsym, Cx{1, 0}) : MPoly(nsym);
            if (i > 0 && j > 0) mag = mag - (k[i] * k[j]) * atS; // P_M entry
            F.v[i * 4 + j] = full - mag;
          }
      } else { // Levi-Civita ε_{a,b,c,d}
        F.ids = {el.a, el.b, el.c, el.d};
        F.v.assign(256, MPoly(nsym));
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
            for (int p = 0; p < 4; ++p)
              for (int q = 0; q < 4; ++q) {
                const double s = levi(i, j, p, q);
                if (s != 0.0) F.v[((i * 4 + j) * 4 + p) * 4 + q] = MPoly::constant(nsym, Cx{s, 0});
              }
      }
      return F;
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
      std::vector<int> unionIds;
      for (const Factor &F : group)
        for (int id : F.ids) {
          bool seen = false;
          for (int u : unionIds)
            if (u == id) {
              seen = true;
              break;
            }
          if (!seen) unionIds.push_back(id);
        }
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
      out.v.assign(outTotal, MPoly(nsym));
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
        MPoly acc(nsym);
        for (int elimVal = 0; elimVal < 4; ++elimVal) {
          idxVal[elimPos] = elimVal;
          MPoly prod = MPoly::constant(nsym, Cx{1, 0});
          bool zero = false;
          for (std::size_t fIdx = 0; fIdx < group.size() && !zero; ++fIdx) {
            int idx = 0;
            for (int p : slotPos[fIdx])
              idx = idx * 4 + idxVal[p];
            const MPoly &e = group[fIdx].v[idx];
            if (e.empty()) {
              zero = true;
              break;
            }
            prod = prod * e;
          }
          if (!zero) acc = acc + prod;
        }
        // Reduce the intermediate IN PLACE: cancel bare-loop atoms (sin²→1-cos² makes k²=l1² a monomial)
        // and collapse sin² powers, so projector `k⊗k·INV` factors don't bloat the running tensor across
        // a long pure-gauge chain (the ZA4 monster nets) — only cleaning at the end is too late.
        if (!units.empty()) acc = divThroughMonomialAtoms(reduce_units(acc, units), atomDen);
        out.v[outFlat] = std::move(acc);
      }
      return out;
    }

    /// Contract a set of dense factors over their shared Lorentz ids into one scalar MPoly, by GREEDY
    /// VARIABLE ELIMINATION: repeatedly pick the id whose incident factors have the smallest combined
    /// id-set (min-degree heuristic), merge those factors, and sum that id out. This keeps every
    /// intermediate tensor treewidth-bounded — a single global 4^G sum (G = all ids) is exponential and
    /// blows up for the quark triangle (G≈14); elimination stays at 4^(few) per step.
    inline MPoly contract_factors(int nsym, std::vector<Factor> facs, const std::vector<MPoly> &atomDen = {},
                                  const std::vector<std::vector<int>> &units = {})
    {
      // Algorithm (greedy min-fill-in elimination), repeated until no shared ids remain:
      //   1. collect the distinct live Lorentz ids;
      //   2. score each id by the size of the combined id-set of its incident factors (its fill-in =
      //      the rank of the intermediate that eliminating it would create) and pick the smallest;
      //   3. partition factors into those touching that id (`group`) and the rest, merge+sum the group
      //      over that id via `eliminate`, push the result back, repeat;
      //   4. once nothing is shared, the leftovers are scalars — multiply their single entries.
      for (;;) {
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
        // choose the id with the smallest incident-factor union (min fill-in).
        int bestId = -1;
        std::size_t bestUnionSize = SIZE_MAX;
        for (int cand : ids) {
          std::vector<int> unionIds;
          for (const Factor &F : facs) {
            bool incident = false;
            for (int id : F.ids)
              if (id == cand) {
                incident = true;
                break;
              }
            if (!incident) continue;
            for (int id : F.ids) {
              bool seen = false;
              for (int u : unionIds)
                if (u == id) {
                  seen = true;
                  break;
                }
              if (!seen) unionIds.push_back(id);
            }
          }
          if (unionIds.size() < bestUnionSize) {
            bestUnionSize = unionIds.size();
            bestId = cand;
          }
        }
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
      // remaining factors are scalars (no ids): multiply their single entries.
      MPoly prod = MPoly::constant(nsym, Cx{1, 0});
      for (const Factor &F : facs)
        prod = prod * (F.v.empty() ? MPoly(nsym) : F.v[0]);
      return prod;
    }

    /// @brief Reduce a Dirac-trace @ref Factor to a scalar @ref MPoly when there is NO surrounding
    ///        Lorentz net. Three cases: (a) no free legs → the scalar entry; (b) the trace is
    ///        structurally zero (an odd γ-count chain — `numeric_dirac` zeroed every entry) → 0; (c)
    ///        every free-leg id is PAIRED (a γ^μ…γ^μ self-contraction, e.g. a dressed numerator whose
    ///        IDENT/SLASH choice left two free legs on the same id) → sum those ids via
    ///        @ref contract_factors. An UNPAIRED free leg with no Lorentz partner is a genuinely
    ///        uncontracted Lorentz index (an ill-formed/mis-split diagram) and still throws.
    ///        This is what lets the dressed-collection combination loop produce odd / self-paired
    ///        chains (which the front-end odd-trace filter never sees) without a spurious abort.
    inline MPoly close_free_legs(int nsym, const Factor &T, const std::vector<MPoly> &atomDen = {},
                                 const std::vector<std::vector<int>> &units = {})
    {
      if (T.ids.empty()) return T.v.empty() ? MPoly::constant(nsym, Cx{1, 0}) : T.v[0];
      bool allZero = true;
      for (const MPoly &v : T.v)
        if (!v.empty()) {
          allZero = false;
          break;
        }
      if (allZero) return MPoly(nsym); // odd γ-count: trace ≡ 0
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
      throw std::runtime_error("numeric_value: Dirac chain has an UNPAIRED free leg but no Lorentz net "
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
      if (fs.empty()) return MPoly::constant(nsym, Cx{1, 0});
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
    inline void fuse_projectors(std::vector<NElem> &e, Cx &coeff)
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

  /// @brief Contract a diagram (Dirac chain ⊗ Lorentz network) to its scalar trace polynomial.
  /// @param nsym     number of user symbols (MPoly variable count)
  /// @param dirac    the closed Dirac chain (may be empty for a pure-gauge diagram)
  /// @param lorentz  the pure-Lorentz network (metrics / vectors / projectors / Levi-Civita)
  /// @param comp     component table `comp[vid]` = 4 MPoly components of fundamental momentum `vid`
  /// @param atomDen  projector denominators `atomDen[aid] = k²` (for monomial cancellation)
  inline MPoly numeric_value(int nsym, const network::DiracNet &dirac, const NNet &lorentz,
                             const std::vector<std::array<MPoly, 4>> &comp, const std::vector<MPoly> &atomDen)
  {
    // a component may hold SEVERAL independent spinor loops (e.g. a quark loop + the projection-closed
    // external line, tied only by gluon propagators) — trace each into its own Lorentz tensor; the
    // gluon legs they share contract via the Lorentz net below. tr(1)=4 rides each loop's trace.
    std::vector<ndetail::Factor> loops = ndetail::dirac_loop_factors(nsym, dirac, comp);
    MPoly result(nsym);
    for (const NTerm &nt : lorentz) {
      // fold same-momentum projector chains (P·P→P, orthogonal→0) before expansion; no-op unless present.
      std::vector<NElem> elems = nt.e;
      Cx co = nt.coeff;
      ndetail::fuse_projectors(elems, co);
      if (co.re == 0 && co.im == 0) continue;
      std::vector<ndetail::Factor> facs = loops;
      facs.reserve(loops.size() + elems.size());
      for (const NElem &el : elems)
        facs.push_back(ndetail::elem_factor(nsym, el, comp));
      // contract_factors consumes its `facs` argument by value — move so the per-term
      // factor list is built once, not copied again into the call.
      MPoly term = ndetail::contract_factors(nsym, std::move(facs));
      term = term * MPoly::constant(nsym, co);
      result = result + term;
    }
    if (lorentz.empty())
      result = ndetail::close_loops(nsym, loops, atomDen);
    return divThroughMonomialAtoms(result, atomDen);
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

  /// @brief Contract a diagram given the Lorentz part as an inv-backend @ref network::NetVal (so the
  ///        generator reuses the existing net-string emission: `proj`/`met`/`vec`/`epsilon` builders).
  ///        Equivalent to @ref numeric_value but reading `network::Elem` instead of @ref NElem.
  inline MPoly numeric_value_netval(int nsym, const network::DiracNet &dirac, const network::NetVal &lor,
                                    const std::vector<std::array<MPoly, 4>> &comp, const std::vector<MPoly> &atomDen,
                                    const std::vector<std::vector<int>> &units = {})
  {
    // trace each independent spinor loop of the component (see numeric_value); the shared gluon legs
    // contract via the Lorentz net below.
    std::vector<ndetail::Factor> loops = ndetail::dirac_loop_factors(nsym, dirac, comp);
    // pre-reduce the atom denominators (idempotent if the caller already did) so monomial-cancellation
    // detection works during the per-step intermediate reduction inside contract_factors.
    std::vector<MPoly> aden = atomDen;
    for (MPoly &a : aden)
      a = reduce_units(a, units);
    MPoly result(nsym);
    for (const network::PTerm &pt : lor) {
      // build the numeric element list, then fold same-momentum projector chains before expansion.
      std::vector<NElem> elems;
      elems.reserve(pt.e.size());
      for (const network::Elem &el : pt.e)
        elems.push_back(elem_to_nelem(el));
      Cx co = pt.coeff;
      ndetail::fuse_projectors(elems, co);
      if (co.re == 0 && co.im == 0) continue;
      std::vector<ndetail::Factor> facs = loops;
      facs.reserve(loops.size() + elems.size());
      for (const NElem &el : elems)
        facs.push_back(ndetail::elem_factor(nsym, el, comp));
      // move the per-term factor list into contract_factors (consumed by value) — avoids
      // a redundant deep copy of every Factor's MPoly entries.
      MPoly term = ndetail::contract_factors(nsym, std::move(facs), aden, units);
      term = term * MPoly::constant(nsym, co);
      result = result + term;
    }
    if (lor.empty())
      result = ndetail::close_loops(nsym, loops, aden, units);
    // sin^2 -> 1 - cos^2 BEFORE cancellation: collapses bare-loop k²=l1²(cos²+sin²) to the monomial l1²
    // (so its atom cancels) and shrinks the polynomial to the FORM angular basis.
    result = reduce_units(result, units);
    return divThroughMonomialAtoms(result, atomDen);
  }

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

  /// @brief One structure option of a dressed numerator slot: a numeric coefficient × a product of
  ///        dressing atoms (`dress`) × a Dirac structure that is either the identity `δ` or a slash
  ///        `γ·p` (`slash` with momentum `vlc`). Higher structures (σ^{μν}, …) are a later extension.
  struct DSlotOpt {
    Cx coeff{1, 0};
    std::vector<int> dress; ///< dressing-atom ids (need not be pre-sorted; merged sorted on use)
    bool slash = false;     ///< true: structure is `γ·p` (use `vlc`); false: identity `δ`
    std::vector<std::pair<double, int>> vlc; ///< slash momentum `Σ coeff·comp(vid)` (when `slash`)
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

  namespace ndetail
  {
    /// Enumerate the Cartesian product of slot-structure choices, contract each concrete Dirac chain
    /// via @p contract (an `MPoly`-returning closure that runs the full Dirac+Lorentz contraction), and
    /// collect the results into one @ref DPoly keyed by the dressing monomial. Each slot is referenced
    /// at most once in @p chain (a numerator occupies one chain position).
    template <class ContractFn>
    inline DPoly dress_collect(int nsym, const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                               ContractFn &&contract)
    {
      DPoly out(nsym);
      const int nSlots = static_cast<int>(slots.size());
      // The dressed chain is always a quark line, so it carries ≥1 closed spinor loop; the loop count is
      // (LoopSep markers)+1. A combination that collapses a WHOLE loop to the identity — every slot on it
      // picks the δ option and the loop carries no fixed γ — leaves that loop's segment empty. split_loops
      // drops empty segments and close_loops returns 1 for an empty factor set, so the loop's tr(1)=4
      // would be silently lost (a 4× undercount, e.g. the all-mass channel of a scalar/Yukawa self-energy
      // quark loop). Count the loops up front so we can restore 4^(#collapsed loops) per combination. This
      // path is reached ONLY for dressed quark lines (≥1 loop): pure-gauge diagrams have no slots and never
      // arrive here, and an existing dressed flow whose every loop keeps a fixed γ (Zq/ZAqbq{1,4,7}, …)
      // never collapses ⇒ nEmpty==0 ⇒ byte-identical.
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
        for (const DChainTok &tok : chain) {
          if (!tok.isSlot) {
            concrete.push_back(tok.fac);
            continue;
          }
          const DSlotOpt &opt = slots[tok.slot][choice[tok.slot]];
          combCoeff = combCoeff * opt.coeff;
          dressMono.insert(dressMono.end(), opt.dress.begin(), opt.dress.end());
          if (opt.slash) concrete.push_back(network::dslash(opt.vlc));
          // identity option: contributes no Dirac token (δ between the surrounding γ's)
        }
        if (!(combCoeff.re == 0 && combCoeff.im == 0)) {
          // restore tr(1)=4 for every spinor loop that collapsed to the identity in this combination
          const int nCollapsed = nloops - static_cast<int>(split_loops(concrete).size());
          for (int c = 0; c < nCollapsed; ++c) combCoeff = combCoeff * Cx{4, 0};
          MPoly mp = contract(concrete);
          if (!mp.empty()) {
            mp = mp * MPoly::constant(nsym, combCoeff);
            out.add(dmono_sorted(std::move(dressMono)), mp);
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
      return out;
    }
  } // namespace ndetail

  /// @brief Dressed analogue of @ref numeric_value_netval: contract a diagram whose Dirac chain carries
  ///        dressed-numerator SLOTS, returning the collected @ref DPoly (one term per dressing monomial).
  ///        Reuses @ref numeric_value_netval per structure combination, so the Dirac/Lorentz contraction
  ///        is the exact same validated code path; only the collection over dressing atoms is new.
  inline DPoly numeric_value_dressed_netval(int nsym, const std::vector<DChainTok> &chain,
                                            const std::vector<DSlot> &slots, const network::NetVal &lor,
                                            const std::vector<std::array<MPoly, 4>> &comp,
                                            const std::vector<MPoly> &atomDen,
                                            const std::vector<std::vector<int>> &units = {})
  {
    return ndetail::dress_collect(nsym, chain, slots, [&](const network::DiracNet &d) {
      return numeric_value_netval(nsym, d, lor, comp, atomDen, units);
    });
  }

  /// @brief Dressed analogue of @ref numeric_value (reading the numeric @ref NNet Lorentz network).
  inline DPoly numeric_value_dressed(int nsym, const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                                     const NNet &lorentz, const std::vector<std::array<MPoly, 4>> &comp,
                                     const std::vector<MPoly> &atomDen)
  {
    return ndetail::dress_collect(nsym, chain, slots,
                                  [&](const network::DiracNet &d) { return numeric_value(nsym, d, lorentz, comp, atomDen); });
  }

  /// @brief Build the projector inverse-atom denominators by scanning the Lorentz nets for every
  ///        projector kind. Each carries its loop momentum's fundamental id (`vid`); a transverse /
  ///        longitudinal / electric projector fills its full atom `atomDen[inv] = k² = Σ_μ comp[μ]²`,
  ///        and an electric / magnetic projector fills its spatial atom
  ///        `atomDen[invS] = |k⃗|² = Σ_{μ=1..3} comp[μ]²` (component 0 = temporal). The result is sized
  ///        to hold every `inv`/`invS` id seen (others are unused all-zero MPolys).
  inline std::vector<MPoly> collect_atom_denoms(int nsym, const std::vector<network::NetVal> &lors,
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
    std::vector<MPoly> atomDen(maxId + 1, MPoly(nsym));
    for (const network::NetVal &nv : lors)
      for (const network::PTerm &pt : nv)
        for (const network::Elem &e : pt.e)
          if (isProj(e)) {
            const auto &cv = comp[e.vid];
            if (e.inv >= 0) { // full k² = Σ_μ comp[μ]²
              MPoly k2(nsym);
              for (int mu = 0; mu < 4; ++mu)
                k2 = k2 + cv[mu] * cv[mu];
              atomDen[e.inv] = std::move(k2);
            }
            if (e.invS >= 0) { // spatial |k⃗|² = Σ_{μ=1..3} comp[μ]² (component 0 = temporal)
              MPoly ks2(nsym);
              for (int mu = 1; mu < 4; ++mu)
                ks2 = ks2 + cv[mu] * cv[mu];
              atomDen[e.invS] = std::move(ks2);
            }
          }
    return atomDen;
  }

  /// @brief Lower a contracted diagram polynomial into the shared env via the existing CSE + Horner
  ///        back-end (`gdetail::best_into`). User symbols intern as env kind 3 (`var`), surviving
  ///        inverse atoms as kind 1 (`inv`). Returns an @ref network::GenProg (real, or complex when the
  ///        polynomial carries an imaginary coefficient).
  /// @param realOnly when true, emit a REAL program even if the polynomial has imaginary coefficients
  ///        — the caller has proven only `Re(this trace)` is consumed (its assembly coefficient is
  ///        real). `Re(Σ c·mono) = Σ Re(c)·mono` for real monomials, so the imaginary half is dead;
  ///        skipping it avoids computing+returning a `std::complex` whose `.imag()` nobody reads.
  inline network::GenProg to_genprog(const MPoly &p, network::GlobalEnv &g, bool realOnly = false)
  {
    // Prune numerically-zero monomials. A NUMERIC frame fixes the external momenta to concrete
    // components, so exact (analytic) cancellations in the trace surface as tiny residual coefficients
    // (~1e-12 … 1e-30 against real coefficients of O(10²)) — pure round-off, not physics. They don't
    // affect the value (the kernel matches FORM to ~1e-13) but each spurious monomial costs runtime
    // arithmetic; on the dense 1/4/7 trace ~half the monomials are such noise. Drop |c| below a RELATIVE
    // tolerance vs the largest coefficient (a clean ~10-order gap separates noise from real terms).
    double maxabs = 0.0;
    for (const auto &[m, c] : p.t)
      maxabs = std::max(maxabs, std::max(std::fabs(c.re), std::fabs(c.im)));
    const double tol = 1e-9 * maxabs;
    auto keep = [&](const Cx &c) { return std::max(std::fabs(c.re), std::fabs(c.im)) >= tol; };
    bool cplx = false;
    if (!realOnly)
      for (const auto &[m, c] : p.t)
        if (keep(c) && (c.im > tol || c.im < -tol)) {
          cplx = true;
          break;
        }
    auto build = [&](bool imag) {
      std::vector<network::LMono> monos;
      monos.reserve(p.t.size());
      for (const auto &[m, c] : p.t) {
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
        network::LMono lm;
        lm.c = imag ? c.im : c.re;
        lm.vp = std::move(vp);
        monos.push_back(std::move(lm));
      }
      return monos;
    };
    expr::rdetail::RBuilder w;
    const int reRoot = network::gdetail::best_into(build(false), w);
    if (!cplx) return network::GenProg{std::move(w.ins), reRoot};
    const int imRoot = network::gdetail::best_into(build(true), w);
    network::GenProg gp{std::move(w.ins), reRoot};
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
  inline network::GenProg to_genprog(const DPoly &p, network::GlobalEnv &g, bool realOnly = false)
  {
    // PER-DRESSING-CHANNEL noise prune (NOT one global tolerance). A DPoly is
    // Σ_d (dressing-product_d)·(kinematic-poly_d); each channel `d` is reweighted at RUNTIME by its
    // dressing product `dr_d` (regulators/Z-factors), which can swing by many orders across the loop
    // domain. A single global `1e-9·max|c|` would set the tolerance from the largest-coefficient channel
    // and DELETE a small-coefficient channel's GENUINE terms at generation time — silently fine while
    // that channel's runtime weight is small, but WRONG the moment a dressing suppresses the large
    // channel (dr→~0) where the small one should dominate. So prune each channel against ITS OWN max,
    // exactly as the @ref MPoly overload does within one trace — which is also what the distributed
    // (collection-off) path does, since there each channel is a separate MPoly trace pruned on its own
    // scale. The "clean ~10-order gap between round-off and physics" holds WITHIN a channel, not across
    // channels at different scales. (`dmonoTol(mp) ≤ global tol`, so this only KEEPS more, never wrongly
    // drops — strictly safer, and identical when all channels share a scale.)
    auto dmonoTol = [](const MPoly &mp) {
      double maxabs = 0.0;
      for (const auto &[m, c] : mp.t)
        maxabs = std::max(maxabs, std::max(std::fabs(c.re), std::fabs(c.im)));
      return 1e-9 * maxabs;
    };
    auto keepAt = [](const Cx &c, double tol) { return std::max(std::fabs(c.re), std::fabs(c.im)) >= tol; };
    bool cplx = false;
    if (!realOnly)
      for (const auto &[d, mp] : p.t) {
        const double tol = dmonoTol(mp);
        for (const auto &[m, c] : mp.t)
          if (keepAt(c, tol) && (c.im > tol || c.im < -tol)) {
            cplx = true;
            break;
          }
        if (cplx) break;
      }
    auto build = [&](bool imag) {
      std::vector<network::LMono> monos;
      for (const auto &[d, mp] : p.t) {
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
        for (const auto &[m, c] : mp.t) {
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
          network::LMono lm;
          lm.c = imag ? c.im : c.re;
          lm.vp = std::move(vp);
          monos.push_back(std::move(lm));
        }
      }
      return monos;
    };
    expr::rdetail::RBuilder w;
    const int reRoot = network::gdetail::best_into(build(false), w);
    if (!cplx) return network::GenProg{std::move(w.ins), reRoot};
    const int imRoot = network::gdetail::best_into(build(true), w);
    network::GenProg gp{std::move(w.ins), reRoot};
    gp.rootIm = imRoot;
    return gp;
  }

} // namespace numtracer::numeric
