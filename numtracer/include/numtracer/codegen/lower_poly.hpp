/// @file lower_poly.hpp
/// @brief Lower a contracted @ref Poly / @ref DressedPoly to the straight-line instruction stream of the
///        emitted kernel.
///
/// Lowering: @ref to_genprog / `lower_into` turn a finished @ref Poly or @ref DressedPoly into straight-line
/// instructions (Horner ordering + common-subexpression elimination, `gdetail::best_into`). Leaves —
/// scalar products, 1/k², dressings — are interned in a @ref network::GlobalEnv shared by all
/// diagrams of a kernel, so each is computed once per kernel call.
///
/// This header is the front half of the lowering: it interns every leaf in the env (raw symbols,
/// 1/k² atoms, dressings, and the two halves of complex leaves), drops round-off-sized coefficients,
/// snaps the rest to @ref kCoeffSnapDigits digits, and splits the polynomial into real and imaginary
/// monomial lists. `gdetail::best_into` (codegen/gen.hpp) then picks a Horner ordering and the Horner
/// core (codegen/lower.hpp) emits the instructions; neither of those knows about @ref Poly.
///
/// Layout: declarations first (always visible), then the bodies, which compile only in the library
/// TU or a header-only build (see core/export.hpp).
#pragma once

#include "numtracer/codegen/gen.hpp"       // network::GlobalEnv / GenProg / LMono / gdetail::best_into
#include "numtracer/core/config.hpp"       // NT_THROW
#include "numtracer/core/envvar.hpp"       // env_int
#include "numtracer/core/export.hpp"       // NUMTRACER_FUNC / NUMTRACER_DEFINE_BODIES
#include "numtracer/poly/dressed_poly.hpp" // DressedPoly / DressingMono
#include "numtracer/poly/poly.hpp"         // Poly
#include "numtracer/poly/tolerances.hpp"   // kNoisePruneRelTol

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <utility>
#include <vector>

namespace numtracer::inline network
{

  // ================================================================================================
  // Declarations: the entry points (always visible)
  // ================================================================================================

  /// @brief Lower a contracted diagram polynomial into the shared env via the CSE + Horner back-end
  ///        (`gdetail::best_into`). User symbols intern as `SymKind::var`, surviving inverse atoms
  ///        as `SymKind::inv`. Returns an @ref network::GenProg (real, or complex when the polynomial
  ///        carries an imaginary coefficient).
  /// @param realOnly when true, emit a REAL program even if the polynomial has imaginary coefficients
  ///        — the caller has proven only `Re(this trace)` is consumed (its assembly coefficient is
  ///        real). `Re(Σ c·mono) = Σ Re(c)·mono` for real monomials, so the imaginary half is dead;
  ///        skipping it avoids computing+returning a `std::complex` whose `.imag()` nobody reads.
  NUMTRACER_FUNC network::GenProg to_genprog(const Poly &p, network::GlobalEnv &g, bool realOnly = false);
  NUMTRACER_FUNC network::GenProg to_genprog(const DressedPoly &p, network::GlobalEnv &g, bool realOnly = false);

  /// @brief Lower ONE polynomial into @p builder; returns `{reRoot, imRoot}` with
  ///        `imRoot == network::kRealProgram` when the polynomial is real. @ref to_genprog is this plus
  ///        a fresh builder.
  NUMTRACER_FUNC std::pair<int, int> lower_into(const Poly &p, network::GlobalEnv &g,
                                                network::rdetail::RBuilder &builder, bool realOnly);
  /// @brief @ref DressedPoly counterpart of the @ref Poly `lower_into`. Each kinematic monomial is emitted
  ///        exactly as in the @ref Poly overload and additionally carries its dressing monomial's
  ///        atoms as `SymKind::dress` leaves (@ref network::GlobalEnv::dr_id). A `DressedPoly` with a single
  ///        empty dressing monomial reduces to exactly the @ref Poly path.
  NUMTRACER_FUNC std::pair<int, int> lower_into(const DressedPoly &p, network::GlobalEnv &g,
                                                network::rdetail::RBuilder &builder, bool realOnly);

  // ================================================================================================
  // Coefficient conditioning and diagnostics (always visible)
  // ================================================================================================

  /// @brief `NT_GEN_POLYSTATS` verbosity: 0 off, 1 the summary counts, 2 the per-monomial key dump.
  ///
  /// Level 2 is a DIFFERENT report, not a superset of level 1: the dump is meant to be piped through
  /// `sort -u` to count distinct monomials across traces, and the summary lines would corrupt that.
  inline int polystats_level()
  {
    static const int lvl = static_cast<int>(env_int("NT_GEN_POLYSTATS", 0));
    return lvl;
  }

  /// @brief Significant decimal digits kept when snapping a lowered coefficient (@ref snap_coeff).
  ///
  /// The numeric backend contracts components in `double`, so one correctly-rounded frame literal
  /// (`Sqrt[3]/2 -> 0.86602540378443865`, from the 120-degree symmetric point) is sprayed by ~1e4
  /// double operations into thousands of near-duplicates of the SAME mathematical constant:
  /// `16.000000000000007` and `16.000000000000004` both mean 16, reached by different association
  /// orders. That matters because the SSA CSE keys constants on RAW BITS
  /// (`codegen/real_cse.hpp` ihash/ieq), so each variant becomes its own RCONST, every downstream
  /// RMUL that consumes it also differs, and one dirty ulp at a leaf duplicates an entire subtree.
  /// The bit-exact folds (`== 1.0` in `rmul`, `== 0.0` in `rconst`) miss for the same reason,
  /// emitting free multiplies.
  ///
  /// The trace sums cancel heavily (a coefficient perturbation is amplified ~1e6), so do NOT pick
  /// this from the literal-collapse count alone: 14 captures the speed win at no accuracy cost, 12
  /// buys a further 2% for a 1000x accuracy loss. Re-measure the table in
  /// docs/NUMTRACER_DESIGN_NOTES.md before changing it.
  inline constexpr int kCoeffSnapDigits = 14;

  /// @brief Round @p v to @ref kCoeffSnapDigits significant decimal digits.
  ///
  /// Applied at the single point where a polynomial coefficient becomes a codegen `LMono::c` — i.e.
  /// AFTER all polynomial arithmetic and after the noise prune (`ldetail::prune_tol`) — so it cannot
  /// perturb a cancellation, only canonicalise what survived. Perturbation is ~1e-12 relative, three
  /// orders below the numeric-vs-FORM correctness gate.
  ///
  /// It is also the chokepoint where a NON-FINITE coefficient is caught. Nothing downstream stops one:
  /// `RBuilder::ieq` (`codegen/real_cse.hpp`) compares constants with `==`, so a NaN never matches
  /// itself and `find_or_add` appends a fresh SSA slot on every call, and the printer would write a
  /// literal `nan` into a committed kernel. A non-finite coefficient is always an upstream bug (a
  /// division by an identically-zero denominator, an uninitialised component), so it aborts here,
  /// where the trace is still identifiable.
  inline double snap_coeff(double v)
  {
    if (!std::isfinite(v))
      NT_THROW(std::runtime_error, "snap_coeff: non-finite (NaN/Inf) polynomial coefficient reached the lowering — it "
                                   "would be emitted as a literal `nan` and would defeat the SSA constant CSE");
    if (v == 0.0) return v;
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.*e", kCoeffSnapDigits - 1, v);
    return std::strtod(buf, nullptr);
  }

#if NUMTRACER_DEFINE_BODIES

  // ================================================================================================
  // 1. Monomial extraction: env interning, noise pruning, real/imaginary split
  // ================================================================================================

  namespace ldetail
  {
    /// @brief One dressing channel of a lowering: the dressing atoms @ref dress shared by every
    ///        monomial of @ref poly. A plain @ref Poly is a single channel with no dressing atoms.
    struct LowerChannel {
      const DressingMono *dress;
      const Poly *poly;
    };

    /// @brief Noise-prune threshold of one channel: @ref kNoisePruneRelTol times its largest
    ///        coefficient component.
    ///
    /// A NUMERIC frame fixes the external momenta to concrete components, so exact (analytic)
    /// cancellations in a trace surface as tiny residual coefficients (~1e-12 … 1e-30 against real
    /// coefficients of O(10²)): round-off, not physics. They do not affect the value, but each costs
    /// runtime arithmetic, and on the dense 1/4/7 trace about half the monomials are such noise.
    inline double prune_tol(const Poly &mp)
    {
      double maxabs = 0.0;
      for (const auto &[m, c] : mp.terms)
        maxabs = std::max(maxabs, std::max(std::fabs(c.re), std::fabs(c.im)));
      return kNoisePruneRelTol * maxabs;
    }
    inline bool keep_coeff(const Cx &c, double tol) { return std::max(std::fabs(c.re), std::fabs(c.im)) >= tol; }

    /// @brief Append `(idOf(a), multiplicity)` for each run of equal ids in the sorted atom list @p atoms.
    template <class Atoms, class IdOf>
    void push_atom_runs(std::vector<std::pair<int, int>> &vp, const Atoms &atoms, IdOf idOf)
    {
      for (int i = 0; i < (int)atoms.size();) {
        int j = i;
        while (j < (int)atoms.size() && atoms[j] == atoms[i])
          ++j;
        vp.push_back({idOf(atoms[i]), j - i});
        i = j;
      }
    }

    /// @brief Append the kept monomials of one channel to the real (and, when @p cplx, imaginary)
    ///        lowering lists.
    ///
    /// Env ids intern in first-seen order, which fixes the emitted kernel: the channel's dressing
    /// atoms (`SymKind::dress`) first, then per monomial its user symbols (`SymKind::var`) and its
    /// surviving inverse atoms (`SymKind::inv`). ONE walk builds both halves of a complex trace; the
    /// im list holds an entry for EVERY kept monomial (`snap(c.im)`, possibly 0.0).
    inline void append_monos(const LowerChannel &ch, network::GlobalEnv &g, bool cplx,
                             std::vector<network::LMono> &monosRe, std::vector<network::LMono> &monosIm)
    {
      const Poly &mp = *ch.poly;
      const double tol = prune_tol(mp);
      std::vector<std::pair<int, int>> drvp;
      push_atom_runs(drvp, *ch.dress, [&](int a) { return g.dr_id(a); });
      for (const auto &[m, c] : mp.terms) {
        if (!keep_coeff(c, tol)) continue;
        std::vector<std::pair<int, int>> vp = drvp;
        for (int k = 0; k < mp.nsym; ++k)
          if (m.e[k] > 0) vp.push_back({g.var_id(k), m.e[k]});
        push_atom_runs(vp, m.atoms, [&](int a) { return g.inv_id(a); });
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

    /// @brief A complex env leaf in one monomial: the env ids of its real and imaginary halves and its
    ///        power (see `network::GlobalEnv::cplxAtom`).
    struct CxLeaf {
      int re, im, n;
    };

    /// @brief Push `(idOf(a), multiplicity)` runs of @p atoms like @ref push_atom_runs, except that a
    ///        run of a complex leaf (@p isCplx) goes to @p cx as its (Re, Im) slot pair instead.
    template <class Atoms, class IsCplx, class IdOf, class ImIdOf>
    void push_atom_runs_cplx(std::vector<std::pair<int, int>> &vp, std::vector<CxLeaf> &cx, const Atoms &atoms,
                             IsCplx isCplx, IdOf idOf, ImIdOf imIdOf)
    {
      for (int i = 0; i < (int)atoms.size();) {
        int j = i;
        while (j < (int)atoms.size() && atoms[j] == atoms[i])
          ++j;
        if (isCplx(atoms[i]))
          cx.push_back({idOf(atoms[i]), imIdOf(atoms[i]), j - i});
        else
          vp.push_back({idOf(atoms[i]), j - i});
        i = j;
      }
    }

    /// @brief @ref append_monos for an env with complex leaves.
    ///
    /// A monomial carrying complex leaves `(re_l + i im_l)^{n_l}` is expanded binomially into real
    /// monomials, each with the complex coefficient `c · C(n,j) i^j`, so every `i` moves into the
    /// coefficients and the trace bodies stay real arithmetic. A monomial with one complex leaf of power
    /// n becomes n+1 monomials. Pieces whose half is an exact zero are dropped from that half.
    inline void append_monos_cplx(const LowerChannel &ch, network::GlobalEnv &g, bool cplx,
                                  std::vector<network::LMono> &monosRe, std::vector<network::LMono> &monosIm)
    {
      const Poly &mp = *ch.poly;
      const double tol = prune_tol(mp);
      std::vector<std::pair<int, int>> drvp;
      std::vector<CxLeaf> drcx;
      push_atom_runs_cplx(
          drvp, drcx, *ch.dress, [&](int a) { return g.dress_cplx(a); }, [&](int a) { return g.dr_id(a); },
          [&](int a) { return g.dr_im_id(a); });
      std::vector<CxLeaf> cx;
      std::vector<std::pair<Cx, std::vector<std::pair<int, int>>>> pieces, next;
      for (const auto &[m, c] : mp.terms) {
        if (!keep_coeff(c, tol)) continue;
        std::vector<std::pair<int, int>> vp = drvp;
        cx = drcx;
        for (int k = 0; k < mp.nsym; ++k)
          if (m.e[k] > 0) vp.push_back({g.var_id(k), m.e[k]});
        push_atom_runs_cplx(
            vp, cx, m.atoms, [&](int a) { return g.atom_cplx(a); }, [&](int a) { return g.inv_id(a); },
            [&](int a) { return g.inv_im_id(a); });
        pieces.assign(1, {c, std::move(vp)});
        for (const CxLeaf &l : cx) {
          next.clear();
          for (const auto &[co, pvp] : pieces) {
            // (re + i im)^n = Σ_j C(n,j) i^j re^(n-j) im^j; i^j cycles 1, i, -1, -i
            double binom = 1.0;
            for (int j = 0; j <= l.n; ++j) {
              static constexpr Cx ipow[4] = {{1.0, 0.0}, {0.0, 1.0}, {-1.0, 0.0}, {0.0, -1.0}};
              const Cx &ij = ipow[j & 3];
              std::vector<std::pair<int, int>> v = pvp;
              if (l.n - j > 0) v.push_back({l.re, l.n - j});
              if (j > 0) v.push_back({l.im, j});
              next.push_back(
                  {Cx{binom * (co.re * ij.re - co.im * ij.im), binom * (co.re * ij.im + co.im * ij.re)}, std::move(v)});
              binom = binom * static_cast<double>(l.n - j) / static_cast<double>(j + 1);
            }
          }
          pieces.swap(next);
        }
        for (auto &[co, pvp] : pieces) {
          std::sort(pvp.begin(), pvp.end());
          const double reC = snap_coeff(co.re), imC = snap_coeff(co.im);
          if (cplx && imC != 0.0) {
            network::LMono lmIm;
            lmIm.c = imC;
            lmIm.vp = pvp; // copy — the re half may move it below
            monosIm.push_back(std::move(lmIm));
          }
          if (reC == 0.0) continue;
          network::LMono lm;
          lm.c = reC;
          lm.vp = std::move(pvp);
          monosRe.push_back(std::move(lm));
        }
      }
    }

    /// @brief Does channel @p ch carry a complex leaf in a kept monomial? Such a monomial makes the
    ///        trace complex even when its coefficient is real.
    inline bool channel_has_cplx_leaf(const LowerChannel &ch, const network::GlobalEnv &g, double tol)
    {
      for (int a : *ch.dress)
        if (g.dress_cplx(a)) return true;
      for (const auto &[m, c] : ch.poly->terms)
        if (keep_coeff(c, tol))
          for (int a : m.atoms)
            if (g.atom_cplx(a)) return true;
      return false;
    }

    /// @brief `NT_GEN_POLYSTATS` report of the monomials handed to the Horner lowering.
    ///
    /// Level 1: the input term count and the kept monomial count. Monomial count is canonical
    /// (independent of any factorisation strategy), so comparing it against the emitted instruction
    /// count says whether a large op count comes from the ALGEBRA or from weak LOWERING.
    /// Level 2: each monomial's key (env id ^ power). Piped through `sort -u` it counts how many
    /// monomials summed over all traces are DISTINCT, i.e. what cross-trace term collection (which
    /// FormTracer gets by summing all diagrams before expanding) could save.
    inline void dump_polystats(std::size_t nTerms, const std::vector<network::LMono> &monos)
    {
      if (polystats_level() == 1) std::fprintf(stderr, "[polystats] terms=%zu kept=%zu\n", nTerms, monos.size());
      if (polystats_level() == 2)
        for (const auto &lm : monos) {
          std::string key;
          for (const auto &[id, pw] : lm.vp)
            key += std::to_string(id) + "^" + std::to_string(pw) + " ";
          std::fprintf(stderr, "[mono] %s\n", key.c_str());
        }
    }

    /// @brief Lower the channels into @p builder as ONE monomial list, so the shared CSE/Horner
    ///        (`gdetail::best_into`) collects the dressing factors across channels.
    ///
    /// The noise prune is PER CHANNEL, not one global tolerance: each channel is reweighted at
    /// runtime by its dressing product, which can swing by many orders across the loop domain, so a
    /// tolerance taken from the largest channel would delete a small channel's genuine terms.
    inline std::pair<int, int> lower_channels(const std::vector<LowerChannel> &chans, network::GlobalEnv &g,
                                              network::rdetail::RBuilder &builder, bool realOnly)
    {
      const bool cplxLeaves = g.has_complex_leaves();
      bool cplx = false;
      if (!realOnly)
        for (const LowerChannel &ch : chans) {
          const double tol = prune_tol(*ch.poly);
          for (const auto &[m, c] : ch.poly->terms)
            if (keep_coeff(c, tol) && (c.im > tol || c.im < -tol)) {
              cplx = true;
              break;
            }
          if (!cplx && cplxLeaves) cplx = channel_has_cplx_leaf(ch, g, tol);
          if (cplx) break;
        }
      std::size_t nTerms = 0;
      for (const LowerChannel &ch : chans)
        nTerms += ch.poly->terms.size();
      std::vector<network::LMono> monosRe, monosIm;
      monosRe.reserve(nTerms);
      if (cplx) monosIm.reserve(nTerms);
      for (const LowerChannel &ch : chans)
        if (cplxLeaves)
          append_monos_cplx(ch, g, cplx, monosRe, monosIm);
        else
          append_monos(ch, g, cplx, monosRe, monosIm);
      dump_polystats(nTerms, monosRe);
      const int reRoot = network::gdetail::best_into(std::move(monosRe), builder);
      if (!cplx) {
        if (polystats_level() == 1) std::fprintf(stderr, "[polystats] ssa instrs=%zu\n", builder.ins.size());
        return {reRoot, network::kRealProgram};
      }
      const int imRoot = network::gdetail::best_into(std::move(monosIm), builder);
      if (polystats_level() == 1) std::fprintf(stderr, "[polystats] ssa instrs=%zu (re+im)\n", builder.ins.size());
      return {reRoot, imRoot};
    }
  } // namespace ldetail

  // ================================================================================================
  // 2. Entry points: lower_into, to_genprog
  // ================================================================================================

  NUMTRACER_FUNC std::pair<int, int> lower_into(const Poly &p, network::GlobalEnv &g,
                                                network::rdetail::RBuilder &builder, bool realOnly)
  {
    static const DressingMono noDress;
    return ldetail::lower_channels({{&noDress, &p}}, g, builder, realOnly);
  }

  NUMTRACER_FUNC network::GenProg to_genprog(const Poly &p, network::GlobalEnv &g, bool realOnly)
  {
    network::rdetail::RBuilder builder;
    const auto [reRoot, imRoot] = lower_into(p, g, builder, realOnly);
    network::GenProg gp{std::move(builder.ins), reRoot};
    gp.rootIm = imRoot;
    return gp;
  }

  NUMTRACER_FUNC std::pair<int, int> lower_into(const DressedPoly &p, network::GlobalEnv &g,
                                                network::rdetail::RBuilder &builder, bool realOnly)
  {
    std::vector<ldetail::LowerChannel> chans;
    chans.reserve(p.terms.size());
    for (const auto &[d, mp] : p.terms)
      chans.push_back({&d, &mp});
    return ldetail::lower_channels(chans, g, builder, realOnly);
  }

  NUMTRACER_FUNC network::GenProg to_genprog(const DressedPoly &p, network::GlobalEnv &g, bool realOnly)
  {
    network::rdetail::RBuilder builder;
    const auto [reRoot, imRoot] = lower_into(p, g, builder, realOnly);
    network::GenProg gp{std::move(builder.ins), reRoot};
    gp.rootIm = imRoot;
    return gp;
  }

#endif // NUMTRACER_DEFINE_BODIES

} // namespace numtracer::inline network
