/// @file trace_table.hpp
/// @brief Trace tables: a trace rewritten as a short polynomial in its momentum scales, whose coefficients
///        depend only on the dimensionless kinematics and are tabulated once per run.
///
/// A trace T is a polynomial in the frame symbols and the surviving inverse atoms `inv_j = 1/Q_j`. It is
/// homogeneous in the momentum scales: the radial loop variable `l` and (optionally) one dimensionful
/// external coordinate `S`. Clearing the atoms gives
///
///     T = N · ∏_j inv_j^{m_j},    N = Σ_{n=n0..n1} a_n · l^n · S^{D-n},
///
/// where `m_j` is the highest power of atom j in T and every `a_n` is a polynomial in the remaining
/// (dimensionless) symbols only — angles and shape variables. Those are the table axes: at run time the
/// kernel reads the `a_n` of its node and evaluates, with `s = max(l, S)` and `r = min(l, S)/s ≤ 1`,
///
///     T = s^d · ∏_j (s^{w_j} inv_j)^{m_j} · r^{n0 or D-n1} · Horner_r(a),     d = D − Σ_j w_j m_j,
///
/// (`w_j` = mass dimension of `Q_j`, 2 or 0). Every factor but `s^d` is bounded and Horner in `r` has no
/// cancellation; the kernel multiplies `s` in last, one factor at a time (`tab_scale`), so a piece overflows
/// only if its value does — `s^d` alone exceeds the float range already at s ~ 1e3, d = 13.
///
/// Everything here is exact polynomial algebra on the generator's own `Poly`: no fitting, and terms that
/// cancel identically vanish (up to the lowering's round-off pruning), so the stored range n0..n1 is
/// the true one. That matters: a coefficient that should vanish but survives as round-off is amplified
/// by `(l/S)^n` at the extremes.
///
/// A dressed trace Σ_α D^α P_α (@ref DressedPoly: dressing monomials D^α, runtime values) and a complex one
/// (Re + i Im) are tabulated piece by piece: each real kinematic polynomial (the real or imaginary part of
/// one P_α) gets its own table form, and the reconstruction multiplies it by D^α at run time.
///
/// Layout: declarations first, bodies in the library TU (see core/export.hpp).
#pragma once

#include "numtracer/codegen/gen.hpp"        // GlobalEnv, FillFormulas, RInstr, edetail emission helpers
#include "numtracer/codegen/lower_poly.hpp" // lower_into
#include "numtracer/codegen/precision.hpp"  // emit_precision
#include "numtracer/core/export.hpp"
#include "numtracer/poly/dressed_poly.hpp"
#include "numtracer/poly/poly.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <map>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace numtracer::inline network
{
  /// @brief The frame symbols that carry mass dimension: the radial loop variable, and optionally one
  ///        external scale (`-1`: the trace has no external scale, e.g. a vacuum diagram at zero momentum).
  struct TableScales {
    int radial = -1;
    int external = -1;
  };

  /// @brief A trace in table form (see the file comment).
  struct TableTrace {
    int d = 0;   ///< homogeneity degree of the trace in the scales
    int D = 0;   ///< homogeneity degree of the cleared numerator N
    int n0 = 0;  ///< lowest live power of the radial variable in N
    int n1 = -1; ///< highest live power; n1 < n0 means T vanishes identically
    struct Inv {
      int atom; ///< inverse-atom id (index into the frame's denominators)
      int m;    ///< its highest power in T
      int w;    ///< mass dimension of its denominator (2, or 0 for a dimensionless one)
    };
    std::vector<Inv> inv;
    std::vector<Poly> coeff; ///< a_{n0} … a_{n1}, polynomials in the non-scale symbols

    int entries() const { return n1 >= n0 ? n1 - n0 + 1 : 0; }
  };

  /// @brief The table form of @p trace, with @p den the frame's denominators (`den[atom] = Q`).
  ///
  /// Returns nullopt when the trace cannot be tabulated: a complex coefficient, or a complex denominator.
  /// Throws when T or a denominator is not homogeneous in @p sc — then the scales are wrong, which is a
  /// caller error, not a property of the trace.
  NUMTRACER_FUNC std::optional<TableTrace> table_form(const Poly &trace, const std::vector<Poly> &den, TableScales sc);

  /// @brief Host evaluation of a table trace from its entries @p a (a_{n0}…a_{n1}), the scales and the inverse
  ///        atoms' values @p invVal (indexed by atom id). The reference the emitted code must reproduce.
  NUMTRACER_FUNC double eval_table_trace(const TableTrace &t, const double *a, double l, double S,
                                         const std::vector<double> &invVal);

  /// @brief One piece of a tabulated trace: the table form of the real or the imaginary part of the kinematic
  ///        polynomial that multiplies one dressing monomial.
  struct TablePiece {
    DressingMono dress; ///< its dressing atoms (empty for an undressed term)
    bool imag = false;  ///< the imaginary part
    TableTrace t;
  };

  /// @brief A trace in table form: T = Σ_pieces D^dress · (i if imag) · T_piece. The dressings stay runtime
  ///        factors (kernel env slots), so a dressed trace costs one table form per dressing monomial.
  struct TabulatedTrace {
    std::vector<TablePiece> pieces;
    bool complex = false; ///< the trace returns `nt_complex_t`

    int entries() const
    {
      int n = 0;
      for (const auto &p : pieces)
        n += p.t.entries();
      return n;
    }
    /// Instructions the reconstruction still executes per call (an estimate; see @ref select_table_traces).
    int recon_ops() const
    {
      int n = 6;
      for (const auto &p : pieces)
        n += p.t.entries() + 5 + 4 * static_cast<int>(p.t.inv.size()) + static_cast<int>(p.dress.size()) +
             std::abs(p.t.d);
      return n;
    }
  };

  /// @brief The table form of a (dressed) trace, piece by piece. @p complex: the trace's lowered program is
  ///        complex; otherwise only its real part is consumed and the imaginary parts are dropped. The real /
  ///        imaginary split uses the lowering's noise threshold, per dressing monomial.
  ///
  /// Returns nullopt for a complex dressing or denominator; throws like @ref table_form on wrong scales.
  NUMTRACER_FUNC std::optional<TabulatedTrace> tabulate(const DressedPoly &trace, const std::vector<Poly> &den,
                                                        TableScales sc, bool complex, const GlobalEnv &g);
  NUMTRACER_FUNC std::optional<TabulatedTrace> tabulate(const Poly &trace, const std::vector<Poly> &den, TableScales sc,
                                                        bool complex, const GlobalEnv &g);

  /// @brief Host evaluation of a tabulated trace from its entries @p a, the scales, the inverse atoms' and
  ///        the dressings' values (indexed by atom / dressing id).
  NUMTRACER_FUNC Cx eval_tabulated_trace(const TabulatedTrace &t, const double *a, double l, double S,
                                         const std::vector<double> &invVal, const std::vector<double> &dressVal);

  /// @brief The balance: keep the table of trace i only if its inline program (@p progs[i]) costs more than
  ///        the reconstruction plus @p balance instructions per byte of table entries read. `balance <= 0`
  ///        keeps every table. Then interns the kept tables' slots in @p g, which must happen before the
  ///        kernel's fill is printed.
  NUMTRACER_FUNC void select_table_traces(std::vector<std::optional<TabulatedTrace>> &tabs,
                                          const std::vector<GenProg> &progs, double balance, TableScales sc,
                                          GlobalEnv &g);

  /// @brief The table fill of several table forms as ONE program (shared CSE): `roots[e]` is entry e, in the
  ///        order table 0's a_{n0}…a_{n1}, then table 1's, …. `env` holds only frame symbols (no atoms).
  struct TableFillProg {
    GlobalEnv env;
    std::vector<RInstr> ins;
    std::vector<int> roots;
  };
  NUMTRACER_FUNC TableFillProg lower_table_fill(const std::vector<const TableTrace *> &traces);

  /// @brief Largest fill emitted as one function, in SSA instructions. nvcc's time per function grows faster
  ///        than its size, while the entries of a fill share little code, so a larger fill is split into
  ///        several functions at the cost of a few percent more work.
  inline constexpr std::size_t kFillChunkInstr = 10000;

  /// @brief Entries [first, last) of a fill, computed by one emitted function from the instructions in keep.
  struct FillChunk {
    std::size_t first = 0, last = 0;
    std::size_t size = 0;   ///< number of instructions kept
    std::vector<char> keep; ///< keep[i]: instruction i is needed by an entry of the chunk
  };

  /// @brief Split @p p into runs of consecutive entries of at most @p budget instructions each (an entry that
  ///        alone needs more is a chunk of its own). Instructions shared between chunks are repeated.
  NUMTRACER_FUNC std::vector<FillChunk> fill_chunks(const TableFillProg &p, std::size_t budget);

  /// @brief Print `decor void name(double *out, const std::size_t stride, <argSig>)`, writing entry e to
  ///        `out[e * stride]`. Always double precision: the fill runs once per run and its result is what
  ///        the reconstruction's accuracy rests on. @p fm supplies the C++ of each frame symbol in terms of
  ///        @p argSig's parameters (a float kernel needs @ref emit_trace_tables' double `powr` for them).
  ///        With @p chunk, only its entries, each from the same statements as in the whole fill.
  NUMTRACER_FUNC void emit_table_fill(std::ostream &out, const TableFillProg &p, const FillFormulas &fm,
                                      const std::string &name, const std::string &argSig, const std::string &decor,
                                      const FillChunk *chunk = nullptr);

  /// @brief Print the reconstruction `decor T name(const T *row, const std::size_t stride, const T *f)` in the
  ///        emitted precision (`nt_complex_t` for a complex trace): `row[e * stride]` is entry `offset + e` of
  ///        the node, `f` the kernel env @p g.
  NUMTRACER_FUNC void emit_table_trace(std::ostream &out, const TabulatedTrace &t, int offset, const std::string &name,
                                       TableScales sc, GlobalEnv &g, const std::string &decor);

  /// @brief Print the table half of a traces header:
  ///        - `tab_fill(out, stride, <argSig>)`, calling one fp64 fill per tabulated trace (`tab_fillN`, so no
  ///          single device function carries every coefficient);
  ///        - `ttrN(row, stride, f)` for every trace N: its reconstruction, a forwarder to `trN(f)` where it
  ///          stays inline, or to `ttr<canon[N]>` where its body duplicates an earlier trace's;
  ///        - `tab_size`, the entries per table node.
  ///        @p argNames are @p argSig's parameter names.
  NUMTRACER_FUNC void emit_trace_tables(std::ostream &out, const std::vector<std::optional<TabulatedTrace>> &tabs,
                                        const std::vector<int> &canon, const FillFormulas &fm,
                                        const std::string &argSig, const std::vector<std::string> &argNames,
                                        TableScales sc, GlobalEnv &g, const std::string &decor);

#if NUMTRACER_DEFINE_BODIES
  namespace tdetail
  {
    /// Scale degree of a monomial: its exponent of the radial plus the external symbol.
    inline int scale_degree(const Mono &m, TableScales sc)
    {
      return m.e[sc.radial] + (sc.external >= 0 ? m.e[sc.external] : 0);
    }

    /// The common scale degree of all terms of @p p; throws (naming @p what) if they differ.
    inline int homogeneous_degree(const Poly &p, TableScales sc, const char *what)
    {
      int deg = -1;
      for (const auto &[m, c] : p.terms) {
        const int dm = scale_degree(m, sc);
        if (deg < 0)
          deg = dm;
        else if (dm != deg)
          NT_THROW(std::runtime_error, (std::string("table_form: ") + what +
                                        " is not homogeneous in the momentum scales (degrees " +
                                        std::to_string(deg) + " and " + std::to_string(dm) +
                                        "): the scale symbols are wrong")
                                           .c_str());
      }
      return deg < 0 ? 0 : deg;
    }

    inline bool is_real(const Poly &p)
    {
      double maxabs = 0, maxim = 0;
      for (const auto &[m, c] : p.terms) {
        maxabs = std::max({maxabs, std::fabs(c.re), std::fabs(c.im)});
        maxim = std::max(maxim, std::fabs(c.im));
      }
      return maxim <= kNoisePruneRelTol * maxabs;
    }

    /// @p p without the terms below the lowering's noise threshold (exact cancellations leave round-off).
    inline Poly pruned(const Poly &p)
    {
      const double tol = ldetail::prune_tol(p);
      PolyScratch s;
      for (const auto &[m, c] : p.terms)
        if (ldetail::keep_coeff(c, tol)) s.push_back({Mono{m.e, m.atoms}, Cx{c.re, 0}});
      return PolyFactory::from_scratch(p.nsym, std::move(s));
    }
  } // namespace tdetail

  NUMTRACER_FUNC std::optional<TableTrace> table_form(const Poly &trace, const std::vector<Poly> &den, TableScales sc)
  {
    if (sc.radial < 0) NT_THROW(std::invalid_argument, "table_form: no radial scale symbol");
    const int ns = trace.nsym;
    TableTrace t;
    if (trace.terms.empty()) return t; // identically zero: n1 < n0
    if (!tdetail::is_real(trace)) return std::nullopt;

    // m_j: the highest power of each atom; the terms grouped by their atom multiset.
    std::map<int, int> mmax;
    std::map<MonoAtoms, PolyScratch> byAtoms;
    for (const auto &[m, c] : trace.terms) {
      for (std::size_t i = 0; i < m.atoms.size();) {
        std::size_t j = i;
        while (j < m.atoms.size() && m.atoms[j] == m.atoms[i]) ++j;
        int &mx = mmax[m.atoms[i]];
        mx = std::max(mx, static_cast<int>(j - i));
        i = j;
      }
      byAtoms[m.atoms].push_back({Mono{m.e, {}}, Cx{c.re, 0}});
    }

    // Denominators: real, homogeneous, and their mass dimension.
    for (const auto &[aid, m] : mmax) {
      if (aid < 0 || aid >= static_cast<int>(den.size()) || den[static_cast<std::size_t>(aid)].empty())
        NT_THROW(std::runtime_error, ("table_form: atom " + std::to_string(aid) + " has no denominator").c_str());
      const Poly &Q = den[static_cast<std::size_t>(aid)];
      if (!tdetail::is_real(Q)) return std::nullopt;
      t.inv.push_back({aid, m, tdetail::homogeneous_degree(Q, sc, "a denominator")});
    }

    // N = Σ_groups G · ∏_j Q_j^{m_j - p_j}, with the powers of Q cached.
    std::map<std::pair<int, int>, Poly> qpow;
    auto power = [&](int aid, int k) -> const Poly & {
      auto it = qpow.find({aid, k});
      if (it != qpow.end()) return it->second;
      Poly r = PolyFactory::constant(ns, Cx{1, 0});
      for (int i = 0; i < k; ++i)
        r = r * den[static_cast<std::size_t>(aid)];
      return qpow.emplace(std::pair{aid, k}, std::move(r)).first->second;
    };
    Poly N = PolyFactory::zero(ns);
    for (auto &[atoms, scratch] : byAtoms) {
      Poly G = PolyFactory::from_scratch(ns, std::move(scratch));
      for (const auto &iv : t.inv) {
        const int p = static_cast<int>(std::count(atoms.begin(), atoms.end(), static_cast<MonoAtomT>(iv.atom)));
        if (iv.m > p) G = G * power(iv.atom, iv.m - p);
      }
      N = N + G;
    }
    N = tdetail::pruned(N);
    if (N.terms.empty()) return t;

    t.D = tdetail::homogeneous_degree(N, sc, "the cleared trace");
    t.d = t.D;
    for (const auto &iv : t.inv)
      t.d -= iv.w * iv.m;

    // Collect by the power of the radial variable; the coefficients keep only the non-scale symbols.
    std::map<int, PolyScratch> byPower;
    for (const auto &[m, c] : N.terms) {
      MonoExp e;
      for (int k = 0; k < ns; ++k)
        if (k != sc.radial && k != sc.external && m.e[k] != 0) e.set(k, m.e[k]);
      byPower[m.e[sc.radial]].push_back({Mono{std::move(e), {}}, c});
    }
    t.n0 = byPower.begin()->first;
    t.n1 = byPower.rbegin()->first;
    for (int n = t.n0; n <= t.n1; ++n) {
      auto it = byPower.find(n);
      t.coeff.push_back(it == byPower.end() ? PolyFactory::zero(ns) : PolyFactory::from_scratch(ns, std::move(it->second)));
    }
    return t;
  }

  NUMTRACER_FUNC double eval_table_trace(const TableTrace &t, const double *a, double l, double S,
                                         const std::vector<double> &invVal)
  {
    if (t.entries() == 0) return 0.0;
    double s, h;
    if (t.n0 == t.n1) { // a single power: N = a l^n S^{D-n}
      s = 1.0;
      h = a[0] * std::pow(l, t.n0) * (t.D > t.n0 ? std::pow(S, t.D - t.n0) : 1.0);
    } else {
      const bool le = l <= S;
      s = le ? S : l;
      const double r = (le ? l : S) / s;
      h = 0;
      if (le) {
        for (int n = t.n1; n >= t.n0; --n)
          h = h * r + a[n - t.n0];
        h *= std::pow(r, t.n0);
      } else {
        for (int n = t.n0; n <= t.n1; ++n)
          h = h * r + a[n - t.n0];
        h *= std::pow(r, t.D - t.n1);
      }
      h *= std::pow(s, t.d);
    }
    for (const auto &iv : t.inv)
      h *= std::pow(std::pow(s, iv.w) * invVal.at(static_cast<std::size_t>(iv.atom)), iv.m);
    return h;
  }

  NUMTRACER_FUNC std::optional<TabulatedTrace> tabulate(const DressedPoly &trace, const std::vector<Poly> &den,
                                                        TableScales sc, bool complex, const GlobalEnv &g)
  {
    if (sc.radial < 0) NT_THROW(std::invalid_argument, "tabulate: no radial scale symbol");
    TabulatedTrace out;
    out.complex = complex;
    for (const auto &[dm, P] : trace.terms) {
      for (int a : dm)
        if (g.dress_cplx(a)) return std::nullopt;
      // the real and imaginary parts, kept by the same rule as the inline lowering (ldetail::lower_channels)
      const double tol = ldetail::prune_tol(P);
      PolyScratch re, im;
      for (const auto &[m, c] : P.terms) {
        if (!ldetail::keep_coeff(c, tol)) continue;
        if (c.re != 0) re.push_back({Mono{m.e, m.atoms}, Cx{c.re, 0}});
        if (complex && (c.im > tol || c.im < -tol)) im.push_back({Mono{m.e, m.atoms}, Cx{c.im, 0}});
      }
      for (const bool imag : {false, true}) {
        PolyScratch &part = imag ? im : re;
        if (part.empty()) continue;
        auto t = table_form(PolyFactory::from_scratch(P.nsym, std::move(part)), den, sc);
        if (!t) return std::nullopt;
        if (t->entries()) out.pieces.push_back({dm, imag, std::move(*t)});
      }
    }
    return out;
  }

  NUMTRACER_FUNC std::optional<TabulatedTrace> tabulate(const Poly &trace, const std::vector<Poly> &den, TableScales sc,
                                                        bool complex, const GlobalEnv &g)
  {
    return tabulate(DressedPolyFactory::fromPoly(trace), den, sc, complex, g);
  }

  NUMTRACER_FUNC Cx eval_tabulated_trace(const TabulatedTrace &t, const double *a, double l, double S,
                                         const std::vector<double> &invVal, const std::vector<double> &dressVal)
  {
    Cx r{0, 0};
    for (const auto &p : t.pieces) {
      double v = eval_table_trace(p.t, a, l, S, invVal);
      for (int d : p.dress)
        v *= dressVal.at(static_cast<std::size_t>(d));
      (p.imag ? r.im : r.re) += v;
      a += p.t.entries();
    }
    return r;
  }

  NUMTRACER_FUNC void select_table_traces(std::vector<std::optional<TabulatedTrace>> &tabs,
                                          const std::vector<GenProg> &progs, double balance, TableScales sc,
                                          GlobalEnv &g)
  {
    if (tabs.size() != progs.size()) NT_THROW(std::invalid_argument, "select_table_traces: one table per trace");
    const double bytes = codegen::emit_single() ? 4. : 8.;
    for (std::size_t i = 0; i < tabs.size(); ++i)
      if (tabs[i] && balance > 0 &&
          !(static_cast<double>(progs[i].ins.size()) > tabs[i]->recon_ops() + balance * bytes * tabs[i]->entries()))
        tabs[i].reset();
    for (const auto &t : tabs) {
      if (!t) continue;
      g.var_id(sc.radial);
      if (sc.external >= 0) g.var_id(sc.external);
      for (const auto &p : t->pieces) {
        for (const auto &iv : p.t.inv)
          g.inv_id(iv.atom);
        for (int d : p.dress)
          g.dr_id(d);
      }
    }
  }

  NUMTRACER_FUNC TableFillProg lower_table_fill(const std::vector<const TableTrace *> &traces)
  {
    TableFillProg p;
    rdetail::RBuilder b;
    for (const TableTrace *t : traces)
      for (const Poly &a : t->coeff)
        p.roots.push_back(lower_into(a, p.env, b, /*realOnly=*/true).first);
    p.ins = std::move(b.ins);
    for (const auto &[kind, a, bb] : p.env.syms)
      if (kind != SymKind::var)
        NT_THROW(std::runtime_error, "lower_table_fill: a table coefficient reads a non-symbol slot");
    return p;
  }

  namespace tdetail
  {
    /// The fill is fp64 whatever the kernel computes in; restored on every exit.
    struct DoubleEmission {
      codegen::EmitPrecision saved = codegen::emit_precision();
      DoubleEmission() { codegen::emit_precision() = codegen::EmitPrecision::Double; }
      ~DoubleEmission() { codegen::emit_precision() = saved; }
    };

    /// `f[i] * f[i] * f[j]` for the dressing atoms @p dress (sorted), or "" for none.
    inline std::string dressing_product(const DressingMono &dress, GlobalEnv &g)
    {
      std::string s;
      for (int d : dress)
        s += (s.empty() ? "f[" : " * f[") + std::to_string(g.dr_id(d)) + "]";
      return s;
    }

    /// The statements computing `h_` = piece @p t without its dressing and without the factor `s^d`, from
    /// entries `row[(off + e) * st]`: every factor printed here is bounded (see the file comment). With an
    /// external scale (@p scaled), `s = s_ = max(l_, e_)`; without one, `s = l_` and `T = a l^d ∏ (l^w inv)^m`.
    inline void emit_piece(std::ostream &out, const TableTrace &t, int off, const std::string &T, bool scaled,
                           GlobalEnv &g)
    {
      auto entry = [&](int n) { return "row[" + std::to_string(off + n - t.n0) + " * st]"; };
      auto rpow = [](int n) { return n ? "powr<" + std::to_string(n) + ">(r_)" : std::string("1"); };
      const std::string s = scaled ? "s_" : "l_";
      if (!scaled)
        out << "    " << T << " h_ = " << entry(t.n0) << ";\n";
      else if (t.n0 == t.n1)
        out << "    " << T << " h_ = " << entry(t.n0) << " * (le_ ? " << rpow(t.n0) << " : " << rpow(t.D - t.n0)
            << ");\n";
      else {
        out << "    " << T << " h_;\n";
        out << "    if (le_) {\n      h_ = " << entry(t.n1) << ";\n";
        for (int n = t.n1 - 1; n >= t.n0; --n)
          out << "      h_ = " << codegen::emit_fma() << "(h_, r_, " << entry(n) << ");\n";
        if (t.n0) out << "      h_ *= powr<" << t.n0 << ">(r_);\n";
        out << "    } else {\n      h_ = " << entry(t.n0) << ";\n";
        for (int n = t.n0 + 1; n <= t.n1; ++n)
          out << "      h_ = " << codegen::emit_fma() << "(h_, r_, " << entry(n) << ");\n";
        if (t.D - t.n1) out << "      h_ *= powr<" << t.D - t.n1 << ">(r_);\n";
        out << "    }\n";
      }
      for (const auto &iv : t.inv) {
        const std::string inv = "f[" + std::to_string(g.inv_id(iv.atom)) + "]";
        const std::string base = iv.w == 2 ? "(" + s + " * " + s + " * " + inv + ")" : inv;
        out << "    h_ *= " << (iv.m == 1 ? base : "powr<" + std::to_string(iv.m) + ">(" + base + ")") << ";\n";
      }
    }
  } // namespace tdetail

  NUMTRACER_FUNC std::vector<FillChunk> fill_chunks(const TableFillProg &p, std::size_t budget)
  {
    const std::size_t n = p.ins.size();
    std::vector<FillChunk> chunks;
    FillChunk cur{0, 0, 0, std::vector<char>(n, 0)};
    std::vector<unsigned> seen(n, 0); // seen[i] == pass: instruction i already collected in this pass
    unsigned pass = 0;
    std::vector<int> stack, fresh;
    // fresh = the instructions entry e needs that the current chunk does not hold yet
    auto collect = [&](std::size_t e) {
      fresh.clear();
      ++pass;
      if (p.roots[e] >= 0) stack.push_back(p.roots[e]);
      while (!stack.empty()) {
        const int i = stack.back();
        stack.pop_back();
        if (cur.keep[i] || seen[i] == pass) continue;
        seen[i] = pass;
        fresh.push_back(i);
        const RInstr &in = p.ins[i];
        if (in.op == RADD || in.op == RMUL) stack.insert(stack.end(), {in.a, in.b});
        if (in.op == RNEG) stack.push_back(in.a);
      }
    };
    for (std::size_t e = 0; e < p.roots.size(); ++e) {
      collect(e);
      if (cur.last > cur.first && cur.size + fresh.size() > budget) {
        chunks.push_back(std::move(cur));
        cur = FillChunk{e, e, 0, std::vector<char>(n, 0)};
        collect(e);
      }
      for (const int i : fresh)
        cur.keep[i] = 1;
      cur.size += fresh.size();
      cur.last = e + 1;
    }
    chunks.push_back(std::move(cur));
    return chunks;
  }

  NUMTRACER_FUNC void emit_table_fill(std::ostream &out, const TableFillProg &p, const FillFormulas &fm,
                                      const std::string &name, const std::string &argSig, const std::string &decor,
                                      const FillChunk *chunk)
  {
    tdetail::DoubleEmission guard;
    out << kLazyPrefix << decor << " void " << name << "(double *out, const std::size_t stride, " << argSig << ") {\n";
    out << std::setprecision(17);
    if (!p.env.syms.empty()) {
      out << fm.prologue;
      out << "  double f[" << p.env.syms.size() << "];\n";
      for (std::size_t i = 0; i < p.env.syms.size(); ++i)
        out << "  f[" << i << "] = " << fm.var(std::get<1>(p.env.syms[i])) << ";\n";
    }
    // The EmitPlan of the whole fill, also for a chunk: its fma folds depend on use counts, and a chunk's own
    // counts would fold differently.
    const edetail::EmitPlan pl = edetail::make_plan(p.ins, p.roots.data(), p.roots.size());
    for (std::size_t i = 0; i < p.ins.size(); ++i)
      if (!chunk || chunk->keep[i]) edetail::emit_stmt(out, p.ins, i, pl);
    for (std::size_t e = chunk ? chunk->first : 0; e < (chunk ? chunk->last : p.roots.size()); ++e) {
      out << "  out[" << e << " * stride] = ";
      if (p.roots[e] < 0)
        out << "0.";
      else
        edetail::emit_operand(out, p.ins, pl, p.roots[e]);
      out << ";\n";
    }
    out << "}\n";
  }

  NUMTRACER_FUNC void emit_table_trace(std::ostream &out, const TabulatedTrace &t, int offset, const std::string &name,
                                       TableScales sc, GlobalEnv &g, const std::string &decor)
  {
    const std::string T = codegen::emit_real_type();
    const std::string zero = edetail::zero_literal();
    out << kLazyPrefix << decor << " " << (t.complex ? std::string("nt_complex_t") : T) << " " << name
        << "([[maybe_unused]] const " << T << " *row, [[maybe_unused]] const std::size_t st, [[maybe_unused]] const "
        << T << " *f) {\n";
    const bool scaled = sc.external >= 0;
    bool part[2] = {false, false};
    for (const auto &p : t.pieces)
      part[p.imag] = true;
    if (t.pieces.empty()) {
      out << "  return " << (t.complex ? "nt_complex_t{" + zero + ", " + zero + "}" : zero) << ";\n}\n";
      return;
    }
    out << "  const " << T << " l_ = f[" << g.var_id(sc.radial) << "];\n";
    if (scaled) {
      out << "  const " << T << " e_ = f[" << g.var_id(sc.external) << "];\n";
      out << "  const bool le_ = l_ <= e_;\n";
      out << "  const " << T << " s_ = le_ ? e_ : l_;\n";
      out << "  const " << T << " r_ = (le_ ? l_ : e_) / s_;\n";
    }
    const char *acc[2] = {"re_", "im_"};
    for (int k = 0; k < 2; ++k)
      if (part[k]) out << "  " << T << " " << acc[k] << ";\n";
    bool assigned[2] = {false, false};
    int off = offset;
    for (const auto &p : t.pieces) {
      out << "  {\n";
      tdetail::emit_piece(out, p.t, off, T, scaled, g);
      // the dressings, then s^d one factor at a time: the running product moves monotonically to the
      // piece's value, so it overflows only if that value does (powr<d>(s) alone overflows float at s ~ 1e3)
      const std::string D = tdetail::dressing_product(p.dress, g);
      if (!D.empty()) out << "    h_ *= " << D << ";\n";
      if (p.t.d) out << "    h_ = tab_scale<" << p.t.d << ">(h_, " << (scaled ? "s_" : "l_") << ");\n";
      const char *a = acc[p.imag];
      out << "    " << a << (assigned[p.imag] ? " += h_;\n  }\n" : " = h_;\n  }\n");
      assigned[p.imag] = true;
      off += p.t.entries();
    }
    out << "  return ";
    if (t.complex)
      out << "nt_complex_t{" << (part[0] ? "re_" : zero) << ", " << (part[1] ? "im_" : zero) << "}";
    else
      out << "re_";
    out << ";\n}\n";
  }

  NUMTRACER_FUNC void emit_trace_tables(std::ostream &out, const std::vector<std::optional<TabulatedTrace>> &tabs,
                                        const std::vector<int> &canon, const FillFormulas &fm,
                                        const std::string &argSig, const std::vector<std::string> &argNames,
                                        TableScales sc, GlobalEnv &g, const std::string &decor)
  {
    if (canon.size() != tabs.size()) NT_THROW(std::invalid_argument, "emit_trace_tables: one canon entry per trace");
    const std::string T = codegen::emit_real_type();
    // The symbol formulas call the kernel's powr<N>, which a float kernel declares for float only: give
    // it a double overload, or the fill would round every shape factor to float.
    if (codegen::emit_single())
      out << "template<int N> " << decor
          << " double powr(double x){ double r=1.; for(int i=0;i<(N<0?-N:N);++i) r*=x; return N<0?1./r:r; }\n";
    // h * x^N as N successive factors (x or 1/x), never forming x^N on its own
    out << "template<int N> " << decor << " " << T << " tab_scale(" << T << " h, const " << T << " x){ const " << T
        << " y = N < 0 ? " << (codegen::emit_single() ? "1.f" : "1.") << " / x : x; for(int i=0;i<(N<0?-N:N);++i) h*=y; return h; }\n";
    std::string args;
    for (const auto &a : argNames)
      args += ", " + a;

    std::vector<int> off(tabs.size(), -1);
    std::ostringstream calls;
    int nEntries = 0, nTab = 0;
    for (std::size_t i = 0; i < tabs.size(); ++i) {
      if (!tabs[i] || canon[i] != static_cast<int>(i) || tabs[i]->entries() == 0) continue;
      std::vector<const TableTrace *> parts;
      for (const auto &p : tabs[i]->pieces)
        parts.push_back(&p.t);
      const TableFillProg prog = lower_table_fill(parts);
      const std::string nm = "tab_fill" + std::to_string(i);
      const std::vector<FillChunk> chunks = fill_chunks(prog, kFillChunkInstr);
      if (chunks.size() == 1)
        emit_table_fill(out, prog, fm, nm, argSig, edetail::eff_decor(decor, prog.ins.size()));
      else {
        for (std::size_t k = 0; k < chunks.size(); ++k)
          emit_table_fill(out, prog, fm, nm + "_c" + std::to_string(k), argSig,
                          edetail::eff_decor(decor, chunks[k].size), &chunks[k]);
        out << kLazyPrefix << decor << " void " << nm << "(double *out, const std::size_t stride, " << argSig
            << ") {\n";
        for (std::size_t k = 0; k < chunks.size(); ++k)
          out << "  " << nm << "_c" << k << "(out, stride" << args << ");\n";
        out << "}\n";
      }
      calls << "  " << nm << "(out + " << nEntries << " * stride, stride" << args << ");\n";
      off[i] = nEntries;
      nEntries += tabs[i]->entries();
      ++nTab;
    }
    out << kLazyPrefix << decor
        << " void tab_fill([[maybe_unused]] double *out, [[maybe_unused]] const std::size_t stride, " << argSig
        << ") {\n"
        << calls.str() << "}\n";
    for (std::size_t i = 0; i < tabs.size(); ++i) {
      const std::string nm = "ttr" + std::to_string(i);
      if (canon[i] != static_cast<int>(i))
        out << kLazyPrefix << decor << " auto " << nm << "(const " << T << " *row, const std::size_t st, const " << T
            << " *f) { return ttr" << canon[i] << "(row, st, f); }\n";
      else if (off[i] >= 0 || (tabs[i] && tabs[i]->entries() == 0))
        emit_table_trace(out, *tabs[i], std::max(off[i], 0), nm, sc, g, decor);
      else
        out << kLazyPrefix << decor << " auto " << nm << "([[maybe_unused]] const " << T
            << " *row, [[maybe_unused]] const std::size_t st, const " << T << " *f) { return tr" << i << "(f); }\n";
    }
    out << "static inline constexpr int tab_size = " << nEntries << ";\n";
    std::fprintf(stderr, "[tables] %d of %zu traces tabulated, %d entries per table node\n", nTab, tabs.size(),
                 nEntries);
  }
#endif // NUMTRACER_DEFINE_BODIES

} // namespace numtracer::inline network
