/// @file frame.hpp
/// @brief @ref numtracer::Frame — the kinematic frame a network is contracted on: the scalar symbols,
///        the momenta's four components in those symbols, and the projector denominators `k²`.
///
/// Everything a contraction needs besides the network itself lives here, so a contraction is one
/// call, `F.trace(chain, net)`, and an evaluation is `F.eval(poly, F.at(...))`:
///
/// ```cpp
/// nt::Frame F;
/// auto P = F.symbol("p"), L = F.symbol("l");
/// auto [C, S] = F.angle("theta");                // S = sqrt(1 - C^2)
/// auto p = F.momentum(P, 0, 0, 0);
/// auto l = F.momentum(L * C, L * S, 0, 0);
/// auto [mu, nu] = F.indices<2>();
/// nt::Poly T = F.trace({nt::slash(p), nt::gamma(mu), nt::slash(l - p), nt::gamma(nu)}, nt::projT(mu, nu, l));
/// double v = F.eval(T, F.at(1.3, 0.8, 0.6)).re;
/// ```
///
/// A frame is built in two stages. First declare its symbols (`symbol`, `angle`); the first time a
/// polynomial is made from them — a momentum, an arithmetic expression — the symbol list is FROZEN,
/// because every polynomial of one frame must share it. Declaring a symbol after that throws.
///
/// Generated code uses the same object: it constructs the frame with its symbol list up front, fills
/// the component table with @ref Frame::set_component, registers the denominators of all its nets
/// once with @ref Frame::add_denominators, and then contracts in parallel through the `const`
/// members.
#pragma once

#include "numtracer/codegen/gen.hpp"              // GlobalEnv / FillFormulas / emit_fill
#include "numtracer/numeric/numeric_contract.hpp" // Poly/DPoly/Mat4, LorentzNet, DiracChain, ndetail::contract*
#include "numtracer/numeric/numeric_driver.hpp"   // poly_to_cpp

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cmath>
#include <initializer_list>
#include <map>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace numtracer::inline numeric
{

  // Forward declarations of the host-only trace-fold templates (`numeric/trace_fold.hpp`). That header
  // spawns threads and is forbidden in the -fno-exceptions net-builder TUs, so `frame.hpp` must NOT
  // include it. The @ref Frame methods that forward here are member templates whose bodies are only
  // instantiated in the generator's main TU, where `trace_fold.hpp` is already included — so these
  // declarations are enough to compile `frame.hpp` everywhere else.
  template <class P, class TraceFn>
  std::vector<P> contract_traces(int nsym, long nCache, unsigned W, TraceFn &&trace);
  template <class P, class TraceFn, class ScaleFn, class Sink>
  void fold_groups_streaming(int nsym, const std::vector<std::vector<int>> &traceIdx,
                             const std::vector<std::vector<Cx>> &subScale,
                             const std::vector<std::vector<int>> &groups, const std::vector<P> &traceTable,
                             long nCache, unsigned W, long window, TraceFn &&trace, ScaleFn &&scale,
                             Sink &&sink);
  template <class TraceFn, class ScaleFn, class Sink>
  void fold_groups_streaming_dressed(int nsym, const std::vector<std::vector<int>> &traceIdx,
                                     const std::vector<std::vector<Cx>> &subScale,
                                     const std::vector<std::vector<DMono>> &subDress,
                                     const std::vector<std::vector<int>> &groups,
                                     const std::vector<Poly> &traceTable, long nCache, unsigned W,
                                     long window, TraceFn &&trace, ScaleFn &&scale, Sink &&sink);

  class Frame;

  /// @brief A scalar symbol of a @ref Frame (a momentum magnitude, an angle cosine, …). Converts to
  ///        the polynomial it stands for, so symbols combine with `+`, `-`, `*` and numbers.
  struct Symbol {
    int id;             ///< position in the frame's symbol list
    const Frame *frame; ///< the frame that declared it
    operator Poly() const;
  };

  /// @brief A point in a frame's symbol space, from @ref Frame::at. Holds a value for EVERY symbol,
  ///        derived ones (an angle's sine) included.
  struct Point {
    std::vector<double> x;
    const Frame *frame = nullptr;
  };

  /// @brief One component of @ref Frame::momentum: a number, a @ref Symbol or a polynomial.
  struct Component {
    Component(double v) : value(v) {}
    Component(int v) : value(v) {}
    Component(Symbol s) : poly(s), isPoly(true), owner(s.frame) {}
    Component(Poly p) : poly(std::move(p)), isPoly(true) {}
    double value = 0;
    Poly poly;
    bool isPoly = false;
    const Frame *owner = nullptr; ///< the frame of a Symbol component (checked by Frame::momentum)
  };

  /// @brief The kinematic frame: symbols, momentum components, and projector denominators.
  ///
  /// Not copyable or movable: @ref Symbol and @ref Point refer back to it.
  class Frame
  {
  public:
    /// An empty frame; declare its symbols with @ref symbol / @ref angle.
    Frame() = default;
    /// A frame with a fixed symbol list (generated code). @p units are the unit-vector groups
    /// `Σ_i x_{g_i}² = 1` among the symbols (e.g. `{cos, sin}`), used to simplify `k²`.
    explicit Frame(std::vector<std::string> names, std::vector<std::vector<int>> units = {}, int nMomenta = 0)
        : names_(std::move(names)), derivedFrom_(names_.size(), -1), units_(std::move(units)), frozen_(true)
    {
      comp_.assign(static_cast<std::size_t>(nMomenta), {zero(), zero(), zero(), zero()});
    }
    Frame(const Frame &) = delete;
    Frame &operator=(const Frame &) = delete;

    // ── symbols ──────────────────────────────────────────────────────────────────────────────────
    /// Declare a scalar symbol. Its value is supplied at evaluation time (@ref at). The name must be
    /// a C++ identifier: it becomes a parameter name in @ref emit_fill.
    Symbol symbol(std::string name)
    {
      require_open("symbol");
      const bool ident = !name.empty() && (std::isalpha(static_cast<unsigned char>(name[0])) || name[0] == '_') &&
                         std::all_of(name.begin(), name.end(),
                                     [](char ch) { return std::isalnum(static_cast<unsigned char>(ch)) || ch == '_'; });
      if (!ident) NT_THROW(std::invalid_argument, ("Frame::symbol: \"" + name + "\" is not a C++ identifier").c_str());
      names_.push_back(std::move(name));
      derivedFrom_.push_back(-1);
      return {static_cast<int>(names_.size()) - 1, this};
    }
    /// Declare an angle θ by name; returns the symbols `{cos_θ, sin_θ}`. Only the cosine is given at
    /// evaluation time — `sin_θ = sqrt(1 - cos_θ²)` is derived — and `cos² + sin² = 1` is used to
    /// simplify (`l²(cos² + sin²) → l²`).
    std::pair<Symbol, Symbol> angle(const std::string &name)
    {
      const Symbol c = symbol("cos_" + name);
      const Symbol s = symbol("sin_" + name);
      derivedFrom_[static_cast<std::size_t>(s.id)] = c.id;
      units_.push_back({c.id, s.id});
      return {c, s};
    }
    /// Number of symbols (the polynomials' variable count).
    int nsym() const { return static_cast<int>(names_.size()); }
    const std::vector<std::string> &symbol_names() const { return names_; }
    /// The unit-vector groups among the symbols.
    const std::vector<std::vector<int>> &units() const { return units_; }

    // ── polynomials (these freeze the symbol list) ─────────────────────────────────────────────────
    Poly zero() const { return Poly(frozen_nsym()); }
    Poly constant(Cx c) const { return Poly::constant(frozen_nsym(), c); }
    Poly constant(double c) const { return constant(Cx{c, 0}); }
    /// The polynomial `x_i` of symbol @p i.
    Poly var(int i) const
    {
      if (i < 0 || i >= frozen_nsym())
        NT_THROW(std::out_of_range, ("Frame::var: symbol " + std::to_string(i) + " out of range").c_str());
      return Poly::var(nsym(), i);
    }
    /// The monomial `c · Π x_i^{e_i}`.
    Poly mono(const std::vector<int> &e, Cx c) const
    {
      if (static_cast<int>(e.size()) != frozen_nsym())
        NT_THROW(std::invalid_argument, "Frame::mono: exponent vector length differs from the symbol count");
      return Poly::mono(nsym(), e, c);
    }
    /// The inverse atom `1/D_aid` as a polynomial (tests and hand-built reductions).
    Poly atom(int aid) const { return Poly::atom(frozen_nsym(), aid); }
    DPoly dzero() const { return DPoly(frozen_nsym()); }
    /// The slashed 4×4 matrix of a momentum given by its components.
    Mat4 slashC(const std::array<Poly, 4> &comp) const { return ::numtracer::numeric::slashC(frozen_nsym(), comp); }

    // ── momenta and indices ──────────────────────────────────────────────────────────────────────
    /// Declare a momentum by its four Euclidean components (component 0 is the heat-bath direction
    /// at finite temperature).
    Momentum momentum(Component c0, Component c1, Component c2, Component c3)
    {
      std::array<Poly, 4> row;
      const Component *cs[4] = {&c0, &c1, &c2, &c3};
      for (int mu = 0; mu < 4; ++mu) {
        const Component &c = *cs[mu];
        if (c.isPoly) {
          if (c.owner && c.owner != this) NT_THROW(std::invalid_argument, "Frame::momentum: a symbol of another frame");
          check_owned(c.poly, "Frame::momentum");
          row[mu] = c.poly.empty() ? zero() : c.poly;
        } else
          row[mu] = c.value == 0.0 ? zero() : constant(c.value);
      }
      comp_.push_back(std::move(row));
      return {{{1.0, static_cast<int>(comp_.size()) - 1}}};
    }
    /// Set component @p mu of momentum @p vid (generated code; the table grows as needed).
    void set_component(int vid, int mu, Poly p)
    {
      if (vid < 0 || mu < 0 || mu > 3) NT_THROW(std::out_of_range, "Frame::set_component: index out of range");
      check_owned(p, "Frame::set_component");
      if (static_cast<std::size_t>(vid) >= comp_.size())
        comp_.resize(static_cast<std::size_t>(vid) + 1, {zero(), zero(), zero(), zero()});
      comp_[static_cast<std::size_t>(vid)][static_cast<std::size_t>(mu)] = std::move(p);
    }
    /// The component table: `components()[vid][mu]`.
    const std::vector<std::array<Poly, 4>> &components() const { return comp_; }
    /// @p N fresh, distinct Lorentz index labels: `auto [mu, nu] = F.indices<2>();`.
    template <int N> std::array<LorentzIndex, N> indices()
    {
      std::array<LorentzIndex, N> r = make_indices<N>(std::make_integer_sequence<int, N>{});
      nextIndex_ += N;
      return r;
    }

    // ── projector denominators ───────────────────────────────────────────────────────────────────
    /// Register the `1/k²` denominators of every projector in @p nets, under the atom ids the nets
    /// carry. Generated code calls this once for all its nets before contracting; hand-built nets do
    /// not need it (@ref trace assigns atoms itself).
    void add_denominators(const std::vector<LorentzNet> &nets)
    {
      std::vector<Poly> den = ndetail::collect_atom_denoms(frozen_nsym(), nets, comp_);
      for (Poly &d : den)
        d = reduce_units(std::move(d), units_);
      if (atomDen_.empty()) {
        atomDen_ = std::move(den);
        return;
      }
      for (std::size_t i = 0; i < den.size(); ++i)
        if (!den[i].empty()) set_denominator(static_cast<int>(i), std::move(den[i]));
    }
    /// The projector denominators: `denominators()[atom] = k²`.
    const std::vector<Poly> &denominators() const { return atomDen_; }
    /// Denominator @p atom as a C++ expression in the symbol names.
    std::string denominator_cpp(int atom) const
    {
      return poly_to_cpp(atomDen_.at(static_cast<std::size_t>(atom)), names_);
    }

    // ── contraction ──────────────────────────────────────────────────────────────────────────────
    /// @brief The trace of the Dirac chain @p chain, contracted with the Lorentz network @p lor, as a
    ///        polynomial in the frame's symbols (projector denominators ride as atoms `1/k²`).
    ///
    /// An empty chain means "no Dirac structure" (the result is just the network's value), and an
    /// empty network is the scalar 1. Every Lorentz index must be contracted.
    ///
    /// This overload assigns atom ids to projectors built without one, which modifies the frame.
    /// A network whose projectors all carry registered atoms leaves the frame untouched, so the call
    /// is then safe from several threads.
    Poly trace(const DiracChain &chain, const LorentzNet &lor = {})
    {
      if (!needs_atoms(lor)) return std::as_const(*this).trace(chain, lor);
      return std::as_const(*this).trace(chain, assign_atoms(lor));
    }
    /// The `const` form: every projector's atoms must already be registered.
    Poly trace(const DiracChain &chain, const LorentzNet &lor = {}) const
    {
      return ndetail::contract(frozen_nsym(), chain, lor, comp_, atomDen_, units_);
    }
    /// A pure Lorentz network contracted to a scalar polynomial:
    /// `F.contract(vec(mu, p) * projT(mu, nu, l) * vec(nu, p))`.
    Poly contract(const LorentzNet &lor) { return trace({}, lor); }
    Poly contract(const LorentzNet &lor) const { return trace({}, lor); }
    /// A chain with dressed-numerator SLOTS, collected per dressing monomial (see @ref DSlotOpt).
    /// Like the plain @ref trace, this overload assigns atoms to projectors built without one.
    DPoly trace(const std::vector<DChainTok> &chain, std::vector<DSlot> slots, LorentzNet lor)
    {
      lor = assign_atoms(std::move(lor));
      for (DSlot &slot : slots)
        for (DSlotOpt &opt : slot)
          for (LorentzFactor &f : opt.netFacs) assign_atoms(f);
      return std::as_const(*this).trace(chain, slots, lor);
    }
    DPoly trace(const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots, const LorentzNet &lor) const
    {
      return ndetail::contract_dressed(frozen_nsym(), chain, slots, lor, comp_, atomDen_, units_);
    }
    /// Generated code: the slots' dressing is carried outside, so the result is a plain polynomial.
    Poly trace_structural(const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                          const LorentzNet &lor) const
    {
      return ndetail::contract_structural(frozen_nsym(), chain, slots, lor, comp_, atomDen_, units_);
    }

    // ── evaluation ───────────────────────────────────────────────────────────────────────────────
    /// A point from the values of the INDEPENDENT symbols, in declaration order (an angle's sine is
    /// derived from its cosine, not given).
    Point at(std::initializer_list<double> values) const { return at(std::vector<double>(values)); }
    /// `F.at(1.3, 0.8, 0.6)`: the same, without the braces (`F.at()` for a frame without symbols).
    template <class... V>
      requires(std::is_arithmetic_v<V> && ...)
    Point at(V... values) const
    {
      return at(std::initializer_list<double>{static_cast<double>(values)...});
    }
    /// The same, from values held in a vector (independent symbols, declaration order).
    Point at(const std::vector<double> &values) const
    {
      if (values.size() != independent_count())
        NT_THROW(std::invalid_argument, ("Frame::at: " + std::to_string(values.size()) + " values given, the frame has " +
                                         std::to_string(independent_count()) + " independent symbols (" +
                                         independent_list() + ")")
                                            .c_str());
      Point pt{std::vector<double>(names_.size(), 0.0), this};
      auto v = values.begin();
      for (std::size_t i = 0; i < names_.size(); ++i)
        if (derivedFrom_[i] < 0) pt.x[i] = *v++;
      fill_derived(pt);
      return pt;
    }
    /// A point from named values: `F.at({{P, 1.3}, {L, 0.8}, {C, 0.6}})`. Every independent symbol
    /// needs one.
    Point at(std::initializer_list<std::pair<Symbol, double>> values) const
    {
      Point pt{std::vector<double>(names_.size(), 0.0), this};
      std::vector<char> seen(names_.size(), 0);
      for (const auto &[s, v] : values) {
        if (s.frame != this) NT_THROW(std::invalid_argument, "Frame::at: a symbol of another frame");
        if (derivedFrom_[static_cast<std::size_t>(s.id)] >= 0)
          NT_THROW(std::invalid_argument, ("Frame::at: " + names_[static_cast<std::size_t>(s.id)] +
                                           " is derived from its cosine; give the cosine instead")
                                              .c_str());
        pt.x[static_cast<std::size_t>(s.id)] = v;
        seen[static_cast<std::size_t>(s.id)] = 1;
      }
      for (std::size_t i = 0; i < names_.size(); ++i)
        if (derivedFrom_[i] < 0 && !seen[i])
          NT_THROW(std::invalid_argument, ("Frame::at: no value for symbol " + names_[i]).c_str());
      fill_derived(pt);
      return pt;
    }
    /// The polynomial's value at @p pt; the projector denominators it carries are evaluated there too.
    Cx eval(const Poly &p, const Point &pt) const { return ndetail::eval(p, checked(pt).x, atom_values(pt, {&p})); }
    /// A dressed polynomial's value; `dressings[id]` is the value of dressing `id`.
    Cx eval(const DPoly &p, const Point &pt, const std::vector<double> &dressings) const
    {
      std::vector<const Poly *> parts;
      for (const auto &[d, mp] : p.terms) parts.push_back(&mp);
      return ndetail::eval(p, checked(pt).x, atom_values(pt, parts), dressings);
    }

    // ── lowering to C++ ──────────────────────────────────────────────────────────────────────────
    /// The `f[]` slot values of a lowered program (@ref to_genprog) at @p pt, for @ref interpret.
    std::vector<double> fill_values(const GlobalEnv &g, const Point &pt) const
    {
      std::vector<double> atoms(atomDen_.size(), std::nan(""));
      for (std::size_t i = 0; i < g.syms.size(); ++i)
        if (std::get<0>(g.syms[i]) == SymKind::inv) atoms.at(static_cast<std::size_t>(std::get<1>(g.syms[i]))) = atom_value(checked(pt), std::get<1>(g.syms[i]));
      std::vector<double> f(g.syms.size());
      for (std::size_t i = 0; i < g.syms.size(); ++i) {
        const auto [kind, a, b] = g.syms[i];
        if (kind == SymKind::var)
          f[i] = pt.x.at(static_cast<std::size_t>(a));
        else if (kind == SymKind::inv)
          f[i] = atoms.at(static_cast<std::size_t>(a));
        else
          NT_THROW(std::invalid_argument, "Frame::fill_values: the program uses a dressing or scalar-product slot");
      }
      return f;
    }
    /// The formulas that compute each `f[]` slot from the independent symbols (see @ref emit_fill).
    FillFormulas fill_formulas() const
    {
      FillFormulas fm;
      fm.var = [this](int id) {
        const int from = derivedFrom_.at(static_cast<std::size_t>(id));
        if (from < 0) return names_[static_cast<std::size_t>(id)];
        const std::string &c = names_[static_cast<std::size_t>(from)];
        return "std::sqrt(1.0 - " + c + "*" + c + ")";
      };
      fm.inv = [this](int id) { return "1.0/(" + denominator_cpp(id) + ")"; };
      fm.dress = [](int) -> std::string {
        NT_THROW(std::invalid_argument, "Frame::fill_formulas: dressings are filled by generated kernels only");
        return {};
      };
      fm.sp = [](int, int) -> std::string {
        NT_THROW(std::invalid_argument, "Frame::fill_formulas: scalar-product slots do not occur here");
        return {};
      };
      return fm;
    }
    /// Print `void name(double* f, double <sym>, …)`: the function that fills the `f[]` array a
    /// lowered program reads, from the independent symbols.
    void emit_fill(std::ostream &out, const GlobalEnv &g, const std::string &name = "fill",
                   const std::string &decor = "static inline") const
    {
      std::string sig;
      for (std::size_t i = 0; i < names_.size(); ++i)
        if (derivedFrom_[i] < 0) sig += (sig.empty() ? "" : ", ") + std::string("double ") + names_[i];
      ::numtracer::network::emit_fill(out, g, name, sig, fill_formulas(), decor);
    }

    // ── trace-fold phases (generated code; member templates over the free templates above) ────────
    template <class P, class TraceFn> std::vector<P> contract_traces(long nCache, unsigned W, TraceFn &&trace) const
    {
      return ::numtracer::numeric::contract_traces<P>(nsym(), nCache, W, std::forward<TraceFn>(trace));
    }
    /// Streaming phase B: folds each group's nets on demand and drains straight to `sink`, so no net
    /// polynomial outlives its group. See `trace_fold.hpp` for the equivalence argument.
    template <class P, class TraceFn, class ScaleFn, class Sink>
    void fold_groups_streaming(const std::vector<std::vector<int>> &traceIdx,
                               const std::vector<std::vector<Cx>> &subScale,
                               const std::vector<std::vector<int>> &groups, const std::vector<P> &traceTable,
                               long nCache, unsigned W, long window, TraceFn &&trace, ScaleFn &&scale,
                               Sink &&sink) const
    {
      ::numtracer::numeric::fold_groups_streaming<P>(nsym(), traceIdx, subScale, groups, traceTable, nCache, W,
                                                     window, std::forward<TraceFn>(trace),
                                                     std::forward<ScaleFn>(scale), std::forward<Sink>(sink));
    }
    /// Streaming phase B, dressed variant: a plain-Poly `traceTable` + the per-sub-term dressing
    /// monomials `subDress` fold into a DPoly per net. See `trace_fold.hpp`.
    template <class TraceFn, class ScaleFn, class Sink>
    void fold_groups_streaming_dressed(const std::vector<std::vector<int>> &traceIdx,
                                       const std::vector<std::vector<Cx>> &subScale,
                                       const std::vector<std::vector<DMono>> &subDress,
                                       const std::vector<std::vector<int>> &groups,
                                       const std::vector<Poly> &traceTable, long nCache, unsigned W,
                                       long window, TraceFn &&trace, ScaleFn &&scale, Sink &&sink) const
    {
      ::numtracer::numeric::fold_groups_streaming_dressed(nsym(), traceIdx, subScale, subDress, groups,
                                                          traceTable, nCache, W, window,
                                                          std::forward<TraceFn>(trace),
                                                          std::forward<ScaleFn>(scale), std::forward<Sink>(sink));
    }

  private:
    std::vector<std::string> names_;               ///< symbol names, in declaration order
    std::vector<int> derivedFrom_;                 ///< per symbol: the cosine an angle's sine derives from, or -1
    std::vector<std::vector<int>> units_;          ///< unit-vector groups (Σ x² = 1)
    /// Set once the first polynomial is made. Atomic because the const members may run in parallel;
    /// they only ever read it once the frame is frozen (generated frames are frozen at construction).
    mutable std::atomic<bool> frozen_{false};
    std::vector<std::array<Poly, 4>> comp_;        ///< comp_[vid][mu]
    std::vector<Poly> atomDen_;                    ///< atomDen_[atom] = k² (empty = unused id)
    std::map<std::pair<bool, Vlc>, int> autoAtom_; ///< (spatial?, momentum up to sign) → assigned atom
    int nextIndex_ = 0;                            ///< next fresh Lorentz label

    template <int N, int... I> std::array<LorentzIndex, N> make_indices(std::integer_sequence<int, I...>) const
    {
      return {LorentzIndex{nextIndex_ + I}...};
    }
    void require_open(const char *what) const
    {
      if (frozen_)
        NT_THROW(std::logic_error, (std::string("Frame::") + what +
                                    ": declare every symbol before making a momentum or polynomial from them")
                                       .c_str());
    }
    int frozen_nsym() const
    {
      if (!frozen_) frozen_ = true; // read-only once frozen: the const contractions run in parallel
      return nsym();
    }
    void check_owned(const Poly &p, const char *what) const
    {
      if (!p.empty() && p.nsym != nsym())
        NT_THROW(std::invalid_argument, (std::string(what) + ": a polynomial of another frame").c_str());
    }
    std::size_t independent_count() const
    {
      std::size_t n = 0;
      for (int d : derivedFrom_) n += d < 0;
      return n;
    }
    std::string independent_list() const
    {
      std::string s;
      for (std::size_t i = 0; i < names_.size(); ++i)
        if (derivedFrom_[i] < 0) s += (s.empty() ? "" : ", ") + names_[i];
      return s;
    }
    void fill_derived(Point &pt) const
    {
      for (std::size_t i = 0; i < names_.size(); ++i)
        if (const int c = derivedFrom_[i]; c >= 0) {
          const double cv = pt.x[static_cast<std::size_t>(c)];
          if (cv < -1.0 || cv > 1.0)
            NT_THROW(std::domain_error, ("Frame::at: cosine " + names_[static_cast<std::size_t>(c)] + " = " +
                                         std::to_string(cv) + " is outside [-1, 1]")
                                            .c_str());
          pt.x[i] = std::sqrt(1.0 - cv * cv);
        }
    }
    const Point &checked(const Point &pt) const
    {
      if (pt.frame != this)
        NT_THROW(std::invalid_argument, "Frame::eval: the point belongs to another frame (use this frame's at())");
      return pt;
    }
    /// `1/k²` of atom @p a at @p pt; refuses an unregistered atom and a vanishing denominator.
    double atom_value(const Point &pt, int a) const
    {
      if (a < 0 || static_cast<std::size_t>(a) >= atomDen_.size() || atomDen_[static_cast<std::size_t>(a)].empty())
        NT_THROW(std::invalid_argument, ("Frame::eval: atom " + std::to_string(a) + " is not registered in this frame").c_str());
      const double d = ndetail::eval(atomDen_[static_cast<std::size_t>(a)], pt.x, {}).re;
      if (d == 0.0)
        NT_THROW(std::domain_error, ("Frame::eval: projector denominator " + std::to_string(a) + " (k^2 = " +
                                     denominator_cpp(a) + ") vanishes at this point")
                                        .c_str());
      return 1.0 / d;
    }
    /// The atom values the polynomials @p parts reference, evaluated at @p pt (others stay NaN).
    std::vector<double> atom_values(const Point &pt, const std::vector<const Poly *> &parts) const
    {
      std::vector<double> v(atomDen_.size(), std::nan(""));
      std::vector<char> done(atomDen_.size(), 0);
      for (const Poly *p : parts)
        for (const auto &[m, c] : p->terms)
          for (int a : m.atoms) {
            if (a >= 0 && static_cast<std::size_t>(a) < done.size() && done[static_cast<std::size_t>(a)]) continue;
            const double val = atom_value(pt, a);
            v[static_cast<std::size_t>(a)] = val;
            done[static_cast<std::size_t>(a)] = 1;
          }
      return v;
    }
    void set_denominator(int id, Poly den)
    {
      if (static_cast<std::size_t>(id) >= atomDen_.size()) atomDen_.resize(static_cast<std::size_t>(id) + 1, zero());
      Poly &slot = atomDen_[static_cast<std::size_t>(id)];
      const auto same = [](const Poly &x, const Poly &y) {
        if (x.terms.size() != y.terms.size()) return false;
        for (std::size_t i = 0; i < x.terms.size(); ++i)
          if (!(x.terms[i].first == y.terms[i].first) || x.terms[i].second.re != y.terms[i].second.re ||
              x.terms[i].second.im != y.terms[i].second.im)
            return false;
        return true;
      };
      if (!slot.empty() && !same(slot, den))
        NT_THROW(std::invalid_argument,
                 ("Frame: atom id " + std::to_string(id) + " already names a different denominator").c_str());
      slot = std::move(den);
    }
    static bool needs_atoms(const LorentzNet &lor)
    {
      for (const LorentzTerm &t : lor)
        for (const LorentzFactor &f : t.e) {
          if (!f.is_projector()) continue;
          if (f.kind != LorentzFactor::ProjM && f.atom < 0) return true;
          if ((f.kind == LorentzFactor::ProjE || f.kind == LorentzFactor::ProjM) && f.atomS < 0) return true;
        }
      return false;
    }
    /// The atom of projector @p f's `k²` (or `|k⃗|²` when @p spatial): one id per momentum up to
    /// sign, its denominator registered on first use.
    int atom_for(const LorentzFactor &f, bool spatial)
    {
      Vlc key = f.vlc.empty() ? Vlc{{1.0, f.vid}} : f.vlc;
      if (key.front().first < 0)
        for (auto &t : key) t.first = -t.first;
      const auto it = autoAtom_.find({spatial, key});
      if (it != autoAtom_.end()) return it->second;
      const std::array<Poly, 4> k = ndetail::factor_momentum(frozen_nsym(), f, comp_);
      Poly den = zero();
      for (int mu = spatial ? 1 : 0; mu < 4; ++mu)
        den = den + k[mu] * k[mu];
      const int id = static_cast<int>(atomDen_.size());
      set_denominator(id, reduce_units(std::move(den), units_));
      autoAtom_.emplace(std::make_pair(spatial, std::move(key)), id);
      return id;
    }
    void assign_atoms(LorentzFactor &f)
    {
      if (!f.is_projector()) return;
      if (f.kind != LorentzFactor::ProjM && f.atom < 0) f.atom = atom_for(f, false);
      if ((f.kind == LorentzFactor::ProjE || f.kind == LorentzFactor::ProjM) && f.atomS < 0)
        f.atomS = atom_for(f, true);
    }
    LorentzNet assign_atoms(LorentzNet lor)
    {
      for (LorentzTerm &t : lor)
        for (LorentzFactor &f : t.e) assign_atoms(f);
      return lor;
    }
  };

  inline Symbol::operator Poly() const { return frame->var(id); }

  // A Symbol converts to its Poly, so `L * C` uses Poly ⊗ Poly (mpoly.hpp); with a number:
  inline Poly operator*(double c, Symbol s) { return c * Poly(s); }
  inline Poly operator*(Symbol s, double c) { return c * Poly(s); }

} // namespace numtracer::numeric
