/// @file sun_net.hpp
/// @brief Numeric SU(N) colour/flavour-factor contraction for the **build-time generator** —
///        contracts a fully-contracted network of structure constants `f^{abc}`, adjoint deltas
///        `δ^{ab}`, fundamental generators `T^a` and fundamental deltas `δ^{ij}` to a single
///        (complex) number, at codegen time, so the kernel never instantiates an SU(N) tensor.
///
/// A single diagram's SU(N) network is several `f`s contracted; for the four-gluon vertex it is
/// too large to carry through a compile-time expression tensor (which would OOM the compiler — even
/// one A4 net), and a sufficiently large *fundamental* (quark-loop) network has the same problem.
/// The group factor is just a number, so the generator computes it here (sparse contraction over
/// the SU(N) generators and structure constants) and emits a `constexpr` literal. This mirrors the
/// Lorentz invariant reduction: the contraction runs in C++ in the generator, and only the result
/// reaches the kernel. (Generator-only — cold path: not part of any consumer kernel.)
///
/// This is a **general** SU(N) contraction: the group rank `N` is a *runtime* value, not a template
/// parameter. The contraction is generator-only (cold), so templating over `N` would buy no runtime
/// speed and would force a fixed enumeration of supported `N`; instead `N` is carried at run time and
/// the SU(N) data (generators + nonzero `f^{abc}`) is built once per `N` and cached (@ref
/// sun_net_detail::sun_data_for). For the tabulated ranks `N ∈ {2,3}` the data is copied verbatim
/// from the validated typed-out @ref numtracer::sun::SUNData tables, so the folded values are byte-identical
/// to the typed path; for any other `N` it is built from the generalized-Gell-Mann construction (the
/// same one the @ref numtracer::sun::SUNBuilder oracle uses).
///
/// A network may carry **two independent groups at once** (e.g. colour SU(3) ⊗ flavour SU(2)).
/// Each @ref SUNFac is tagged with its group's rank `g`; the two groups have disjoint label spaces,
/// so the value factorises — @ref sun_value_cx contracts each group separately and multiplies.
#pragma once

#include "numtracer/core/export.hpp" // NUMTRACER_FUNC / NUMTRACER_DEFINE_BODIES (compiled vs header-only)
#include "numtracer/sun/sun_data.hpp" // SUNData<N> tables, Mat<N>, matmul, trace, Cx, sun_detail::FEntry
#include "numtracer/core/config.hpp"  // NT_THROW (exception-optional guard for -fno-exceptions builds)

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string> // the open-leg guard's diagnostic (assert_no_open_labels)
#include <vector>

namespace numtracer::network {

/// @brief An SU(N) colour/flavour-network factor.
///
/// `kind` selects the structure; `g` is the group rank (e.g. 3 = SU(3) colour, 2 = SU(2) flavour),
/// so one network can mix groups. The index fields carry adjoint (extent `N²−1`) or fundamental
/// (extent `N`) labels depending on `kind`:
///   - kind 0 `delta_adj(a,b)`     : `a,b` adjoint
///   - kind 1 `f(a,b,c)`           : `a,b,c` adjoint
///   - kind 2 `T^a_{ij}` generator : `a` adjoint, `b = i` fund (row), `c = j` fund (col)
///   - kind 3 `delta_fund(i,j)`    : `a = i`, `b = j` fundamental
///   - kind 4 `diag_fund(i,j)`     : `a = i`, `b = j` fundamental — a δ^{ij} carrying a *per-component*
///                                   dressing, i.e. the diagonal insertion `diag(D_0..D_{N-1})`
///   - kind 5 `diag_adj(a,b)`      : `a,b` adjoint — a δ^{ab} carrying a per-component adjoint
///                                   dressing, i.e. `diag(D_0..D_{N²−2})`
///
/// Kinds 4/5 are the *group-diagonal dressing* insertions: instead of folding to a constant, a net
/// carrying them folds (via @ref sun_value_dressed) to a polynomial `Σ_a c_a D_a` over named runtime
/// dressing symbols. The per-component dressing is carried by @ref SUNFac::comp2dr: entry `v` is the
/// dressing-id for component `v`, or `-1` to *drop* that component (it contributes nothing). This lets
/// the front-end dress only selected components (e.g. the Cartan directions of a condensate) and drop
/// the rest — no dead terms. `comp2dr` is empty for the plain (kinds 0–3) factors.
/// @brief The structure a @ref SUNFac selects. The explicit 0–5 values are only for self-documentation;
///        no code outside this header depends on the specific numbers.
enum class SUNFacKind : int {
  DeltaAdj = 0,  ///< `δ^{ab}` (adjoint Kronecker)
  F = 1,         ///< `f^{abc}` (structure constant)
  T = 2,         ///< `(T^a)_{ij}` (fundamental generator)
  DeltaFund = 3, ///< `δ^{ij}` (fundamental Kronecker)
  DiagFund = 4,  ///< `diag_fund(i,j;dr)` (fundamental group-diagonal dressing)
  DiagAdj = 5,   ///< `diag_adj(a,b;dr)` (adjoint group-diagonal dressing)
};

struct SUNFac {
  SUNFacKind kind; ///< which structure (see @ref SUNFacKind and the table above)
  int g;        ///< group rank `N` (3 colour, 2 flavour, …)
  int a, b, c;  ///< adjoint/fundamental index labels (see kind table above; `c` unused for kinds 0/3/4/5)
  std::vector<int> comp2dr; ///< kinds 4/5 only: per-component dressing-id (`-1` = drop); empty otherwise
};
/// @brief A fully-contracted SU(N) colour/flavour network (every label appears exactly twice).
using SUNNet = std::vector<SUNFac>;

/// @brief One monomial of a dressed SU(N) value: a constant times a product of named dressing symbols.
///
/// `dress` is a sorted list of dressing-ids — each names one runtime dressing `D^{dr}` (a kind-4/5
/// component that survived the fold). A repeated id is a power. An empty `dress` is a plain constant.
/// The codegen maps each id to a scalar dressing symbol and emits the standard `name(scale)` token.
struct SUNTerm {
  Cx coeff{1.0, 0.0};
  std::vector<int> dress; ///< sorted dressing-ids (repetition = power)
};
/// @brief A dressed SU(N) value: `Σ_t coeff_t · Π D^{dr}` (a small polynomial over the named runtime
///        dressing symbols). With no kind-4/5 factor it is a single constant term.
using SUNPoly = std::vector<SUNTerm>;

/// @brief Builders for SU(N) colour/flavour-network factors.
///
/// Each returns one @ref SUNFac tagged with its group rank `g` (e.g. `3` colour, `2` flavour);
/// assemble them into a @ref SUNNet and contract with @ref sun_value_cx. The names mirror the
/// algebra: `f` (structure constant), `deltaAdj`/`deltaFund` (adjoint/fundamental Kronecker),
/// `T` (fundamental generator).
namespace SUN {
inline SUNFac f(int g, int a, int b, int c) { return {SUNFacKind::F, g, a, b, c, {}}; }   ///< structure constant `f^{abc}`
inline SUNFac deltaAdj(int g, int a, int b) { return {SUNFacKind::DeltaAdj, g, a, b, -1, {}}; }  ///< adjoint Kronecker `δ^{ab}`
inline SUNFac T(int g, int a, int i, int j) { return {SUNFacKind::T, g, a, i, j, {}}; }   ///< fundamental generator `(T^a)_{ij}`
inline SUNFac deltaFund(int g, int i, int j) { return {SUNFacKind::DeltaFund, g, i, j, -1, {}}; } ///< fundamental Kronecker `δ^{ij}`
/// fundamental diagonal dressing: `δ^{ij}` carrying a per-component dressing map (component → id, `-1` = drop).
inline SUNFac diagFund(int g, int i, int j, std::vector<int> comp2dr) { return {SUNFacKind::DiagFund, g, i, j, -1, std::move(comp2dr)}; }
/// adjoint diagonal dressing: `δ^{ab}` carrying a per-component dressing map (component → id, `-1` = drop).
inline SUNFac diagAdj(int g, int a, int b, std::vector<int> comp2dr) { return {SUNFacKind::DiagAdj, g, a, b, -1, std::move(comp2dr)}; }
} // namespace SUN

/// @brief A factory that binds the SU(N) group rank `g` once and mints @ref SUNFac factors without
///        repeating it at every call — the colour/flavour analogue of @ref numtracer::numeric::LorentzEnv
///        (which binds `nsym`). Named `SUNEnv` rather than "ColourEnv" because SU(N) also carries flavour;
///        a single kernel may mix ranks (e.g. colour SU(3) ⊗ flavour SU(2)), so instantiate one `SUNEnv`
///        per rank. Thin ergonomic wrapper over the `SUN::` builders (which stay public), so it needs no
///        friendship. `sun_value_cx`/`sun_value_dressed` take a fully-built @ref SUNNet and are unchanged.
class SUNEnv {
public:
  explicit SUNEnv(int g) : g_(g) {}
  int rank() const { return g_; }

  SUNFac f(int a, int b, int c) const { return SUN::f(g_, a, b, c); }
  SUNFac deltaAdj(int a, int b) const { return SUN::deltaAdj(g_, a, b); }
  SUNFac T(int a, int i, int j) const { return SUN::T(g_, a, i, j); }
  SUNFac deltaFund(int i, int j) const { return SUN::deltaFund(g_, i, j); }
  SUNFac diagFund(int i, int j, std::vector<int> comp2dr) const { return SUN::diagFund(g_, i, j, std::move(comp2dr)); }
  SUNFac diagAdj(int a, int b, std::vector<int> comp2dr) const { return SUN::diagAdj(g_, a, b, std::move(comp2dr)); }

private:
  int g_; ///< group rank `N` (3 colour, 2 flavour, …)
};

// The whole SU(N) contraction machinery (dynamic matrices, generator/f-table build+cache, group
// contraction) is internal to the generator: its bodies compile only in the library TU (or a
// header-only build). A normal consumer sees just the SUNFac/SUNNet types + SUN builders above and
// links the sun_value* entry points declared below. See core/export.hpp.
#if NUMTRACER_DEFINE_BODIES
namespace sun_net_detail {

/// @brief A runtime-sized N×N complex matrix (row-major) for the cold generator path.
///
/// The hot SU(N) oracle uses the stack-allocated, compile-time-sized @ref numtracer::Mat; here `N`
/// is a runtime value, so a heap-backed dynamic matrix is used instead. Generator-only, so the heap
/// allocation and lack of unrolling are irrelevant.
struct DynMat {
  int n = 0;                                  ///< Dimension `N`.
  std::vector<std::complex<double>> data;     ///< Row-major element storage (`n*n`).
  DynMat() = default;
  explicit DynMat(int N) : n(N), data(N * N) {}
  std::complex<double> &operator()(int i, int j) { return data[i * n + j]; }
  std::complex<double> operator()(int i, int j) const { return data[i * n + j]; }
};

/// @brief Matrix product, skipping structurally zero entries of the left operand (generators are sparse).
inline DynMat dmatmul(const DynMat &x, const DynMat &y) {
  DynMat r(x.n);
  for (int i = 0; i < x.n; ++i)
    for (int k = 0; k < x.n; ++k) {
      const std::complex<double> xik = x(i, k);
      if (xik == std::complex<double>{0, 0}) continue;
      for (int j = 0; j < x.n; ++j) r(i, j) += xik * y(k, j);
    }
  return r;
}

/// @brief Matrix trace.
inline std::complex<double> dtrace(const DynMat &x) {
  std::complex<double> t{0, 0};
  for (int i = 0; i < x.n; ++i) t += x(i, i);
  return t;
}

/// @brief Disjoint-set (union-find) over index labels, with path-halving on `find`.
///
/// A δ (or a diagonal-dressing δ) identifies its two indices; `find` returns a label's class
/// representative, so every Lorentz/colour label sharing a class is one contracted index. Adjoint and
/// fundamental labels are disjoint integers, so they never merge across sectors. Shared by both
/// @ref contract_group and @ref contract_group_dressed (identical index-identification step).
struct UnionFind {
  std::vector<int> par; ///< parent pointers; a root is its own parent
  explicit UnionFind(int maxLabel) : par(maxLabel + 1) {
    for (int i = 0; i <= maxLabel; ++i) par[i] = i;
  }
  int find(int x) {
    while (par[x] != x) {
      par[x] = par[par[x]];
      x = par[x];
    }
    return x;
  }
  void unite(int x, int y) { par[find(x)] = find(y); }
};

/// @brief The SU(N) data the contraction reads: the fundamental generators and the nonzero `f^{abc}`.
struct SUNDyn {
  int N = 0;                                  ///< Group rank.
  std::vector<DynMat> gens;                   ///< The `N²−1` fundamental generators `T^a`.
  std::vector<sun::sun_detail::FEntry> f_nz;  ///< The nonzero structure constants `f^{abc}`.
};

/// @brief Copy the validated typed-out @ref numtracer::sun::SUNData<N> tables into a runtime @ref SUNDyn.
///
/// Used for the tabulated ranks so the folded values are byte-identical to the typed path.
template <int N> SUNDyn seed_from_table() {
  static_assert(sun::kHasSUNData<N>, "seed_from_table needs a typed-out SUNData<N>");
  SUNDyn d;
  d.N = N;
  constexpr int adj = N * N - 1;
  d.gens.assign(adj, DynMat(N));
  for (int a = 0; a < adj; ++a)
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        d.gens[a](i, j) = sun::SUNData<N>::generators[a](i, j);
  for (const auto &e : sun::SUNData<N>::f_nonzeros) d.f_nz.push_back(e);
  return d;
}

/// @brief Build SU(N) data for an arbitrary rank from the generalized-Gell-Mann construction.
///
/// Mirrors @ref numtracer::sun::SUNBuilder::build_generators / build_structure_constants exactly, with
/// `f^{abc} = -2 i\,\mathrm{tr}([T^a,T^b]T^c)`, using the dynamic matrix.
inline SUNDyn build_oracle(int N) {
  if (N < 1) NT_THROW(std::runtime_error, "sun_net: build_oracle requires group rank N >= 1");
  SUNDyn d;
  d.N = N;
  const int adj = N * N - 1; // ≥ 0 for N ≥ 1 (guarded above)
  d.gens.assign(adj, DynMat(N));
  const std::complex<double> I{0, 1};
  int idx = 0;
  for (int k = 0; k < N; ++k)
    for (int l = k + 1; l < N; ++l) {
      DynMat m(N); m(k, l) = {0.5, 0}; m(l, k) = {0.5, 0};
      d.gens[idx++] = m;
      DynMat n(N); n(k, l) = -I * 0.5; n(l, k) = I * 0.5;
      d.gens[idx++] = n;
    }
  for (int l = 1; l <= N - 1; ++l) {
    DynMat m(N);
    const double norm = std::sqrt(1.0 / (2.0 * l * (l + 1)));
    for (int j = 0; j < l; ++j) m(j, j) = {norm, 0};
    m(l, l) = {-norm * l, 0};
    d.gens[idx++] = m;
  }
  if (idx != adj) NT_THROW(std::logic_error, "sun_net: wrong generator count");
  // structure constants below this magnitude are floating-point round-off in the commutator trace,
  // not genuine nonzero entries, so they are dropped.
  constexpr double kStructConstTol = 1e-12;
  for (int a = 0; a < adj; ++a)
    for (int b = 0; b < adj; ++b) {
      const DynMat ab = dmatmul(d.gens[a], d.gens[b]);
      const DynMat ba = dmatmul(d.gens[b], d.gens[a]);
      DynMat comm(N);
      for (std::size_t i = 0; i < comm.data.size(); ++i) comm.data[i] = ab.data[i] - ba.data[i];
      for (int c = 0; c < adj; ++c) {
        const std::complex<double> tr = dtrace(dmatmul(comm, d.gens[c]));
        const double val = (std::complex<double>{0, -2} * tr).real();
        if (std::fabs(val) > kStructConstTol) d.f_nz.push_back({a, b, c, val});
      }
    }
  return d;
}

/// @brief The SU(N) data for rank `N`, built once and cached.
///
/// Tabulated ranks (`N ∈ {2,3}`) are seeded byte-identically from the typed-out tables; any other
/// `N` is built from the generalized-Gell-Mann oracle. Thread-safe: the generator contracts nets in
/// parallel. `std::map` node addresses are stable across inserts and a built @ref SUNDyn is never
/// mutated, so the returned reference stays valid for concurrent read-only use.
inline const SUNDyn &sun_data_for(int N) {
  static std::mutex cacheMutex;
  static std::map<int, SUNDyn> cache;
  std::lock_guard<std::mutex> lock(cacheMutex);
  auto it = cache.find(N);
  if (it == cache.end()) {
    if (N < 1) NT_THROW(std::runtime_error, "sun_net: group rank N must be >= 1");
    SUNDyn d = (N == 2) ? seed_from_table<2>() : (N == 3) ? seed_from_table<3>() : build_oracle(N);
    it = cache.emplace(N, std::move(d)).first;
  }
  return it->second;
}

/// @brief Reject an OPEN (once-occurring) colour/flavour label before it is silently closed away.
///
/// The adjoint sector has the same hazard as the Lorentz one (see
/// @ref numtracer::numeric::ndetail::assert_no_open_ids), by two different mechanisms:
///   * an open adjoint index on an `f` or a generator is DENSE-SUMMED over `0..N²−2` — `f^{abc}f^{abd}`
///     with `c`,`d` open returned 24 instead of the rank-2 `N δ^{cd}`;
///   * an open index on a `δ` is worse: union-find identifies the δ's two labels, the class is then
///     "touched only by δ" and is counted as a CLOSED LOOP — a lone `δ^{ab}` returned `N²−1 = 8`.
/// Only the open FUNDAMENTAL chain was caught (by @ref extract_cycles), and only when a generator sits
/// on it; a fundamental `δ` line with a free end closed to `N` just like the adjoint case.
///
/// Counting occurrences here — before the union-find, which is what destroys the evidence — catches all
/// of them uniformly. Labels are counted per group (@ref sun_value_cx already partitions by rank, and
/// each group owns a disjoint label space), a label used twice within one factor (`δ^{aa}`, a legal
/// closed loop) counts 2 and passes, and counts ≥ 3 are left alone so no valid net changes value.
inline void assert_no_open_labels(int N, const std::vector<const SUNFac *> &net) {
  std::map<int, int> cnt; // label -> occurrences (nets are small; the map keeps the report ordered)
  for (const SUNFac *f : net) {
    ++cnt[f->a];
    ++cnt[f->b];
    // only `f^{abc}` and `T^a_{ij}` carry a third label; every other kind leaves `c` at -1
    if (f->kind == SUNFacKind::F || f->kind == SUNFacKind::T) ++cnt[f->c];
  }
  for (const auto &[lbl, n] : cnt)
    if (n == 1) {
      const std::string msg =
          "sun_net: SU(" + std::to_string(N) + ") index label " + std::to_string(lbl) +
          " is OPEN (occurs once) — the colour/flavour net does not close to a scalar. Tie the leg off "
          "(a δ, an f, a generator trace, or a colour projector). (Contracting an open leg would "
          "silently sum or δ-close it and return a meaningless number.)";
      NT_THROW(std::runtime_error, msg.c_str());
    }
}

/// @brief Walk the generator chain into its disjoint cycles.
///
/// Each generator `T^a_{ij}` is a directed edge rowClass(i) → colClass(j); a closed contraction
/// makes these edges disjoint cycles, one `tr(T^{a_1}…T^{a_m})` per cycle.
///
/// `field` picks WHAT each visited generator contributes to its cycle, and is the only difference
/// between the two things callers want:
///   - `0` → the generator's ADJOINT class. All the undressed trace needs; a dense matrix product
///           sums the fundamental indices implicitly. Used by @ref contract_group.
///   - `2` → the generator's COLUMN class, i.e. the FUNDAMENTAL index sitting immediately AFTER it,
///           so the last entry is the one closing the trace. A `diagFund` weights one such index per
///           component, so the dressed contraction has to know which class sits where. Used by
///           @ref contract_group_dressed.
///
/// One walk, so the two results are aligned entry-for-entry by construction. They were two copies of
/// this function differing in that single subscript, with a comment asking the reader to keep the
/// orders in step by hand.
inline std::vector<std::vector<int>> extract_cycles(const std::vector<std::array<int, 3>> &gens,
                                                    std::size_t field)
{
  std::map<int, std::size_t> rowToGen; // rowClass(i) -> generator index
  for (std::size_t gi = 0; gi < gens.size(); ++gi) rowToGen[gens[gi][1]] = gi;
  std::vector<char> seen(gens.size(), 0);
  std::vector<std::vector<int>> cycles;
  for (std::size_t start = 0; start < gens.size(); ++start) {
    if (seen[start]) continue;
    std::vector<int> cycle;
    std::size_t curGen = start;
    do {
      seen[curGen] = 1;
      cycle.push_back(gens[curGen][field]);
      auto it = rowToGen.find(gens[curGen][2]); // next generator: its row == this col
      if (it == rowToGen.end())
        NT_THROW(std::runtime_error, "sun_net: open fundamental chain (only closed traces supported)");
      curGen = it->second;
    } while (curGen != start);
    cycles.push_back(std::move(cycle));
  }
  return cycles;
}

/// Adjoint classes per cycle — see @ref extract_cycles.
inline std::vector<std::vector<int>> extract_cycles_adj(const std::vector<std::array<int, 3>> &gens)
{
  return extract_cycles(gens, 0);
}

/// Fundamental (column) classes per cycle, aligned with @ref extract_cycles_adj — see @ref extract_cycles.
inline std::vector<std::vector<int>> extract_cycles_fund(const std::vector<std::array<int, 3>> &gens)
{
  return extract_cycles(gens, 2);
}

/// @brief The product of generator traces `∏_cycles tr(T^{a_1}…T^{a_m})` for a fully-pinned adjoint
///        assignment @p classVal (class → component `0..Adim-1`). Shared by both contractions.
inline Cx loop_prod(const SUNDyn &dat, const std::vector<std::vector<int>> &cycles,
                    const std::map<int, int> &classVal) {
  Cx prod{1.0, 0.0};
  for (const auto &cycle : cycles) {
    DynMat m = dat.gens[classVal.at(cycle[0])];
    for (std::size_t k = 1; k < cycle.size(); ++k)
      m = dmatmul(m, dat.gens[classVal.at(cycle[k])]);
    const std::complex<double> tr = dtrace(m);
    prod = prod * Cx{tr.real(), tr.imag()};
  }
  return prod;
}

/// @brief Contract a single-group network (all factors share rank `N`) to its scalar value.
///
/// Generalises the adjoint-only union-find + sparse-`f`-backtracking to also handle fundamental
/// generator traces: the fundamental indices form closed cycles (one per quark loop), each cycle a
/// `tr(T^{a_1}…T^{a_m})`; the adjoint indices on those generators are pinned by the `f`-backtracking
/// when shared with an `f`, else summed densely. Adjoint indices touched only by `δ^{ab}` close an
/// adjoint loop (`δ^{aa} = N²−1`); fundamental indices touched only by `δ^{ij}` close a fundamental
/// loop (`δ^{ii} = N`).
inline Cx contract_group(int N, const std::vector<const SUNFac *> &net) {
  // Algorithm:
  //   1. union-find: each δ identifies its two indices, so every label collapses to an index CLASS.
  //   2. classify classes by sector (adjoint / fundamental) and by what touches them (δ only, an `f`,
  //      a generator). A class touched only by δ is a CLOSED LOOP: δ^{aa}=N²−1 adjoint / δ^{ii}=N fund.
  //   3. the generators' fundamental indices form closed directed cycles → one tr(T^{a_1}…T^{a_m}) each.
  //   4. sum the value: backtrack over the nonzero `f^{abc}` table to pin the f-shared adjoint classes
  //      (`sumF`), dense-sum the remaining generator-only adjoint classes (`sumGen`), and for each
  //      consistent assignment fold in the product of generator traces (`loopProd`). × the closed-loop
  //      scalar. Snap √3-residual zeros from the generator-table arithmetic.
  // Step 1 (union-find) destroys the evidence an open leg leaves, so the closure check comes first.
  assert_no_open_labels(N, net);
  const SUNDyn &dat = sun_data_for(N);
  const int Adim = N * N - 1; // adjoint dimension

  // ---- union-find over labels: a delta identifies its two indices (adjoint *or* fundamental) ----
  int maxlbl = -1;
  for (const SUNFac *f : net) {
    maxlbl = std::max({maxlbl, f->a, f->b});
    if (f->kind == SUNFacKind::F || f->kind == SUNFacKind::T) maxlbl = std::max(maxlbl, f->c);
  }
  if (maxlbl < 0) return Cx{1.0, 0.0}; // empty colour net: no factor ⇒ identity (avoids UnionFind(-1))
  UnionFind uf(maxlbl);
  auto find = [&](int x) { return uf.find(x); };
  for (const SUNFac *f : net)
    if (f->kind == SUNFacKind::DeltaAdj || f->kind == SUNFacKind::DeltaFund)
      uf.unite(f->a, f->b); // delta_adj / delta_fund identify their two indices

  // ---- classify the index classes touched in each sector --------------------------------------
  std::set<int> adjClasses, fClasses, genAdjClasses;       // adjoint
  std::set<int> fundClasses, genFundClasses;               // fundamental
  std::vector<std::array<int, 3>> fTriples;                // f^{abc} index triples (by class)
  std::vector<std::array<int, 3>> gens;                    // {adjClass, rowClass(i), colClass(j)} per generator
  for (const SUNFac *f : net) {
    switch (f->kind) {
    case SUNFacKind::DeltaAdj:
      adjClasses.insert(find(f->a)); adjClasses.insert(find(f->b));
      break;
    case SUNFacKind::F: { // f^{abc}
      const int a = find(f->a), b = find(f->b), c = find(f->c);
      adjClasses.insert(a); adjClasses.insert(b); adjClasses.insert(c);
      fClasses.insert(a); fClasses.insert(b); fClasses.insert(c);
      fTriples.push_back({a, b, c});
      break;
    }
    case SUNFacKind::T: { // T^a_{ij}
      const int a = find(f->a), i = find(f->b), j = find(f->c);
      adjClasses.insert(a); genAdjClasses.insert(a);
      fundClasses.insert(i); fundClasses.insert(j);
      genFundClasses.insert(i); genFundClasses.insert(j);
      gens.push_back({a, i, j});
      break;
    }
    case SUNFacKind::DeltaFund:
      fundClasses.insert(find(f->a)); fundClasses.insert(find(f->b));
      break;
    default: NT_THROW(std::runtime_error, "sun_net: unknown SUNFac kind");
    }
  }

  // ---- closed-loop scalar factors ----
  // closed adjoint loops (δ^{aa} = N²−1): classes touched only by delta_adj (not by f, not a generator).
  double factorScalar = 1.0;
  for (int c : adjClasses)
    if (!fClasses.count(c) && !genAdjClasses.count(c)) factorScalar *= static_cast<double>(Adim);
  // closed fundamental loops (δ^{ii} = N): fundamental classes touched only by delta_fund.
  for (int c : fundClasses)
    if (!genFundClasses.count(c)) factorScalar *= static_cast<double>(N);

  // generator adjoint classes not pinned by an f are summed densely below.
  std::vector<int> genOnly;
  for (int c : genAdjClasses)
    if (!fClasses.count(c)) genOnly.push_back(c);

  // ---- fundamental-cycle extraction (generator traces) ----
  const std::vector<std::vector<int>> cycles = extract_cycles_adj(gens);
  // the product of generator traces for a fully-pinned adjoint assignment (`classVal`: class -> 0..Adim-1).
  auto loopProd = [&](const std::map<int, int> &classVal) { return loop_prod(dat, cycles, classVal); };

  // dense sum over the gen-only adjoint classes (after the f-classes are pinned).
  Cx total{0.0, 0.0};
  std::map<int, int> classVal; // adjoint class -> its pinned/summed component value
  auto sumGen = [&](auto &&self, std::size_t gi, Cx fProd) -> void {
    if (gi == genOnly.size()) { total = total + fProd * loopProd(classVal); return; }
    const int cls = genOnly[gi];
    for (int v = 0; v < Adim; ++v) { classVal[cls] = v; self(self, gi + 1, fProd); }
    classVal.erase(cls);
  };

  // ---- sparse sum ∏ f over consistent assignments of the f-classes (then fold in generator traces) ----
  // Backtrack over the nonzero f-table: each f-triple's three classes must agree with whatever earlier
  // triples already pinned them (`consistent`); newly-pinned classes (`newA/newB/newC`) are un-pinned on return.
  auto sumF = [&](auto &&self, std::size_t fi, Cx fProd) -> void {
    if (fi == fTriples.size()) { sumGen(sumGen, 0, fProd); return; }
    const auto [ca, cb, cc] = fTriples[fi];
    for (const auto &e : dat.f_nz) {
      auto consistent = [&](int cl, int v) { auto it = classVal.find(cl); return it == classVal.end() || it->second == v; };
      if (!consistent(ca, e.a) || !consistent(cb, e.b) || !consistent(cc, e.c)) continue;
      const bool newA = !classVal.count(ca), newB = !classVal.count(cb), newC = !classVal.count(cc);
      classVal[ca] = e.a; classVal[cb] = e.b; classVal[cc] = e.c;
      self(self, fi + 1, fProd * Cx{e.v, 0.0});
      if (newA) classVal.erase(ca);
      if (newB) classVal.erase(cb);
      if (newC) classVal.erase(cc);
    }
  };
  sumF(sumF, 0, Cx{1.0, 0.0});

  Cx r = total * Cx{factorScalar, 0.0};
  // colour/flavour factors are exact rationals (× generator traces); snap the √3-residual zeros that
  // the generator-table arithmetic leaves just above 0 (well below any genuine rational magnitude).
  constexpr double kZeroSnapTol = 1e-9;
  if (std::fabs(r.re) < kZeroSnapTol) r.re = 0.0;
  if (std::fabs(r.im) < kZeroSnapTol) r.im = 0.0;
  return r;
}

// ---- SUNPoly helpers (dressed contraction) --------------------------------------------------
/// @brief Add `c · Π D_{dress}` into `p`, merging into an existing term with the same dress key.
inline void poly_add_term(SUNPoly &p, Cx c, std::vector<int> key) {
  if (c.re == 0.0 && c.im == 0.0) return; // adding 0 changes nothing
  std::sort(key.begin(), key.end());
  for (auto &t : p)
    if (t.dress == key) { t.coeff = t.coeff + c; return; }
  p.push_back(SUNTerm{c, std::move(key)});
}
/// @brief Polynomial product: distribute over terms, concatenating (then sorting) the dress keys.
inline SUNPoly poly_mul(const SUNPoly &a, const SUNPoly &b) {
  SUNPoly r;
  for (const auto &s : a)
    for (const auto &t : b) {
      std::vector<int> key;
      key.reserve(s.dress.size() + t.dress.size());
      key.insert(key.end(), s.dress.begin(), s.dress.end());
      key.insert(key.end(), t.dress.begin(), t.dress.end());
      poly_add_term(r, s.coeff * t.coeff, std::move(key));
    }
  return r;
}

/// @brief Generator-trace product for a pinned adjoint assignment, with per-component FUNDAMENTAL
///        dressings inserted on the cycles. The dressed analogue of @ref loop_prod.
///
/// A `diagFund` sitting on a generator cycle weights one summed fundamental index per component, so
/// the trace is no longer a plain matrix product: its value depends on which component that index
/// takes. Pinning a position to a value `v` turns the diagonal factor into `D_v · e_v e_vᵀ`, which
/// splits the cyclic product into runs between consecutive pinned positions:
///
///     tr(M_0 … M_{m-1} with pins at p_1…p_r)
///        = Σ_{v_1…v_r} ∏_j (M_{p_j+1} … M_{p_{j+1}})_{v_j v_{j+1}} · ∏_j D^{(c_{p_j})}_{v_j}
///
/// so a cycle carrying `r` distinct dressed classes costs `N^r` terms — N=3 and r≤2 in the quark
/// flows this exists for. Positions whose classes were united by the union-find share one value, so
/// the enumeration runs over DISTINCT classes, not positions. Cycles with no dressed class take the
/// dense path and reproduce @ref loop_prod exactly.
inline SUNPoly loop_poly_dressed(const SUNDyn &dat, int N, const std::vector<std::vector<int>> &cycles,
                                 const std::vector<std::vector<int>> &cyclesFund,
                                 const std::map<int, int> &classVal,
                                 const std::map<int, std::vector<const std::vector<int> *>> &fundDiag) {
  SUNPoly prod{SUNTerm{Cx{1.0, 0.0}, {}}};
  for (std::size_t ci = 0; ci < cycles.size(); ++ci) {
    const std::vector<int> &cyc = cycles[ci], &fnd = cyclesFund[ci];
    const std::size_t m = cyc.size();
    std::vector<std::size_t> pinned; // positions carrying a fundamental dressing
    for (std::size_t k = 0; k < m; ++k)
      if (fundDiag.count(fnd[k])) pinned.push_back(k);

    if (pinned.empty()) { // undressed cycle: dense trace, identical to loop_prod
      DynMat mm = dat.gens[classVal.at(cyc[0])];
      for (std::size_t k = 1; k < m; ++k) mm = dmatmul(mm, dat.gens[classVal.at(cyc[k])]);
      const std::complex<double> tr = dtrace(mm);
      prod = poly_mul(prod, SUNPoly{SUNTerm{Cx{tr.real(), tr.imag()}, {}}});
      continue;
    }

    // run j = M_{p_j+1} … M_{p_{j+1}} (cyclic; a single matrix when the pins are adjacent, and the
    // whole cycle when there is only one pin).
    const std::size_t r = pinned.size();
    std::vector<DynMat> runs;
    runs.reserve(r);
    for (std::size_t j = 0; j < r; ++j) {
      const std::size_t from = (pinned[j] + 1) % m, to = pinned[(j + 1) % r];
      DynMat mm(N);
      for (int i = 0; i < N; ++i) mm(i, i) = {1.0, 0.0};
      for (std::size_t s = 0; s < m; ++s) {
        const std::size_t idx = (from + s) % m;
        mm = dmatmul(mm, dat.gens[classVal.at(cyc[idx])]);
        if (idx == to) break;
      }
      runs.push_back(std::move(mm));
    }

    // distinct dressed classes on this cycle; positions sharing a class share a value
    std::vector<int> dclasses;
    for (std::size_t j = 0; j < r; ++j)
      if (std::find(dclasses.begin(), dclasses.end(), fnd[pinned[j]]) == dclasses.end())
        dclasses.push_back(fnd[pinned[j]]);
    const std::size_t nd = dclasses.size();
    auto slotOf = [&](std::size_t j) {
      for (std::size_t d = 0; d < nd; ++d)
        if (dclasses[d] == fnd[pinned[j]]) return d;
      return nd; // unreachable: every pinned position's class is in dclasses
    };

    long long tot = 1;
    for (std::size_t d = 0; d < nd; ++d) tot *= N;
    SUNPoly cyclePoly;
    std::vector<int> val(nd, 0);
    for (long long t = 0; t < tot; ++t) {
      long long x = t;
      for (std::size_t d = 0; d < nd; ++d) { val[d] = static_cast<int>(x % N); x /= N; }
      std::vector<int> key; // every diag factor on every dressed class of this cycle
      bool drop = false;
      for (std::size_t d = 0; d < nd && !drop; ++d)
        for (const std::vector<int> *mp : fundDiag.at(dclasses[d])) {
          const int dr = (*mp)[val[d]];
          if (dr < 0) { drop = true; break; } // a dropped component kills this assignment
          key.push_back(dr);
        }
      if (drop) continue;
      Cx c{1.0, 0.0};
      for (std::size_t j = 0; j < r; ++j) {
        const std::complex<double> e = runs[j](val[slotOf(j)], val[slotOf((j + 1) % r)]);
        c = c * Cx{e.real(), e.imag()};
      }
      if (c.re == 0.0 && c.im == 0.0) continue;
      poly_add_term(cyclePoly, c, std::move(key));
    }
    prod = poly_mul(prod, cyclePoly);
  }
  return prod;
}

/// @brief Dressed single-group contraction: like @ref contract_group but folds **group-diagonal
///        dressing** factors (kinds 4/5) into a @ref SUNPoly `Σ_a c_a D_a` instead of one number.
///
/// A `diag_adj`/`diag_fund` factor is a δ that also pins a per-component dressing map (@ref
/// SUNFac::comp2dr, component → dressing-id, `-1` = drop) on its (closed) index class: a diag-dressed
/// *closed loop* becomes `Σ_a Π D_a` over the components its factors keep (replacing `δ^{aa}=N²−1` /
/// `δ^{ii}=N`); a diag-dressed *adjoint index that is summed/pinned* by the `f`/generator algebra
/// attaches each factor's `comp2dr[value]` to that assignment's contribution (and drops the assignment
/// when any factor drops that value). The plain (kinds 0–3) algebra is identical to @ref contract_group,
/// so a net with no diag factor yields a single constant term equal to `contract_group` (the @ref
/// sun_value_dressed gate routes such nets to the fast path).
///
/// A per-component *fundamental* dressing on a **generator** line (a `diag_fund` whose index is also
/// a generator row/col) is handled by @ref loop_poly_dressed, which expands the affected generator
/// trace over that index's components instead of matrix-multiplying through it. This is the quark
/// use case: a colour-diagonal propagator always sits between two `T^a` vertices.
inline SUNPoly contract_group_dressed(int N, const std::vector<const SUNFac *> &net) {
  assert_no_open_labels(N, net); // before the union-find, as in contract_group
  const SUNDyn &dat = sun_data_for(N);
  const int Adim = N * N - 1;

  // ---- union-find: a delta OR a diagonal-dressing identifies its two indices ----
  int maxlbl = -1;
  for (const SUNFac *f : net) {
    maxlbl = std::max({maxlbl, f->a, f->b});
    if (f->kind == SUNFacKind::F || f->kind == SUNFacKind::T) maxlbl = std::max(maxlbl, f->c);
  }
  // empty net: no factor ⇒ constant-1 polynomial (avoids UnionFind(-1)). In practice unreachable —
  // sun_value_dressed routes non-diagonal nets to sun_value_cx — but kept symmetric with contract_group.
  if (maxlbl < 0) return SUNPoly{SUNTerm{Cx{1.0, 0.0}, {}}};
  UnionFind uf(maxlbl);
  auto find = [&](int x) { return uf.find(x); };
  for (const SUNFac *f : net)
    if (f->kind == SUNFacKind::DeltaAdj || f->kind == SUNFacKind::DeltaFund ||
        f->kind == SUNFacKind::DiagFund || f->kind == SUNFacKind::DiagAdj)
      uf.unite(f->a, f->b);

  // ---- classify, recording the per-component dressings sitting on each class ----
  std::set<int> adjClasses, fClasses, genAdjClasses, fundClasses, genFundClasses;
  std::vector<std::array<int, 3>> fTriples, gens;
  // class -> the per-component dressing maps of every diag factor sitting on it. Each map is
  // component → dressing-id (`-1` = drop that component); the class value is the product over factors.
  std::map<int, std::vector<const std::vector<int> *>> adjDiag, fundDiag;
  for (const SUNFac *f : net) {
    switch (f->kind) {
    case SUNFacKind::DeltaAdj:
      adjClasses.insert(find(f->a)); adjClasses.insert(find(f->b));
      break;
    case SUNFacKind::F: {
      const int a = find(f->a), b = find(f->b), c = find(f->c);
      adjClasses.insert(a); adjClasses.insert(b); adjClasses.insert(c);
      fClasses.insert(a); fClasses.insert(b); fClasses.insert(c);
      fTriples.push_back({a, b, c});
      break;
    }
    case SUNFacKind::T: {
      const int a = find(f->a), i = find(f->b), j = find(f->c);
      adjClasses.insert(a); genAdjClasses.insert(a);
      fundClasses.insert(i); fundClasses.insert(j);
      genFundClasses.insert(i); genFundClasses.insert(j);
      gens.push_back({a, i, j});
      break;
    }
    case SUNFacKind::DeltaFund:
      fundClasses.insert(find(f->a)); fundClasses.insert(find(f->b));
      break;
    case SUNFacKind::DiagFund: { // diag_fund: δ^{ij} with a per-component fundamental dressing
      const int i = find(f->a);
      fundClasses.insert(i); fundClasses.insert(find(f->b));
      fundDiag[i].push_back(&f->comp2dr);
      break;
    }
    case SUNFacKind::DiagAdj: { // diag_adj: δ^{ab} with a per-component adjoint dressing
      const int a = find(f->a);
      adjClasses.insert(a); adjClasses.insert(find(f->b));
      adjDiag[a].push_back(&f->comp2dr);
      break;
    }
    default: NT_THROW(std::runtime_error, "sun_net: unknown SUNFac kind");
    }
  }
  // A per-component fundamental dressing on a generator line is NOT handled here: it stays on the
  // cycle and is folded by loop_poly_dressed below. The closed-loop pass just has to leave it alone
  // (it already does — it skips genFundClasses), so there is nothing to reject.

  // ---- closed-loop factors: plain loops fold to a scalar; diag-dressed loops to a Σ_a SUNPoly ----
  double factorScalar = 1.0;
  SUNPoly closedPoly{SUNTerm{Cx{1.0, 0.0}, {}}};
  // A diag-dressed closed loop `Σ_v Π_factor D^{factor}_v`. Component `v` contributes only if *every*
  // factor on the loop keeps it (`comp2dr[v] != -1`); a component any factor drops vanishes entirely.
  auto diagLoop = [](int dim, const std::vector<const std::vector<int> *> &factors) {
    SUNPoly fac;
    for (int v = 0; v < dim; ++v) {
      std::vector<int> key;
      bool drop = false;
      for (const std::vector<int> *m : factors) {
        const int dr = (*m)[v];
        if (dr < 0) { drop = true; break; }
        key.push_back(dr);
      }
      if (drop) continue;
      poly_add_term(fac, Cx{1.0, 0.0}, std::move(key));
    }
    return fac;
  };
  for (int c : adjClasses)
    if (!fClasses.count(c) && !genAdjClasses.count(c)) {
      auto it = adjDiag.find(c);
      if (it == adjDiag.end()) factorScalar *= static_cast<double>(Adim);
      else closedPoly = poly_mul(closedPoly, diagLoop(Adim, it->second));
    }
  for (int c : fundClasses)
    if (!genFundClasses.count(c)) {
      auto it = fundDiag.find(c);
      if (it == fundDiag.end()) factorScalar *= static_cast<double>(N);
      else closedPoly = poly_mul(closedPoly, diagLoop(N, it->second));
    }

  // generator adjoint classes not pinned by an f are summed densely; those also diag-dressed tag.
  std::vector<int> genOnly;
  for (int c : genAdjClasses)
    if (!fClasses.count(c)) genOnly.push_back(c);
  std::map<int, std::vector<const std::vector<int> *>> asgDiag; // diag-dressed adjoint classes pinned in the assignment
  for (const auto &kv : adjDiag)
    if (fClasses.count(kv.first) || genAdjClasses.count(kv.first)) asgDiag[kv.first] = kv.second;

  // ---- fundamental-cycle extraction (generator traces) — shared with contract_group ----
  const std::vector<std::vector<int>> cycles = extract_cycles_adj(gens);
  // the fundamental classes riding those same cycles, needed only when a diagFund sits on one
  const std::vector<std::vector<int>> cyclesFund = extract_cycles_fund(gens);

  // ---- assignment sum (sparse f-backtracking + dense gen-only), tagging diag-dressed values ----
  SUNPoly total; // 0
  std::map<int, int> classVal; // adjoint class -> its pinned/summed component value
  auto emit = [&](Cx fProd) {
    std::vector<int> key;
    for (const auto &kv : asgDiag) {
      const int val = classVal.at(kv.first);
      for (const std::vector<int> *m : kv.second) {
        const int dr = (*m)[val];
        if (dr < 0) return; // this assignment lands on a dropped component ⇒ contributes 0
        key.push_back(dr);
      }
    }
    // The generator traces are a POLYNOMIAL now, not a scalar: a diagFund on a cycle tags each
    // fundamental component with its own dressing id. Every cycle term multiplies the adjoint-side
    // key built above.
    const SUNPoly lp = loop_poly_dressed(dat, N, cycles, cyclesFund, classVal, fundDiag);
    for (const auto &t : lp) {
      std::vector<int> k2 = key;
      k2.insert(k2.end(), t.dress.begin(), t.dress.end());
      poly_add_term(total, fProd * t.coeff, std::move(k2));
    }
  };
  auto sumGen = [&](auto &&self, std::size_t gi, Cx fProd) -> void {
    if (gi == genOnly.size()) { emit(fProd); return; }
    const int cls = genOnly[gi];
    for (int v = 0; v < Adim; ++v) { classVal[cls] = v; self(self, gi + 1, fProd); }
    classVal.erase(cls);
  };
  auto sumF = [&](auto &&self, std::size_t fi, Cx fProd) -> void {
    if (fi == fTriples.size()) { sumGen(sumGen, 0, fProd); return; }
    const auto [ca, cb, cc] = fTriples[fi];
    for (const auto &e : dat.f_nz) {
      auto consistent = [&](int cl, int v) { auto it = classVal.find(cl); return it == classVal.end() || it->second == v; };
      if (!consistent(ca, e.a) || !consistent(cb, e.b) || !consistent(cc, e.c)) continue;
      const bool newA = !classVal.count(ca), newB = !classVal.count(cb), newC = !classVal.count(cc);
      classVal[ca] = e.a; classVal[cb] = e.b; classVal[cc] = e.c;
      self(self, fi + 1, fProd * Cx{e.v, 0.0});
      if (newA) classVal.erase(ca);
      if (newB) classVal.erase(cb);
      if (newC) classVal.erase(cc);
    }
  };
  sumF(sumF, 0, Cx{1.0, 0.0});

  // ---- combine: (assignment sum) × (closed-loop diag polys) × factorScalar; snap residual zeros ----
  SUNPoly r = poly_mul(total, closedPoly);
  constexpr double kZeroSnapTol = 1e-9;
  for (auto &t : r) {
    t.coeff = t.coeff * Cx{factorScalar, 0.0};
    if (std::fabs(t.coeff.re) < kZeroSnapTol) t.coeff.re = 0.0;
    if (std::fabs(t.coeff.im) < kZeroSnapTol) t.coeff.im = 0.0;
  }
  r.erase(std::remove_if(r.begin(), r.end(), [](const SUNTerm &t) { return t.coeff.re == 0.0 && t.coeff.im == 0.0; }),
          r.end());
  return r;
}

} // namespace sun_net_detail
#endif // NUMTRACER_DEFINE_BODIES

// Public SU(N) entry points: declared always, defined once (library TU / header-only build).
NUMTRACER_FUNC Cx sun_value_cx(const SUNNet &net);
NUMTRACER_FUNC double sun_value(const SUNNet &net);
NUMTRACER_FUNC SUNPoly sun_value_dressed(const SUNNet &net);

#if NUMTRACER_DEFINE_BODIES
/// @brief Contract a (possibly two-group) SU(N) colour/flavour network to its complex scalar value.
///
/// The factors are partitioned by their group rank `g`; each group has a disjoint label space, so
/// the value is the product of the per-group contractions. An empty network is the identity (value 1).
/// Works for any rank `g ≥ 1` (@ref sun_net_detail::contract_group builds the SU(g) data on demand).
NUMTRACER_FUNC Cx sun_value_cx(const SUNNet &net) {
  std::set<int> groups;
  for (const SUNFac &f : net) groups.insert(f.g);
  Cx r{1.0, 0.0};
  for (int g : groups) {
    std::vector<const SUNFac *> sub;
    for (const SUNFac &f : net)
      if (f.g == g) sub.push_back(&f);
    r = r * sun_net_detail::contract_group(g, sub);
  }
  return r;
}

/// @brief Contract an SU(N) network to its real scalar value (backward-compatible entry point).
///
/// Kept returning `double` so existing generator call sites are untouched; adjoint SU(3) colour
/// nets are real, so taking the real part is exact.
NUMTRACER_FUNC double sun_value(const SUNNet &net) { return sun_value_cx(net).re; }

/// @brief Contract a (possibly two-group) SU(N) network carrying **group-diagonal dressings** to a
///        @ref SUNPoly — `Σ_t coeff_t · Π D^{dr}`.
///
/// If the net carries no kind-4/5 (diagonal-dressing) factor this returns the single constant term
/// `{sun_value_cx(net), {}}` — *byte-identical* to the undressed fold, so existing flows are
/// unaffected. Otherwise each group is contracted with @ref sun_net_detail::contract_group_dressed
/// and the per-group polynomials are multiplied (disjoint label spaces, so the value factorises).
NUMTRACER_FUNC SUNPoly sun_value_dressed(const SUNNet &net) {
  bool hasDiag = false;
  for (const SUNFac &f : net)
    if (f.kind == SUNFacKind::DiagFund || f.kind == SUNFacKind::DiagAdj) { hasDiag = true; break; }
  if (!hasDiag) return SUNPoly{SUNTerm{sun_value_cx(net), {}}};

  std::set<int> groups;
  for (const SUNFac &f : net) groups.insert(f.g);
  SUNPoly r{SUNTerm{Cx{1.0, 0.0}, {}}};
  for (int g : groups) {
    std::vector<const SUNFac *> sub;
    for (const SUNFac &f : net)
      if (f.g == g) sub.push_back(&f);
    r = sun_net_detail::poly_mul(r, sun_net_detail::contract_group_dressed(g, sub));
  }
  return r;
}
#endif // NUMTRACER_DEFINE_BODIES

} // namespace numtracer::network
