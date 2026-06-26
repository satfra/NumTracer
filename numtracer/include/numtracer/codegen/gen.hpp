/// @file gen.hpp
/// @brief Build-time generator support for the numeric kernels: lower several diagrams
///        into a **shared** fundamental-symbol environment, and emit each as a straight-line C++
///        function. Runtime-only (uses `<ostream>`) — this is the codegen side, not the kernel.
///
/// Contracting a large trace at compile time is RAM/time-prohibitive in GCC (it retains every
/// intermediate), whereas the same contraction runs in a fraction of the memory and time at
/// runtime. So the contraction is a codegen step: a small C++ program builds the network (the
/// contraction stays in C++), contracts each diagram numerically → lowers it (see
/// @ref numtracer::numeric::to_genprog), and **prints** a committed header of straight-line
/// `trN(const double* f)` functions — the same flat form as a FORM reference kernel. The kernel
/// computes the shared `f[]` once per call and invokes them.
#pragma once

#include "numtracer/codegen/lower.hpp"
#include "numtracer/network/network.hpp" // splitmix64_finalise / hash_combine hash helpers + <cstdint>

#include <climits>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

namespace numtracer::network
{

  /// @brief A shared fundamental-symbol environment: assigns one global env id (`f[]` index) per
  ///        distinct symbol (a scalar product `sp(a,b)` or an inverse `inv(id)`), across all the
  ///        diagrams of a kernel — so the ~few fundamental symbols are computed once per call.
  struct GlobalEnv {
    std::vector<std::tuple<int, int, int>> syms; ///< env id i → (kind 0=sp/1=inv, a, b); inv is (1,id,0).
    std::vector<int> bucket;                     ///< open-addressed index into @ref syms (or -1).
    std::size_t mask = 0;                        ///< bucket.size()-1 (power of two); 0 while empty.

    static std::uint64_t hsym(int k, int a, int b)
    {
      return hash_combine(
          hash_combine(hash_combine(splitmix64_finalise(2), static_cast<std::uint64_t>(static_cast<unsigned>(k))),
                       static_cast<std::uint64_t>(static_cast<unsigned>(a))),
          static_cast<std::uint64_t>(static_cast<unsigned>(b)));
    }
    void rehash(std::size_t cap)
    {
      bucket.assign(cap, -1);
      mask = cap - 1;
      for (std::size_t s = 0; s < syms.size(); ++s) {
        auto [k, a, b] = syms[s];
        std::size_t p = hsym(k, a, b) & mask;
        while (bucket[p] != -1)
          p = (p + 1) & mask;
        bucket[p] = static_cast<int>(s);
      }
    }
    /// First-seen lookup of symbol `(k,a,b)`; appends on miss so env ids stay in first-seen order.
    int intern(int k, int a, int b)
    {
      const std::uint64_t h = hsym(k, a, b);
      if (mask) {
        std::size_t p = h & mask;
        while (bucket[p] != -1) {
          if (syms[bucket[p]] == std::tuple<int, int, int>{k, a, b}) return bucket[p];
          p = (p + 1) & mask;
        }
      }
      if ((syms.size() + 1) * 10 >= (mask + 1) * 7) rehash(mask == 0 ? 16 : (mask + 1) * 2);
      const int s = static_cast<int>(syms.size());
      syms.push_back({k, a, b});
      std::size_t p = h & mask;
      while (bucket[p] != -1)
        p = (p + 1) & mask;
      bucket[p] = s;
      return s;
    }
    int inv_id(int v) { return intern(1, v, 0); }
    /// A raw USER-SYMBOL leaf (kind 3): a kernel argument (a momentum component / angle) the numeric
    /// backend interns directly. It fills its `f[]` slot from the argument verbatim (see
    /// @ref FillFormulas::var) and rides the polynomial as a monomial variable — the inv backend never
    /// emits kind 3, so its env layout / fill are byte-identical.
    int var_id(int v) { return intern(3, v, 0); }
    /// A DRESSING symbol (kind 2): an opaque runtime leaf (a propagator dressing / regulator call)
    /// the kernel evaluates once and stores in `f[]`. It enters the polynomial only as a monomial
    /// factor, so it lowers exactly like an `inv` symbol — the difference is purely how the kernel
    /// fills its slot (a dressing C++ expression vs `1/q²`); see @ref FillFormulas::dress.
    int dr_id(int id) { return intern(2, id, 0); }
  };

  /// @brief Sentinel for @ref GenProg::rootIm: the program has no imaginary part (it is purely real).
  inline constexpr int kRealProgram = INT_MIN;

  /// @brief A lowered diagram program over the shared env (its `RVAR` ids are global `f[]` indices).
  ///
  /// `rootIm == kRealProgram` marks a purely REAL program (`root` is the result slot; emitted as
  /// `double`). Otherwise the program is COMPLEX: the real and imaginary parts share one instruction
  /// stream `ins`, `root` is the real-part slot and `rootIm` the imaginary-part slot (either `-1` if
  /// that part is structurally zero); emitted as `std::complex<double>`. A complex `tr_i` arises when a
  /// complex colour factor (an imaginary non-abelian-vertex T-trace) was folded into the polynomial —
  /// the kernel multiplies it by the (also complex) runtime dressing coefficient and the consumer takes
  /// Re, so the imaginary part must survive to the kernel rather than being dropped at this lowering.
  struct GenProg {
    std::vector<expr::RInstr> ins;
    int root = -1;
    int rootIm = kRealProgram;
  };

  /// @brief Lower a polynomial's monomials into the shared env via CSE (`RBuilder`) + Horner.
  ///
  /// The input is first sorted into a **canonical order** so the emitted program (env-id layout, Horner
  /// pivot choices, hence the operation count) depends only on the polynomial as a *set*, not on the
  /// order the upstream reduction happened to produce it. This makes the kernel reproducible across the
  /// serial and (term-reordering) sharded-parallel rebase paths — the same op count either way.
  namespace gdetail
  {
    /// The deterministic monomial orderings the greedy (order-sensitive) Horner is tried on; the caller
    /// keeps the factorization with the fewest emitted ops — reproducible and never worse than any single
    /// order.
    inline std::vector<std::vector<LMono>> make_orderings(const std::vector<LMono> &monos)
    {
      auto compareAscending = [](const LMono &x, const LMono &y) { return x.vp < y.vp; };
      auto degree = [](const LMono &m) {
        int d = 0;
        for (auto [id, e] : m.vp)
          d += e;
        return d;
      };
      std::vector<std::vector<LMono>> orders;
      orders.push_back(monos); // as-built (canonical)
      {
        auto v = monos;
        std::sort(v.begin(), v.end(), compareAscending);
        orders.push_back(std::move(v));
      }
      {
        auto v = monos;
        std::sort(v.begin(), v.end(), compareAscending);
        std::reverse(v.begin(), v.end());
        orders.push_back(std::move(v));
      }
      {
        auto v = monos;
        std::sort(v.begin(), v.end(), [](const LMono &x, const LMono &y) {
          if (x.vp.size() != y.vp.size()) return x.vp.size() > y.vp.size(); // most-factors first
          return x.vp < y.vp;
        });
        orders.push_back(std::move(v));
      }
      {
        auto v = monos;
        std::sort(v.begin(), v.end(), [](const LMono &x, const LMono &y) {
          if (x.vp.size() != y.vp.size()) return x.vp.size() < y.vp.size(); // fewest-factors first
          return x.vp < y.vp;
        });
        orders.push_back(std::move(v));
      }
      {
        auto v = monos;
        std::sort(v.begin(), v.end(), [&](const LMono &x, const LMono &y) {
          const int dx = degree(x), dy = degree(y);
          if (dx != dy) return dx > dy; // highest total degree first
          return x.vp < y.vp;
        });
        orders.push_back(std::move(v));
      }
      {
        auto v = monos;
        std::sort(v.begin(), v.end(), [&](const LMono &x, const LMono &y) {
          const int dx = degree(x), dy = degree(y);
          if (dx != dy) return dx < dy; // lowest total degree first
          return x.vp < y.vp;
        });
        orders.push_back(std::move(v));
      }
      {
        auto v = monos;
        std::sort(v.begin(), v.end(), [](const LMono &x, const LMono &y) { return x.vp > y.vp; });
        orders.push_back(std::move(v));
      } // reverse-lex
      return orders;
    }
  } // namespace gdetail

  namespace gdetail
  {
    /// Pick the cheapest Horner ordering of @p monos (costed on scratch builders), then replay it into
    /// @p w (so several parts — e.g. a trace's real and imaginary halves — share one CSE stream). Returns
    /// the result slot in @p w.
    inline int best_into(const std::vector<LMono> &monos, expr::rdetail::RBuilder &w)
    {
      // The greedy Horner is order-sensitive, so we cost several deterministic orderings and keep the
      // cheapest. But each trial is a FULL horner pass: for the dense 1/4/7 traces (tens of thousands of
      // monomials) the 8-way sweep dominates GENERATION while changing the op count only ~1% (the op count
      // tracks the monomial count, not the pivot order — see the noise-prune notes). So scale the sweep
      // down with the polynomial size; NT_GEN_HORNER_ORDERS=<n> overrides.
      std::size_t numOrderings = 8;
      if (const char *e = std::getenv("NT_GEN_HORNER_ORDERS")) {
        int v = std::atoi(e);
        if (v > 0) numOrderings = (std::size_t)v;
      } else if (monos.size() > 2000)
        numOrderings = 1;
      else if (monos.size() > 500)
        numOrderings = 3;
      if (numOrderings <= 1) return horner(w, monos); // canonical (as-built) order only — no sweep
      auto orders = make_orderings(monos);
      if (numOrderings > orders.size()) numOrderings = orders.size();
      std::size_t bestIdx = 0, bestOps = 0;
      bool have = false;
      for (std::size_t i = 0; i < numOrderings; ++i) {
        expr::rdetail::RBuilder scratch; // cost this ordering on a throwaway builder
        horner(scratch, orders[i]);
        if (!have || scratch.ins.size() < bestOps) {
          bestOps = scratch.ins.size();
          bestIdx = i;
          have = true;
        }
      }
      return horner(w, std::move(orders[bestIdx]));
    }
  } // namespace gdetail

  /// @brief Print a lowered program as a straight-line C++ function `name(const double* f)`.
  ///        `decor` is the function decorator/prefix (e.g. `"static __host__ __device__ inline"`
  ///        for a CUDA-callable kernel); the default keeps the emitted bytes unchanged.
  inline void emit_cpp(std::ostream &out, const GenProg &p, const std::string &name,
                       const std::string &decor = "static inline")
  {
    // Compile-time escape hatch (NT_GEN_NOINLINE_TRACES): force the per-diagram trace functions
    // OUT-OF-LINE. The default all-inlined emission is the runtime-optimal one — the register
    // spill it can cause once many large traces fuse into one kernel (≈four-gluon scale) is benign
    // (the loop is fp64-pipe-bound, not occupancy/memory-bound, and inlining keeps cross-trace CSE),
    // so out-of-lining is measurably SLOWER (≈7–19%). But a 100k+-line all-inlined kernel is
    // expensive to compile: out-of-lining the traces roughly halves the nvcc compile time and RAM
    // of the 147k-line ZA4_147 kernel. `__attribute__((noinline))` is honoured by both g++ (host)
    // and nvcc (device); fill()/powr stay inline. See tests/gpu/README.md for the measurements.
    std::string effDecor = decor;
    if (std::getenv("NT_GEN_NOINLINE_TRACES")) {
      const std::string kw = " inline";
      if (effDecor.size() >= kw.size() && effDecor.compare(effDecor.size() - kw.size(), kw.size(), kw) == 0)
        effDecor.replace(effDecor.size() - kw.size(), kw.size(), " __attribute__((noinline))");
      else
        effDecor += " __attribute__((noinline))";
    }
    // COMPLEX trace (folded imaginary colour): real + imaginary halves share the instruction stream;
    // return std::complex<double>{re, im}. The kernel multiplies it by the (complex) dressing
    // coefficient and the consumer takes std::real — so the imaginary part reaches the kernel.
    if (p.rootIm != kRealProgram) {
      auto slot = [](int r) { return r < 0 ? std::string("0.0") : "s" + std::to_string(r); };
      out << effDecor << " std::complex<double> " << name << "(const double *f) {\n";
      out << std::setprecision(17);
      for (std::size_t i = 0; i < p.ins.size(); ++i) {
        const expr::RInstr &in = p.ins[i];
        out << "  const double s" << i << " = ";
        switch (in.op) {
        case expr::RCONST:
          out << in.k;
          break;
        case expr::RVAR:
          out << "f[" << in.a << "]";
          break;
        case expr::RADD:
          out << "s" << in.a << "+s" << in.b;
          break;
        case expr::RSUB:
          out << "s" << in.a << "-s" << in.b;
          break;
        case expr::RMUL:
          out << "s" << in.a << "*s" << in.b;
          break;
        default:
          out << "-s" << in.a;
          break; // RNEG
        }
        out << ";\n";
      }
      out << "  return std::complex<double>{" << slot(p.root) << ", " << slot(p.rootIm) << "};\n}\n";
      return;
    }
    out << effDecor << " double " << name << "(const double *f) {\n";
    if (p.root < 0) {
      out << "  return 0.0;\n}\n";
      return;
    }
    out << std::setprecision(17);
    for (std::size_t i = 0; i < p.ins.size(); ++i) {
      const expr::RInstr &in = p.ins[i];
      out << "  const double s" << i << " = ";
      switch (in.op) {
      case expr::RCONST:
        out << in.k;
        break;
      case expr::RVAR:
        out << "f[" << in.a << "]";
        break;
      case expr::RADD:
        out << "s" << in.a << "+s" << in.b;
        break;
      case expr::RSUB:
        out << "s" << in.a << "-s" << in.b;
        break;
      case expr::RMUL:
        out << "s" << in.a << "*s" << in.b;
        break;
      default:
        out << "-s" << in.a;
        break; // RNEG
      }
      out << ";\n";
    }
    out << "  return s" << p.root << ";\n}\n";
  }

  /// @brief Print the shared env layout as a comment (which `f[i]` is which symbol) so the codegen /
  ///        kernel knows the formula to fill each slot with.
  inline void emit_env_layout(std::ostream &out, const GlobalEnv &g)
  {
    out << "// fundamental-symbol env layout (fill f[i] per call):\n";
    for (std::size_t i = 0; i < g.syms.size(); ++i) {
      auto [kind, a, b] = g.syms[i];
      if (kind == 0)
        out << "//   f[" << i << "] = sp(" << a << "," << b << ")\n";
      else if (kind == 1)
        out << "//   f[" << i << "] = inv(" << a << ")\n";
      else if (kind == 2)
        out << "//   f[" << i << "] = dress(" << a << ")\n";
      else
        out << "//   f[" << i << "] = var(" << a << ")\n";
    }
  }

  /// @brief Formula providers for the fundamental symbols: given the symbol the reduction assigned to
  ///        an `f[]` slot, return the C++ expression that computes it from the kernel scalars.
  ///        `sp(a,b)` is a scalar product of fundamental momenta `a·b`; `inv(id)` is `1/q_id²`.
  struct FillFormulas {
    std::function<std::string(int /*a*/, int /*b*/)> sp; ///< C++ for the scalar product `sp(a,b)`.
    std::function<std::string(int /*invId*/)> inv;       ///< C++ for the inverse `inv(invId)`.
    std::function<std::string(int /*drId*/)> dress;      ///< C++ for the dressing/regulator call `dress(drId)`.
    std::function<std::string(int /*varId*/)> var;       ///< C++ for a raw user symbol (numeric backend, kind 3).
  };

  /// @brief Print a `fill(double* f, <args>)` that fills every shared `f[]` slot from the kernel
  ///        scalars, using the supplied formula providers. Emitting this **from the generator**
  ///        (which alone knows the reduced layout) decouples the kernel from the `f[]` ordering: the
  ///        kernel just calls `fill(...)` then the `trN(f)`. `argSig` is the parameter list (e.g.
  ///        `"double l1, double cos1, double cos2, double p"`).
  inline void emit_fill(std::ostream &out, const GlobalEnv &g, const std::string &name, const std::string &argSig,
                        const FillFormulas &fm, const std::string &decor = "static inline")
  {
    out << std::setprecision(17);
    out << decor << " void " << name << "(double *f, " << argSig << ") {\n";
    for (std::size_t i = 0; i < g.syms.size(); ++i) {
      auto [kind, a, b] = g.syms[i];
      out << "  f[" << i << "] = "
          << (kind == 0   ? fm.sp(a, b)
              : kind == 1 ? fm.inv(a)
              : kind == 2 ? fm.dress(a)
                          : fm.var(a))
          << ";\n";
    }
    out << "}\n";
  }

} // namespace numtracer::network
