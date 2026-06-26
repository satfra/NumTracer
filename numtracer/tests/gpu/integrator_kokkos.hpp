// NumTracer GPU tests — Kokkos twin of numtracer/cuda/integrator.cuh: the same two-phase
// map/reduce (phase 1: one parallel_for writes every weighted integrand value to a device
// View; phase 2: a pure per-grid-point team reduction, DiFfRG QuadratureIntegrator::map
// style) expressed in Kokkos instead of raw CUDA + cub. Built against the DiFfRG-bundled
// Kokkos (CUDA backend) for an apples-to-apples comparison with the raw-CUDA integrator:
// same generated __host__ __device__ kernels, same quadrature, same buffer layout.
//
// Differences from the CUDA twin, by construction:
//   * the integrand is passed as a trivially-copyable functor EVAL (operator()(l1,cos1,
//     cos2[,phi],p)) instead of a KERNEL class + dressing parameter pack — extended
//     __device__ lambdas and parameter packs are a fragile combination under nvcc;
//   * block size / launch configuration are chosen by Kokkos (RangePolicy occupancy
//     heuristics), not pinned to 128 threads with __launch_bounds__;
//   * phases are timed with fence-bounded std::chrono (Kokkos has no cudaEvent analogue).
#pragma once

#include <Kokkos_Core.hpp>

#include <chrono>
#include <vector>

namespace numtracer::kokkos {

// kernels return complex<double> or plain double depending on the flow's diagram content
KOKKOS_INLINE_FUNCTION double real_of(double x) { return x; }
template <class C> KOKKOS_INLINE_FUNCTION double real_of(const C &z) { return z.real(); }

// arity dispatch OUTSIDE the lambda: nvcc forbids first-capturing a variable of an extended
// __host__ __device__ lambda inside an if-constexpr context.
template <int NV, class EVAL>
KOKKOS_INLINE_FUNCTION double eval_at(const EVAL &ev, const double (&x)[4], double p) {
  if constexpr (NV == 3)
    return ev(x[0], x[1], x[2], p);
  else
    return ev(x[0], x[1], x[2], x[3], p);
}

struct MapResult {
  std::vector<double> integral;
  double phase1_ms = 0; ///< integrand evaluation (fence-bounded wall time)
  double phase2_ms = 0; ///< segmented team reduction (fence-bounded wall time)
  double wall_ms = 0;   ///< whole map() incl. H2D/D2H
  long long evals = 0;
};

/// @brief Two-phase quadrature integrator over a grid of external momenta. EVAL is a
///        trivially-copyable device-callable functor; NV its number of integration
///        variables (3 or 4). Views are grow-only across map() calls.
template <class EVAL, int NV> class SegmentedQuadIntegrator {
  using exec = Kokkos::DefaultExecutionSpace;
  using mem = typename exec::memory_space;
  using view = Kokkos::View<double *, mem>;

public:
  void set_axis(int d, const std::vector<double> &nodes, const std::vector<double> &weights) {
    nodes_[d] = to_device(nodes, "nt_nodes");
    weights_[d] = to_device(weights, "nt_weights");
    n_[d] = nodes.size();
  }

  MapResult map(const std::vector<double> &p_grid, double prefactor, const EVAL &ev) {
    const auto t0 = std::chrono::steady_clock::now();
    const int n_p = (int)p_grid.size();
    long long Q = 1;
    Kokkos::Array<view, NV> nd, wt;
    Kokkos::Array<long long, NV> nn;
    for (int d = 0; d < NV; ++d) {
      nd[d] = nodes_[d];
      wt[d] = weights_[d];
      nn[d] = n_[d];
      Q *= n_[d];
    }
    const long long total = Q * n_p;
    if ((size_t)total > buf_.extent(0)) buf_ = view("nt_buf", total);
    if ((size_t)n_p > out_.extent(0)) out_ = view("nt_out", n_p);
    p_ = to_device(p_grid, "nt_p");
    const view buf = buf_, out = out_, pv = p_;

    // phase 1: integrand to memory (flat index decode identical to the CUDA twin)
    Kokkos::fence();
    const auto t1 = std::chrono::steady_clock::now();
    // LaunchBounds<128>: same codegen pinning as the CUDA twin's __launch_bounds__(128) —
    // measured +14% on ZA4 over letting ptxas assume unbounded block sizes.
    Kokkos::parallel_for(
        "nt_phase1_kokkos",
        Kokkos::RangePolicy<exec, Kokkos::LaunchBounds<128>>(0, total),
        KOKKOS_LAMBDA(const long long idx) {
          long long r = idx;
          int i[NV];
          for (int d = NV - 1; d >= 0; --d) {
            i[d] = (int)(r % nn[d]);
            r /= nn[d];
          }
          const int ip = (int)r;
          double w = prefactor;
          double x[4] = {};
          for (int d = 0; d < NV; ++d) {
            w *= wt[d](i[d]);
            x[d] = nd[d](i[d]);
          }
          const double v = eval_at<NV>(ev, x, pv(ip));
          buf(idx) = w * x[0] * x[0] * x[0] * v;
        });
    Kokkos::fence();
    const auto t2 = std::chrono::steady_clock::now();

    // phase 2: pure team reduction, one team per grid point (DiFfRG map() style)
    Kokkos::parallel_for(
        "nt_phase2_kokkos", Kokkos::TeamPolicy<exec>(n_p, Kokkos::AUTO, 32),
        KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<exec>::member_type &team) {
          const long long k = team.league_rank();
          double r = 0;
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(team, Q),
              [&](const long long j, double &u) { u += buf(k * Q + j); }, r);
          Kokkos::single(Kokkos::PerTeam(team), [&]() { out(k) = r; });
        });
    Kokkos::fence();
    const auto t3 = std::chrono::steady_clock::now();

    MapResult res;
    res.integral.resize(n_p);
    auto h = Kokkos::create_mirror_view(out);
    Kokkos::deep_copy(h, out);
    for (int i = 0; i < n_p; ++i) res.integral[i] = h(i);
    res.phase1_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
    res.phase2_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();
    res.wall_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    res.evals = total;
    return res;
  }

  size_t buffer_bytes() const { return buf_.extent(0) * sizeof(double); }

private:
  static view to_device(const std::vector<double> &v, const char *name) {
    view d(name, v.size());
    auto h = Kokkos::create_mirror_view(d);
    for (size_t i = 0; i < v.size(); ++i) h(i) = v[i];
    Kokkos::deep_copy(d, h);
    return d;
  }

  view nodes_[NV], weights_[NV], buf_, out_, p_;
  long long n_[NV] = {};
};

} // namespace numtracer::kokkos
