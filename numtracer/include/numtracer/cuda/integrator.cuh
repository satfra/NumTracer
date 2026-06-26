/// @file integrator.cuh
/// @brief Two-phase CUDA quadrature integrator over a grid of external momenta, modeled on
///        DiFfRG's QuadratureIntegrator::map(): phase 1 evaluates the weighted integrand at
///        every (grid point, quadrature node) pair into one device buffer, phase 2 is a PURE
///        segmented reduction (one segment per grid point). The integrand is a generated
///        NumTracer kernel class whose `kernel`/`constant` are `__host__ __device__`
///        (generation option `"Decorator" -> "static __host__ __device__ inline"`).
///
/// Header-only and dependency-free beyond the CUDA toolkit (CUB ships with it): quadrature
/// nodes/weights are CONSUMED as host arrays (the tests obtain them from GSL), never computed
/// here. The cos1 axis is expected to carry Chebyshev-2 weights, i.e. the sqrt(1-cos1^2)
/// angular jacobian is baked into the weights and must NOT appear in the integrand — the
/// same convention as DiFfRG's Integrator_p2_4D_{2,3}ang.
#pragma once

#include <cub/cub.cuh>
#include <cuda_runtime.h>

#include <chrono>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace numtracer::cuda {

/// @brief Abort with file/line on any CUDA error (a test utility — failing loudly is correct).
#define NT_CUDA_CHECK(call)                                                                        \
  do {                                                                                             \
    const cudaError_t nt_err_ = (call);                                                            \
    if (nt_err_ != cudaSuccess) {                                                                  \
      std::fprintf(stderr, "CUDA error %s at %s:%d: %s\n", cudaGetErrorName(nt_err_), __FILE__,    \
                   __LINE__, cudaGetErrorString(nt_err_));                                         \
      std::exit(2);                                                                                \
    }                                                                                              \
  } while (0)

__host__ __device__ inline double real_of(double x) { return x; }
// generic: any complex-like type with .real() (std::complex on host, cuda::std::complex in CUDA TUs)
template <class C> __host__ __device__ inline double real_of(const C &z) { return z.real(); }

/// @brief Device view of the NV quadrature axes (pointers into device memory).
template <int NV> struct QuadDev {
  const double *nodes[NV];
  const double *weights[NV];
  int n[NV];
};

/// @brief Phase-1 CUDA thread-block size (tuned for warp efficiency / occupancy).
inline constexpr int kPhase1BlockSize = 128;

/// @brief Phase 1: one thread per (grid point, quadrature node) pair, grid-stride.
///        Mixed-radix decode of the flat index (innermost axis fastest, grid point slowest, so
///        each grid point owns one CONTIGUOUS segment of `out` for the segmented reduction).
///        Writes  prefactor * l1^3 * (prod_d w_d) * Re kernel(l1, cos1, cos2[, phi], p, k, ...).
template <class KERNEL, int NV, class... Dress>
__global__ void __launch_bounds__(kPhase1BlockSize) nt_phase1(QuadDev<NV> quad_dev, const double *p_grid, long long total,
                                                 double k, double prefactor, double *out, Dress... dress) {
  for (long long idx = blockIdx.x * (long long)blockDim.x + threadIdx.x; idx < total;
       idx += (long long)gridDim.x * blockDim.x) {
    long long r = idx;
    int i[NV];
#pragma unroll
    for (int d = NV - 1; d >= 0; --d) {
      i[d] = (int)(r % quad_dev.n[d]);
      r /= quad_dev.n[d];
    }
    const int ip = (int)r;
    double w = prefactor;
#pragma unroll
    for (int d = 0; d < NV; ++d) w *= quad_dev.weights[d][i[d]];
    const double l1 = quad_dev.nodes[0][i[0]];
    double v;
    if constexpr (NV == 3)
      v = real_of(KERNEL::kernel(l1, quad_dev.nodes[1][i[1]], quad_dev.nodes[2][i[2]], p_grid[ip], k, dress...));
    else
      v = real_of(KERNEL::kernel(l1, quad_dev.nodes[1][i[1]], quad_dev.nodes[2][i[2]], quad_dev.nodes[3][i[3]], p_grid[ip], k, dress...));
    out[idx] = w * l1 * l1 * l1 * v; // l1^3: the 4D radial integration element
  }
}

/// @brief Result of one map(): per-grid-point integrals + phase timings.
struct MapResult {
  std::vector<double> integral; ///< one integral per external-momentum grid point
  float phase1_ms = 0;          ///< integrand evaluation (cudaEvent)
  float phase2_ms = 0;          ///< segmented reduction (cudaEvent)
  double wall_ms = 0;           ///< whole map() incl. H2D/D2H (std::chrono)
  long long evals = 0;          ///< kernel evaluations = n_p * prod(orders)
};

/// @brief The two-phase integrator. KERNEL is a generated NumTracer kernel class
///        (e.g. `numtracer_kernels::my_gpu_kernel<MyRegulator>`), NV its number of
///        integration variables (3: l1,cos1,cos2 — 4: +phi). All device allocations are
///        grow-only so repeated map() calls (timing loops) do not reallocate.
template <class KERNEL, int NV> class SegmentedQuadIntegrator {
public:
  SegmentedQuadIntegrator() {
    for (int d = 0; d < NV; ++d) {
      d_nodes_[d] = d_weights_[d] = nullptr;
      axis_sizes_[d] = 0;
    }
    NT_CUDA_CHECK(cudaEventCreate(&e0_));
    NT_CUDA_CHECK(cudaEventCreate(&e1_));
    NT_CUDA_CHECK(cudaEventCreate(&e2_));
  }
  ~SegmentedQuadIntegrator() {
    for (int d = 0; d < NV; ++d) {
      cudaFree(d_nodes_[d]);
      cudaFree(d_weights_[d]);
    }
    cudaFree(d_buf_);
    cudaFree(d_offsets_);
    cudaFree(d_out_);
    cudaFree(d_p_);
    cudaFree(d_tmp_);
    cudaEventDestroy(e0_);
    cudaEventDestroy(e1_);
    cudaEventDestroy(e2_);
  }
  SegmentedQuadIntegrator(const SegmentedQuadIntegrator &) = delete;
  SegmentedQuadIntegrator &operator=(const SegmentedQuadIntegrator &) = delete;

  /// @brief Upload one quadrature axis (host nodes/weights, already mapped to the target
  ///        interval — GSL's fixed tables include the interval scaling in the weights).
  void set_axis(int d, const std::vector<double> &nodes, const std::vector<double> &weights) {
    const size_t n = nodes.size();
    if ((int)n != axis_sizes_[d]) {
      cudaFree(d_nodes_[d]);
      cudaFree(d_weights_[d]);
      NT_CUDA_CHECK(cudaMalloc(&d_nodes_[d], n * sizeof(double)));
      NT_CUDA_CHECK(cudaMalloc(&d_weights_[d], n * sizeof(double)));
      axis_sizes_[d] = (int)n;
    }
    NT_CUDA_CHECK(cudaMemcpy(d_nodes_[d], nodes.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    NT_CUDA_CHECK(cudaMemcpy(d_weights_[d], weights.data(), n * sizeof(double), cudaMemcpyHostToDevice));
  }

  /// @brief Integrate for every external momentum in `p_grid` in one phase-1 launch + one
  ///        segmented reduction. `dress...` are the (small, trivially-copyable) dressing
  ///        stubs forwarded by value to the device kernel.
  template <class... Dress>
  MapResult map(const std::vector<double> &p_grid, double k, double prefactor, const Dress &...dress) {
    const auto t0 = std::chrono::steady_clock::now();
    const int n_p = (int)p_grid.size();
    long long total_quad_points = 1;
    QuadDev<NV> quad_dev;
    for (int d = 0; d < NV; ++d) {
      quad_dev.nodes[d] = d_nodes_[d];
      quad_dev.weights[d] = d_weights_[d];
      quad_dev.n[d] = axis_sizes_[d];
      total_quad_points *= axis_sizes_[d];
    }
    const long long total = total_quad_points * n_p;

    grow(d_buf_, buf_cap_, (size_t)total);
    grow(d_out_, out_cap_, (size_t)n_p);
    grow(d_p_, p_cap_, (size_t)n_p);
    NT_CUDA_CHECK(cudaMemcpy(d_p_, p_grid.data(), n_p * sizeof(double), cudaMemcpyHostToDevice));
    if ((size_t)(n_p + 1) > off_cap_) {
      cudaFree(d_offsets_);
      NT_CUDA_CHECK(cudaMalloc(&d_offsets_, (n_p + 1) * sizeof(long long)));
      off_cap_ = n_p + 1;
    }
    std::vector<long long> offsets(n_p + 1);
    for (int i = 0; i <= n_p; ++i) offsets[i] = (long long)i * total_quad_points;
    NT_CUDA_CHECK(cudaMemcpy(d_offsets_, offsets.data(), (n_p + 1) * sizeof(long long), cudaMemcpyHostToDevice));

    // phase 1: integrand to memory
    const int block = kPhase1BlockSize;
    const long long want = (total + block - 1) / block;
    // cap the launch to CUDA's max grid x-dimension (2^20 blocks); nt_phase1 grid-strides over any tail.
    const int grid = (int)(want < (1 << 20) ? want : (1 << 20));
    NT_CUDA_CHECK(cudaEventRecord(e0_));
    nt_phase1<KERNEL, NV><<<grid, block>>>(quad_dev, d_p_, total, k, prefactor, d_buf_, dress...);
    NT_CUDA_CHECK(cudaGetLastError());
    NT_CUDA_CHECK(cudaEventRecord(e1_));

    // phase 2: pure segmented reduction, one segment per grid point
    size_t tmp_bytes = 0;
    NT_CUDA_CHECK(cub::DeviceSegmentedReduce::Sum(nullptr, tmp_bytes, d_buf_, d_out_, n_p, d_offsets_,
                                                  d_offsets_ + 1));
    if (tmp_bytes > tmp_cap_) {
      cudaFree(d_tmp_);
      NT_CUDA_CHECK(cudaMalloc(&d_tmp_, tmp_bytes));
      tmp_cap_ = tmp_bytes;
    }
    NT_CUDA_CHECK(cub::DeviceSegmentedReduce::Sum(d_tmp_, tmp_bytes, d_buf_, d_out_, n_p, d_offsets_,
                                                  d_offsets_ + 1));
    NT_CUDA_CHECK(cudaEventRecord(e2_));
    NT_CUDA_CHECK(cudaEventSynchronize(e2_));

    MapResult res;
    res.integral.resize(n_p);
    NT_CUDA_CHECK(cudaMemcpy(res.integral.data(), d_out_, n_p * sizeof(double), cudaMemcpyDeviceToHost));
    // the loop-independent constant() once per grid point (host call — the kernels are __host__ too)
    for (int i = 0; i < n_p; ++i) res.integral[i] += real_of(KERNEL::constant(p_grid[i], k, dress...));
    NT_CUDA_CHECK(cudaEventElapsedTime(&res.phase1_ms, e0_, e1_));
    NT_CUDA_CHECK(cudaEventElapsedTime(&res.phase2_ms, e1_, e2_));
    res.wall_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    res.evals = total;
    return res;
  }

  /// @brief Bytes currently held by the phase-1 value buffer (reporting).
  size_t buffer_bytes() const { return buf_cap_ * sizeof(double); }

private:
  /// Grow-only device buffer: realloc to `need` doubles only when it exceeds the current `cap`, so
  /// repeated map() calls of similar size reuse the allocation (never shrinks; freed in the dtor).
  static void grow(double *&ptr, size_t &cap, size_t need) {
    if (need <= cap) return;
    cudaFree(ptr);
    NT_CUDA_CHECK(cudaMalloc(&ptr, need * sizeof(double)));
    cap = need;
  }

  double *d_nodes_[NV], *d_weights_[NV];
  int axis_sizes_[NV];
  double *d_buf_ = nullptr;
  size_t buf_cap_ = 0;
  long long *d_offsets_ = nullptr;
  size_t off_cap_ = 0;
  double *d_out_ = nullptr;
  size_t out_cap_ = 0;
  double *d_p_ = nullptr;
  size_t p_cap_ = 0;
  void *d_tmp_ = nullptr;
  size_t tmp_cap_ = 0;
  cudaEvent_t e0_, e1_, e2_;
};

} // namespace numtracer::cuda
