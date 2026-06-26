// NumTracer GPU tests — CPU same-quadrature reference integrators. Implemented in separate
// plain-C++ TUs (ref_za3.cpp / ref_za4.cpp) that evaluate the committed HOST kernels
// (ZA3_num_kernel.hh / ZA4_num_kernel.hh): the .cu main TU never has to compile the host
// kernel headers, and the GPU result is checked against an independently compiled path.
// Accumulation is long double so the reference's summation error stays far below the
// validation tolerance (cub's pairwise device sums are more accurate than a naive serial sum).
#pragma once

#include "quadrature.hpp"

#include <vector>

// Integrals over the (l1, cos1, cos2[, phi]) quadrature grid, one per p in p_grid:
//   Σ prefactor · l1³ · Πw · Re kernel(...)  + constant(p, k, ...).
std::vector<double> cpu_reference_za3(const Axis ax[3], const std::vector<double> &p_grid, double k,
                                      double prefactor);
std::vector<double> cpu_reference_za4(const Axis ax[4], const std::vector<double> &p_grid, double k,
                                      double prefactor);
