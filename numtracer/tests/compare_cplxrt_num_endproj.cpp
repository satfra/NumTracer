// The "ComplexEndProjection" half of the finite-density gate (see compare_cplxrt_num.cpp).
//
// gen_cplxrt_numeric.wls emits the same traced kernel a third time with "ComplexEndProjection" ->
// True: no symbolic Pure/RePart projection and no probe, one body that assembles the full complex
// integrand and returns ntRe(...) of it. That is real(complex body) by construction, so it must match
// the complex oracle pointwise; the two bodies are CSE'd differently, so the comparison is graded
// against the oracle's overall scale.
#include "Cplxrtend_num_kernel.hh"

double cplxrt_endproj(double l0, double l1, double cos1, double p0, double p, double muq, double Ep)
{
  return numtracer_kernels::Cplxrtend_num_kernel::kernel(l0, l1, cos1, p0, p, muq, Ep);
}
