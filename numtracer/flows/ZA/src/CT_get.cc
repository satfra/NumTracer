#include "../kernel.hh"

#include "../ZA.hh"

void ZA_integrator::get(double& dest, const double& p, const Spline& ZA)
{
  integrator.get(dest, p,  ZA);
}