#include "../kernel.hh"

#include "../ZA.hh"

DiFfRG::GPU_exec ZA_integrator::map(double* dest, const LogarithmicCoordinates1D<double>& coordinates, const Spline& ZA)
{
  return integrator.map(dest, coordinates, ZA);
}