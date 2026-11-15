/// @file
/// @brief HHT Cartesian/finite-rotation factory; formulation H2-H4 and explicit
/// lifetimes.
#ifndef FESA_ANALYSIS_HHT_INTEGRATOR_H_
#define FESA_ANALYSIS_HHT_INTEGRATOR_H_

#include "fesa/analysis/time_integrator.h"

namespace fesa {
/// Step integration settings portion; no increment controller or input parser.
/// Future DynamicStepSettings owns this portion and outlives its integrator.
struct HhtParameters {
  double alpha = -.05;  ///< Dimensionless, finite in [-1/3,0].
};
/// Creates a Cartesian/finite-rotation HHT integrator borrowing immutable
/// parameters. Caller owns *out and deletes through TimeIntegrator. Parameters
/// must outlive it and remain unchanged. Failure nulls *out; inputs unchanged.
/// beta=(1-alpha)^2/4 and gamma=1/2-alpha are derived, not independent inputs.
/// Each call rebuilds coefficients for dt; no state or numerical buffer cache.
/// See docs/HHT/formulation.md and rotational_integrator.md
/// for units, validity, chart guards and per-node R2 block layout.
StatusCode CreateHhtIntegrator(const HhtParameters* parameters,
                               TimeIntegrator** out, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_HHT_INTEGRATOR_H_
