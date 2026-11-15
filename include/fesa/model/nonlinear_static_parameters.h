/// @file
/// @brief Immutable nonlinear STATIC input parameters, consistent model time.
#ifndef FESA_MODEL_NONLINEAR_STATIC_PARAMETERS_H_
#define FESA_MODEL_NONLINEAR_STATIC_PARAMETERS_H_
#include "fesa/base/diagnostic.h"
namespace fesa {
/// Model-layer scalar definition, to be owned by StepDefinition. No solver or
/// vendor data. first-step NLGEOM is an input/preflight responsibility.
/// @see docs/MITC4/increment_execution_contract.md and
/// abaqus_increment_control.md.
struct NonlinearStaticParameters {
  double period = 1;
  double initial_increment = 1;
  double minimum_increment = 1e-5;
  bool has_maximum = false;
  double maximum_increment = 0;  ///< Unused unless has_maximum.
  int increment_limit = 100;     ///< STEP INC, committed increments only.
  bool extrapolate = true;       ///< STEP EXTRAPOLATION=LINEAR, false for NO.
};
/// Required finite positive period/initial/min; min<=initial. Explicit max is
/// finite positive and >=initial. INC positive. No clamp or time epsilon.
/// Borrowed inputs; failure changes only optional diagnostic, no retained data.
StatusCode ValidateNonlinearStaticParameters(
    const NonlinearStaticParameters* parameters, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_MODEL_NONLINEAR_STATIC_PARAMETERS_H_
