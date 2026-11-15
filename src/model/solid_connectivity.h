/// @file
/// @brief ID-based conformity of typed solid interpolation traces.
#ifndef FESA_MODEL_SOLID_CONNECTIVITY_H_
#define FESA_MODEL_SOLID_CONNECTIVITY_H_
#include "fesa/model/domain.h"
namespace fesa {
/// Checks shared triangle3/triangle6/quad4 traces after node/DOF membership.
/// Borrows nonnull Domain for this call; changes no definitions or mappings.
/// Failure records an unsupported interface/element in optional d. Local
/// allocation exceptions propagate to the caller's allocation-failure boundary.
/// No geometry, results, numeric buffers, or cross-call cache is used.
StatusCode ValidateSolidConnectivity(const Domain* domain, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_MODEL_SOLID_CONNECTIVITY_H_
