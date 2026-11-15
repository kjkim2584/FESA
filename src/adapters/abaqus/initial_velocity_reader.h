/// @file
/// @brief Resolves approved translational/angular initial velocity input.
#ifndef FESA_ADAPTERS_ABAQUS_INITIAL_VELOCITY_READER_H_
#define FESA_ADAPTERS_ABAQUS_INITIAL_VELOCITY_READER_H_

#include "fesa/analysis/dof_manager.h"
#include "fesa/io/input_deck.h"
#include "fesa/math/vector.h"

namespace fesa {
namespace abaqus {
/// Resolves validated targets into global Cartesian velocity (L/T, rad/T).
/// Borrows a matching Domain and first dynamic step's 3N/6N mapping.
/// Components must exist on each target node; 3N rejects angular components.
/// Unspecified entries are zero; duplicate expanded DOFs are input errors,
/// including equal values. Fixed nonzero velocity is rejected. Participation,
/// mass and start-route validation belong to the dynamic analysis caller.
/// Caller owns *out on success; creation-only out is null on every failure.
/// Inputs remain unchanged. Diagnostic source_file borrows deck's path.
StatusCode ResolveInitialVelocity(const InputDeck* deck, const Domain* domain,
                                  const DofManager* dofs, Vector** out,
                                  Diagnostic* d);
}  // namespace abaqus
}  // namespace fesa
#endif  // FESA_ADAPTERS_ABAQUS_INITIAL_VELOCITY_READER_H_
