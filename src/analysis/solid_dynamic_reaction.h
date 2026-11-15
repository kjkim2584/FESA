/// @file
/// @brief End-state reaction estimates for linear Solid/MITC4.
#ifndef FESA_ANALYSIS_SOLID_DYNAMIC_REACTION_H_
#define FESA_ANALYSIS_SOLID_DYNAMIC_REACTION_H_

#include "fesa/analysis/analysis_model.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/elements/element_workspace.h"

namespace fesa {
/// Borrowed HHT history and current individual external contributions.
struct HhtSolidForceInput {
  const Vector* previous_displacement =
      nullptr;  ///< Global committed u [L, rad].
  const Vector* previous_velocity = nullptr; ///< Required for damped history.
  const Vector* external_origins =
      nullptr;          ///< Current forces/moments [F, FL].
  double alpha = -.05;  ///< Dimensionless [-1/3,0].
};
/// Evaluates H1 with both current and previous structural/damping resistance,
/// and current algorithmic inertia. previous_velocity is required when a
/// damped element contributes history with alpha!=0.
/// Same physical layout/support/worker contract as reaction below. Includes
/// fixed rows; caller checks free rows only. weighted_external is already
/// (1+alpha)*p_next-alpha*p_previous. No recovery/decision/state mutation.
/// All inputs borrowed through join; residual must be distinct from inputs.
/// Creates *origins in original element order: structural n,inertia n for each
/// element, then hht.external_origins in supplied order (null means none).
/// A material with nonzero Rayleigh adds damping n after inertia n, preserving
/// the previous layout/count for undamped materials. All are current values.
/// n is each element's actual local DOF count, without stride padding.
/// No weighting, assembly, constraint filter or zero removal in origins.
/// Caller owns successful origins; failure nulls output and frees partials.
/// Failure preserves residual until final CopyFrom; backend failure there
/// invalidates residual. hht/origins required, output creation-only.
StatusCode ComputeSolidHhtResidual(
    const AnalysisModel* model, const DofManager* dofs, MassMode mode,
    const Vector* displacement, const Vector* velocity,
    const Vector* acceleration, const Vector* weighted_external,
    const HhtSolidForceInput* hht, ElementWorkspace** workers, int worker_count,
    Vector* residual, Vector** origins, Diagnostic* d);

/// Evaluates A27/O1 estimate: f_int+M*A+C*v-external on fixed DOFs, zero elsewhere.
/// Inputs are borrowed global node-ID ordered 3N/6N Cartesian u [L, rad],
/// v [L/T, rad/T],
/// a [L/T^2, rad/T^2], external [F, FL]. Fixed v/a must be zero. Supports
/// supported linear Solid (including different formulations/local sizes) or
/// homogeneous small-rotation MITC4; kDefault mass is resolved per element.
/// No HHT
/// weighting, convergence test, recovery or commit. Workers are distinct
/// borrowed exclusive contexts through joined dispatch; MITC4 permits null
/// contexts/array. Writes reaction [F, FL], disjoint from every input. Failure
/// preserves output until final CopyFrom; a failure of that backend invalidates
/// output. Never changes official node/element results or analysis state. No
/// cached buffers.
/// @see docs/HHT/formulation.md O1 and implementation_contract.md.
StatusCode ComputeSolidDynamicReaction(
    const AnalysisModel* model, const DofManager* dofs, MassMode mode,
    const Vector* displacement, const Vector* velocity,
    const Vector* acceleration, const Vector* external,
    ElementWorkspace** workers, int worker_count, Vector* reaction,
    Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_SOLID_DYNAMIC_REACTION_H_
