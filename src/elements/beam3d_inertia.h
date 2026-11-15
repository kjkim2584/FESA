/// @file
/// @brief Exact physical14 beam inertia and independent spatial q/V/A partials.
#ifndef FESA_ELEMENTS_BEAM3D_INERTIA_H_
#define FESA_ELEMENTS_BEAM3D_INERTIA_H_

#include "elements/beam3d_kinematics.h"
#include "elements/beam3d_time_kinematics.h"
#include "fesa/base/diagnostic.h"
#include "fesa/math/matrix.h"
#include "fesa/math/vector.h"
#include "fesa/model/beam_section.h"
#include "fesa/model/element_dynamic.h"
#include "fesa/model/element_trial.h"

namespace fesa {

/// Integrates physical inertia on the reference length with four Gauss points.
/// All inputs/workspaces are borrowed for this call and never retained. pose
/// uses global spatial [xa3,wa3,xb3,wb3,eta2], m/rad/m; V and A use the same
/// order per s and s^2. pose.direction is required but does not seed GI.
/// section requires positive finite rho*A and rho*(I11+I22); J is unrelated.
/// force14 is required, overwritten with [F,FL,F]. In kForces all matrices
/// must be null; kForcesAndLinearization requires disjoint14x14 GI/CI/MI.
/// GI fixes V/A, CI fixes pose/A, MI fixes pose/V. No HHT weights are applied.
/// Outputs cannot alias inputs or workspace getters. Preflight argument and
/// material failures preserve outputs. Later failure invalidates all requested
/// outputs; caller discards them. Geometry/branch typed failures and original
/// resource/backend diagnostics propagate. Workspaces are overwritten and may
/// be reused after failure; getters expire at the next corresponding Evaluate.
/// No official model/state or history is modified. All14 physical/internal
/// rows/columns are retained, using ordered mixed variations without half skew.
/// @see docs/Linear static/b33_inertia_kernel_contract.md.
StatusCode EvaluateBeamInertia(
    const BeamKinematicInput* pose, const BeamSection* section,
    const Vector* velocity14, const Vector* acceleration14,
    BeamKinematicWorkspace* current, BeamTimeKinematicWorkspace* time,
    DynamicEvaluationRequest request, Vector* inertia_force14, Matrix* gi14,
    Matrix* ci14, Matrix* mi14, TrialEvaluationFailure* failure, Diagnostic* d);

}  // namespace fesa

#endif  // FESA_ELEMENTS_BEAM3D_INERTIA_H_
