/// @file
/// @brief Fixed-step dead loads, reaction release and finite prescriptions.
#ifndef FESA_ANALYSIS_NONLINEAR_STEP_CONDITIONS_H_
#define FESA_ANALYSIS_NONLINEAR_STEP_CONDITIONS_H_

#include "fesa/analysis/analysis_state.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/analysis/nonlinear_problem.h"

namespace fesa {
/// Borrows same-step immutable conditions/mapping and copies of the actual
/// STEP-start state, never the latest increment's changing committed buffers.
/// Physical node-major u/reaction; unit wxyz orientations[N,4] only when
/// rotations exist (otherwise null). Global L/rad,
/// F/FL. These snapshots stay fixed across all attempts/cutbacks of this step.
struct NonlinearStepConditionInput {
  const EffectiveConditions* conditions = nullptr;
  const DofManager* dofs = nullptr;
  const Vector* start_displacement = nullptr;
  const Matrix* start_orientations = nullptr;
  const Vector* start_reaction = nullptr;
  double lambda = 0, committed_lambda = 0;
};
/// Creates global dead load p(lambda) and individually nonzero force/moment
/// origins, including +(1-lambda)R_start for released physical constraints.
/// All outputs creation-only, distinct pointer locations. Caller deletes each
/// successful vector; empty origin lists are null. Failure nulls all outputs,
/// releases partial allocations and preserves every input/physical state.
/// No stiffness, density, numerical history or result output is evaluated.
StatusCode EvaluateNonlinearStepLoads(const NonlinearStepConditionInput* input,
                                      Vector** external, Vector** force_origins,
                                      Vector** moment_origins, Diagnostic* d);
/// Applies only prescribed translations and full rotation groups to prepared
/// candidate, preserving all unconstrained components. Uses actual start u/q,
/// previous INPUT rotation commands for continuation; first addition uses
/// principal Log(q_start). Any failure discards candidate, preserves committed
/// and fixed snapshots. Required failure cleared; only chart/path failure
/// returns kRotationBranch/kEquilibriumFailure. No input retained or output
/// projection/commit performed. All objects borrowed for call.
StatusCode ApplyNonlinearStepPrescriptions(
    const NonlinearStepConditionInput* input, AnalysisState* state,
    NonlinearFailure* failure, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_NONLINEAR_STEP_CONDITIONS_H_
