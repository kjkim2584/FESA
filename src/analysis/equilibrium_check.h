/// @file
/// @brief Recomputed first-step equilibrium with S1 field convergence rules.
#ifndef FESA_ANALYSIS_EQUILIBRIUM_CHECK_H_
#define FESA_ANALYSIS_EQUILIBRIUM_CHECK_H_
#include "analysis/static_assembler.h"
#include "fesa/analysis/analysis_state.h"
#include "fesa/analysis/residual_report.h"
#include "fesa/io/control_settings.h"
namespace fesa {
/// Deterministic contribution/byte bounds, without allocation or computation.
/// Caller uses before HDF5 truncate; evaluator reuses identical count checks.
StatusCode ValidateEquilibriumCounts(const AnalysisModel* model,
                                     const EffectiveConditions* conditions,
                                     Diagnostic* d);
/// S1 scalar decision seam: field contains measured means/norms/corrections.
/// Does not change measurements. Sets branch/limits/passed; rejects invalid
/// scales. displacement selects the small-increment/Lc exception.
StatusCode EvaluateFieldCriteria(const FieldControls* controls,
                                 bool displacement,
                                 double characteristic_length,
                                 FieldResidual* field, Diagnostic* d);
/// Evaluates force/displacement and moment/rotation together using the S1
/// all-active-field linear condition in abaqus_controls.md. Borrows all inputs
/// for this call; controls and both distinct fields must be non-null. d may be
/// null. Changes only both fields' phi, branch, limits and passed flag. An
/// equilibrium failure retains complete field diagnostics; invalid/numerical
/// failure invalidates the decision. No-DOF fields are excluded from the gate.
StatusCode EvaluateJointFieldCriteria(const ControlSettings* controls,
                                      double characteristic_length,
                                      FieldResidual* force,
                                      FieldResidual* moment, Diagnostic* d);
struct FieldConvergenceOptions;
/// Nonlinear extension of the same field criteria. Options are borrowed and
/// define actual/estimated correction validity and the Rn/Rp selection.
/// The small-increment exception requires an available actual correction;
/// an unchanged initial/predictor pose is not a measured negligible response.
/// Legacy callers pass through the unchanged options-free entry point above.
StatusCode EvaluateJointFieldCriteriaWithOptions(
    const ControlSettings* controls, double characteristic_length,
    const FieldConvergenceOptions* force_options,
    const FieldConvergenceOptions* moment_options, FieldResidual* force,
    FieldResidual* moment, Diagnostic* d);
/// Only step 1, initial committed U/UR=0: actual increment and correction are
/// candidate_u. Later steps are unsupported until committed history is passed.
/// Candidate prescribed components must equal their effective input values.
/// Recomputes element forces from borrowed candidate; never uses stored
/// results. Overwrites workspace/report; failure invalidates scratch, preserves
/// model. Uses actual physical 3/6 component mappings and element 3n force
/// contributions. Three-component models exclude the rotation field and its
/// conversion. Each element must supply a supported characteristic length; no
/// fallback length is inferred here. kEquilibriumFailure retains full report,
/// including reaction=fint-p on constrained components and zero elsewhere. No
/// state/history commit occurs.
StatusCode CheckEquilibrium(const AnalysisModel* model, const DofManager* dofs,
                            const EffectiveConditions* conditions,
                            const ControlSettings* controls,
                            const Vector* candidate_u,
                            StaticWorkspace* workspace, ResidualReport* report,
                            Diagnostic* d);
/// History-aware evaluator. The borrowed state supplies committed total U/UR,
/// raw means and current maxima; never mutates it. Optional workers are
/// borrowed until all disjoint contribution tasks join; null owns serial
/// scratch. Candidate minus committed is the actual increment/correction.
/// Incomplete numerical reports are discarded.
StatusCode CheckEquilibrium(const AnalysisModel* model, const DofManager* dofs,
                            const EffectiveConditions* conditions,
                            const ControlSettings* controls,
                            const Vector* candidate_u,
                            const AnalysisState* state,
                            StaticWorkspace* workspace, ResidualReport* report,
                            Diagnostic* d, ElementWorkspace** workers = nullptr,
                            int worker_count = 1);
}  // namespace fesa
#endif  // FESA_ANALYSIS_EQUILIBRIUM_CHECK_H_
