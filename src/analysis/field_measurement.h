/// @file
/// @brief Shared physical field measurements, without convergence decisions.
#ifndef FESA_ANALYSIS_FIELD_MEASUREMENT_H_
#define FESA_ANALYSIS_FIELD_MEASUREMENT_H_

#include "fesa/analysis/dof_manager.h"
#include "fesa/analysis/residual_report.h"
#include "fesa/io/control_settings.h"

namespace fesa {
class AnalysisModel;
class AnalysisState;
class EffectiveConditions;
/// Creates field0 force[F] or field1 moment[FL] origins in original element
/// ID, block, connectivity/component order, then appends already typed external
/// origins. Each element has block_count physical blocks of its actual DOFs.
/// With material_damping=true, adds one final block only for a Solid/MITC4
/// material with nonzero Rayleigh scalars. Default false preserves old callers.
/// With beam_internal=true, requires Beam12+2 and packs14 per block, selecting
/// eta entries into force origins. Default false preserves linear Beam12.
/// Retains signs, zeros and multiplicity, including fixed entries. Supports
/// three/six DOFs per node; exact packed shape required. Null packed is valid
/// only for no elements. Borrows inputs for call; no state, weighting or sums.
/// Required creation output is null on empty success/failure; caller deletes
/// success. Inputs unchanged; scratch released before return.
StatusCode CreateElementFieldFlux(const AnalysisModel* model,
                                  const Vector* element_values,
                                  int block_count, int field,
                                  const Vector* external_origins,
                                  Vector** out, Diagnostic* d,
                                  bool material_damping = false,
                                  bool beam_internal = false);
/// Splits raw linear dynamic origins (element structural n,inertia n, optional
/// damping n for nonzero material Rayleigh, then
/// CLOADs in EffectiveConditions order) into force[F] and moment[FL] fields.
/// Supports actual-size Solid 3DOF and MITC4 6DOF layouts; retains within-field
/// order, signs, zeros and fixed contributions. No weighting, numerical
/// summation, history or decision. Requires matching bound conditions and exact
/// shape. Borrows inputs only for call. Distinct creation-only outputs
/// required; caller owns success. Absent field and all failure outputs are
/// null; inputs remain unchanged. Uses library gather, call-local scratch and
/// no cache.
StatusCode SplitLinearDynamicFluxOrigins(const AnalysisModel* model,
                                         const EffectiveConditions* conditions,
                                         const Vector* origins, Vector** force,
                                         Vector** moment, Diagnostic* d);
/// Borrows values/controls/output for the call. Origin values retain signs,
/// zeros and multiplicity, in force [F] or moment [F*L]; no temporal weights.
/// previous is the nonnegative active baseline in the same field units.
/// Uses abs(Q)>=epsilon_l*previous unless max(abs(Q))<previous/10, when all
/// values participate. Zero threshold includes zeros. Empty selection fails;
/// no invented denominator/scale. Null values means inactive field (no-op,
/// controls not consumed). Output required, field is 0 force or 1 moment.
/// Writes only contribution_count/maximum_flux/spatial_mean; other residual
/// data remain untouched. Preflight preserves output; later failure invalidates
/// those summaries. No time averaging, field conversion, decision or commit.
/// @see docs/HHT/convergence_contract.md R11L and MITC4 field history
/// contract.
StatusCode ComputeSpatialFluxMean(const Vector* values,
                                  const FieldControls* controls,
                                  double previous, int field,
                                  FieldResidual* output, Diagnostic* d);
/// Measures only free equations of field0 force or field1 moment, in physical
/// global residual order. Borrows inputs/output; sets dof_count, maximum,
/// original node/component at maximum, l2 and rms. Empty field writes zeros.
/// Other summaries unchanged. Invalid arguments preserve output; later failure
/// invalidates measured fields. No global resultant diagnostic is computed.
StatusCode MeasureResidualField(const Vector* residual, const DofManager* dofs,
                                int field, FieldResidual* output,
                                Diagnostic* d);
/// Measures nonlinear Beam force residuals including 2E internal entries [F].
/// Physical free translations precede sorted internal entries when maxima tie.
/// Reports node/component or element/end, clearing the other location. Inputs
/// borrowed; preflight preserves output, numerical failure invalidates it.
StatusCode MeasureBeamResidualField(const Vector* physical_residual,
                                    const Vector* internal_residual,
                                    const DofManager* dofs,
                                    FieldResidual* output, Diagnostic* d);
/// Adds internal eta increment/correction maxima [L] to an already measured
/// force field. Requires matching prepared Beam state/mapping. Null correction
/// leaves its maximum unchanged. Borrows all inputs; no state/history mutation.
/// Preflight preserves output; numerical failure invalidates measurements.
StatusCode MeasureBeamInternalMotion(const DofManager* dofs,
                                     const AnalysisState* state,
                                     const Vector* actual_correction,
                                     FieldResidual* force, Diagnostic* d);
/// Linear 3N dynamic measurements in model units/global physical order.
/// Uses free residuals, participating (including fixed) candidate-committed
/// increments and candidate-predictor actual direct-solve corrections [L].
/// Borrows all pointers; never changes vectors/state. Updates residual and
/// increment/correction measurements only. Preflight preserves output; later
/// failure invalidates those fields. No fake zero correction or time history.
StatusCode MeasureLinearDynamicField(const DofManager* dofs,
                                     const Vector* residual,
                                     const Vector* committed,
                                     const Vector* candidate,
                                     const Vector* predictor,
                                     FieldResidual* output, Diagnostic* d);
/// Field-explicit linear variant for 3/6DOF physical mappings. Field0 selects
/// translations (force residual F, increments/corrections L); field1 selects
/// rotations (moment residual FL, increments/corrections rad). Absent fields
/// write zero measurements. Free residuals and participating physical motion
/// follow the contract above; fixed motion participates, excluded nodes do not.
/// Inputs may alias each other and are never modified. Invalid field/shape
/// preserves output; backend failure invalidates measurements. No convergence
/// decision, global resultant, scale conversion or history changes.
StatusCode MeasureLinearDynamicField(const DofManager* dofs, int field,
                                     const Vector* residual,
                                     const Vector* committed,
                                     const Vector* candidate,
                                     const Vector* predictor,
                                     FieldResidual* output, Diagnostic* d);
/// Shared nonlinear participating motion measurements, in physical node order.
/// Borrows matching mapping/prepared nonlinear state and optional actual
/// global[du,dw]. Measures u_trial-u_committed and spatial
/// Log(R_trial*R_committed^T), including prescribed DOFs, excluding nonparticipating nodes.
/// Writes only increment_max/correction_max of participating force/moment
/// fields, in L/rad. Null correction leaves correction_max unchanged (NA is
/// tracked by the caller); absent fields and other summary members untouched.
/// Distinct field outputs and required rotation_branch; clears branch on entry
/// and sets true only for a recoverable Log chart failure. No state/history
/// mutation, convergence/scale decision, global resultant or retained scratch.
/// Preflight preserves fields; after numerical failure discard measurements.
StatusCode MeasureNonlinearMotionFields(const DofManager* dofs,
                                        const AnalysisState* state,
                                        const Vector* actual_correction,
                                        FieldResidual* force,
                                        FieldResidual* moment,
                                        bool* rotation_branch, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_FIELD_MEASUREMENT_H_
