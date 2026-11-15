/// @file
/// @brief Same-candidate Cauchy S/current-F LE and independent tensor recovery.
#ifndef FESA_ELEMENTS_INCREMENTAL_SOLID_RECOVERY_H_
#define FESA_ELEMENTS_INCREMENTAL_SOLID_RECOVERY_H_
#include "elements/incremental_elastic_workspace.h"
#include "fesa/elements/solid_reference_data.h"
#include "fesa/model/element_trial.h"
namespace fesa {
/// Test-only thread-local observer at official recovery entry; caller resets
/// before borrowed callback state expires. No independent validation is run.
void SetIncrementalRecoveryObserver(void (*observer)());
/// Borrows ng current F3×3 pointers and same-candidate ng×6 Cauchy stress,
/// Ordinary types use physical F; HEX8 requires its worker's modified_f_data
/// (Fbar) array, as frozen by the coupled HEX8 contract. Neither is history.
/// tensor order [11,22,33,12,13,23]. Writes IP ng×6 and element-node n×6 S/LE;
/// stress F/L², LE dimensionless. Preflight preserves outputs; subsequent
/// failure invalidates all outputs. No material history/commit/file mutation.
/// SPD-log domain failure reports geometry but is fatal in the recovery phase.
StatusCode RecoverIncrementalSolidPoints(
    const SolidReferenceData* data, const Matrix* const* candidate_f,
    const Matrix* candidate_stress, IncrementalElasticWorkspace* workspace,
    Matrix* ip_stress, Matrix* ip_strain, Matrix* nodal_stress,
    Matrix* nodal_strain, TrialEvaluationFailure* failure, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_INCREMENTAL_SOLID_RECOVERY_H_
