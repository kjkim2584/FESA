/// @file
/// @brief Call-local finite S/LE sampling, no force/tangent/mass or
/// publication.
#ifndef FESA_ELEMENTS_MITC4_NONLINEAR_RECOVERY_H_
#define FESA_ELEMENTS_MITC4_NONLINEAR_RECOVERY_H_
#include "fesa/math/matrix.h"
#include "fesa/model/element_trial.h"
#include "fesa/model/node.h"
#include "fesa/model/shell_section.h"
namespace fesa {
/// Uses only translations3n and unit wxyz orientations nx4, never body spin or
/// principal UR. Four initial flat nodes, linear isotropic shell with no rho
/// requirement. Outputs are distinct caller-owned16x6 S/LE and3x3 row axes.
/// Rows are4plane Gauss points for each zeta[-a,+a,-1,+1]. S tensor shear,
/// LE engineering shear, order[11,22,33,12,13,23], one current center basis.
/// Required failure flag cleared on entry. Invalid initial definition/input,
/// allocation/library are fatal; valid-initial current geometry or positive
/// stretch domain failures typed geometry. d optional, includes failed point.
/// No input/state/official result mutation, cache or retained resources.
/// Preflight preserves outputs; later failure invalidates them all.
/// @see docs/MITC4/nonlinear_recovery_execution_contract.md.
StatusCode RecoverMitc4FiniteSamples(
    const Node* const* nodes, int count, const ShellSection* section,
    const Vector* translations, const Matrix* orientations,
    Matrix* stress_samples, Matrix* strain_samples, Matrix* output_axes,
    TrialEvaluationFailure* failure, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_MITC4_NONLINEAR_RECOVERY_H_
