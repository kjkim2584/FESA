/// @file
/// @brief Nonlinear solid material-point output and F11 tensor recovery, M9.
#ifndef FESA_ELEMENTS_SOLID_NONLINEAR_RECOVERY_H_
#define FESA_ELEMENTS_SOLID_NONLINEAR_RECOVERY_H_

#include "fesa/elements/solid_reference_data.h"
#include "fesa/model/elastic_material.h"
#include "fesa/model/element_trial.h"

namespace fesa {
/// Borrows finite dimensionless F3x3 and material. Caller explicitly supplies
/// actual F for ordinary TL or Fbar for C3D8. Required caller-owned distinct
/// 3x3 outputs: spatial Cauchy stress [force/length^2] and LE=log(FF^T)/2
/// [dimensionless], both full tensors in the global basis (tensor shear).
/// Outputs cannot alias F or each other. Preflight preserves numerical outputs;
/// computation failure invalidates them. Required failure is cleared on entry;
/// nonpositive/singular J or SPD-log domain failure sets kGeometry. Backend,
/// allocation and arithmetic failures remain fatal. No state is committed.
/// A recovery caller must terminate on failure per HHT P05-P07, not retry DH.
StatusCode ComputeSolidPointOutput(const Matrix* deformation_gradient,
    const ElasticMaterial* material, Matrix* cauchy, Matrix* logarithmic_strain,
    TrialEvaluationFailure* failure, Diagnostic* d);

/// Borrows the immutable full-formulation stiffness table and its ng borrowed
/// material-point F pointers, each 3x3, with the same F/Fbar choice as above.
/// Caller owns four distinct outputs: IP S/LE ng x6, element-node S/LE n x6,
/// ordered [11,22,33,12,13,23], tensor shear, in global force/length^2 and
/// dimensionless units respectively. Converts IPs before F11 extrapolation.
/// All inputs/outputs must be disjoint. Preflight preserves numerical outputs;
/// computation failure invalidates requested payload. Required failure is
/// cleared on entry. No validity flag, element state or writer is changed;
/// caller keeps results unconfirmed until the common recovery/commit succeeds.
StatusCode RecoverSolidPointOutputs(const SolidReferenceData* data,
    const Matrix* const* material_gradients, const ElasticMaterial* material,
    Matrix* ip_stress, Matrix* ip_strain, Matrix* nodal_stress,
    Matrix* nodal_strain, TrialEvaluationFailure* failure, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_SOLID_NONLINEAR_RECOVERY_H_
