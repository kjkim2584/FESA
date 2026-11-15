/// @file
/// @brief Private selected-D structural residual and exact spatial tangent.
#ifndef FESA_ELEMENTS_BEAM3D_NONLINEAR_OPS_H_
#define FESA_ELEMENTS_BEAM3D_NONLINEAR_OPS_H_

#include "elements/beam3d_kinematics.h"
#include "fesa/elements/beam3d.h"

namespace fesa {
/// Prepares current global positions and material t/n1 rows from immutable
/// model coordinates, translations6 and unit wxyz orientations2x4. Borrows
/// inputs for the call; owns/reclaims all temporary storage. Outputs are
/// distinct Matrix(2,3), separate from pose/initial axes. Validation preserves
/// outputs; later numerical failure invalidates them. No material/eta/state
/// evaluation or mutation. Shared by structural, geometry and motion adapters.
/// Optional separate2x3 relative_positions receives (X-X_a)+u, computed
/// before X+u rounding. Same failure/lifetime contract; null skips this output.
StatusCode PrepareBeamNonlinearPose(const Beam3D* beam,
                                    const ElementTrialState* state,
                                    Matrix* positions, Matrix* tangents,
                                    Matrix* normals, Diagnostic* d,
                                    Matrix* relative_positions = nullptr);
/// Borrows all inputs for this call; no ownership or model/result mutation.
/// Current translations6 and unit wxyz orientations2x4 define the pose;
/// eta2 is in m. Force14 and tangent14x14 use spatial [ua,wa,ub,wb,eta2].
/// The four nodal 3-blocks are global force N/moment Nm; eta residuals are N.
/// Tangent columns differentiate m/rad coordinates without symmetrization.
/// Axial Green strain uses Gauss3; bending/twist use Gauss4. Both material
/// and geometric tangent terms differentiate this same selected energy.
/// kMaterialTangent requests only integral B^T C B, with the same coordinates
/// and quadrature; this approximate Jacobian omits stress-weighted mixed
/// derivatives. It never changes an internal-force evaluation or model state.
/// Only requested outputs may be nonnull and must be separate from inputs.
/// Preflight errors preserve outputs; later failures invalidate requested
/// outputs. K1 geometry/branch and original backend diagnostics propagate.
/// Optional history inputs are validated but never used to define the pose.
/// Workspace fields expire at its next evaluation; no scratch is retained.
/// All pose/eta/workspace/failure pointers are required; d is optional.
/// @see docs/Linear static/b33_element_response_contract.md.
StatusCode EvaluateBeamNonlinearResponse(
    const Beam3D* beam, const ElementTrialState* state,
    const Vector* internal_position, BeamKinematicWorkspace* workspace,
    ElementEvaluationRequest request, Vector* force, Matrix* tangent,
    TrialEvaluationFailure* failure, Diagnostic* d);
}  // namespace fesa

#endif  // FESA_ELEMENTS_BEAM3D_NONLINEAR_OPS_H_
