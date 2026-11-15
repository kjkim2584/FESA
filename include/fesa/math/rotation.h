/// @file
/// @brief Unit-quaternion SO(3) operations for spatial finite corrections.
#ifndef FESA_MATH_ROTATION_H_
#define FESA_MATH_ROTATION_H_

#include "fesa/math/matrix.h"

namespace fesa {
/// All objects borrowed for call; inputs/outputs distinct. Quaternions are
/// unit wxyz, map initial global axes to current global axes, retain cover
/// sign. Vectors are radians. Invalid shape/alias/nonfinite/unit data preserve
/// out; after backend failure discard out. d optional. No retained data or
/// cache. Unit tolerance4096*epsilon permits rounding normalization only.
StatusCode NormalizeQuaternion(const Vector* q, Vector* out, Diagnostic* d);
/// Exp(v), no sign canonicalization even for angles exceeding pi. Callers guard
/// actual increment paths separately; output alone cannot detect a4pi path.
StatusCode RotationExp(const Vector* v, Vector* out, Diagnostic* d);
/// out=left*right (Hamilton product), normalized without changing cover sign.
StatusCode ComposeRotations(const Vector* left, const Vector* right,
                            Vector* out, Diagnostic* d);
/// Unit q to initial-to-current global3x3 rotation, rows current components.
StatusCode RotationMatrix(const Vector* q, Matrix* out, Diagnostic* d);
/// Principal global axis-angle output norm<=pi; q/-q agree. At exact pi the
/// largest absolute axis component is positive, ties x,y,z. Output only;
/// never feed canonical sign back into increment history.
StatusCode PrincipalRotationVector(const Vector* q, Vector* out, Diagnostic* d);
/// b=Log(ref^-1*current) in reference/body axes, without sign canonicalization.
/// Requires relative w>256*epsilon. Branch failure sets *branch_exceeded=true,
/// returns kEquilibriumFailure and preserves out. Other failures set false.
/// branch_exceeded required; callers also guard each applied correction path.
StatusCode BodyRotationIncrement(const Vector* reference, const Vector* current,
                                 Vector* out, bool* branch_exceeded,
                                 Diagnostic* d);
/// Spatial counterpart Log(current*ref^-1), same chart/failure contract.
StatusCode SpatialRotationIncrement(const Vector* reference,
                                    const Vector* current, Vector* out,
                                    bool* branch_exceeded, Diagnostic* d);
/// Inverse left Jacobian I-.5[b]x+c[b]x^2, norm(b)<pi; maps spatial local
/// perturbation of Exp(b) to db. Solà appendixB145-146; small-angle series
/// c=1/12+theta^2/720+theta^4/30240 for theta<=.01.
StatusCode InverseLeftRotationJacobian(const Vector* b, Matrix* out,
                                       Diagnostic* d);
}  // namespace fesa
#endif  // FESA_MATH_ROTATION_H_
