/// @file
/// @brief Read-only finite-state element input and requested numerical outputs.
#ifndef FESA_MODEL_ELEMENT_TRIAL_H_
#define FESA_MODEL_ELEMENT_TRIAL_H_
#include "fesa/math/matrix.h"
namespace fesa {
struct IncrementalSolidTrialContext;
/// Compute only the requested output; unrequested output pointers must be null.
enum class ElementEvaluationRequest {
  kInternalForce,
  kTangent,
  kBoth,
  /// B33 only: tangent-only integral B^T C B, excluding geometric stiffness.
  /// An approximate static Newton Jacobian, not a consistent derivative.
  kMaterialTangent
};
/// Recoverable numerical failure, separate from invalid input/resource/library
/// status. kNone does not turn a non-OK StatusCode into success.
enum class TrialEvaluationFailure { kNone, kGeometry, kRotationBranch };
/// Caller owns input buffers, borrowed for one evaluation only. No
/// AnalysisState or element-owned official output pointer is retained or
/// mutated.
struct ElementTrialState {
  const IncrementalSolidTrialContext* incremental_solid = nullptr;
  ///< Optional call-lifetime baseline and separate candidate stress output.
  ///< StVK/Shell callers leave null; no official/committed mutation is allowed.
  const Vector* translations =
      nullptr;  ///< Global3n, node-major, length units.
  const Matrix* orientations = nullptr;  ///< Unit wxyz n×4, initial→current.
  const Matrix* committed_orientations = nullptr;  ///< Last accepted q, n×4.
  const Vector* body_spin = nullptr;  ///< Current a_c+bodyLog(q_c^-1 q), rad3n.
  /// Optional current Beam eta_a/eta_b in connectivity order, m. Caller owns
  /// this call-lifetime borrow; never mutated or retained. Other families
  /// must leave null. Current/committed pose callers supply matching-time eta.
  const Vector* beam_internal_position = nullptr;
};
}  // namespace fesa
#endif  // FESA_MODEL_ELEMENT_TRIAL_H_
