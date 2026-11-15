/// @file
/// @brief Exact time Taylor coefficients and independent q/V/A Jacobians.
#ifndef FESA_ELEMENTS_BEAM3D_TIME_KINEMATICS_H_
#define FESA_ELEMENTS_BEAM3D_TIME_KINEMATICS_H_

#include "elements/beam3d_kinematics.h"

namespace fesa {

/// Independent spatial pose, physical velocity, or acceleration derivative.
enum class BeamTimeDerivative {
  kNone,
  kDisplacement,
  kVelocity,
  kAcceleration
};

/// Thread-confined owner of reusable buffers; no model or input is retained.
/// Time coefficient2 equals one half of the instantaneous second derivative.
/// q/V/A are independent spatial coordinates; eta uses m, m/s and m/s².
/// @see docs/Linear static/b33_inertia_kernel_contract.md.
class BeamTimeKinematicWorkspace final {
 public:
  /// Creates fixed storage. Caller deletes *out; failure sets *out=nullptr.
  /// out is required, must not already own storage; d is optional.
  static StatusCode Create(BeamTimeKinematicWorkspace** out, Diagnostic* d);
  /// Frees owned storage and invalidates all borrowed coefficient pointers.
  ~BeamTimeKinematicWorkspace();
  /// Workspace identity and borrowed pointers cannot be copied or moved.
  BeamTimeKinematicWorkspace(const BeamTimeKinematicWorkspace&) = delete;
  BeamTimeKinematicWorkspace& operator=(const BeamTimeKinematicWorkspace&) =
      delete;
  BeamTimeKinematicWorkspace(BeamTimeKinematicWorkspace&&) = delete;
  BeamTimeKinematicWorkspace& operator=(BeamTimeKinematicWorkspace&&) = delete;

  /// Overwrites coefficients0..2 at xi in [0,1] without allocating.
  /// pose uses the K1 units, triad tolerance, domain and branch checks.
  /// pose.direction is required and validated but does not seed this kernel.
  /// velocity14 and acceleration14 are finite global spatial [xa3,wa3,xb3,
  /// wb3,eta2] data, in m/s, rad/s and m/s², rad/s² respectively.
  /// All inputs are required, read-only and borrowed for the call. failure is
  /// a required output pointer, overwritten with the typed evaluation failure.
  /// derivative selects all14 independent identity seeds; kNone skips them.
  /// Invalid inputs return kInvalidArgument; geometry/branch errors return
  /// kEquilibriumFailure with the typed failure. Backend/resource diagnostics
  /// are preserved. failure is cleared on entry. Any failure invalidates all
  /// getters; the next call recomputes every coefficient from scratch.
  StatusCode Evaluate(const BeamKinematicInput* pose, const Vector* velocity14,
                      const Vector* acceleration14,
                      BeamTimeDerivative derivative, double xi,
                      TrialEvaluationFailure* failure, Diagnostic* d);

  /// Borrows coefficient order0..2, shape3x1 (strain4x1, stretch1x1).
  /// Fields retain K1 units divided by s^order. Null before success, after any
  /// failure, or for invalid field/order. Pointers expire at next Evaluate.
  const Matrix* Coefficient(BeamKinematicField field, int order) const;
  /// Borrows coefficient Jacobian3x14 (strain4x14, stretch1x14), same lifetime.
  /// Column units follow the selected independent q/V/A coordinates.
  /// Null also when the last successful request selected kNone.
  const Matrix* Derivative(BeamKinematicField field, int order) const;

 private:
  /// Sole-owned fixed numerical buffers, released by the destructor.
  struct Storage;
  /// Factory-only constructor; Create performs fallible storage allocation.
  BeamTimeKinematicWorkspace();
  Storage* storage_ = nullptr;  ///< Sole owner; no borrowed model state.
};

}  // namespace fesa

#endif  // FESA_ELEMENTS_BEAM3D_TIME_KINEMATICS_H_
