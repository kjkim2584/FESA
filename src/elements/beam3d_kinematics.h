/// @file
/// @brief Exact selected-D beam point fields and ordered spatial variations.
#ifndef FESA_ELEMENTS_BEAM3D_KINEMATICS_H_
#define FESA_ELEMENTS_BEAM3D_KINEMATICS_H_

#include "fesa/base/diagnostic.h"
#include "fesa/math/matrix.h"
#include "fesa/math/vector.h"
#include "fesa/model/element_trial.h"

namespace fesa {

/// Global xyz position (m), unitless triad, and [epsilon,K1,K2,tau] packet.
/// epsilon=(r' dot r'-1)/2 is dimensionless; K1=-n2 dot r'', K2=n1 dot r'',
/// tau=n2 dot n1' have units 1/m. Primes use initial arc length S.
/// Stretch lambda=norm(r') is dimensionless and independent of Green rounding.
enum class BeamKinematicField {
  kPosition,
  kTangent,
  kNormal1,
  kNormal2,
  kStrain,
  kStretch
};

/// Requested spatial blocks; full preserves the original ordered mixed graph.
enum class BeamKinematicDerivative { kValue, kFirst, kFull };

/// Inputs borrowed only during Evaluate; no input is changed or retained.
struct BeamKinematicInput {
  double length = 0;  ///< Initial straight length, m, strictly positive.
  const Matrix* positions = nullptr;  ///< 2x3 global positions, a/b rows, m.
  /// Optional2x3 positions relative to the SAME constant origin, (X-O)+u.
  /// Borrowed; same spatial/time variations as positions. Used only for S
  /// derivatives, never absolute kPosition. Null uses positions. Provider
  /// guarantees the common fixed origin; shape/finiteness are validated.
  const Matrix* relative_positions = nullptr;
  const Matrix* tangents = nullptr;  ///< 2x3 current unit tangents, global xyz.
  const Matrix* normals =
      nullptr;  ///< 2x3 first unit normals, orthogonal to T.
  const Vector* internal_position = nullptr;  ///< eta_a/eta_b, length2, m.
  /// Spatial direction [xa3,wa3,xb3,wb3,eta2], length14; x/eta in m, w in rad.
  const Vector* direction = nullptr;
};

/// Thread-confined owner of reusable math buffers, independent of model state.
/// Endpoint rotations use Exp(virtual) Exp(direction) R, with fixed global
/// virtual basis. Mixed is the directional derivative of that spatial B;
/// no coordinate-Hessian half-skew correction is applied.
class BeamKinematicWorkspace final {
 public:
  /// Creates zeroed storage; caller deletes *out. Failure sets *out=nullptr.
  /// out is required and must not already own an allocation; d is optional.
  static StatusCode Create(BeamKinematicWorkspace** out, Diagnostic* d);
  /// Releases every owned buffer; invalidates all borrowed field pointers.
  ~BeamKinematicWorkspace();
  /// Owning workspace identity and borrowed pointers cannot be copied/moved.
  BeamKinematicWorkspace(const BeamKinematicWorkspace&) = delete;
  BeamKinematicWorkspace& operator=(const BeamKinematicWorkspace&) = delete;
  BeamKinematicWorkspace(BeamKinematicWorkspace&&) = delete;
  BeamKinematicWorkspace& operator=(BeamKinematicWorkspace&&) = delete;

  /// Overwrites all point fields at finite xi in [0,1]. All inputs and failure
  /// are required. Input shapes, finiteness, and triad unit/orthogonality are
  /// validated; inputs are never normalized to repair them. Unit/orthogonality
  /// tolerance is 64*epsilon, a roundoff allowance for three-term dot products.
  /// Transport denominators <=64*epsilon and |pi-|phi||<=64*epsilon identify
  /// numerically unresolved rotation branches, not a physical rotation cap.
  /// Nonpositive endpoint stretches/zero centerline derivative return geometry
  /// failure; singular transport/pi twist return rotation-branch failure.
  /// Invalid arguments return kInvalidArgument. Resource/backend errors retain
  /// original diagnostics. failure is cleared on entry when nonnull. Any
  /// failure invalidates all getters; the next call fully recomputes fields.
  /// Value/first requests omit later blocks; their getters return null. Full
  /// is the legacy default. Storage is overwritten and inputs are preserved.
  StatusCode Evaluate(const BeamKinematicInput* input, double xi,
                      TrialEvaluationFailure* failure, Diagnostic* d,
                      BeamKinematicDerivative derivative =
                          BeamKinematicDerivative::kFull);

  /// Invalidates all borrowed getters and exact-key reuse without changing
  /// model/state or reallocating owned buffers. Only the owning worker calls.
  void Invalidate();

  /// Borrows value (3x1, strain4x1, stretch1x1).
  /// Null before success or after any failure, including stretch.
  /// All getter pointers expire at the next Evaluate or destruction.
  const Matrix* Value(BeamKinematicField field) const;
  /// Borrows spatial first variation (3x14, strain4x14, stretch1x14).
  /// Same lifetime; column units follow the spatial input coordinates.
  const Matrix* Variation(BeamKinematicField field) const;
  /// Borrows directional variation (3x1, strain4x1, stretch1x1), same lifetime.
  const Matrix* Directional(BeamKinematicField field) const;
  /// Borrows ordered directional/virtual (3x14, strain4x14, stretch1x14).
  /// Same lifetime; column units follow the spatial input coordinates.
  const Matrix* Mixed(BeamKinematicField field) const;

 private:
  /// Private owned numerical buffers, released on destruction.
  struct Storage;
  /// Factory-only constructor; fallible buffer creation happens in Create.
  BeamKinematicWorkspace();
  Storage* storage_ = nullptr;  ///< Sole owner; no model pointers retained.
};

}  // namespace fesa

#endif  // FESA_ELEMENTS_BEAM3D_KINEMATICS_H_
