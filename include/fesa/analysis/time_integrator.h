/// @file
/// @brief Borrowed translational kinematics; no state ownership or commit.
#ifndef FESA_ANALYSIS_TIME_INTEGRATOR_H_
#define FESA_ANALYSIS_TIME_INTEGRATOR_H_

#include "fesa/math/matrix.h"
#include "fesa/math/vector.h"

namespace fesa {
/// Call-scoped read-only Cartesian translations, identical component order.
/// The caller owns all three distinct n>0 vectors, in L, L/T and L/T^2.
/// Caller may apply the same Cartesian algebra to linear shell small-rotation
/// arrays with consistent per-component units. Finite rotation uses the
/// separate rotational role below.
/// HHT acceleration is the algorithmic state used by H3, not a separately
/// reconstructed physical end-equilibrium acceleration. Initial physical
/// acceleration seeds this state; continuous steps inherit it unchanged.
struct TranslationalStateView {
  const Vector* displacement = nullptr;  ///< Borrowed committed u.
  const Vector* velocity = nullptr;      ///< Borrowed committed v.
  const Vector* acceleration = nullptr;  ///< Borrowed committed a.
};
/// Call-scoped writable candidate; caller owns three same-sized vectors.
/// All six state/candidate vectors must be distinct. valid is kinematic
/// validity only, never convergence, acceptance, recovery or commit status.
struct TranslationalCandidate {
  Vector* displacement = nullptr;  ///< Borrowed trial u, L.
  Vector* velocity = nullptr;      ///< Borrowed trial v, L/T.
  Vector* acceleration = nullptr;  ///< Borrowed trial a, L/T^2.
  bool valid = false;              ///< Set true only on complete success.
};
/// Scalar identity multipliers for Cartesian translation linearization.
struct TranslationalLinearization {
  double acceleration_per_displacement = 0;  ///< da/du, 1/T^2.
  double velocity_per_displacement = 0;      ///< dv/du, 1/T.
};
/// Call-scoped finite rotation state, borrowed n-by-4 unit wxyz orientations
/// and distinct3n global spatial omega/eta, rad/T and rad/T^2. n>0.
/// HHT eta is the R1 algorithmic angular acceleration. Substituting it into
/// physical inertia does not define a separate physical acceleration recovery.
struct RotationalStateView {
  const Matrix* orientations = nullptr;  ///< Borrowed committed pose.
  const Vector* velocity = nullptr;      ///< Borrowed committed omega.
  const Vector* acceleration = nullptr;  ///< Borrowed committed eta.
};
/// Caller-owned candidate with same shapes; all state/candidate objects are
/// distinct. Quaternion cover sign is retained; no body-spin history here.
struct RotationalCandidate {
  Matrix* orientations =
      nullptr;                 ///< Borrowed predicted/caller-corrected pose.
  Vector* velocity = nullptr;  ///< Borrowed trial omega.
  Vector* acceleration = nullptr;  ///< Borrowed trial eta.
  bool valid = false;              ///< Kinematics only, never acceptance.
};
/// Borrowed disjoint n-by-9 matrices: each node row holds a row-major3x3 block
/// with global spatial correction columns. Units1/T and1/T^2. Caller owns.
struct RotationalLinearization {
  Matrix* velocity_per_rotation = nullptr;      ///< d(omega)/d(delta_w).
  Matrix* acceleration_per_rotation = nullptr;  ///< d(eta)/d(delta_w).
};
/// Replaceable Cartesian/finite-rotation integration role; no state ownership.
/// Numerical scratch is not retained across calls. Inputs must stay alive
/// through each call; committed vectors remain unchanged on all outcomes.
class TimeIntegrator {
 public:
  virtual ~TimeIntegrator() = default;
  TimeIntegrator(const TimeIntegrator&) = delete;
  TimeIntegrator& operator=(const TimeIntegrator&) = delete;
  TimeIntegrator(TimeIntegrator&&) = delete;
  TimeIntegrator& operator=(TimeIntegrator&&) = delete;
  /// Writes predictor u/v and zero a using current dt [T]. Sets valid=false
  /// on entry. Argument errors preserve vectors; numerical/backend failures
  /// may partially write candidate vectors, which must then be discarded.
  /// Success sets valid=true. Null Diagnostic is allowed on every method.
  virtual StatusCode PrepareKinematics(const TranslationalStateView* committed,
                                       double dt,
                                       TranslationalCandidate* candidate,
                                       Diagnostic* d) const = 0;
  /// Uses caller-solved end displacement, preserving it; overwrites v/a.
  /// Same lifetime/failure contract as Prepare. No prior valid flag required.
  virtual StatusCode UpdateKinematics(const TranslationalStateView* committed,
                                      double dt,
                                      TranslationalCandidate* candidate,
                                      Diagnostic* d) const = 0;
  /// Recomputes current-dt derivatives; failure preserves both output scalars.
  virtual StatusCode EvaluateKinematicLinearization(
      double dt, TranslationalLinearization* out, Diagnostic* d) const = 0;
  /// R1 predictor Exp(theta_p)*q_n plus spatial omega/eta. All buffers
  /// borrowed; state unchanged. Sets valid=false and required branch=false on
  /// entry. Shape/alias errors preserve data; subsequent failures invalidate
  /// candidate. Predictors crossing the existing pi chart fail with
  /// branch=true. No commit.
  virtual StatusCode PrepareRotationalKinematics(
      const RotationalStateView* committed, double dt,
      RotationalCandidate* candidate, bool* branch, Diagnostic* d) const = 0;
  /// Preserves caller-corrected pose and overwrites omega/eta by R1. Same
  /// lifetime/failure contract as predictor. Caller guards actual correction
  /// paths through RotationState; no principal-UR/history reconstruction.
  virtual StatusCode UpdateRotationalKinematics(
      const RotationalStateView* committed, double dt,
      RotationalCandidate* candidate, bool* branch, Diagnostic* d) const = 0;
  /// R2 derivatives in per-node row-major3x3 blocks; no global dense matrix.
  /// Borrows views/out; preserves every state/candidate value including valid.
  /// Clears required branch; chart failure sets it true. Validation preserves
  /// out, later failure invalidates both blocks. No retained scratch.
  virtual StatusCode EvaluateRotationalLinearization(
      const RotationalStateView* committed,
      const RotationalCandidate* candidate, double dt,
      RotationalLinearization* out, bool* branch, Diagnostic* d) const = 0;

 protected:
  TimeIntegrator() = default;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_TIME_INTEGRATOR_H_
