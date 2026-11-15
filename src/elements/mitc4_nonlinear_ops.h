/// @file
/// @brief Shared call-local finite MITC4 basis checks and covariant transforms.
#ifndef FESA_ELEMENTS_MITC4_NONLINEAR_OPS_H_
#define FESA_ELEMENTS_MITC4_NONLINEAR_OPS_H_
#include "fesa/math/matrix.h"
#include "fesa/model/element_trial.h"
namespace fesa {
/// Owns exclusive scratch for repeated volume checks within one Trial or
/// Rayleigh evaluation. No input/state pointers or geometry values retained.
/// Destroy on return, including failure; never share between workers.
class Mitc4VolumeWorkspace {
 public:
  Mitc4VolumeWorkspace() = default;
  ~Mitc4VolumeWorkspace();
  Mitc4VolumeWorkspace(const Mitc4VolumeWorkspace&) = delete;
  Mitc4VolumeWorkspace& operator=(const Mitc4VolumeWorkspace&) = delete;
  Mitc4VolumeWorkspace(Mitc4VolumeWorkspace&&) = delete;
  Mitc4VolumeWorkspace& operator=(Mitc4VolumeWorkspace&&) = delete;
  /// Allocates four Vector3 and one Matrix3x3 once (168-byte payload).
  /// Failure frees partial storage and leaves an empty, retryable workspace.
  StatusCode Initialize(Diagnostic* d);
  /// Same numerical and failure contract as ValidateMitc4VolumeBasis below.
  /// Requires successful Initialize. Borrows inputs/scalar outputs only for
  /// this call; preserves them on failure. Scratch is invalid after failure,
  /// fully overwritten on retry. Performs no allocations or input mutation.
  StatusCode Validate(const Vector* const* initial,
                      const Vector* const* current,
                      double* reference_determinant, double* jacobian,
                      TrialEvaluationFailure* failure, Diagnostic* d);
 private:
  /// Releases all owned scratch; leaves an empty, retryable workspace.
  void Clear();
  Vector* cross_ = nullptr;
  Vector* normalized_[3] = {};
  Matrix* skew_ = nullptr;
};
/// Borrows three initial/current global basis vectors of size3, ordered
/// r,s,zeta. Reference determinant has L^3 units, J is dimensionless. Outputs
/// required, distinct and preserved on failure. Required failure cleared on
/// entry. Initial nonpositive volume is fatal input; current normalized
/// triple<=256eps or nonpositive/nonfinite J is typed geometry. Backend
/// failures stay fatal. Inputs may alias each other; no mutation, allocation
/// retained or cache.
StatusCode ValidateMitc4VolumeBasis(const Vector* const* initial,
                                    const Vector* const* current,
                                    double* reference_determinant,
                                    double* jacobian,
                                    TrialEvaluationFailure* failure,
                                    Diagnostic* d);
/// External-math cross product of finite global vectors3; out3 distinct from
/// inputs. Validation preserves out; later failure invalidates. No retention.
StatusCode CrossMitc4Vectors(const Vector* a, const Vector* b, Vector* out,
                             Diagnostic* d);
/// Writes the same P directly from an already solved H3x3, avoiding a second
/// initial-basis LU when recovery also needs F. Borrowed inputs, full out5x5.
/// Invalid input preserves out; later failure invalidates. No retention.
StatusCode BuildMitc4NonlinearTransformFromH(const Matrix* h, Matrix* out,
                                             Diagnostic* d);
/// Writes full5x5 P mapping [Err,Ess,Ers,Erz,Esz] to local engineering5
/// [aa,bb,ab,bn,na]. initial_basis3x3 stores global G vectors as rows;
/// frame_transpose3x3 is Q0^T. Uses H=G_columns^-1 Q0^T, excludes Ezz.
/// Borrowed finite inputs, separate out; validation preserves out, backend
/// failure invalidates out. No state retained; initial LU singular is fatal.
StatusCode BuildMitc4NonlinearCartesianTransform(const Matrix* initial_basis,
                                                 const Matrix* frame_transpose,
                                                 Matrix* out, Diagnostic* d);
/// Fully overwrites selector5x7 at finite natural r,s in[-1,1]. Raw order:
/// rr,ss,rs,rz(0,-1),rz(0,1),sz(-1,0),sz(1,0), all same zeta.
/// Borrowed caller-owned out; invalid arguments preserve, later failure
/// invalidates. No numerical buffers retained.
StatusCode BuildMitc4NonlinearTyingSelector(double r, double s, Matrix* out,
                                            Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_MITC4_NONLINEAR_OPS_H_
