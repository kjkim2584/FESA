/// @file
/// @brief Private current kinetic and spatial Cauchy Rayleigh beam metrics.
#ifndef FESA_ELEMENTS_BEAM3D_DYNAMIC_METRICS_H_
#define FESA_ELEMENTS_BEAM3D_DYNAMIC_METRICS_H_

#include "elements/beam3d_kinematics.h"
#include "fesa/base/diagnostic.h"
#include "fesa/math/matrix.h"
#include "fesa/model/beam_section.h"
#include "fesa/model/element_trial.h"

namespace fesa {

/// Complete preserves the legacy mandatory mass output. DirectionOnly
/// requires mass=null and a directional output; unconsumed M is not formed.
enum class BeamKineticMetricRequest {
  kMassAndOptionalDirection,
  kDirectionOnly
};

/// Thread-confined owner of fixed spin row scratch; no point/model retention.
/// Factory allocates once, then point evaluation does not allocate. This
/// object can be reused by the integrated metric and inertia consumers.
class BeamKineticMetricWorkspace final {
 public:
  /// Caller deletes success; failure sets *out=nullptr and releases partial
  /// allocation. Required out must not hold an existing allocation.
  static StatusCode Create(BeamKineticMetricWorkspace** out, Diagnostic* d);
  /// Releases owned rows; no borrowed point pointer is retained.
  ~BeamKineticMetricWorkspace();
  BeamKineticMetricWorkspace(const BeamKineticMetricWorkspace&) = delete;
  BeamKineticMetricWorkspace& operator=(const BeamKineticMetricWorkspace&) =
      delete;
  BeamKineticMetricWorkspace(BeamKineticMetricWorkspace&&) = delete;
  BeamKineticMetricWorkspace& operator=(BeamKineticMetricWorkspace&&) = delete;

 private:
  /// Factory-only construction; caller never supplies scratch dimensions.
  BeamKineticMetricWorkspace() = default;
  friend StatusCode ComputeBeamPointKineticMetric(const BeamKinematicWorkspace*,
                                                  double, double,
                                                  BeamKineticMetricWorkspace*,
                                                  Matrix*, Matrix*, Diagnostic*,
                                                  BeamKineticMetricRequest);
  Matrix* spin_ = nullptr;  ///< Owned1x14 n2^T Jn1.
  Matrix* spin_directional_ = nullptr;  ///< Owned1x14 exact Dspin.
};

/// Overwrites mass14x14 and optional directional14x14 at a successful point.
/// All pointers are borrowed for the call; point is never changed. rho_area
/// is positive finite mass/length; polar_line_inertia=rho*(I11+I22) is positive
/// finite mass*length, independent of the stiffness torsion constant J.
/// Outputs must be distinct and cannot alias any point getter. Current
/// getters, valid coefficients, and matching shapes are required. Argument
/// failures preserve outputs; numerical failures invalidate partial outputs.
/// Ordered mixed spatial variations are used without Hessian conversion.
/// Required scratch is borrowed, overwritten, and reusable after failure;
/// point evaluation performs no allocation. Backend diagnostics propagate.
/// @see docs/Linear static/b33_dynamic_metric_contract.md.
StatusCode ComputeBeamPointKineticMetric(
    const BeamKinematicWorkspace* point, double rho_area,
    double polar_line_inertia, BeamKineticMetricWorkspace* scratch,
    Matrix* mass, Matrix* mass_directional, Diagnostic* d,
    BeamKineticMetricRequest request =
        BeamKineticMetricRequest::kMassAndOptionalDirection);

/// Integrates current mass and Rayleigh damping on initial length with Gauss4.
/// Required pose/section/workspace/failure are borrowed; no model is retained
/// or changed. pose.direction is required length14 even when derivatives are
/// omitted. Coordinates are global spatial [xa3,wa3,xb3,wb3,eta2] in m/rad.
/// M/C are required14x14; DM/DC are either both null or both14x14. All output
/// matrices are distinct and disjoint from pose and workspace getter storage.
/// Values overwrite outputs. Directional values follow pose.direction
/// with the physical spatial basis fixed. Damping uses the explicit kStretch
/// field and Cauchy rate metric; Rayleigh alpha is unrelated to HHT alpha.
/// Preflight argument/material errors preserve outputs. Later failures make
/// partial outputs invalid; a K1 failure invalidates the current point; callers
/// discard these buffers. K1 geometry/rotation typed failures and backend
/// diagnostics propagate. Missing density is kInputError. No HHT weighting,
/// force, global assembly, state commit, or reference acceptance is performed.
/// Caller owns workspace and outputs; next point evaluation expires getters.
/// @see docs/Linear static/b33_dynamic_metric_contract.md.
StatusCode EvaluateBeamDynamicMetrics(
    const BeamKinematicInput* pose, const BeamSection* section,
    BeamKinematicWorkspace* workspace, Matrix* mass, Matrix* mass_directional,
    Matrix* damping, Matrix* damping_directional,
    TrialEvaluationFailure* failure, Diagnostic* d);

/// Worker-owned Rayleigh buffers and four exact point workspaces. Base starts
/// a candidate; scalar material coefficients are captured, model pointers are
/// not retained. Only pose.direction may change through direction requests;
/// all other pose/section data must remain immutable until a fresh base.
class BeamRayleighWorkspace final {
 public:
  /// Caller deletes success; failure nulls required out and frees partial data.
  static StatusCode Create(BeamRayleighWorkspace** out, Diagnostic* d);
  /// Releases every owned point/work buffer; no model/state ownership.
  ~BeamRayleighWorkspace();
  BeamRayleighWorkspace(const BeamRayleighWorkspace&) = delete;
  BeamRayleighWorkspace& operator=(const BeamRayleighWorkspace&) = delete;
  BeamRayleighWorkspace(BeamRayleighWorkspace&&) = delete;
  BeamRayleighWorkspace& operator=(BeamRayleighWorkspace&&) = delete;
  /// Starts a fresh candidate and overwrites damping14x14. Optional current
  /// consistent mass is borrowed ONLY from the same ElementDynamic call, pose,
  /// section and unweighted physical14 coordinates. It is never retained.
  /// Null mass integrates M only when alpha!=0. Preflight shape/alias/material
  /// errors preserve output; later failures invalidate output and all four
  /// point caches. Required pose/section/damping/failure are borrowed; d
  /// optional. Typed failures pass.
  StatusCode EvaluateBase(const BeamKinematicInput* pose,
                          const BeamSection* section,
                          const Matrix* current_mass, Matrix* damping,
                          TrialEvaluationFailure* failure, Diagnostic* d);
  /// Overwrites DC14x14 for pose.direction, in existing point0..3 sum order.
  /// Requires successful base for the same immutable candidate. Only consumed
  /// DM/section terms are formed; C/M are not recomputed. Any failure ends the
  /// candidate and invalidates all four point caches; discard partial output.
  StatusCode EvaluateDirection(const BeamKinematicInput* pose, Matrix* dc,
                               TrialEvaluationFailure* failure, Diagnostic* d);

 private:
  struct Storage;
  BeamRayleighWorkspace() = default;
  Storage* storage_ = nullptr;  ///< Owned numeric-only candidate/work buffers.
};

}  // namespace fesa

#endif  // FESA_ELEMENTS_BEAM3D_DYNAMIC_METRICS_H_
