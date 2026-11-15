/// @file
/// @brief First-step trial equilibrium diagnostics; never committed history.
#ifndef FESA_ANALYSIS_RESIDUAL_REPORT_H_
#define FESA_ANALYSIS_RESIDUAL_REPORT_H_
#include "fesa/math/vector.h"
namespace fesa {
/// S1 decision branch; no-target fields do not take part in convergence.
enum class ResidualBranch { kNoDofs, kZeroFlux, kLinear, kGeneral };
/// Force or moment field summary in global coordinates and model units.
struct FieldResidual {
  int dof_count = 0;  ///< Free components checked, not contribution count.
  int contribution_count = 0;  ///< Six per beam plus individual active CLOADs.
  int maximum_node = 0;        ///< Original ID at largest free residual.
  int maximum_dof = 0;         ///< Global DOF 1..6 at maximum.
  double maximum = 0;          ///< Absolute free residual, force or moment.
  double l2 = 0;               ///< Field-only Euclidean residual norm.
  double rms = 0;              ///< l2/sqrt(dof_count), zero with no targets.
  double spatial_mean = 0;     ///< Sum absolute unassembled flux/count.
  double raw_mean = 0;         ///< Trial raw single-increment mean.
  double used_mean = 0;        ///< Simultaneously converted decision scale.
  double maximum_flux = 0;     ///< Trial maximum absolute flux contribution.
  double phi = 0;              ///< maximum/used_mean.
  double increment_max =
      0;  ///< All participating field components, including
          ///< prescribed DOFs; length/radians. Free nodes excluded.
  double correction_max = 0;    ///< Actual first correction, not a fake zero.
  double residual_limit = 0;    ///< Applied normalized residual criterion.
  double correction_limit = 0;  ///< Applied absolute correction criterion.
  bool converted = false;       ///< Used mean differs due to field conversion.
  bool passed = true;           ///< No-target fields pass vacuously.
  ResidualBranch branch = ResidualBranch::kNoDofs;
  int maximum_element = 0;  ///< Beam ID for an internal maximum; otherwise 0.
  int maximum_internal_end = -1;  ///< Connectivity end 0/1; otherwise -1.
};
/// Caller owns this report, its two summaries and one full reaction vector.
/// A failed computation invalidates report data; an equilibrium failure leaves
/// complete diagnostics available. No model/state/history mutation occurs.
class ResidualReport {
 public:
  /// Allocates for a positive global count of physical 3/6 DOF layouts; failure
  /// nulls *out.
  static StatusCode Create(int global_count, ResidualReport** out,
                           Diagnostic* d);
  ~ResidualReport();
  ResidualReport(const ResidualReport&) = delete;
  ResidualReport& operator=(const ResidualReport&) = delete;
  ResidualReport(ResidualReport&&) = delete;
  ResidualReport& operator=(ResidualReport&&) = delete;
  /// Borrowed summaries/reactions valid until this is destroyed or reevaluated.
  const FieldResidual* Displacement() const;
  const FieldResidual* Rotation() const;
  const Vector* Reaction() const;
  /// Initial arithmetic mean beam length; zero for an element-free model.
  double CharacteristicLength() const;

 private:
  friend class EquilibriumEvaluator;
  ResidualReport() = default;
  FieldResidual* displacement_ = nullptr;  ///< Owned.
  FieldResidual* rotation_ = nullptr;      ///< Owned.
  Vector* reaction_ = nullptr;  ///< Owned, [RF,RM] for all original node rows.
  double length_ = 0;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_RESIDUAL_REPORT_H_
