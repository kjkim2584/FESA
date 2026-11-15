/// @file
/// @brief Call-local read-only physical element pose in connectivity order.
#ifndef FESA_ANALYSIS_NONLINEAR_ELEMENT_POSE_H_
#define FESA_ANALYSIS_NONLINEAR_ELEMENT_POSE_H_

#include "fesa/analysis/analysis_state.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/model/element.h"
#include "fesa/model/element_trial.h"
#include "fesa/model/incremental_solid_trial_context.h"

namespace fesa {
/// Owns one element's copied pose/view and actual physical global indices.
/// Incremental trials borrow committed IP stress from AnalysisState until this
/// pose is deleted; state must outlive the pose and not commit concurrently.
/// No model/mapping pointer is retained. Create/destroy inside
/// one element evaluation, never cache across iterations/increments.
class NonlinearElementPose final {
 public:
  /// Borrows matching nonlinear state/mapping/element during creation only.
  /// committed=false requires a prepared candidate and copies current
  /// u/q/body-spin; its trial reference q is the last committed q.
  /// committed=true reads the committed u/q/body-spin even without a candidate,
  /// using that same q as trial reference (zero extra body increment). No
  /// outputUR reconstruction. Supports actual3/6DOF elements; threeDOF has no
  /// rotational buffers. Caller deletes success. Failure nulls *out, cleans
  /// partial storage, and preserves all inputs, state/result validity and tags.
  /// No normalization, geometry, element force, tangent or official result
  /// computation.
  /// Beam trials additionally own a copied eta2 selected by original element
  /// ID from the matching current/committed state. GlobalDofs/DofCount stay
  /// physical12. Missing/mismatched internal mapping/state fails before
  /// publish. Optional geometric-stress selector is copied into an incremental
  /// trial context only; no state mutation or cross-call retention. Default
  /// false.
  static StatusCode Create(const Element* element, const DofManager* dofs,
                           const AnalysisState* state, bool committed,
                           NonlinearElementPose** out, Diagnostic* d,
                           IncrementalElasticWorkspace* scratch = nullptr,
                           bool use_committed_geometric_stress = false);
  ~NonlinearElementPose();
  NonlinearElementPose(const NonlinearElementPose&) = delete;
  NonlinearElementPose& operator=(const NonlinearElementPose&) = delete;
  NonlinearElementPose(NonlinearElementPose&&) = delete;
  NonlinearElementPose& operator=(NonlinearElementPose&&) = delete;
  /// Owned view borrowed until deletion; translations3n/global L, quaternion
  /// n-by-4 initial-to-current unit wxyz, body-spin3n/rad. No mutable access.
  const ElementTrialState* Trial() const;
  /// Owned connectivity/component-ordered physical mapping, length DofCount.
  /// Borrowed read-only until deletion; never free it or treat as equations.
  const int* GlobalDofs() const;
  int DofCount() const;  ///< Number of actual local physical components.

 private:
  NonlinearElementPose() = default;
  /// Fills unpublished owner; caller cleans this on every failure.
  StatusCode Initialize(const Element* element, const DofManager* dofs,
                        const AnalysisState* state, bool committed,
                        Diagnostic* d, IncrementalElasticWorkspace* scratch,
                        bool use_committed_geometric_stress);
  int dof_count_ = 0;
  int* global_dofs_ = nullptr;          ///< Owned physical indices.
  ElementTrialState* trial_ = nullptr;  ///< Owned view of the stores below.
  Vector* translations_ = nullptr;      ///< Owned3n, global L.
  Vector* internal_position_ = nullptr;  ///< Owned optional Beam eta2, L.
  Vector* spin_ = nullptr;  ///< Owned optional3n, material body-spin rad.
  Matrix* orientations_ = nullptr;  ///< Owned optional n-by-4 current wxyz.
  Matrix* committed_ = nullptr;     ///< Owned optional n-by-4 reference wxyz.
  Vector* committed_translations_ = nullptr;  ///< Owned optional global3n L.
  Matrix* candidate_stress_ = nullptr;  ///< Owned optional ng×6 stress output.
  IncrementalSolidTrialContext* incremental_ = nullptr;  ///< Owned call view.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_NONLINEAR_ELEMENT_POSE_H_
