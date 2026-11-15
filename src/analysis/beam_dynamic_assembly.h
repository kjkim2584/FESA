/// @file
/// @brief Call-local finite Beam force assembly and full-mass initialization.
#ifndef FESA_ANALYSIS_BEAM_DYNAMIC_ASSEMBLY_H_
#define FESA_ANALYSIS_BEAM_DYNAMIC_ASSEMBLY_H_
#include "fesa/analysis/analysis_model.h"
#include "fesa/analysis/analysis_state.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/analysis/linear_solver.h"
#include "fesa/analysis/nonlinear_problem.h"
#include "fesa/elements/element_workspace.h"
namespace fesa {
/// Sole owner of six force stores; retains no model, state or mapping.
class BeamDynamicForces final {
 public:
  /// Borrows matching Beam model/mapping/state through joined execution.
  /// Consistent mass only. Current reads complete prepared joint motion;
  /// committed reads joint committed pose/motion without requiring candidate.
  /// No mutation of state/results. Caller deletes success; failure nulls out
  /// and selects typed element failure after join. Units/order follow
  /// docs/Linear static/b33_dynamic_assembly_contract.md.
  /// Optional worker-exclusive contexts are borrowed until join; count must
  /// cover the configured team. Null/zero preserves call-local allocation.
  static StatusCode Create(const AnalysisModel* model, const DofManager* dofs,
                           const AnalysisState* state, MassMode mode,
                           bool committed, BeamDynamicForces** out,
                           NonlinearFailure* failure, Diagnostic* d,
                           const ElementWorkspace* const* workspaces = nullptr,
                           int workspace_count = 0);
  ~BeamDynamicForces();
  BeamDynamicForces(const BeamDynamicForces&) = delete;
  BeamDynamicForces& operator=(const BeamDynamicForces&) = delete;
  BeamDynamicForces(BeamDynamicForces&&) = delete;
  BeamDynamicForces& operator=(BeamDynamicForces&&) = delete;
  /// Borrowed physical6N global [F1,F2,F3,M1,M2,M3], including fixed rows.
  const Vector* Inertia() const;
  const Vector* Damping() const;
  /// Borrowed internal2E force [F], increasing original ID, ends a/b.
  const Vector* InternalInertia() const;
  const Vector* InternalDamping() const;
  /// Borrowed unassembled14E physical12/internal2, original element ID order.
  /// Signs, zeros and multiplicity preserved; all getters share owner lifetime.
  const Vector* InertiaOrigins() const;
  const Vector* DampingOrigins() const;

 private:
  BeamDynamicForces() = default;
  Vector* inertia_ = nullptr;           ///< Owned physical6N.
  Vector* damping_ = nullptr;           ///< Owned physical6N.
  Vector* internal_inertia_ = nullptr;  ///< Owned internal2E.
  Vector* internal_damping_ = nullptr;  ///< Owned internal2E.
  Vector* inertia_origins_ = nullptr;   ///< Owned origins14E.
  Vector* damping_origins_ = nullptr;   ///< Owned origins14E.
};
/// Solves full instantaneous MI_ff over free physical+internal equations.
/// Borrows prepared pose and explicit physical6N/internal2E velocities,
/// external6N, and solver. Fixed/excluded physical velocity must be zero.
/// Outputs physical6N/internal2E are disjoint from each other, inputs and all
/// state stores. On success fixed/excluded acceleration is zero. No state,
/// official-result, time/history or element mutation. Preflight preserves
/// outputs; later failure invalidates outputs. Uses staging until solve ends.
StatusCode ComputeInitialBeamAcceleration(
    const AnalysisModel* model, const DofManager* dofs,
    const AnalysisState* state, const Vector* physical_velocity,
    const Vector* internal_velocity, const Vector* external,
    LinearSolver* solver, Vector* physical_acceleration,
    Vector* internal_acceleration, NonlinearFailure* failure, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_BEAM_DYNAMIC_ASSEMBLY_H_
