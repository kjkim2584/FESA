/// @file
/// @brief Call-local current/committed finite-Shell physical force assembly.
#ifndef FESA_ANALYSIS_FINITE_SHELL_DYNAMIC_FORCES_H_
#define FESA_ANALYSIS_FINITE_SHELL_DYNAMIC_FORCES_H_

#include "fesa/analysis/analysis_model.h"
#include "fesa/analysis/analysis_state.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/analysis/nonlinear_problem.h"
#include "fesa/elements/element_workspace.h"

namespace fesa {
/// Owns global and unassembled physical E0 forces for one evaluation only.
/// Retains no model/mapping/state pointer; never cache across iterations/steps.
class FiniteShellDynamicForces final {
 public:
  /// Coordinator-only. Borrows a matching Domain-layout6N state/model/mapping,
  /// supported positive-density MITC4 elements and mass mode through joined
  /// execution. Uses configured thread count and exclusive call-local scratch.
  /// committed=false needs prepared q/body-spin and staged motion; true reads
  /// committed pose/motion and does not require a candidate. No time weighting,
  /// external/structural force, convergence, recovery, output or state
  /// mutation. Caller deletes success. Required failure is cleared/typed via
  /// common post-join selection. Every error nulls *out and cleans partial
  /// storage; input/physical state/result tags remain unchanged. @see
  /// docs/HHT/implementation_contract.md for physical units/order/ownership.
  static StatusCode Create(const AnalysisModel* model, const DofManager* dofs,
                           const AnalysisState* state, MassMode mode,
                           bool committed, FiniteShellDynamicForces** out,
                           NonlinearFailure* failure, Diagnostic* d);
  ~FiniteShellDynamicForces();
  FiniteShellDynamicForces(const FiniteShellDynamicForces&) = delete;
  FiniteShellDynamicForces& operator=(const FiniteShellDynamicForces&) = delete;
  FiniteShellDynamicForces(FiniteShellDynamicForces&&) = delete;
  FiniteShellDynamicForces& operator=(FiniteShellDynamicForces&&) = delete;
  /// Borrowed global6N [F1,F2,F3,M1,M2,M3] in node-ID/global axes, F/FL.
  /// Includes fixed components; excluded components are zero. Owner lifetime.
  const Vector* Inertia() const;
  const Vector* Damping() const;
  /// Borrowed [element ID,connectivity6n] force/moment entries, same units and
  /// source time. Preserves signs, zeros and multiplicity before global sum.
  const Vector* InertiaOrigins() const;
  const Vector* DampingOrigins() const;

 private:
  FiniteShellDynamicForces() = default;
  /// Creates unpublished stores; caller cleans this owner on any failure.
  StatusCode Initialize(const AnalysisModel* model, const DofManager* dofs,
                        const AnalysisState* state, MassMode mode,
                        bool committed, NonlinearFailure* failure,
                        Diagnostic* d);
  Vector* inertia_ = nullptr;          ///< Owned global physical fI.
  Vector* damping_ = nullptr;          ///< Owned global physical fD.
  Vector* inertia_origins_ = nullptr;  ///< Owned unassembled fI.
  Vector* damping_origins_ = nullptr;  ///< Owned unassembled fD.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_FINITE_SHELL_DYNAMIC_FORCES_H_
