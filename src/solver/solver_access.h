/// @file
/// @brief Internal read-only execution inspection for integration verification.
#ifndef FESA_SOLVER_SOLVER_ACCESS_H_
#define FESA_SOLVER_SOLVER_ACCESS_H_
#include "fesa/solver/fesa_solver.h"
namespace fesa {
/// Internal orchestration observation, no public solver options or vendor API.
enum class SolverStage { kBeforePreflight, kModelWritten, kSolved, kCommitted };
/// Thread-local borrowed nonthrowing callback; tests reset before ending scope.
/// Hooks run only outside the WriteFrame-success/Commit atomic memory gate.
void SetSolverObserver(void (*callback)(SolverStage));
/// Internal controller notification; callback may configure adapter fault
/// seams.
void NotifySolverStage(SolverStage stage);
/// Returned objects are borrowed until next Run/destruction; no mutation seam.
class SolverAccess {
 public:
  /// Current execution state, possibly null before preparation.
  static const AnalysisState* State(const FESASolver* solver) {
    return solver->state_;
  }
  /// Thermal transactional state; null mechanically or before execution.
  static const ThermalState* ThermalExecutionState(const FESASolver* solver) {
    return solver->thermal_state_;
  }
  /// Last fully committed thermal step snapshot, null before a complete step.
  static const ThermalConditions* ThermalConditionsSnapshot(
      const FESASolver* solver) {
    return solver->thermal_conditions_;
  }
  /// Current definitions/results, possibly null on input failure.
  static const Domain* Model(const FESASolver* solver) {
    return solver->domain_;
  }
  /// Last completed step's conditions, borrowed until next Run/destruction.
  static const EffectiveConditions* Conditions(const FESASolver* solver) {
    return solver->conditions_;
  }
  /// Whole-step completion, distinct from the latest committed increment.
  static int LastCompletedStep(const FESASolver* solver) {
    return solver->last_completed_step_;
  }
  /// Report object survives until next Run; its fields obey Analysis lifetime.
  static const AnalysisReport* Report(const FESASolver* solver) {
    return solver->analysis_report_;
  }
};
}  // namespace fesa
#endif  // FESA_SOLVER_SOLVER_ACCESS_H_
