/// @file
/// @brief Common force-first full Newton loop, independent of element families.
#ifndef FESA_ANALYSIS_NEWTON_RAPHSON_H_
#define FESA_ANALYSIS_NEWTON_RAPHSON_H_
#include "fesa/analysis/linear_solver.h"
#include "fesa/analysis/newton_controls.h"
#include "fesa/analysis/newton_report.h"
#include "fesa/analysis/nonlinear_problem.h"
#include "fesa/io/control_settings.h"
namespace fesa {
/// Automatic retains IR/IC prediction and ISc protection. DIRECT omits slow
/// prediction and uses IC as actual correction cap; convergence is tested
/// first at the cap. Static callers retain the default automatic behavior.
enum class NewtonExecution { kAutomatic, kDirect };
/// Analysis owns this stateless loop. No borrowed input survives Solve.
/// @see docs/MITC4/newton_execution_contract.md for equations and sequencing.
class NewtonRaphson final {
 public:
  /// Caller owns success; failure nulls output. No numerical buffer allocation.
  static StatusCode Create(NewtonRaphson** out, Diagnostic* d);
  ~NewtonRaphson() = default;
  NewtonRaphson(const NewtonRaphson&) = delete;
  NewtonRaphson& operator=(const NewtonRaphson&) = delete;
  NewtonRaphson(NewtonRaphson&&) = delete;
  NewtonRaphson& operator=(NewtonRaphson&&) = delete;
  /// Borrowed required problem/solver/controls/counts/report, optional
  /// progress/d. Resets report/history per attempt. Evaluates i0 without fake
  /// correction. Active equations honor RequiresInitialCorrection before
  /// accepting i0; scalar criteria remain unchanged. Solver called directly.
  /// Converged trial is retained. No increment size, recovery,
  /// output/flush/commit/state discard. Recoverable outcomes return
  /// kEquilibriumFailure; fatal status preserved. Logging failure forces abort;
  /// if numeric failure already exists it remains primary and
  /// Info.logging_failed records the independent abort reason.
  StatusCode Solve(NonlinearProblem* problem, LinearSolver* solver,
                   const ControlSettings* controls,
                   const NewtonControls* counts, NewtonReport* report,
                   NewtonProgress* progress, Diagnostic* d,
                   NewtonExecution execution = NewtonExecution::kAutomatic);

 private:
  NewtonRaphson() = default;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_NEWTON_RAPHSON_H_
